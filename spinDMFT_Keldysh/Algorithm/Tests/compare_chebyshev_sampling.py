#!/usr/bin/env python3
"""Compare dense/Chebyshev observables using HDF5 jackknife standard errors.

Retains every run and its log in the requested output directory. Distinct
seeds make the dense/Chebyshev comparisons independent. Grid entries are
correlated, so the report lists pointwise sigma distances without treating
them as independent observations or requiring every entry inside one sigma.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np


def scalar(value):
    return np.asarray(value).reshape(-1)[0].item()


def read(path):
    with h5py.File(path) as f:
        values = {key: np.asarray(f["results"][key]) for key in
                  ("Re_correlation", "Im_correlation", "Re_magnetization", "Im_magnetization")}
        errors = {key: np.asarray(f["runtimedata"][key + "_sample_stds"]) for key in values}
        attrs = f["runtimedata"].attrs
        diagnostics = {key: np.asarray(attrs[key]).reshape(-1).tolist() for key in
                       ("gaussian_polynomial_degrees", "gaussian_spectral_bounds",
                        "gaussian_covariance_approximation_errors", "gaussian_largest_factorization_dimensions")}
        assert all(np.isfinite(x).all() for x in (*values.values(), *errors.values()))
        return values, errors, diagnostics


def compare(a, b):
    result = {}
    for key in a[0]:
        difference = np.abs(a[0][key] - b[0][key])
        se = np.hypot(a[1][key], b[1][key])
        # Algebraic identities can have only roundoff-sized estimated errors.
        mask = se > 1e-10
        assert np.max(difference[~mask], initial=0) < 1e-9, (key, "deterministic identity differs")
        z = difference[mask] / se[mask]
        result[key] = {"points": int(z.size), "max_combined_sigma": float(z.max(initial=0)),
                       "rms_combined_sigma": float(np.sqrt(np.mean(z*z))) if z.size else 0.,
                       "fraction_within_1sigma": float(np.mean(z <= 1)) if z.size else 1.,
                       "fraction_within_2sigma": float(np.mean(z <= 2)) if z.size else 1.,
                       "fraction_within_3sigma": float(np.mean(z <= 3)) if z.size else 1.,
                       "max_absolute_difference": float(difference.max(initial=0))}
        # Conservative family-wise regression guard over many correlated
        # entries; this is not an assertion of one-sigma simultaneous coverage.
        assert result[key]["max_combined_sigma"] < 5., (key, result[key])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--samples", type=int, default=8192)
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base = ["--numImagTimeSteps=6", "--numRealTimeSteps=6", "--beta=0.8", "--Tmax=0.5",
            f"--numSamplesPerCore={args.samples}", "--numBlocks=64", "--cstype=D",
            "--Bname=z", "--Babs=0.7", "--dontrescaleB", "--iterlimit=1", "--JQ=0.4",
            "--samplingStrategy=independent"]
    cases = {
        "initial_bath": [],
        "feedback": ["--iterlimit=2", "--JL=-0.4"],
        "harmonic_response": ["--bath=harmonic", "--bathCoupling=0.3", "--bathOmega=1.1",
                              "--bathComponent=x", "--beta=1.2", "--Tmax=0.8"],
        "refined_response": ["--bath=harmonic", "--bathCoupling=0.25", "--bathOmega=1.1",
                             "--bathComponent=x", "--realTimeSubsteps=3"],
    }
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    report = {"samples_per_run": args.samples, "cases": {}}

    def run(case, mode, options, seed):
        work = output / f"{case}_{mode}"
        work.mkdir(exist_ok=False)
        keys = {x.split("=", 1)[0] for x in options}
        command = [str(executable), *[x for x in base if x.split("=", 1)[0] not in keys],
                   *options, f"--gaussianFactorization={'chebyshev' if mode.startswith('chebyshev') else 'dense'}",
                   f"--seed={seed}"]
        started = time.monotonic()
        with (work / "run.log").open("w") as log:
            subprocess.run(command, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=args.timeout)
        files = list(work.rglob("*.hdf5"))
        assert len(files) == 1, files
        data = read(files[0])
        metadata = {"file": str(files[0]), "seconds": time.monotonic()-started,
                    "command": command, "diagnostics": data[2]}
        if mode.startswith("chebyshev"):
            d = data[2]
            assert max(d["gaussian_largest_factorization_dimensions"]) == 0
            tolerance = next((float(x.split("=", 1)[1]) for x in options
                              if x.startswith("--gaussianChebyshevTolerance=")), 1e-3)
            assert all(e <= tolerance*s*(1+1e-10) for e, s in zip(
                d["gaussian_covariance_approximation_errors"], d["gaussian_spectral_bounds"]))
        return data, metadata

    for index, (case, options) in enumerate(cases.items()):
        dense, dense_meta = run(case, "dense", options, 734 + 10*index)
        new, new_meta = run(case, "chebyshev", options, 735 + 10*index)
        entry = {"dense": dense_meta, "chebyshev": new_meta, "comparison": compare(dense, new)}
        if case == "harmonic_response":
            tighter, tight_meta = run(case, "chebyshev_tighter",
                                      [*options, "--gaussianChebyshevTolerance=0.00025"], 735 + 10*index)
            entry["tighter"] = tight_meta
            entry["dense_vs_tighter"] = compare(dense, tighter)
            # Same latent draws isolate the polynomial approximation change;
            # use combined errors here as a conservative scale, not a paired test.
            entry["tolerance_change_over_combined_se"] = compare(new, tighter)
        report["cases"][case] = entry
        (output / "report.json").write_text(json.dumps(report, indent=2)+"\n")
        worst = max(x["max_combined_sigma"] for x in entry["comparison"].values())
        print(f"{case}: largest dense/Chebyshev difference {worst:.3f} combined standard errors; "
              f"degree {new_meta['diagnostics']['gaussian_polynomial_degrees']}", flush=True)
    print(f"Report: {output / 'report.json'}")


if __name__ == "__main__":
    main()
