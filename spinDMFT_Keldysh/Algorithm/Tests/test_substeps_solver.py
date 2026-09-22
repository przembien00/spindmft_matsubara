#!/usr/bin/env python3
"""FFT defaults, dense/FFT substep CLI, and solver/HDF5 regression tests.

Usage: python3 test_substeps_solver.py /absolute/path/executable_DOUBLE.out
Requires numpy, h5py and working MPI; output is isolated in temporary directories.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def scalar(value):
    return np.asarray(value).reshape(-1)[0]


def text(value):
    item = scalar(value)
    return item.decode() if isinstance(item, bytes) else str(item)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--endpoint-reference", type=Path, help="previous executable for exact endpoint regression")
    parser.add_argument("--reference", type=Path, help="previous executable for exact q=0 and q=1 regression")
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    common = ["--numImagTimeSteps=4", "--numRealTimeSteps=4", "--beta=0.8", "--Tmax=0.3",
              "--numSamplesPerCore=32", "--numBlocks=4", "--seed=734", "--cstype=D",
              "--Bname=z", "--Babs=0.7", "--dontrescaleB", "--iterlimit=1", "--JQ=0.1"]
    checked = []
    with tempfile.TemporaryDirectory(prefix="substeps-solver-tests-") as directory:
        root = Path(directory)

        def run(name, options, error=None):
            work = root / name
            work.mkdir()
            option_keys = {x.split("=", 1)[0] for x in options}
            base_options = [x for x in common if x.split("=", 1)[0] not in option_keys]
            result = subprocess.run([str(executable), *base_options, *options], cwd=work, env=env,
                                    text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
            files = list(work.rglob("*.hdf5"))
            if error:
                assert result.returncode != 0 and error in result.stdout, result.stdout
                assert not files, "invalid options wrote simulation data"
                checked.append(name)
                return
            assert result.returncode == 0 and len(files) == 1, result.stdout
            with h5py.File(files[0]) as f:
                p = f["parameters"].attrs
                factor = next((x.split("=", 1)[1] for x in options if x.startswith("--gaussianFactorization=")), "fft")
                assert text(p["gaussian_factorization"]) == factor
                q = int(scalar(p["real_time_substeps"]))
                assert float(scalar(p["delta_real_propagation_t"])) == float(scalar(p["delta_real_t"]))/max(1, q)
                assert text(p["propagator"]) == ("gauss-cf4" if q else "endpoint-cfet4")
                assert ("__prop=cf4" in files[0].name) == (q > 0)
                if q and factor == "fft":
                    assert "signed-frequency" in text(p["real_time_field_interpolation"])
                elif q == 0:
                    assert text(p["real_time_field_interpolation"]) == "native endpoints"
                else:
                    assert text(p["real_time_field_interpolation"]) == "four-point cubic interpolation"
                assert ("__substeps=" in files[0].name) == (q > 1)
                if not any(x.startswith("--fftCrossFrequencyCutoff=") for x in options):
                    assert scalar(p["fft_cross_frequency_cutoff"]) == -1
                    np.testing.assert_array_equal(f["runtimedata"].attrs["gaussian_covariance_approximation_errors"], 0.)
                arrays = {k: f["results"][k][()] for k in
                          ("Re_correlation", "Im_correlation", "Re_magnetization", "Im_magnetization")}
                nr = 1 + int(next((x.split("=", 1)[1] for x in options if x.startswith("--numRealTimeSteps=")), "4"))
                nm = 1 + int(next((x.split("=", 1)[1] for x in options if x.startswith("--numImagTimeSteps=")), "4"))
                for key, array in arrays.items():
                    assert array.shape[0] == nr and np.isfinite(array).all(), (key, array)
                assert arrays["Re_correlation"].shape[-1] == nm
                rank = f["runtimedata"].attrs["gaussian_factor_latent_dimensions"]
                if q <= 1 and args.reference:
                    previous_work = root / (name+"_previous")
                    previous_work.mkdir()
                    previous = subprocess.run([str(args.reference.resolve()), *base_options, *options],
                        cwd=previous_work, env=env, text=True, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, timeout=60)
                    assert previous.returncode == 0, previous.stdout
                    previous_files = list(previous_work.rglob("*.hdf5"))
                    assert len(previous_files) == 1
                    with h5py.File(previous_files[0]) as old:
                        for key in f["results"]:
                            np.testing.assert_array_equal(f["results"][key][()], old["results"][key][()],
                                                          err_msg=name+"/"+key)
                        np.testing.assert_array_equal(rank, old["runtimedata"].attrs["gaussian_factor_latent_dimensions"])
                    checked.append(name+"_matches_previous")
            checked.append(name)
            return arrays, rank

        for name, value in (("negative", "-2"), ("fractional", "1.5"),
                            ("overflow", "184467440737095516160"), ("grid_overflow", "18446744073709551615")):
            run(name, ["--realTimeSubsteps="+value], "realTimeSubsteps")
        run("svd_cf4_rejected", ["--realTimeSubsteps=3", "--gaussianFactorization=svd"], "use realTimeSubsteps=0 for svd")
        independent = ["--samplingStrategy=independent"]
        baseline, rank = run("defaults", independent)
        explicit, _ = run("explicit_defaults", independent + ["--gaussianFactorization=fft",
                           "--fftCrossFrequencyCutoff=-1", "--realTimeSubsteps=1"])
        for key in baseline:
            np.testing.assert_array_equal(baseline[key], explicit[key])
        for q in (0, 1, 3):
            for strategy in ("independent", "pcn"):
                for normalization in ("partition-function", "closed-contour"):
                    name = f"q{q}_{strategy}_{normalization}"
                    options = ["--realTimeSubsteps="+str(q), "--samplingStrategy="+strategy,
                        "--mhBurnIn=2", "--correlationNormalization="+normalization]
                    result, refined_rank = run(name, options)
                    if q == 0 and args.endpoint_reference:
                        work = root / (name+"_reference")
                        work.mkdir()
                        legacy_options = [x for x in options if not x.startswith("--realTimeSubsteps=")]
                        legacy = subprocess.run([str(args.endpoint_reference.resolve()), *common, *legacy_options],
                            cwd=work, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
                        assert legacy.returncode == 0, legacy.stdout
                        files = list(work.rglob("*.hdf5"))
                        assert len(files) == 1
                        with h5py.File(files[0]) as f:
                            for key, values in result.items():
                                np.testing.assert_array_equal(values, f["results"][key][()])
                        checked.append(name+"_matches_previous_endpoint")
                    np.testing.assert_array_equal(refined_rank, rank)
        for factor in ("dense", "weighted-dense"):
            factor_rank = None
            strategies = ("independent", "pcn") if factor == "dense" else ("independent",)
            for q in (0, 1, 2, 3, 4):
                for strategy in strategies:
                    for normalization in ("partition-function", "closed-contour"):
                        options = [f"--gaussianFactorization={factor}", f"--realTimeSubsteps={q}",
                            f"--samplingStrategy={strategy}", "--mhBurnIn=2",
                            f"--correlationNormalization={normalization}"]
                        if factor == "weighted-dense":
                            options += ["--gaussianWeightEta=2", "--gaussianWeightKappa=8"]
                        _, refined_rank = run(f"{factor}_q{q}_{strategy}_{normalization}", options)
                        if factor_rank is None:
                            factor_rank = refined_rank
                        np.testing.assert_array_equal(refined_rank, factor_rank)
        run("q0_short_grid", independent + ["--realTimeSubsteps=0", "--numImagTimeSteps=1", "--numRealTimeSteps=1"])
        for factor in ("dense", "svd"):
            run("q0_"+factor, independent + ["--realTimeSubsteps=0", "--gaussianFactorization="+factor])
        run("removed_cf4_option", ["--cf4Propagator"], "unrecognised option")
        run("cutoff_opt_in", independent + ["--realTimeSubsteps=2", "--fftCrossFrequencyCutoff=3"])
    print(f"Passed {len(checked)} dense/FFT substep CLI and HDF5 checks.")


if __name__ == "__main__":
    main()
