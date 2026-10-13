#!/usr/bin/env python3
"""Executable/HDF5 regression checks. Requires numpy, h5py, and working MPI.

Usage: python test_weighted_solver.py /absolute/path/executable_DOUBLE.out
All solver output is isolated in temporary directories.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def text(attribute):
    value = np.asarray(attribute).reshape(-1)[0]
    return value.decode() if isinstance(value, bytes) else str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    common = ["--numImagTimeSteps=4", "--numRealTimeSteps=4", "--beta=0.8",
              "--Tmax=0.3", "--numSamplesPerCore=64", "--numBlocks=8", "--seed=734",
              "--cstype=D", "--Bname=z", "--Babs=0.7", "--dontrescaleB", "--iterlimit=2"]
    independent = ["--samplingStrategy=independent"]
    weighted = ["--gaussianFactorization=weighted-dense"]
    records = []
    with tempfile.TemporaryDirectory(prefix="weighted-solver-tests-") as directory:
        root = Path(directory)

        def run(name, options, error=None):
            work = root / name
            work.mkdir()
            result = subprocess.run([str(executable), *common, *options], cwd=work,
                                    env=env, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=60)
            if error is not None:
                expected = (error,) if isinstance(error, str) else error
                assert result.returncode != 0 and any(message in result.stdout for message in expected), result.stdout
                assert not list(work.rglob("*.hdf5")), "invalid options wrote simulation data"
                records.append({"case": name, "rejected": True})
                return
            assert result.returncode == 0, result.stdout
            files = list(work.rglob("*.hdf5"))
            assert len(files) == 1, (files, result.stdout)
            return files[0]

        invalid = [
            ("pcn", weighted, "requires samplingStrategy=independent"),
            ("ignored_dense", independent + ["--gaussianWeightEta=2"], "require gaussianFactorization"),
            ("ignored_fft", independent + ["--gaussianFactorization=fft", "--gaussianWeightEta=2"], "require gaussianFactorization"),
            ("zero", independent + weighted + ["--gaussianWeightM=0"], "strictly positive"),
            ("negative", independent + weighted + ["--gaussianWeightEta=-1"], "strictly positive"),
            ("nan", independent + weighted + ["--gaussianWeightKappa=nan"], ("gaussianWeightKappa", "strictly positive")),
            ("infinity", independent + weighted + ["--gaussianWeightM=inf"], ("gaussianWeightM", "strictly positive")),
        ]
        for name, options, error in invalid:
            run(name, options, error)

        cases = [
            ("zero_dense", ["--gaussianFactorization=dense", "--JQ=0", "--noselfcons"], None, False),
            ("zero_weighted", weighted + ["--JQ=0", "--noselfcons"], (1., 2., 8.), False),
            ("zero_weighted_cf4", weighted + ["--JQ=0", "--noselfcons", "--realTimeSubsteps=1"], (1., 2., 8.), True),
            ("harmonic_weighted", weighted + ["--bath=harmonic", "--bathCoupling=0.2", "--bathComponent=z"], (1., 2., 8.), False),
            ("harmonic_custom_cf4", weighted + ["--bath=harmonic", "--bathCoupling=0.2", "--bathComponent=z",
                "--gaussianWeightM=0.1", "--gaussianWeightEta=8", "--gaussianWeightKappa=2", "--realTimeSubsteps=1",
                "--correlationNormalization=closed-contour"], (0.1, 8., 2.), True),
            ("selfconsistent_weighted", weighted + ["--JQ=0.2"], (1., 2., 8.), False),
        ]
        dense_result = None
        for name, options, weights, cf4 in cases:
            path = run(name, independent + options + ([] if cf4 else ["--realTimeSubsteps=0"]))
            with h5py.File(path) as file:
                p, result = file["parameters"].attrs, file["results"]
                assert text(p["sampling_strategy"]) == "independent"
                assert text(p["propagator"]) == ("gauss-cf4" if cf4 else "endpoint-cfet4")
                if weights:
                    assert text(p["gaussian_factorization"]) == "weighted-dense"
                    actual = tuple(float(np.asarray(p[key]).reshape(-1)[0]) for key in
                                   ("gaussian_weight_M", "gaussian_weight_eta", "gaussian_weight_kappa"))
                    assert actual == weights
                    assert "kappa=(V_+-V_-)/2" in text(p["gaussian_weight_basis"])
                    assert "original physical" in text(p["gaussian_reconstruction_basis"])
                    assert "__noiseW=" in path.name
                else:
                    assert text(p["gaussian_factorization"]) == "dense"
                    assert "gaussian_weight_M" not in p and "__noiseW=" not in path.name
                arrays = {key: result[key][()] for key in
                          ("Re_correlation", "Im_correlation", "Re_magnetization", "Im_magnetization")}
                assert all(np.isfinite(value).all() for value in arrays.values())
                error = float(np.max(file["runtimedata"].attrs["gaussian_factor_reconstruction_errors"]))
                assert error < 1e-10
                if name.startswith("zero"):
                    magnetization = arrays["Re_magnetization"]
                    assert magnetization.shape == (5, 3)
                    np.testing.assert_allclose(magnetization[:, 2], -0.5 * np.tanh(0.8 * 0.7 / 2), atol=1e-12)
                    np.testing.assert_allclose(arrays["Im_magnetization"], 0, atol=1e-12)
                    if dense_result is None:
                        dense_result = arrays
                    else:
                        for key in arrays:
                            np.testing.assert_allclose(arrays[key], dense_result[key], atol=1e-12)
                records.append({"case": name, "physical_covariance_residual": error})
    print(json.dumps({"passed": len(records), "cases": records}, indent=2))


if __name__ == "__main__":
    main()
