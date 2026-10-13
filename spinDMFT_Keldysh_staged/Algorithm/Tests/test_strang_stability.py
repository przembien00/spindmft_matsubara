#!/usr/bin/env python3
"""Optional long-time regression for the reported polarized Strang instability.

Uses one rank and 60000 samples per stage, with beta=2.5, JL=-2, h_z=0.5,
Tmax=15, 100 real intervals, 30 imaginary intervals and seed 734.
This takes a few minutes and is deliberately outside the default CTest suite.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def run(executable, work, options):
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    result = subprocess.run([str(executable), *options], cwd=work, env=env,
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            timeout=600)
    (work / "run.log").write_text(result.stdout)
    return result, list(work.rglob("*.hdf5"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--baseline", type=Path, help="optional saved executable from before the fix")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    options = ["--realTimeEvolution=staggered-strang", "--spinmodel=ISO", "--cstype=C",
               "--JL=-2", "--beta=2.5", "--Bname=z", "--Babs=0.5", "--dontrescaleB",
               "--numImagTimeSteps=30", "--numRealTimeSteps=100", "--Tmax=15",
               "--numSamplesPerCore=60000", "--preparationSamplesPerCore=60000",
               "--numBlocks=24", "--iterlimit=12", "--preparationIterlimit=20", "--seed=734"]
    with tempfile.TemporaryDirectory(prefix="strang-stability-") as temporary:
        root = args.output_dir.resolve() if args.output_dir else Path(temporary)
        # Unique run directories preserve any previous comparison results.
        root.mkdir(parents=True, exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix="balanced_", dir=root))
        result, files = run(args.executable.resolve(strict=True), work, options)
        assert result.returncode == 0 and len(files) == 1, result.stdout
        with h5py.File(files[0]) as f:
            re = f["results/Re_correlation"][()]
            im = f["results/Im_correlation"][()]
            maximum = float(np.max(np.abs(re + 1j*im)))
            assert np.isfinite(re).all() and np.isfinite(im).all()
            assert maximum < 0.5, "the reported long-time correlation explosion returned"
            assert np.max(f["preparation"].attrs["frozen_MM_residuals"]) < 1e-12
            assert b"sqrt(2)*eta,nu/sqrt(2)" in f["parameters"].attrs["staggered_gaussian_basis"]
            assert f["runtimedata"].attrs["absolute_iteration_errors"][-1] < 0.5
            report = {"fixed_output": str(files[0]), "max_abs_correlation": maximum,
                      "real_iterations": int(f["runtimedata"].attrs["num_Iterations"])}
        if args.baseline:
            before = Path(tempfile.mkdtemp(prefix="unbalanced_", dir=root))
            old, old_files = run(args.baseline.resolve(strict=True), before, options)
            if old_files:
                with h5py.File(old_files[0]) as f:
                    a = f["results/Re_correlation"][()] + 1j*f["results/Im_correlation"][()]
                    report["baseline_max_abs_correlation"] = float(np.nanmax(np.abs(a)))
                    assert report["baseline_max_abs_correlation"] > 1., "baseline no longer reproduces the failure"
            else:
                assert old.returncode != 0 and "non-finite sample weight" in old.stdout, old.stdout
            report["baseline_exit_code"] = old.returncode
        (root / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
