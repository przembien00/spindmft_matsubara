#!/usr/bin/env python3
"""Staggered Strang CLI, solver and HDF5 checks in isolated temporary folders."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import h5py
import numpy as np


def scalar(value):
    return np.asarray(value).reshape(-1)[0]


def text(value):
    item = scalar(value)
    return item.decode() if isinstance(item, bytes) else str(item)


def main():
    executable = Path(sys.argv[1]).resolve(strict=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    base = ["--numImagTimeSteps=4", "--numRealTimeSteps=4", "--beta=0.8", "--Tmax=0.3",
            "--numSamplesPerCore=32", "--numBlocks=4", "--seed=734", "--cstype=D",
            "--Bname=z", "--Babs=0.7", "--dontrescaleB", "--iterlimit=1", "--JQ=0.1"]
    checked = []
    with tempfile.TemporaryDirectory(prefix="staggered-strang-tests-") as directory:
        def run(name, options, error=None):
            work = Path(directory) / name
            work.mkdir()
            keys = {x.split("=", 1)[0] for x in options}
            command = [str(executable), *[x for x in base if x.split("=", 1)[0] not in keys], *options]
            result = subprocess.run(command, cwd=work, env=env, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
            files = list(work.rglob("*.hdf5"))
            if error:
                assert result.returncode != 0 and error in result.stdout, result.stdout
                assert not files, "invalid parameters produced output"
            else:
                assert result.returncode == 0 and len(files) == 1, result.stdout
                with h5py.File(files[0]) as f:
                    p = f["parameters"].attrs
                    assert "fft_cross_frequency_cutoff" not in p
                    staggered = True
                    assert text(p["real_time_evolution"]) == "staggered-strang"
                    assert text(p["propagator"]) == "staggered-strang"
                    factor = next((x.split("=", 1)[1] for x in options if x.startswith("--gaussianFactorization=")),
                                  "dense")
                    assert text(p["gaussian_factorization"]) == factor
                    assert "__prop=" not in files[0].name
                    q = int(scalar(p["real_time_substeps"]))
                    if staggered:
                        assert "direct joint" in text(p["real_time_field_interpolation"])
                        assert "exactly zero" in text(p["staggered_covariance"])
                        assert ("__substeps=" in files[0].name) == (q > 1)
                    assert float(scalar(p["delta_real_propagation_t"])) == float(scalar(p["delta_real_t"])) / q
                    nr = 1 + int(next((x.split("=", 1)[1] for x in options if x.startswith("--numRealTimeSteps=")), "4"))
                    for key in ("Re_correlation", "Im_correlation", "Re_magnetization", "Im_magnetization"):
                        array = f["results"][key][()]
                        assert array.shape[0] == nr and np.isfinite(array).all(), key
                    assert np.max(f["runtimedata"].attrs["gaussian_factor_reconstruction_errors"]) < 1e-8
                    if factor == "chebyshev":
                        r = f["runtimedata"].attrs
                        assert np.max(r["gaussian_largest_factorization_dimensions"]) == 0
                        assert np.min(r["gaussian_polynomial_degrees"]) >= 16
                        assert np.all(r["gaussian_covariance_approximation_errors"] <=
                                      float(scalar(p["gaussian_chebyshev_tolerance"])) *
                                      r["gaussian_spectral_bounds"] * (1+1e-10))
            checked.append(name)

        run("default", [])
        option = ["--realTimeEvolution=staggered-strang"]
        run("dense", option)
        run("feedback", [*option, "--iterlimit=2", "--JL=-0.4"])
        run("substeps", [*option, "--realTimeSubsteps=3", "--spinInsertionStrategy=prefix"])
        run("short_grid", [*option, "--numRealTimeSteps=1"])
        run("weighted", [*option, "--gaussianFactorization=weighted-dense", "--gaussianWeightM=0.1"])
        run("chebyshev", ["--gaussianFactorization=chebyshev"])
        run("chebyshev_symmetric", ["--gaussianFactorization=chebyshev", "--cstype=A", "--Babs=0"])
        run("chebyshev_closed", ["--gaussianFactorization=chebyshev", "--correlationNormalization=closed-contour"])
        run("chebyshev_fft", ["--gaussianFactorization=chebyshev", "--realTimeSubsteps=5"])
        run("chebyshev_pcn", ["--gaussianFactorization=chebyshev", "--samplingStrategy=pcn"],
            "chebyshev requires samplingStrategy=independent")
        run("chebyshev_bad_tolerance", ["--gaussianFactorization=chebyshev", "--gaussianChebyshevTolerance=0"],
            "Chebyshev tolerance must lie")
        run("chebyshev_bad_degree", ["--gaussianFactorization=chebyshev", "--gaussianChebyshevMaxDegree=15"],
            "Chebyshev tolerance must lie")
        run("chebyshev_insufficient_degree", ["--gaussianFactorization=chebyshev", "--gaussianChebyshevMaxDegree=16"],
            "Chebyshev covariance error bound exceeds")
        run("chebyshev_option_dense", ["--gaussianChebyshevTolerance=0.01"],
            "gaussianChebyshev options require gaussianFactorization=chebyshev")
        run("pcn", [*option, "--samplingStrategy=pcn", "--mhBurnIn=2", "--realTimeSubsteps=2"])
        run("pcn_closed", [*option, "--samplingStrategy=pcn", "--mhBurnIn=2", "--correlationNormalization=closed-contour"])
        run("harmonic", [*option, "--bath=harmonic", "--bathCoupling=0.05"])
        run("spin_one", [*option, "--spin=1"])
        run("removed_cfet", ["--realTimeEvolution=cfet"], "other evolution methods have been removed")
        run("bad_evolution", ["--realTimeEvolution=unknown"], "realTimeEvolution must be")
        run("zero_substeps", [*option, "--realTimeSubsteps=0"], "requires realTimeSubsteps >= 1")
        run("short_imaginary", [*option, "--numImagTimeSteps=2"], "requires numImagTimeSteps >= 3")
        run("fft", [*option, "--gaussianFactorization=fft"], "gaussianFactorization must be dense or weighted-dense")
        run("removed_cutoff_default", [*option, "--fftCrossFrequencyCutoff=-1"], "unrecognised option")
        run("counterterm", [*option, "--contactCounterterm"], "unrecognised option")
        run("cutoff", [*option, "--fftCrossFrequencyCutoff=1"], "unrecognised option")
        run("weighted_pcn", [*option, "--gaussianFactorization=weighted-dense", "--samplingStrategy=pcn"],
            "weighted-dense requires samplingStrategy=independent")
    print(f"Passed {len(checked)} solver checks: {', '.join(checked)}")


if __name__ == "__main__":
    main()
