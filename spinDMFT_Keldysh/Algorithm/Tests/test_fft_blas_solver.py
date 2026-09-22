#!/usr/bin/env python3
"""Compare current FFT sampling against a saved pre-BLAS solver executable.

Usage: python3 test_fft_blas_solver.py CURRENT_EXECUTABLE REFERENCE_EXECUTABLE
Uses small fixed-seed runs at fixed input covariance, checks observables/errors and pCN decisions, and
keeps all HDF5 output in a temporary directory. Requires numpy and h5py.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('reference', type=Path)
    args = parser.parse_args()
    executables = [p.resolve(strict=True) for p in (args.executable, args.reference)]
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    common = ['--numImagTimeSteps=6', '--numRealTimeSteps=8', '--beta=0.8', '--Tmax=1.2',
              '--numSamplesPerCore=65', '--numBlocks=4', '--seed=734', '--cstype=C',
              '--Bname=z', '--Babs=0.7', '--dontrescaleB', '--iterlimit=1', '--JQ=0.1',
              '--mhBurnIn=3', '--gaussianFactorization=fft']
    cases = 0
    maximum_difference = 0.
    with tempfile.TemporaryDirectory(prefix='fft-blas-regression-') as directory:
        root = Path(directory)
        for q in (0, 1, 3):
            for cutoff in (-1, 3):
                for strategy in ('independent', 'pcn'):
                    for normalization in ('partition-function', 'closed-contour'):
                        options = common + [f'--realTimeSubsteps={q}', f'--fftCrossFrequencyCutoff={cutoff}',
                            '--samplingStrategy='+strategy, '--correlationNormalization='+normalization]
                        files = []
                        for version, executable in enumerate(executables):
                            work = root / f'case{cases}_{version}'
                            work.mkdir()
                            proc = subprocess.run([str(executable), *options], cwd=work, env=env,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
                            assert proc.returncode == 0, proc.stdout
                            outputs = list(work.rglob('*.hdf5'))
                            assert len(outputs) == 1, (outputs, proc.stdout)
                            files.append(outputs[0])
                        with h5py.File(files[0]) as current, h5py.File(files[1]) as reference:
                            for group in ('results', 'runtimedata'):
                                assert set(current[group]) == set(reference[group])
                                for key, dataset in current[group].items():
                                    actual, expected = dataset[()], reference[group][key][()]
                                    if np.issubdtype(np.asarray(actual).dtype, np.number):
                                        np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-12,
                                            err_msg=f'q={q}, cutoff={cutoff}, {strategy}, {normalization}, {group}/{key}')
                                        if np.size(actual):
                                            maximum_difference = max(maximum_difference, float(np.max(np.abs(actual-expected))))
                                    else:
                                        np.testing.assert_array_equal(actual, expected)
                            # Sampling statistics are stored as attributes; pCN
                            # must retain exactly the same accept/reject decisions.
                            for key in current['runtimedata'].attrs:
                                if key in ('mh_acceptance_rates', 'mh_nonpositive_rejection_rates') or any(
                                        x in key for x in ('latent_dimensions', 'largest_factorization')):
                                    np.testing.assert_array_equal(current['runtimedata'].attrs[key], reference['runtimedata'].attrs[key])
                        cases += 1
    print(f'Passed {cases} pre/post-BLAS FFT solver comparisons; maximum dataset difference {maximum_difference:.3g}.')


if __name__ == '__main__':
    main()
