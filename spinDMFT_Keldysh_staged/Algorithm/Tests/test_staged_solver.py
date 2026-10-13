#!/usr/bin/env python3
"""Executable/HDF5 tests for staged thermal preparation (numpy, h5py required)."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def read_complex(group, name):
    return group['Re_' + name][...] + 1j * group['Im_' + name][...]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--mpi-ranks', type=int, default=1)
    parser.add_argument('--evolution', choices=['cfet', 'staggered-strang'], default='cfet')
    args = parser.parse_args()
    executable = args.executable.resolve()
    records = []
    common = ['--beta=1', '--Bname=z', '--Babs=0.5', '--dontrescaleB',
              '--cstype=C', '--numImagTimeSteps=6', '--numRealTimeSteps=6',
              '--Tmax=0.5', '--numSamplesPerCore=64', '--numBlocks=8',
              '--iterlimit=5', '--seed=17']
    common.append('--realTimeEvolution=' + args.evolution)
    launcher = ['mpirun', '-n', str(args.mpi_ranks)] if args.mpi_ranks > 1 else []
    with tempfile.TemporaryDirectory(prefix='staged_solver_') as temporary:
        root = args.output_dir.resolve() if args.output_dir else Path(temporary)
        root.mkdir(parents=True, exist_ok=True)

        def run(name, options, status=0, message=None):
            work = Path(tempfile.mkdtemp(prefix=name + '_', dir=root))
            overridden = {option.split("=", 1)[0] for option in options}
            defaults = [option for option in common if option.split("=", 1)[0] not in overridden]
            command = [*launcher, str(executable), *defaults, *options]
            result = subprocess.run(command, cwd=work, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    env=dict(os.environ), timeout=90)
            (work / 'run.log').write_text(result.stdout)
            assert result.returncode == status, (name, result.returncode, result.stdout)
            if message:
                assert message in result.stdout, result.stdout
            files = list(work.rglob('*.hdf5'))
            if status == 1:
                assert not files, files
                records.append({'case': name, 'rejected': True})
                return None
            assert len(files) == 1, (name, files, result.stdout)
            assert files[0].name.startswith('staged__')
            with h5py.File(files[0]) as f:
                assert text(f['parameters'].attrs['joint_matsubara_pcn']) == 'no'
                p = f['parameters'].attrs
                evolution = next((v.split('=', 1)[1] for v in options
                                  if v.startswith('--realTimeEvolution=')), args.evolution)
                assert text(p['real_time_evolution']) == evolution
                assert ('__prop=staggered-strang' in files[0].name) == (evolution == 'staggered-strang')
                if evolution == 'staggered-strang' and p['num_RealTimeSteps'] > 0:
                    q = int(p['real_time_substeps'])
                    assert q >= 1
                    assert text(p['propagator']) == 'staggered-strang'
                    assert 'direct joint' in text(p['real_time_field_interpolation'])
                    assert 'exactly zero' in text(p['staggered_covariance'])
                    assert 'sqrt(2)*eta,nu/sqrt(2)' in text(p['staggered_gaussian_basis'])
                    assert ('__substeps=' in files[0].name) == (q > 1)
                    assert p['delta_real_propagation_t'] == p['delta_real_t'] / q
            records.append({'case': name, 'exit_code': status, 'output': str(files[0])})
            return files[0]

        spin = {
            'x': np.array([[0, 1], [1, 0]], complex) / 2,
            'y': np.array([[0, -1j], [1j, 0]], complex) / 2,
            'z': np.diag([0.5, -0.5]).astype(complex),
        }
        energies = np.array([0.25, -0.25])
        partition = np.exp(-energies).sum()
        expected_m = -0.5 * np.tanh(0.25)
        for substeps in ([1, 3] if args.evolution == 'staggered-strang' else [0, 1]):
            # Exercise both a partial final batch and fewer than 32 samples.
            zero_samples = 37 if substeps in (0, 3) else 17
            options = ['--JQ=0', '--realTimeSubsteps=' + str(substeps),
                       '--numSamplesPerCore=' + str(zero_samples)]
            path = run('zero_bath_cf' + str(substeps), options)
            with h5py.File(path) as f:
                assert text(f['preparation'].attrs['termination']) == 'by convergence'
                assert text(f['preparation'].attrs['real_time_started']) == 'yes'
                threshold = 5 * 0.25 / np.sqrt(3 * zero_samples * args.mpi_ranks)
                np.testing.assert_allclose(
                    f['preparation'].attrs['absolute_iteration_error_threshold'], threshold,
                    atol=0, rtol=1e-14)
                assert text(f['runtimedata'].attrs['termination']) == 'by convergence'
                assert np.all(f['preparation'].attrs['frozen_MM_residuals'] == 0)
                assert f['preparation'].attrs['preparation_latent_dimension'] == 0
                assert np.all(read_complex(f['preparation'], 'covariance') == 0)
                m = read_complex(f['results'], 'magnetization')
                np.testing.assert_allclose(m[:, 2], expected_m, atol=2e-13, rtol=0)
                labels = [text(v) for v in f['results/correlation_direction_labels'][...]]
                actual = read_complex(f['results'], 'correlation')
                tau = f['results'].attrs['imaginary_time_edges']
                times = f['results'].attrs['real_times']
                expected = np.empty_like(actual)
                for ti, t in enumerate(times):
                    u = np.diag(np.exp(-1j * t * energies))
                    for ci, (a, b) in enumerate(labels):
                        sa = u.conj().T @ spin[a] @ u
                        for ki, s in enumerate(tau):
                            insertion = np.diag(np.exp(-(1-s)*energies)) @ spin[b] @ np.diag(np.exp(-s*energies))
                            expected[ti, ci, ki] = np.trace(insertion @ sa) / partition
                np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=0)
                records[-1]['max_exact_correlation_error'] = float(np.max(np.abs(actual-expected)))

        failed = run('preparation_limit', ['--JQ=0', '--preparationIterlimit=1'], 2,
                     'real-time stage was not started')
        with h5py.File(failed) as f:
            assert text(f['preparation'].attrs['real_time_started']) == 'no'
            assert text(f['preparation'].attrs['termination']) == 'by iteration limit'
            assert f['results/Re_correlation'].shape[0] == 1
            assert 'NOT_CONVERGED' in str(failed)

        failed = run('realtime_limit', ['--JQ=0', '--iterlimit=1', '--preparationIterlimit=5'], 3)
        with h5py.File(failed) as f:
            assert text(f['preparation'].attrs['termination']) == 'by convergence'
            assert text(f['runtimedata'].attrs['termination']) == 'by iteration limit'
            assert 'NOT_CONVERGED' in str(failed)

        for removed in ['--conditionalSampling', '--matsubaraGaugeScale=1']:
            run('reject_removed_' + removed[2:].split('=')[0], [removed], 1, 'unrecognised option')
        for option in ['--gaussianFactorization=fft',
                       '--gaussianFactorization=weighted-dense',
                       '--spinInsertionStrategy=closed-contour',
                       '--correlationNormalization=closed-contour']:
            run('reject_' + option[2:].replace('=', '_'), [option], 1, 'staged solver requires')
        run('reject_noselfcons', ['--noselfcons'], 1, 'cannot prepare a converged')
        run('reject_cutoff', ['--fftCrossFrequencyCutoff=2'], 1, 'does not support an FFT cutoff')
        run('reject_evolution', ['--realTimeEvolution=unknown'], 1, 'realTimeEvolution must be')
        if args.evolution == 'staggered-strang':
            run('reject_zero_substeps', ['--realTimeSubsteps=0'], 1, 'requires realTimeSubsteps >= 1')
            run('reject_short_imaginary', ['--numImagTimeSteps=2'], 1, 'requires numImagTimeSteps >= 3')
            run('short_real_grid', ['--JQ=0', '--numRealTimeSteps=1', '--realTimeSubsteps=2'])
        harmonic_options = ['--beta=0.5', '--cstype=D', '--bath=harmonic',
                   '--bathComponent=z', '--bathCoupling=0.25', '--numSamplesPerCore=8192',
                   '--preparationSamplesPerCore=4096', '--numBlocks=32']
        path = run('harmonic_dephasing', harmonic_options)
        if args.evolution == 'staggered-strang':
            reference = run('harmonic_preparation_reference', [*harmonic_options, '--realTimeEvolution=cfet'])
            with h5py.File(path) as new, h5py.File(reference) as old:
                for name in new['preparation']:
                    np.testing.assert_array_equal(new['preparation'][name][...], old['preparation'][name][...],
                                                  err_msg='Matsubara preparation changed: ' + name)
                for name in ('termination', 'num_Iterations', 'samples_per_core', 'generated_seed'):
                    assert new['preparation'].attrs[name] == old['preparation'].attrs[name]
        with h5py.File(path) as f:
            g = f['preparation']
            np.testing.assert_allclose(
                g.attrs['absolute_iteration_error_threshold'],
                5 * 0.25 / np.sqrt(3 * 4096 * args.mpi_ranks), atol=0, rtol=1e-14)
            A = read_complex(g, 'covariance')
            F = read_complex(g, 'factor')
            np.testing.assert_allclose(F @ F.T, A, atol=1e-11, rtol=1e-10)
            assert text(f['parameters'].attrs['field_sampling']) == (
                'full-contour dense Takagi; fixed Matsubara pseudo-covariance and mean')
            assert 'freshly factorized full contour' in text(g.attrs['factor_scope'])
            assert 'auxiliary marginal may change' in text(g.attrs['frozen_sampling_scope'])
            assert 'matsubara_gauge_scale' not in g.attrs
            assert 'conditional_sampling' not in f['parameters'].attrs
            assert g.attrs['preparation_latent_dimension'] <= A.shape[0]
            assert np.max(g.attrs['frozen_MM_residuals']) < 1e-12
            assert g.attrs['samples_per_core'] == 4096
            assert text(g.attrs['termination']) == 'prescribed harmonic bath'
            m = read_complex(f['results'], 'magnetization')
            errors = f['runtimedata/Re_magnetization_sample_stds'][..., 2]
            exact = -0.5*np.tanh(0.125)
            assert np.all(np.abs(m[:, 2].real-exact) < 6*errors+3e-4), (m[:, 2], errors)
            labels = [text(v) for v in f['results/correlation_direction_labels'][...]]
            zz = labels.index('zz')
            value = f['results/Re_correlation'][:, zz, 0]
            error = f['runtimedata/Re_correlation_sample_stds'][:, zz, 0]
            assert np.all(np.abs(value-0.25) < 6*error+3e-4), (value, error)
            # Raw real-time origin data remain independently measured and are
            # not replaced with the smaller-budget preparation output.
            assert not np.array_equal(f['results/Re_correlation'][0], g['Re_correlation'][0])
            records[-1]['max_MM_residual'] = float(np.max(g.attrs['frozen_MM_residuals']))
            records[-1]['max_stationarity_error'] = float(np.max(np.abs(m[:, 2].real-exact)))

        path = run('selfconsistent_finite_field', ['--beta=0.5', '--JQ=1',
                   '--numSamplesPerCore=2048', '--preparationSamplesPerCore=4096',
                   '--numBlocks=16', '--iterlimit=15', '--preparationIterlimit=15',
                   '--mixingAlpha=0.7'])
        with h5py.File(path) as f:
            assert text(f['preparation'].attrs['termination']) == 'by convergence'
            assert text(f['runtimedata'].attrs['termination']) == 'by convergence'
            assert np.max(f['preparation'].attrs['frozen_MM_residuals']) < 1e-12
            assert np.max(f['runtimedata'].attrs['gaussian_factor_reconstruction_errors']) < 1e-8
            assert f['results/Re_correlation'].shape == (7, 4, 7)
            assert np.all(np.isfinite(read_complex(f['results'], 'correlation')))
            records[-1]['preparation_iterations'] = int(f['preparation'].attrs['num_Iterations'])
            records[-1]['realtime_iterations'] = int(f['runtimedata'].attrs['num_Iterations'])
            records[-1]['max_MM_residual'] = float(np.max(f['preparation'].attrs['frozen_MM_residuals']))

            m = read_complex(f['results'], 'magnetization')[:, 2]
            records[-1]['max_origin_relative_magnetization_drift'] = float(np.max(np.abs(m-m[0])))

        report = json.dumps({'mpi_ranks': args.mpi_ranks, 
                             'real_time_evolution': args.evolution,
                             'cases': records}, indent=2)
        if args.output_dir:
            (root / 'validation.json').write_text(report + '\n')
        print(report)


if __name__ == '__main__':
    main()
