#!/usr/bin/env python3
"""Staged Matsubara-pCN regression against exact and direct thermal sampling."""
import argparse
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
    parser.add_argument('--mpi-ranks', type=int, default=1)
    args = parser.parse_args()
    executable = args.executable.resolve()
    launcher = ['mpirun', '-n', str(args.mpi_ranks)] if args.mpi_ranks > 1 else []
    common = ['--jointMatsubaraPCN', '--mhStepSize=0.7',
              '--mhBurnIn=512', '--beta=0.5', '--Bname=z', '--Babs=0.5', '--dontrescaleB',
              '--cstype=C', '--numImagTimeSteps=6', '--numRealTimeSteps=3',
              '--Tmax=0.3', '--numSamplesPerCore=128', '--numBlocks=16',
              '--iterlimit=15', '--preparationIterlimit=15', '--seed=31415']
    with tempfile.TemporaryDirectory(prefix='staged_pcn_') as root:
        def run(name, options, status=0, use_flag=True):
            work = Path(root) / name
            work.mkdir()
            overridden = {v.split('=', 1)[0] for v in options}
            defaults = [v for v in common if v.split('=', 1)[0] not in overridden]
            if not use_flag:
                defaults = [v for v in defaults if v != '--jointMatsubaraPCN']
            result = subprocess.run([*launcher, str(executable), *defaults, *options],
                                    cwd=work, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, env=dict(os.environ), timeout=90)
            # Invalid inputs now propagate as uncaught exceptions; the abort
            # status depends on whether the executable is launched through MPI.
            if status == 1:
                assert result.returncode != 0, (name, result.returncode, result.stdout)
            else:
                assert result.returncode == status, (name, result.returncode, result.stdout)
            files = list(work.rglob('*.hdf5'))
            if status == 1:
                assert not files
                return result.stdout
            assert len(files) == 1, (name, files)
            print(name, 'passed')
            return files[0]

        assert 'conflicts' in run('conflict', ['--samplingStrategy=independent'], 1)
        assert 'conflicts' in run('disabled_conflict', ['--jointMatsubaraPCN=false', '--samplingStrategy=pcn'], 1)
        for removed in ['--conditionalSampling', '--matsubaraGaugeScale=1']:
            assert 'unrecognised option' in run('removed_' + removed[2:].split('=')[0], [removed], 1)
        for step in ['nan', 'inf', '0']:
            assert 'mhStepSize must lie' in run('bad_step_' + step, ['--mhStepSize=' + step], 1)

        def check_metadata(f, path):
            p = f['parameters'].attrs
            for removed in ('convergence_rule', 'factor_scope', 'frozen_state_kind', 'frozen_sampling_scope', 'uncertainty_scope'):
                assert removed not in f['preparation'].attrs
            for part in ('Re', 'Im'):
                values = f['preparation'][part + '_correlation']
                assert values.ndim == 2
                for suffix in ('_sample_stds', '_tau_int'):
                    assert f['preparation'][part + '_correlation' + suffix].shape == values.shape
            assert f['results/Re_correlation'].ndim == 3
            assert 'jointMpcn_absZ__' in path.name
            assert text(p['sampling_strategy']) == 'pcn'
            assert text(p['joint_matsubara_pcn']) == 'yes'
            assert 'conditional_sampling' not in p
            assert 'null space' in text(p['real_residual_sampling'])
            assert '|Z_M|' in text(p['pcn_sampling_weight'])
            assert 'Matsubara' in text(p['pcn_updated_coordinates'])
            assert 'including rejection' in text(p['real_residual_sampling'])
            assert 'normalization_strategy' not in p
            assert text(f['preparation'].attrs['sampling_strategy']) == 'pcn'
            for g in [f['preparation'], f['runtimedata']]:
                rates = g.attrs['mh_acceptance_rates']
                assert len(rates) == g.attrs['num_Iterations']
                assert np.all((rates >= 0) & (rates <= 1))
                assert np.all(g.attrs['mh_zero_weight_rejection_rates'] == 0)
                assert np.all(np.isfinite(g.attrs['blocking_curve_max_error']))
                assert np.all(np.isfinite(g.attrs['average_phase_magnitudes']))
                assert 'Re_correlation_tau_int' in g
            assert np.max(f['preparation'].attrs['frozen_MM_residuals']) < 1e-12

        alias = run('strategy_alias', ['--samplingStrategy=pcn', '--JQ=0'], use_flag=False)
        with h5py.File(alias) as f:
            check_metadata(f, alias)

        for evolution, substeps, spin in [('staggered-strang', 3, '1/2'),
                                          ('staggered-strang', 2, '1/2'),
                                          ('staggered-strang', 1, '1')]:
            path = run('zero_' + evolution + '_' + str(substeps) + '_' + spin.replace('/', 'o'),
                       ['--JQ=0', '--realTimeEvolution=' + evolution,
                        '--realTimeSubsteps=' + str(substeps), '--spin=' + spin])
            with h5py.File(path) as f:
                check_metadata(f, path)
                s = 0.5 if spin == '1/2' else 1.
                eigenvalues = np.arange(-s, s + 1)
                weights = np.exp(-0.25 * eigenvalues)
                exact_m = np.sum(weights * eigenvalues) / weights.sum()
                exact_zz = np.sum(weights * eigenvalues**2) / weights.sum()
                m = read_complex(f['results'], 'magnetization')
                np.testing.assert_allclose(m[:, 2], exact_m, atol=2e-12, rtol=0)
                labels = [text(v) for v in f['results/correlation_direction_labels'][...]]
                zz = labels.index('zz')
                np.testing.assert_allclose(read_complex(f['results'], 'correlation')[:, zz], exact_zz,
                                           atol=2e-12, rtol=0)
                assert np.all(f['runtimedata'].attrs['mh_acceptance_rates'] == 1)
                assert np.all(f['preparation'].attrs['mh_acceptance_rates'] == 1)

        harmonic = ['--bath=harmonic', '--bathComponent=z', '--bathCoupling=0.25',
                    '--cstype=D', '--numSamplesPerCore=8192',
                    '--preparationSamplesPerCore=4096', '--numBlocks=32']
        for evolution in ['staggered-strang']:
            opts = [*harmonic, '--realTimeEvolution=' + evolution]
            pcn = run('harmonic_pcn_' + evolution, opts)
            direct = run('harmonic_direct_' + evolution, [*opts, '--jointMatsubaraPCN=false', '--samplingStrategy=independent'])
            with h5py.File(pcn) as new, h5py.File(direct) as old:
                check_metadata(new, pcn)
                assert np.any(new['runtimedata'].attrs['mh_acceptance_rates'] < 1)
                assert np.any(new['preparation'].attrs['mh_acceptance_rates'] < 1)
                for name in ['correlation', 'magnetization']:
                    for part in ['Re', 'Im']:
                        key = part + '_' + name
                        delta = np.abs(new['results'][key][...] - old['results'][key][...])
                        errors = np.hypot(new['runtimedata'][key + '_sample_stds'][...],
                                          old['runtimedata'][key + '_sample_stds'][...])
                        assert np.all(delta < 6 * errors + 2e-3), (evolution, key, delta.max(), errors.max())
                exact_m = -0.5 * np.tanh(0.125)
                m = read_complex(new['results'], 'magnetization')[:, 2]
                errors = new['runtimedata/Re_magnetization_sample_stds'][..., 2]
                assert np.all(np.abs(m.real - exact_m) < 6 * errors + 2e-3)
                assert not any(key in new['preparation'] for key in
                               ('Re_factor', 'Im_factor', 'Re_covariance', 'Im_covariance'))

        path = run('selfconsistent', ['--beta=0.3', '--JQ=1',
                   '--numSamplesPerCore=2048', '--preparationSamplesPerCore=4096', '--numBlocks=32'])
        with h5py.File(path) as f:
            check_metadata(f, path)
            assert text(f['preparation'].attrs['termination']) == 'by convergence'
            assert text(f['runtimedata'].attrs['termination']) == 'by convergence'
            assert np.all(np.isfinite(read_complex(f['results'], 'correlation')))


if __name__ == '__main__':
    main()
