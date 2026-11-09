#!/usr/bin/env python3
"""Two-time executable/HDF5 regression, including MPI and both samplers."""
import argparse
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def complex_data(group, name):
    return group['Re_' + name][...] + 1j * group['Im_' + name][...]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--mpi-ranks', type=int, default=1)
    args = parser.parse_args()
    executable = args.executable.resolve()
    launcher = ['mpirun', '-n', str(args.mpi_ranks)] if args.mpi_ranks > 1 else []
    common = ['--beta=1', '--Bname=z', '--Babs=0.5', '--dontrescaleB',
              '--cstype=D', '--numImagTimeSteps=4', '--numRealTimeSteps=4',
              '--Tmax=0.6', '--numSamplesPerCore=32', '--numBlocks=8',
              '--iterlimit=4', '--seed=17']
    with tempfile.TemporaryDirectory(prefix='staged_two_time_') as directory:
        root = Path(directory)

        def run(name, options, status=0):
            work = root / name
            work.mkdir()
            overridden = {x.split('=', 1)[0] for x in options}
            defaults = [x for x in common if x.split('=', 1)[0] not in overridden]
            result = subprocess.run([*launcher, str(executable), *defaults, *options],
                                    cwd=work, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=90)
            allowed = status if isinstance(status, tuple) else (status,)
            assert result.returncode in allowed, (name, result.returncode, result.stdout)
            files = list(work.rglob('*.hdf5'))
            assert len(files) == 1, (name, files, result.stdout)
            return files[0]

        # Default and explicit true retain the same sampling and estimators.
        default = run('default', ['--JQ=0'])
        explicit = run('explicit', ['--JQ=0', '--timeTranslationInvariant=true'])
        with h5py.File(default) as a, h5py.File(explicit) as b:
            for group in ('results', 'preparation'):
                for name in a[group]:
                    np.testing.assert_array_equal(a[group][name][...], b[group][name][...])
            assert 'Re_correlation_greater' not in a['results']
            assert text(a['parameters'].attrs['time_translation_invariant']) == 'yes'

        for pcn in (False, True):
            for spin in (0.5, 1.):
                options = ['--JQ=0', '--timeTranslationInvariant=false',
                           '--spin=' + ('1/2' if spin == 0.5 else '1'),
                           '--realTimeSubsteps=2']
                if pcn:
                    options += ['--jointMatsubaraPCN', '--mhBurnIn=4']
                path = run(f'full_{pcn}_{spin}', options)
                assert 'two_time__' in path.name
                with h5py.File(path) as f, h5py.File(default) as baseline:
                    assert text(f['parameters'].attrs['time_translation_invariant']) == 'no'
                    result = f['results']
                    times = result.attrs['real_times']
                    labels = [text(x) for x in result['correlation_direction_labels'][...]]
                    greater = complex_data(result, 'correlation_greater')
                    lesser = complex_data(result, 'correlation_lesser')
                    assert greater.shape == lesser.shape == (5, 9, 5)
                    assert result['Re_correlation'].shape == (5, 9, 5)
                    assert f['preparation/Re_correlation'].shape == (9, 5)
                    assert np.max(f['preparation'].attrs['frozen_MM_residuals']) < 1e-10
                    m = np.arange(spin, -spin - 0.5, -1)
                    raising = np.zeros((len(m), len(m)), complex)
                    for j in range(1, len(m)):
                        raising[j-1, j] = np.sqrt(spin*(spin+1)-m[j]*(m[j]+1))
                    ops = {'x': (raising + raising.T) / 2,
                           'y': (raising - raising.T) / (2j), 'z': np.diag(m)}
                    energies = 0.5*m
                    rho = np.diag(np.exp(-energies) / np.exp(-energies).sum())
                    heisenberg = []
                    for t in times:
                        u = np.diag(np.exp(-1j*t*energies))
                        heisenberg.append({a: u.conj().T @ op @ u for a, op in ops.items()})
                    for t in range(len(times)):
                        for s in range(len(times)):
                            for p, (a, b) in enumerate(labels):
                                expected = np.trace(rho @ heisenberg[t][a] @ heisenberg[s][b])
                                expected_lesser = np.trace(rho @ heisenberg[s][b] @ heisenberg[t][a])
                                np.testing.assert_allclose(greater[t, p, s], expected, atol=3e-12)
                                np.testing.assert_allclose(lesser[t, p, s], expected_lesser, atol=3e-12)
                    mixed = complex_data(result, 'correlation')
                    np.testing.assert_allclose(greater[:, :, 0], mixed[:, :, -1], atol=3e-12)
                    np.testing.assert_allclose(lesser[:, :, 0], mixed[:, :, 0], atol=3e-12)
                    for ordering in ('greater', 'lesser'):
                        for part in ('Re', 'Im'):
                            errors = f[f'runtimedata/{part}_correlation_{ordering}_sample_stds'][...]
                            assert errors.shape == (5, 9, 5)
                            assert np.all(np.isfinite(errors)) and np.all(errors >= 0)
                            if pcn:
                                assert f[f'runtimedata/{part}_correlation_{ordering}_tau_int'].shape == errors.shape
                    if spin == 0.5 and not pcn:
                        for name in f['preparation']:
                            np.testing.assert_array_equal(f['preparation'][name][...], baseline['preparation'][name][...])

        # Exercise actual feedback for a nonzero bath on two iterations.
        for pcn in (False, True):
            options = ['--timeTranslationInvariant=false', '--JQ=0.1',
                       '--reliterror=1', '--iterlimit=2', '--preparationIterlimit=10']
            if pcn:
                options = [x for x in options if not x.startswith('--reliterror=')]
                options += ['--jointMatsubaraPCN', '--mhBurnIn=16', '--mhStepSize=1',
                            '--numSamplesPerCore=256', '--reliterror=5']
            path = run('feedback_' + str(pcn), options, status=(0, 3))
            with h5py.File(path) as f:
                assert f['runtimedata'].attrs['num_Iterations'] == 2
                assert np.max(f['preparation'].attrs['frozen_MM_residuals']) < 1e-10
                assert np.all(np.isfinite(f['results/Re_correlation_greater'][...]))

        # The frozen preparation also has a valid output if it fails.
        path = run('preparation_failure', ['--timeTranslationInvariant=false',
                                          '--preparationIterlimit=1', '--JQ=0', '--reliterror=0.1'], status=2)
        with h5py.File(path) as f:
            assert text(f['preparation'].attrs['real_time_started']) == 'no'
            assert 'Re_correlation_greater' not in f['results']
        work = root / 'invalid'
        work.mkdir()
        result = subprocess.run([*launcher, str(executable), *common,
                                 '--timeTranslationInvariant=false', '--constantMagnetization'],
                                cwd=work, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, timeout=30)
        assert result.returncode != 0 and 'incompatible' in result.stdout
        assert not list(work.rglob('*.hdf5'))
    print('Two-time solver/HDF5 checks passed')


if __name__ == '__main__':
    main()
