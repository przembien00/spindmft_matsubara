#!/usr/bin/env python3
"""Pulse-file validation and exact driven-spin executable/HDF5 tests."""
import argparse
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np


def text(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def complex_data(group, name):
    return group['Re_' + name][...] + 1j * group['Im_' + name][...]


def spin_operators(spin):
    m = np.arange(spin, -spin - 0.5, -1)
    raising = np.zeros((len(m), len(m)), complex)
    for j in range(1, len(m)):
        raising[j-1, j] = np.sqrt(spin*(spin+1)-m[j]*(m[j]+1))
    return [(raising+raising.T)/2, (raising-raising.T)/(2j), np.diag(m)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--mpi-ranks', type=int, default=1)
    args = parser.parse_args()
    executable = args.executable.resolve()
    launcher = ['mpirun', '-n', str(args.mpi_ranks)] if args.mpi_ranks > 1 else []
    common = ['--JQ=0', '--Bname=z', '--Babs=0.5', '--dontrescaleB', '--beta=1',
              '--cstype=D', '--numImagTimeSteps=6', '--numRealTimeSteps=5', '--Tmax=1',
              '--numSamplesPerCore=64', '--numBlocks=8', '--iterlimit=5', '--seed=17']
    with tempfile.TemporaryDirectory(prefix='staged_pulses_') as temporary:
        root = Path(temporary)

        def write(name, times, axes, angles):
            path = root / (name + '.hdf5')
            with h5py.File(path, 'w') as file:
                group = file.create_group('pulses')
                group['times'], group['axes'], group['angles'] = times, axes, angles
            return path

        def run(name, options, status=(0,)):
            work = root / name
            work.mkdir()
            overridden = {x.split('=', 1)[0] for x in options}
            defaults = [x for x in common if x.split('=', 1)[0] not in overridden]
            result = subprocess.run([*launcher, str(executable), *defaults, *options],
                                    cwd=work, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, timeout=90)
            if status is None:
                assert result.returncode != 0, result.stdout
                assert not list(work.rglob('*.hdf5'))
                return result.stdout
            assert result.returncode in status, (name, result.returncode, result.stdout)
            files = list(work.rglob('*.hdf5'))
            assert len(files) == 1, (files, result.stdout)
            return files[0]

        # Distinct noncommuting axes, an origin pulse, interior and edge
        # pulses, two coincident pulses, and a pulse at the final time.
        times = np.array([0., 0.07, 0.2, 0.3, 0.3, 0.4, 1.])
        axes = np.array([[0, 2, 0], [1, 1, 0], [0, 0, 1], [1, 0, 0],
                         [0, 1, 0], [0, 0, 1], [1, 0, 0]], float)
        angles = np.array([np.pi/2, -0.8, np.pi, 0.3, 0.6, -0.2, np.pi/2])
        path = write('sequence', times, axes, angles)
        normalized_axes = axes / np.linalg.norm(axes, axis=1)[:, None]
        replay = None
        for spin in (0.5, 1.):
            ops = spin_operators(spin)
            energy = 0.5*np.diag(ops[2])
            weights = np.exp(-energy)
            rho = np.diag(weights/weights.sum())
            rotations = []
            for axis, angle in zip(normalized_axes, angles):
                values, vectors = np.linalg.eigh(sum(a*s for a, s in zip(axis, ops)))
                rotations.append((vectors*np.exp(1j*angle*values)) @ vectors.conj().T)
            baseline = run('baseline_' + str(spin), ['--spin=' + str(spin)])
            for pcn in (False, True):
                for q in (1, 2, 3):
                    options = ['--pulseSequence=' + str(path), '--spin=' + str(spin),
                               '--realTimeSubsteps=' + str(q)]
                    if pcn:
                        options += ['--jointMatsubaraPCN', '--mhBurnIn=4']
                    output = run(f'exact_{spin}_{pcn}_{q}', options)
                    replay = output
                    assert 'pulsed__' in output.name
                    with h5py.File(output) as file, h5py.File(baseline) as base:
                        assert text(file['parameters'].attrs['time_translation_invariant']) == 'no'
                        assert file['parameters'].attrs['number_of_pulses'] == len(times)
                        np.testing.assert_allclose(file['pulses/times'][...], times, atol=1e-15)
                        np.testing.assert_allclose(file['pulses/axes'][...], normalized_axes, atol=1e-15)
                        np.testing.assert_array_equal(file['pulses/angles'][...], angles)
                        for name in ('Re_correlation', 'Im_correlation', 'Re_magnetization', 'Im_magnetization'):
                            np.testing.assert_allclose(file['preparation'][name][...], base['preparation'][name][...], atol=4e-13)
                        assert np.max(file['preparation'].attrs['frozen_MM_residuals']) < 1e-12
                        result = file['results']
                        grid = result.attrs['real_times']
                        tau = result.attrs['imaginary_time_edges']
                        labels = [text(x) for x in result['correlation_direction_labels'][...]]
                        heisenberg = []
                        evolution = np.eye(len(energy), dtype=complex)
                        cursor = 0
                        while cursor < len(times) and times[cursor] == 0:
                            evolution = rotations[cursor] @ evolution
                            cursor += 1
                        heisenberg.append([evolution.conj().T @ s @ evolution for s in ops])
                        h = (grid[1]-grid[0])/q
                        for t in range(1, len(grid)):
                            for j in range(q):
                                start = grid[t-1] + j*h
                                stop = grid[t-1] + (j+1)*h
                                events = []
                                while cursor < len(times) and times[cursor] <= stop+1e-14:
                                    events.append(cursor)
                                    cursor += 1
                                generator = np.diag(energy).astype(complex)
                                explicit = []
                                for event in events:
                                    if j+1 == q and abs(times[event]-stop) <= 1e-14:
                                        generator -= angles[event]/h * sum(a*s for a, s in zip(normalized_axes[event], ops))
                                    else:
                                        explicit.append(event)
                                eigenvalues, eigenvectors = np.linalg.eigh(generator)
                                def evolve(duration):
                                    return (eigenvectors*np.exp(-1j*duration*eigenvalues)) @ eigenvectors.conj().T
                                last = start
                                for event in explicit:
                                    evolution = rotations[event] @ evolve(times[event]-last) @ evolution
                                    last = times[event]
                                evolution = evolve(stop-last) @ evolution
                            heisenberg.append([evolution.conj().T @ s @ evolution for s in ops])
                        magnetization = complex_data(result, 'magnetization')
                        greater = complex_data(result, 'correlation_greater')
                        lesser = complex_data(result, 'correlation_lesser')
                        mixed = complex_data(result, 'correlation')
                        for t in range(len(grid)):
                            np.testing.assert_allclose(magnetization[t], [np.trace(rho@s) for s in heisenberg[t]], atol=4e-12)
                            for p, (a, b) in enumerate(labels):
                                a, b = 'xyz'.index(a), 'xyz'.index(b)
                                for s in range(len(grid)):
                                    np.testing.assert_allclose(greater[t, p, s], np.trace(rho @ heisenberg[t][a] @ heisenberg[s][b]), atol=4e-12)
                                    np.testing.assert_allclose(lesser[t, p, s], np.trace(rho @ heisenberg[s][b] @ heisenberg[t][a]), atol=4e-12)
                                for k, imaginary in enumerate(tau):
                                    insertion = np.diag(np.exp(-(1-imaginary)*energy)) @ ops[b] @ np.diag(np.exp(-imaginary*energy))
                                    np.testing.assert_allclose(mixed[t, p, k], np.trace(insertion @ heisenberg[t][a])/weights.sum(), atol=4e-12)
                        for ordering in ('greater', 'lesser'):
                            for part in ('Re', 'Im'):
                                errors = file[f'runtimedata/{part}_correlation_{ordering}_sample_stds'][...]
                                assert np.all(np.isfinite(errors)) and np.all(errors >= 0)

        run('replay', ['--pulseSequence=' + str(replay)])
        empty = write('empty', np.empty(0), np.empty((0, 3)), np.empty(0))
        default = run('default', [])
        empty_output = run('empty_output', ['--pulseSequence=' + str(empty)])
        with h5py.File(default) as a, h5py.File(empty_output) as b:
            for name in a['results']:
                np.testing.assert_array_equal(a['results'][name][...], b['results'][name][...])

        for pcn in (False, True):
            options = ['--pulseSequence=' + str(path), '--JQ=0.1', '--JL=0.1',
                       '--numSamplesPerCore=256', '--preparationIterlimit=30', '--iterlimit=3']
            if pcn:
                options += ['--jointMatsubaraPCN', '--mhStepSize=1', '--mhBurnIn=16']
            output = run('feedback_' + str(pcn), options, status=(0, 3))
            with h5py.File(output) as file:
                assert np.max(file['preparation'].attrs['frozen_MM_residuals']) < 1e-12
                assert np.all(np.isfinite(file['results/Re_correlation_greater'][...]))

        invalid = [([-0.1], [[1, 0, 0]], [1.]), ([1.1], [[1, 0, 0]], [1.]),
                   ([np.nan], [[1, 0, 0]], [1.]), ([0.2], [[0, 0, 0]], [1.]),
                   ([0.2], [[np.inf, 0, 0]], [1.]), ([0.2], [[1, 0, 0]], [np.nan]),
                   ([0.4, 0.2], [[1, 0, 0], [0, 1, 0]], [1., 2.]),
                   ([0.2], [[1, 0]], [1.]), ([0.2], [[1, 0, 0]], [1., 2.]),
                   ([[0.2]], [[1, 0, 0]], [1.]), ([0.2], [[1, 0, 0]], [b'x'])]
        for i, (t, n, phi) in enumerate(invalid):
            bad = write('bad_' + str(i), t, n, phi)
            assert 'pulseSequence' in run('invalid_' + str(i), ['--pulseSequence=' + str(bad)], status=None)
        for i, option in enumerate(['--timeTranslationInvariant=true', '--constantMagnetization', '--cstype=C']):
            run('incompatible_' + str(i), ['--pulseSequence=' + str(path), option], status=None)
        missing = root / 'missing.hdf5'
        run('missing_file', ['--pulseSequence=' + str(missing)], status=None)
        incomplete = write('incomplete', [0.2], [[1, 0, 0]], [1.])
        with h5py.File(incomplete, 'a') as file:
            del file['pulses/axes']
        run('missing_dataset', ['--pulseSequence=' + str(incomplete)], status=None)
    print('Pulse solver/HDF5 checks passed')


if __name__ == '__main__':
    main()
