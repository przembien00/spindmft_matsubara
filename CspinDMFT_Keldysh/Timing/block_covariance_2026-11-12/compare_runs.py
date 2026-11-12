"""Complete matched old/new runs, including stored errors and diagnostic histories.

Usage: python compare_runs.py [--ranks 2]
The baseline executable before.out is captured before editing the solver.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import time

import h5py
import numpy as np

out = Path(__file__).resolve().parent
solver = out.parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--ranks', type=int, default=1)
parser.add_argument('--check', action='store_true')
args = parser.parse_args()
env = dict(os.environ)
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'BLAZE_NUM_THREADS']:
    env[key] = '1'
rows = []


def contents(path):
    data = {}
    with h5py.File(path) as f:
        for stage in ['preparation', 'results', 'runtimedata']:
            group = f[stage]
            for name, value in group.attrs.items():
                data[stage + '/@' + name] = value

            def visit(name, value):
                if isinstance(value, h5py.Dataset):
                    data[stage + '/' + name] = value[()]

            group.visititems(visit)
    return data


for symmetry in 'ACD':
    for strategy in ['independent', 'pcn']:
        for init in ['default', 'isolated']:
            tag = f'{symmetry}_{strategy}_{init}_r{args.ranks}_'+('checked' if args.check else 'default')
            row = dict(tag=tag, runs=[])
            outputs = []
            for version, exe in [('before', out / 'before.out'),
                                 ('after', solver / 'executable_DOUBLE.out')]:
                destination = out / tag / version
                destination.mkdir(parents=True, exist_ok=True)
                config = solver / 'Configuration_Data/CubicLattice/Cubic_3D_N=2_NN_J=-0.40825.hdf5'
                command = [str(exe), f'--config={config}', '--project=CubicLattice',
                    '--dstproject=Results', '--beta=.5', '--numImagTimeSteps=6',
                    '--numRealTimeSteps=5', '--realTimeSubsteps=2', '--Tmax=.3',
                    '--numSamplesPerCore=32', '--preparationSamplesPerCore=32',
                    '--numBlocks=8', '--numSamplesPerSet=3', '--seed=613',
                    f'--cstype={symmetry}', f'--samplingStrategy={strategy}',
                    '--unifield=0' if symmetry == 'A' else '--unifield=.2',
                    '--preparationIterlimit=12', '--iterlimit=6', '--itermode=absolute',
                    '--absiterror=.08', '--mhBurnIn=8']
                if version == 'after' and args.check:
                    command += ['--checkGaussianReconstruction']
                if init == 'isolated':
                    command += ['--initdcorr=isolated', '--initndcorr=isolated']
                if args.ranks > 1:
                    command = ['mpirun', '-n', str(args.ranks)] + command
                start = time.monotonic()
                completed = subprocess.run(command, cwd=destination, env=env,
                    capture_output=True, text=True, timeout=180)
                (destination / 'run.log').write_text(completed.stdout + completed.stderr)
                row['runs'].append(dict(version=version, command=command,
                    returncode=completed.returncode, seconds=time.monotonic()-start,
                    sha256=hashlib.sha256(exe.read_bytes()).hexdigest()))
                assert completed.returncode in (0, 3), (tag, version, completed.stdout, completed.stderr)
                files = list(destination.rglob('*.hdf5'))
                outputs.append(contents(max(files, key=lambda p: p.stat().st_mtime_ns)))
            assert row['runs'][0]['returncode'] == row['runs'][1]['returncode']
            assert outputs[0].keys() == outputs[1].keys()
            maximum = 0.
            differences = []
            for key in outputs[0]:
                before, after = np.asarray(outputs[0][key]), np.asarray(outputs[1][key])
                if not args.check and key.endswith('/@gaussian_factor_reconstruction_errors'):
                    assert np.all(np.isnan(after))
                    continue
                if before.dtype.kind in 'fc':
                    np.testing.assert_allclose(after, before, rtol=2e-10, atol=2e-12, err_msg=tag+': '+key)
                    finite = np.isfinite(before) & np.isfinite(after)
                    delta = float(np.max(np.abs(after[finite]-before[finite]))) if finite.any() else 0.
                    maximum = max(maximum, delta)
                    if delta:
                        differences.append([key, delta])
                else:
                    np.testing.assert_array_equal(after, before, err_msg=tag+': '+key)
            row.update(compared=len(outputs[0])-(0 if args.check else 2),
                       skipped_reconstruction_histories=0 if args.check else 2, max_absolute_difference=maximum,
                       differences=differences)
            rows.append(row)
            (out / (f'comparison_r{args.ranks}_'+('checked' if args.check else 'default')+'.json')).write_text(json.dumps(rows, indent=2)+'\n')
            print(tag, row['compared'], maximum, flush=True)
