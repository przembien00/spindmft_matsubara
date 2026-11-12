"""Fresh-process matched audit-grid timings and peak RSS; no probe overhead."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import time

out = Path(__file__).resolve().parent
solver = out.parents[1]
env = dict(os.environ)
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'BLAZE_NUM_THREADS']:
    env[key] = '1'
parser = argparse.ArgumentParser()
parser.add_argument('--after-only', action='store_true')
args = parser.parse_args()
rows = ([row for row in json.loads((out/'benchmarks.json').read_text()) if row['version']=='before']
        if args.after_only else [])
for n, symmetry in [(2, 'A'), (2, 'C'), (4, 'A'), (4, 'C')]:
    for version, exe in [('before', out/'before.out'), ('after', solver/'executable_DOUBLE.out')]:
        if args.after_only and version == 'before':
            continue
        destination = out/'benchmarks'/f'N{n}_{symmetry}'/version
        destination.mkdir(parents=True, exist_ok=True)
        command = [str(exe),
            f'--config={solver}/Configuration_Data/CubicLattice/Cubic_3D_N={n}_NN_J=-0.40825.hdf5',
            '--project=CubicLattice', '--dstproject=Results', f'--cstype={symmetry}',
            '--unifield=0' if symmetry == 'A' else '--unifield=.5', '--beta=5',
            '--numImagTimeSteps=50', '--numRealTimeSteps=150', '--Tmax=10',
            '--realTimeSubsteps=1', '--numSamplesPerCore=32', '--preparationSamplesPerCore=32',
            '--numBlocks=8', '--numSamplesPerSet=1', '--noselfcons', '--seed=613']
        print('Starting', n, symmetry, version, flush=True)
        executable_hash = hashlib.sha256(exe.read_bytes()).hexdigest()
        start = time.monotonic()
        result = subprocess.run(['/usr/bin/time', '-l']+command, cwd=destination,
            env=env, capture_output=True, text=True, timeout=900)
        (destination/'run.log').write_text(result.stdout+result.stderr)
        assert result.returncode == 0, (command, result.returncode, result.stderr)
        row = dict(sites=n, symmetry=symmetry, version=version,
            returncode=result.returncode, wall_seconds=time.monotonic()-start,
            executable_sha256=executable_hash, command=command)
        match = re.search(r'(\d+)\s+maximum resident set size', result.stderr)
        row['peak_rss_bytes'] = int(match[1]) if match else None
        for label in ['initialization', 'Matsubara field distribution', 'Matsubara sampling',
                      'joint field distribution', 'mixed real-time sampling']:
            match = re.search(re.escape(label)+r' done \(0\).*?it took me ([\deE.+-]+) seconds', result.stdout)
            row[label] = float(match[1]) if match else None
        rows.append(row)
        (out/'benchmarks.json').write_text(json.dumps(rows, indent=2)+'\n')
        print(row['wall_seconds'], row['peak_rss_bytes'], flush=True)
