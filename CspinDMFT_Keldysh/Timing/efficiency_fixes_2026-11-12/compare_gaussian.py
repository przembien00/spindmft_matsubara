"""Compare the old/new real-linear maps, both covariances, and fixed-seed draws."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess

import numpy as np

out = Path(__file__).resolve().parent
repo = out.parents[2]
parser = argparse.ArgumentParser()
parser.add_argument('--skip-build', action='store_true')
args = parser.parse_args()
env = dict(os.environ, VECLIB_MAXIMUM_THREADS='1', OMP_NUM_THREADS='1')
results = []
for version, source in [('before', out/'shared_before.cpp'),
                        ('after', repo/'cpp_libs/Complex_Gaussian/Complex_Gaussian.cpp')]:
    exe = out/f'law_{version}.out'
    if not args.skip_build:
        subprocess.run(['/usr/bin/c++', '-O3', '-DNDEBUG', '-std=c++17', '-DUSE_DOUBLE',
            '-I', str(repo/'cpp_libs'), '-I', str(repo/'cpp_libs/Complex_Gaussian'),
            '-isystem', '/opt/homebrew/include', str(out/'gaussian_law.cpp'),
            str(source), '-framework', 'Accelerate', '-o', str(exe)], check=True)
    path = out/f'gaussian_{version}.txt'
    with path.open('w') as handle:
        subprocess.run([str(exe)], stdout=handle, env=env, check=True)
    results.append(np.loadtxt(path))
np.testing.assert_allclose(results[1], results[0], rtol=1e-12, atol=1e-13)
report = dict(objectives=['unweighted', 'explicit weights', 'contour weights'],
    fixtures='dense complex, repeated blocks, global cutoff, zero modes, tolerated asymmetry',
    compared=['real-linear factor', 'E[z z^T]', 'E[z z^dagger]', '19 fixed-seed draws'],
    cases=18, underflow_rejection_checked=True, exactly_equal=bool(np.array_equal(*results)),
    max_absolute_difference=float(np.max(np.abs(results[1]-results[0]))),
    source_sha256={version: hashlib.sha256(source.read_bytes()).hexdigest()
        for version, source in [('before', out/'shared_before.cpp'),
            ('after', repo/'cpp_libs/Complex_Gaussian/Complex_Gaussian.cpp')]})
(out/'gaussian_comparison.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
