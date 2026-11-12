"""Verify benchmark outputs and populate the implementation report/manifest."""
from pathlib import Path
import hashlib
import json

import h5py
import numpy as np

out = Path(__file__).resolve().parent
repo = out.parents[2]
rows = json.loads((out/'benchmarks.json').read_text())
assert len(rows) == 8
table = ['| Sites / symmetry | Setup s, old → new | Sampling s, old → new | Total wall s, old → new | Peak RSS GB, old → new |',
         '|---|---:|---:|---:|---:|']
comparisons = []


def contents(path):
    data = {}
    with h5py.File(path) as f:
        for stage in ['preparation', 'results', 'runtimedata']:
            group = f[stage]
            data.update({stage+'/@'+name: value for name, value in group.attrs.items()})

            def visit(name, value):
                if isinstance(value, h5py.Dataset):
                    data[stage+'/'+name] = value[()]

            group.visititems(visit)
    return data


for n, symmetry in [(2, 'A'), (2, 'C'), (4, 'A'), (4, 'C')]:
    before, after = [next(r for r in rows if (r['sites'],r['symmetry'],r['version']) == (n,symmetry,v))
                     for v in ['before', 'after']]
    assert before['returncode'] == after['returncode'] == 0
    table.append(f"| {n} / {symmetry} | {before['joint field distribution']:.3f} → {after['joint field distribution']:.3f} | "
        f"{before['mixed real-time sampling']:.4f} → {after['mixed real-time sampling']:.4f} | "
        f"{before['wall_seconds']:.3f} → {after['wall_seconds']:.3f} | "
        f"{before['peak_rss_bytes']/1e9:.3f} → {after['peak_rss_bytes']/1e9:.3f} |")
    outputs = []
    for version in ['before', 'after']:
        files = (out/'benchmarks'/f'N{n}_{symmetry}'/version).rglob('*.hdf5')
        outputs.append(contents(max(files, key=lambda p: p.stat().st_mtime_ns)))
    assert outputs[0].keys() == outputs[1].keys()
    maximum = 0.
    for key in outputs[0]:
        a, b = np.asarray(outputs[0][key]), np.asarray(outputs[1][key])
        if a.dtype.kind in 'fc':
            np.testing.assert_allclose(a, b, rtol=2e-10, atol=2e-12, err_msg=f'{n}/{symmetry}: {key}')
            finite = np.isfinite(a) & np.isfinite(b)
            if finite.any():
                maximum = max(maximum, float(np.max(np.abs(a[finite]-b[finite]))))
        else:
            np.testing.assert_array_equal(a,b)
    comparisons.append(dict(sites=n,symmetry=symmetry,compared=len(outputs[0]),max_absolute_difference=maximum))
(out/'benchmark_comparison.json').write_text(json.dumps(comparisons,indent=2)+'\n')
report = out/'REPORT.md'
text = report.read_text()
start = text.index('<!-- BENCHMARK_TABLE -->')
end = text.index('\nSampling is ',start)
text = text[:start]+'<!-- BENCHMARK_TABLE -->\n\n'+'\n'.join(table)+'\n'+text[end:]
report.write_text(text)
paths = [repo/'cpp_libs/Complex_Gaussian/Complex_Gaussian.cpp',
         repo/'CspinDMFT_Keldysh/README.md']
paths += sorted((repo/'CspinDMFT_Keldysh/Algorithm').rglob('*.cpp'))
paths += sorted((repo/'CspinDMFT_Keldysh/Algorithm').rglob('*.h'))
paths = [p for p in paths if 'build' not in p.parts]
paths += [repo/'CspinDMFT_Keldysh/Algorithm/CMakeLists.txt_']
manifest = dict(date='2026-11-12',
    source_sha256={str(p.relative_to(repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    executable_sha256=hashlib.sha256((repo/'CspinDMFT_Keldysh/executable_DOUBLE.out').read_bytes()).hexdigest(),
    baseline_sha256=hashlib.sha256((out/'before.out').read_bytes()).hexdigest())
assert all(r['executable_sha256']==manifest['executable_sha256'] for r in rows if r['version']=='after')
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('\n'.join(table))
print(json.dumps(comparisons,indent=2))
