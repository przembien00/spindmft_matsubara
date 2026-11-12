from pathlib import Path
import re,json
out=Path(__file__).resolve().parent
rows=[]
for n,s in [(2,'A'),(2,'C'),(4,'A'),(4,'C')]:
 text=(out/f'N{n}_{s}.log').read_text()
 audit={}
 for line in text.splitlines():
  if line.startswith('AUDIT '):
   d={k:float(v) for k,v in re.findall(r'(\w+)=([\deE.+-]+)',line)}
   if 'batch' in d:audit.setdefault('batches',{})[int(d.pop('batch'))]=d
   else:audit.update(d)
 audit.update(sites=n,symmetry=s,setup_seconds=float(re.search(r'joint field distribution done .*?it took me ([\deE.+-]+) seconds',text)[1]))
 rows.append(audit)
(out/'measurements.json').write_text(json.dumps(rows,indent=2)+'\n')
table=['## Current measurements','', '| N / symmetry | Joint setup s | Draw ms/sample, batch 1 | Draw ms/sample, batch 32 | Trajectory ms/sample | Peak RSS GB |','|---|---:|---:|---:|---:|---:|']
for r in rows:
 table.append(f'| {r["sites"]} / {r["symmetry"]} | {r["setup_seconds"]:.3f} | {1000*r["batches"][1]["draw_seconds_per_sample"]:.4f} | {1000*r["batches"][32]["draw_seconds_per_sample"]:.4f} | {1000*r["measure_seconds_per_sample"]:.3f} | {r["peak_rss_platform_units"]/1e9:.3f} |')
table+=['','The four-spin C trajectory breaks down into 3.77 ms imaginary evolution, 4.54 ms Matsubara insertions, 22.31 ms real propagation, and **57.77 ms observable evaluation**. Thus optimize the measurement/propagation path before expecting large overall gains from field batching.','', 'For N=2 C, batching 32 reduces the draw cost about 7.65x; for N=4 C, about 9.53x. Combining measured draw and trajectory costs predicts only approximately **16.5%** and **1.9%** reductions, respectively, in those two per-sample components. These exclude setup, statistics and I/O, and are not measured whole-run speedups.','', 'Across all batches the maximum fixed-factor field difference was **6.01e-15**. The fresh C probe RSS was 2.20 GB for four spins; an earlier cold-cache probe reached 2.45 GB. Treat this as observed allocation/RSS variation, not a measured optimization. The older large-grid 4.74 GB result is not a current-grid measurement.','', '| N / symmetry | Field dimension | Retained latent dimension | Largest complex block |','|---|---:|---:|---:|']
for r in rows:table.append(f'| {r["sites"]} / {r["symmetry"]} | {int(r["dimension"])} | {int(r["latent"])} | {int(r["largest_block"])} |')
p=out/'REPORT.md';p.write_text(p.read_text().replace('<!-- LIVE_RESULTS -->','\n'.join(table)))
