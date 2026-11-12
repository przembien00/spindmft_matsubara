from pathlib import Path
import os,subprocess,time,json,hashlib
out=Path(__file__).resolve().parent
root=out.parents[1]
env=dict(os.environ)
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','BLAZE_NUM_THREADS']: env[key]='1'
rows=[]
for n,sym in [(2,'A'),(2,'C'),(4,'A'),(4,'C')]:
 tag=f'N{n}_{sym}'
 cmd=[str(out/'probe.out'),f'--config={root}/Configuration_Data/CubicLattice/Cubic_3D_N={n}_NN_J=-0.40825.hdf5','--project=CubicLattice','--dstproject=Efficiency_audit_2026-11-12',f'--cstype={sym}',f'--unifield={0 if sym=="A" else .5}','--beta=5','--numImagTimeSteps=50','--numRealTimeSteps=150','--Tmax=10','--realTimeSubsteps=1','--numSamplesPerCore=4','--preparationSamplesPerCore=4','--numBlocks=2','--noselfcons','--seed=613',f'--fileext={tag}']
 print('Starting',tag,flush=True)
 start=time.monotonic()
 with (out/f'{tag}.log').open('w') as log:
  result=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
 rows.append(dict(tag=tag,command=cmd,returncode=result.returncode,elapsed_seconds=time.monotonic()-start,executable_sha256=hashlib.sha256((out/'probe.out').read_bytes()).hexdigest()))
 (out/'runs.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(tag,result.returncode,rows[-1]['elapsed_seconds'],flush=True)
 if result.returncode: raise SystemExit(result.returncode)
