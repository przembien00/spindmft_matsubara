"""Exact physics and native CLI/storage regressions; optionally pass mpirun -n 2."""
from pathlib import Path
import subprocess
import sys
import re
import tempfile
import h5py
import numpy as np

exe=str(Path(sys.argv[1]).resolve())
launcher=sys.argv[2:]

def config(path,bath=False):
    path.parent.mkdir(parents=True,exist_ok=True)
    with h5py.File(path,'w') as f:
        f.attrs['num_Spins']=2
        f['spin-spin couplings']=[[0.,.7],[.7,0.]]
        f['expectation weights']=np.zeros((2,2))
        w=np.zeros((2,2,2,2))
        if bath:
            w[0,0,0,0]=w[1,1,1,1]=.03
            w[0,1,0,1]=w[1,0,0,1]=.005
        f['correlation weights']=w
        f['correlation_categories']=np.array([[0,0],[0,1],[1,1]],dtype=np.uint32)

def run(root,tag,extra=(),code=0):
    options=dict(config='dimer',project='TestProject',dstproject='Results',fileext=tag,
        spinspinmodel='ISO',beta='1.7',numImagTimeSteps='4',numRealTimeSteps='3',
        realTimeSubsteps='2',Tmax='.09',numSamplesPerCore='24',numBlocks='4',
        preparationIterlimit='3',iterlimit='3',itermode='absolute',absiterror='1e-9',seed='12345',initndcorr='isolated',initdcorr='isolated')
    flags=[]
    for item in extra:
        key,sep,value=item.removeprefix('--').partition('=') if sys.version_info>=(3,9) else item[2:].partition('=')
        if sep: options[key]=value
        else: flags.append(item)
    if not any(key in options for key in ('unifield','stagfield','chemshift')): options['unifield']='.4'
    args=launcher+[exe]+['--'+key+'='+str(value) for key,value in options.items()]+flags
    before=set(root.rglob('*.hdf5'))
    completed=subprocess.run(args,cwd=root,capture_output=True,text=True,timeout=120)
    assert (completed.returncode!=0 if code==1 else completed.returncode==code),(args,completed.returncode,completed.stdout,completed.stderr)
    created=set(root.rglob('*.hdf5'))-before
    if code==1:
        assert not created
        return None
    assert len(created)==1,(created,completed.stdout)
    output=created.pop()
    with h5py.File(output) as f:
        parameters=f['parameters'].attrs
        assert not any('(original)' in name for name in parameters)
        for removed in ['algorithm','**algorithm','normalization_strategy','real_time_evolution']:
            assert removed not in parameters
        assert text(parameters['filename'])==output.name
        configuration=Path(options['config']).name
        if Path(configuration).suffix not in ['.hdf5','.h5']: configuration+='.hdf5'
        assert text(parameters['configuration_src_file'])==configuration
        assert text(parameters['project_name'])==options['project']
        assert text(parameters['destination project_name'])==options['dstproject']
        assert text(parameters['filename_extension'])==options['fileext']
        assert text(parameters['correlation_symmetry_type'])==options.get('cstype',options.get('stype','D'))
        assert f['parameters'].attrs.get_id('correlation_symmetry_type').dtype.kind=='S'
        assert 'imaginary_time_edges' not in f['results'].attrs
        assert 'real_times' not in f['results'].attrs
        assert 'num_samples_per_core_history' not in f['preparation'].attrs
        assert 'real_time_started' not in f['preparation'].attrs
        terminal=re.sub(r'\x1b\[[0-9;]*m','',completed.stdout)
        category_key='correlation_categories, stored according to the hierarchy n, i-j'
        pairs=','.join(f'{int(i)+1}-{int(j)+1}' for i,j in f['parameters'][category_key][:])
        assert terminal.count('The following correlations enter the self consistency (numbering starts at 1): '+pairs+'\n')==1
        assert terminal.count('ID of this Run: ')==1
        for label in ['mf_model','config_file','num_Spins','num_HilbertSpaceDimension','uncoupled_spins',
                      'spin_list','spinspin_cmodel info','spinmf_cmodel info','chemical_shift info',
                      'rescale_meanfield','rescale_spinspin','num_ImagTimeSteps','delta_real_t',
                      'num_PreparationSamplesPerCore','sampling_strategy','staged_field_sampling','delta_imag_t']:
            assert len(re.findall(r'^'+re.escape(label)+r'\s*:',terminal,re.M))==1,(label,terminal)
        prep_iterations=int(f['preparation'].attrs['num_Iterations'])
        real_iterations=int(f['runtimedata'].attrs['num_Iterations'])
        iterations=prep_iterations+real_iterations
        assert not re.search(r'^delta_t\s*:',terminal,re.M)
        assert not re.search(r'^propagator\s*:',terminal,re.M)
        assert 'staggered-strang' not in terminal
        assert terminal.count('['+'='*50+'] 100%')==iterations,terminal
        estimate='estimated duration of '+('pCN simulation' if options.get('samplingStrategy')=='pcn' else 'Monte-Carlo simulation')
        assert terminal.count(estimate)==iterations,terminal
        for label in ['Iteration step finished:','largest Gaussian factor block',
                      'closed-contour trace residual','MPI communication done (0)',
                      'iteration-step finalization done (0)']:
            assert terminal.count(label)==iterations,(label,terminal)
        assert terminal.count('--- Matsubara iteration ')==prep_iterations
        assert terminal.count('--- Mixed/real-time iteration ')==real_iterations
        for label in ['initialization done (0)','Matsubara self consistency done (0)',
                      'storing done (0)','everything done (0)']:
            assert terminal.count(label)==1,(label,terminal)
        if real_iterations: assert terminal.count('real-time self consistency done (0)')==1
        if options.get('samplingStrategy')=='pcn':
            assert terminal.count('pcn_target')==1
            assert terminal.count('MH acceptance rate')==iterations
            assert terminal.count('zero-modulus rejection rate')==iterations
        else:
            assert 'pcn_target' not in terminal and 'MH acceptance rate' not in terminal
        assert ('J(0,1)' in terminal)==(int(parameters['num_Spins'])>=5)
    return output

def correlations(group):
    # Dense view for comparisons only; unmeasured pairs stay zero and are
    # never asserted to be reconstructed from a stochastic trajectory.
    out=np.zeros((4,6,6,5),complex)
    for pair in group['Re_correlation']:
        i,j=(int(site)-1 for site in pair.split('-'))
        re=group['Re_correlation'][pair][:]
        im=group['Im_correlation'][pair][:]
        out[:,3*i:3*i+3,3*j:3*j+3,:]=(re+1j*im).reshape(4,3,3,5)
    return out

def text(value):
    return value.decode() if isinstance(value,bytes) else value

with tempfile.TemporaryDirectory(prefix='cluster-native-') as directory:
    root=Path(directory)
    conf=root/'Configuration_Data/TestProject/dimer.hdf5'
    config(conf)
    exact=run(root,'exact')
    assert exact.parent==root/'Data/Results'
    assert exact.name.startswith('spinmodel=ISO__config=dimer__beta=1.7__unifield=0.4')
    assert '__prop=' not in exact.name
    assert '__substeps=2_exact.hdf5' in exact.name
    s=[np.array([[0,1],[1,0]],complex)/2,np.array([[0,-1j],[1j,0]],complex)/2,np.diag([.5,-.5])]
    ops=[np.kron(v,np.eye(2)) for v in s]+[np.kron(np.eye(2),v) for v in s]
    H=.7*sum(ops[a]@ops[3+a] for a in range(3))+.4*(ops[2]+ops[5])
    e,v=np.linalg.eigh(H)
    def exp(factor): return (v*np.exp(factor*e))@v.conj().T
    rho=exp(-1.7);z=np.trace(rho)
    with h5py.File(exact) as f:
        assert set(f)=={'parameters','runtimedata','results','preparation','timedata'}
        p=f['parameters'].attrs
        assert p['num_Cores']==(int(launcher[-1]) if launcher else 1)
        assert p['num_ImagTimePoints']==5 and p['num_RealTimePoints']==4
        assert text(p['adaptive_sample_size'])=='no'
        assert p['unifield']==.4 and 'magnetic_field' not in p
        np.testing.assert_array_equal(np.fromstring(text(p['chemical_shift']),sep=','),[.4,.4])
        np.testing.assert_array_equal(f['runtimedata'].attrs['num_samples_per_core_history'],24.)
        assert text(f['runtimedata'].attrs['termination'])=='by convergence'
        assert list(f['results/correlation_direction_labels'][:])==[b'xx',b'xy',b'xz',b'yx',b'yy',b'yz',b'zx',b'zy',b'zz']
        expected_pairs={'1-1','1-2','2-2'}
        assert set(f['results/Re_correlation'])==expected_pairs
        assert set(f['preparation/Re_correlation'])==expected_pairs
        assert set(f['runtimedata/Re_correlation_sample_stds'])==expected_pairs
        assert 'configuration_correlation_categories' not in f['parameters']
        np.testing.assert_array_equal(f['parameters/correlation_categories, stored according to the hierarchy n, i-j'][:],[[0,0],[0,1],[1,1]])
        assert f['preparation/Re_correlation/1-2'].shape==(9,5)
        assert f['results/Re_correlation/1-2'].shape==(4,9,5)
        values=correlations(f['results'])
        for it,t in enumerate(np.arange(int(p['num_RealTimePoints']))*p['delta_real_t']):
            for a in range(6):
                sa=exp(1j*t)@ops[a]@exp(-1j*t)
                for b in range(6):
                    if f'{a//3+1}-{b//3+1}' not in f['results/Re_correlation']: continue
                    for k,tau in enumerate(np.arange(int(p['num_ImagTimePoints']))*p['delta_imag_t']):
                        reference=np.trace(exp(-(1.7-tau))@ops[b]@exp(-tau)@sa)/z
                        np.testing.assert_allclose(values[it,a,b,k],reference,atol=3e-12)
        assert np.max(np.abs(values[...,0]-values[...,-1]))>.01
        for data in f['runtimedata/Re_correlation_sample_stds'].values():
            if isinstance(data,h5py.Dataset): assert np.max(data[:])<1e-12
        np.testing.assert_allclose(f['runtimedata'].attrs['closed_contour_ratio_Re'],1.,atol=1e-12)
        original=exact.read_bytes()
    again=run(root,'exact')
    assert again.stem==exact.stem+'X'
    assert exact.read_bytes()==original
    pcn_exact=run(root,'pcn_exact',['--samplingStrategy=pcn','--mhBurnIn=4'])
    with h5py.File(pcn_exact) as f:
        np.testing.assert_allclose(correlations(f['results']),values,atol=3e-12)
        assert f['runtimedata'].attrs['mh_acceptance_rates'][-1]==1.
    # Both field helpers reproduce explicit site fields and exact thermal moments.
    for option,fields in [('unifield',[.4,.4]),('unifield',[-.4,-.4]),
                          ('unifield',[0.,0.]),('stagfield',[.4,-.4]),('stagfield',[-.4,.4])]:
        tag=option+str(fields[0])
        field_run=run(root,tag,['--'+option+'='+str(fields[0])])
        explicit=run(root,tag+'_explicit',['--chemshift='+','.join(map(str,fields))])
        h=.7*sum(ops[a]@ops[3+a] for a in range(3))+fields[0]*ops[2]+fields[1]*ops[5]
        eigen,basis=np.linalg.eigh(h)
        density=(basis*np.exp(-1.7*eigen))@basis.conj().T; density/=np.trace(density)
        reference=np.array([np.trace(density@op) for op in ops]).reshape(2,3)
        with h5py.File(field_run) as f,h5py.File(explicit) as g:
            np.testing.assert_allclose(correlations(f['results']),correlations(g['results']),atol=3e-12)
            actual=f['results/Re_magnetization'][:]+1j*f['results/Im_magnetization'][:]
            np.testing.assert_allclose(actual,np.repeat(reference[None],4,axis=0),atol=3e-12)
            assert f['parameters'].attrs[option]==fields[0]
            np.testing.assert_array_equal(np.fromstring(text(f['parameters'].attrs['chemical_shift']),sep=','),fields)
            assert '__'+option+'=' in field_run.name
    for extra in [
        ['--unifield=.1','--chemshift=.1,-.2'],
        ['--stagfield=.1','--chemshift=.1,-.2'],
        ['--unifield=0','--stagfield=0'],
        ['--unifield=nan'],['--stagfield=inf'],
        ['--Bname=z'],['--Babs=.4'],['--Bphi=0'],['--Btheta=0']]:
        run(root,'invalid_field',extra,code=1)
    # Match the finite-T graph convention: color from site 1, independent of rescss.
    for name,couplings in [
        ('chain',[[0.,.7,0.],[.7,0.,-.2],[0.,-.2,0.]]),
        ('odd_cycle',[[0.,.7,.3],[.7,0.,-.2],[.3,-.2,0.]]),
        ('disconnected',[[0.,.7,0.],[.7,0.,0.],[0.,0.,0.]])]:
        path=conf.with_name(name+'.hdf5')
        with h5py.File(path,'w') as f:
            f.attrs['num_Spins']=3
            f['spin-spin couplings']=couplings
            f['expectation weights']=np.zeros((3,3))
            f['correlation weights']=np.zeros((3,3,3,3))
        extra=['--config='+name,'--stagfield=.2','--rescss=0']
        if name!='chain':
            run(root,name,extra,code=1)
        else:
            chain=run(root,name,extra)
            with h5py.File(chain) as f:
                np.testing.assert_array_equal(np.fromstring(text(f['parameters'].attrs['chemical_shift']),sep=','),[.2,-.2,.2])
                moments=f['results/Re_magnetization'][:, :, 2]
                expected=-.5*np.tanh(1.7*np.array([.2,-.2,.2])/2.)
                np.testing.assert_allclose(moments,np.repeat(expected[None],4,axis=0),atol=3e-12)
    # Initialization can round-trip the native per-site-pair preparation tensors.
    imported=run(root,'imported',['--loadinit','--initcorrfile='+str(exact)])
    with h5py.File(imported) as f: np.testing.assert_allclose(correlations(f['results']),values,atol=3e-12)
    run(root,'invalid',['--numImagTimeSteps=2'],code=1)
    run(root,'invalid',['--timeTranslationInvariant=false'],code=1)
    run(root,'invalid',['--samplingStrategy=pcn','--jointMatsubaraPCN=false'],code=1)
    # Cluster XXZ/XY and local fields/chemical shifts are actual Hamiltonian inputs.
    anisotropic=run(root,'anisotropic',['--spinspinmodel=XY','--lambda=.6','--rho=.2',
        '--spinmfmodel=XXZ','--smflambda=1.3','--chemshift=.1,-.2','--extraint=quadrupolar','--intstrength=.08'])
    with h5py.File(anisotropic) as f:
        actual=f['results/Re_magnetization'][:]+1j*f['results/Im_magnetization'][:]
        h=.7*(ops[0]@ops[3]+ops[1]@ops[4]+.6*ops[2]@ops[5]+.2*(ops[0]@ops[4]+ops[1]@ops[3]))
        h+=.1*ops[2]-.2*ops[5]+3*.08*(ops[2]@ops[2]+ops[5]@ops[5])
        eigen,basis=np.linalg.eigh(h); density=(basis*np.exp(-1.7*eigen))@basis.conj().T; density/=np.trace(density)
        reference=np.array([np.trace(density@op) for op in ops]).reshape(2,3)
        np.testing.assert_allclose(actual,np.repeat(reference[None],4,axis=0),atol=3e-12)
    # A/B/C containers and the saved tensors use the same reduced basis.
    for symmetry,labels in [('A',[b'xx']),('B',[b'xx',b'zz']),('C',[b'xx',b'xy',b'yx',b'zz'])]:
        options=['--cstype='+symmetry]
        if symmetry=='A': options.append('--unifield=0')
        reduced=run(root,'reduced_'+symmetry,options)
        with h5py.File(reduced) as f:
            assert list(f['results/correlation_direction_labels'][:])==labels
            assert f['preparation/Re_correlation/1-1'].shape==(len(labels),5)
            assert f['results/Re_correlation/1-1'].shape==(4,len(labels),5)
            magnetization=f['results/Re_magnetization'][:]
            np.testing.assert_allclose(magnetization[...,:2],0.,atol=0.)
            if symmetry=='A': np.testing.assert_allclose(magnetization,0.,atol=0.)
            else: assert np.max(np.abs(magnetization[...,2]))>.01
    config(conf,bath=True)
    for strategy in ['independent','pcn']:
        bath=run(root,'bath_'+strategy,['--samplingStrategy='+strategy,'--reliterror=100','--absiterror=1',
            '--numSamplesPerCore=48','--mhBurnIn=8','--preparationSamplesPerCore=24','--numSamplesPerSet=3'])
        with h5py.File(bath) as f:
            assert np.all(np.isfinite(correlations(f['results'])))
            assert np.max(f['runtimedata/Re_correlation_sample_stds/1-1'][:])>1e-8
            assert text(f['parameters'].attrs['check_gaussian_reconstruction'])=='no'
            assert np.all(np.isnan(f['runtimedata'].attrs['gaussian_factor_reconstruction_errors']))
            assert 0<f['runtimedata'].attrs['average_phase_magnitudes'][-1]<=1+1e-12
            assert 0<=f['runtimedata'].attrs['mh_acceptance_rates'][-1]<=1
            assert f['preparation'].attrs['samples_per_block']==6
            assert f['runtimedata'].attrs['samples_per_block']==12
            # Saving must retain the measured origin, while the preparation
            # statistics remain separate from the last real-time statistics.
            assert np.max(np.abs(f['results/Re_correlation/1-1'][0]-f['preparation/Re_correlation/1-1'][:]))>1e-10
            assert np.max(f['preparation/Re_correlation_sample_stds/1-1'][:])>1e-8
            for i in range(2):
                diagonal=f['preparation/Re_correlation'][f'{i+1}-{i+1}'][:]
                np.testing.assert_allclose(diagonal[[0,4,8],0],.25,atol=1e-12)
    for strategy in ['independent','pcn']:
        base=['--samplingStrategy='+strategy,'--reliterror=100','--absiterror=1',
            '--numSamplesPerCore=48','--mhBurnIn=8','--preparationSamplesPerCore=24','--numSamplesPerSet=3']
        unchecked=run(root,'unchecked_'+strategy,base)
        checked=run(root,'checked_'+strategy,base+['--checkGaussianReconstruction'])
        with h5py.File(unchecked) as f,h5py.File(checked) as g:
            assert text(g['parameters'].attrs['check_gaussian_reconstruction'])=='yes'
            for stage in ['preparation','runtimedata']:
                assert np.all(np.isnan(f[stage].attrs['gaussian_factor_reconstruction_errors']))
                assert np.all(g[stage].attrs['gaussian_factor_reconstruction_errors']<1e-8)
                assert f[stage].attrs['termination']==g[stage].attrs['termination']
                assert f[stage].attrs['num_Iterations']==g[stage].attrs['num_Iterations']
            for stage in ['preparation','results']:
                for part in ['Re','Im']:
                    for pair in f[stage+'/'+part+'_correlation']:
                        np.testing.assert_array_equal(f[stage+'/'+part+'_correlation/'+pair][:],g[stage+'/'+part+'_correlation/'+pair][:])
    reduced_config=conf.with_name('representatives.hdf5')
    with h5py.File(conf) as source,h5py.File(reduced_config,'w') as target:
        target.attrs['num_Spins']=2
        for key in ['spin-spin couplings','expectation weights']:
            target[key]=source[key][:]
        weights=source['correlation weights'][:]
        weights[1,1,0,0]=weights[1,1,1,1]; weights[1,1,1,1]=0.
        target['correlation weights']=weights
        target['correlation_categories']=np.array([[0,1],[0,0]],dtype=np.uint32)
    for strategy in ['independent','pcn']:
        sparse=run(root,'representatives_'+strategy,['--config=representatives','--samplingStrategy='+strategy,
            '--absiterror=1','--mhBurnIn=8'])
        with h5py.File(sparse) as f:
            for stage in ['preparation','results']:
                for part in ['Re','Im']:
                    assert set(f[stage+'/'+part+'_correlation'])=={'1-1','1-2'}
            assert set(f['runtimedata/Re_correlation_sample_stds'])=={'1-1','1-2'}
            np.testing.assert_array_equal(f['parameters/correlation_categories, stored according to the hierarchy n, i-j'][:],[[0,1],[0,0]])
            assert np.all(np.isfinite(correlations(f['results'])))
        run(root,'sparse_import_'+strategy,['--config=representatives','--loadinit','--initcorrfile='+str(sparse),'--absiterror=1'])
    with h5py.File(reduced_config,'r+') as f:
        weights=f['correlation weights'][:]
        weights[:,:,1,0]=weights[:,:,0,1]
        f['correlation weights'][:]=weights
    mirrored=run(root,'symmetric_weights',['--config=representatives','--samplingStrategy=pcn',
        '--absiterror=1','--mhBurnIn=8'])
    with h5py.File(sparse) as reference,h5py.File(mirrored) as actual:
        for stage in ['preparation','results']:
            for part in ['Re','Im']:
                group=stage+'/'+part+'_correlation'
                for pair in reference[group]:
                    np.testing.assert_array_equal(reference[group+'/'+pair][:],actual[group+'/'+pair][:])
    with h5py.File(reduced_config,'r+') as f:
        weights=f['correlation weights'][:]; weights[0,0,1,1]=.01
        f['correlation weights'][:]=weights
    run(root,'unlisted_category_weight',['--config=representatives'],code=1)
    stopped=run(root,'stopped',['--absiterror=1e-15','--preparationIterlimit=1'],code=2)
    assert stopped.parent.name=='NOT_CONVERGED' and '__preparation_failed' in stopped.name
    with h5py.File(stopped) as f:
        assert f['runtimedata'].attrs['num_Iterations']==0
        assert text(f['preparation'].attrs['termination'])=='by iteration limit'
    harmonic=run(root,'harmonic',['--bath=harmonic','--bathCoupling=.1','--bathComponent=z'])
    with h5py.File(harmonic) as f:
        assert text(f['runtimedata'].attrs['termination'])=='prescribed harmonic bath'
        assert np.all(np.isfinite(correlations(f['results'])))
    run(root,'removed_adaptive',['--adaptive'],code=1)
    run(root,'removed_staterrtolerance',['--staterrtolerance=.1'],code=1)
    refined=run(root,'refined',['--loadinit','--initcorrfile='+str(exact),'--extrapolate','--numImagTimeSteps=8','--noselfcons'])
    with h5py.File(refined) as f: assert f['preparation/Re_correlation/1-1'].shape==(9,9)
    (root/'terminate_987').touch()
    manual=run(root,'manual',['--runID=987','--absiterror=1e-15'],code=2)
    with h5py.File(manual) as f: assert text(f['preparation'].attrs['termination'])=='by hand'
    # Explicit prescribed/no-feedback run is one sampling iteration per stage.
    fixed=run(root,'fixed',['--noselfcons'])
    with h5py.File(fixed) as f: assert f['runtimedata'].attrs['num_Iterations']==1
print('native options/storage, exact dynamics, pCN, initialization and convergence checks passed')
