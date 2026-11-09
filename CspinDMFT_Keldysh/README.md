# CspinDMFT_Keldysh

Thermal cluster spinDMFT in real time, using the correlation-replica embedding
of `CspinDMFT`/`CspinDMFT_finite_T` and the equilibrium staged Keldysh algorithm
of `spinDMFT_Keldysh_staged`.

## Build

```sh
cd CspinDMFT_Keldysh/Algorithm
cp CMakeLists.txt_ CMakeLists.txt
cmake -S . -B build
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

Requires C++17, MPI, HDF5, LAPACK/BLAS, Boost program_options and the repository's
Blaze headers. The executable is `CspinDMFT_Keldysh/executable_DOUBLE.out`.
The Gaussian and joint-Matsubara row-space implementations are compiled from
`spinDMFT_Keldysh_staged`; fixes to those implementations are shared.

## Two-spin AFM cluster in a cubic lattice

From `CspinDMFT_Keldysh/`:

```sh
mpirun -n 4 ./executable_DOUBLE.out \
  --project=CubicLattice \
  --config=Cubic_3D_N=2_NN_J=0.5 \
  --spinspinmodel=ISO --beta=1 \
  --numImagTimeSteps=32 --numRealTimeSteps=50 --Tmax=1 \
  --realTimeSubsteps=2 \
  --numSamplesPerCore=5000 --numBlocks=20 \
  --preparationIterlimit=50 --iterlimit=50 --reliterror=5 \
  --fileext=AFM
```

The local Hamiltonian uses `+J S_i.S_j`, so positive `J=0.5` is AFM.
The configuration is read from
`Configuration_Data/CubicLattice/Cubic_3D_N=2_NN_J=0.5.hdf5`.
A direct configuration path is also accepted. Configuration files are read
without modification. `--project` selects the source configuration subfolder
and the default data destination; `--dstproject` overrides the destination.

For a converged run the example saves

```text
Data/CubicLattice/spinmodel=ISO__config=Cubic_3D_N=2_NN_J=0.5__beta=1__substeps=2_AFM.hdf5
```

Runs stopped by the iteration limit or a termination file go under
`Data/<destination>/NOT_CONVERGED/`. Failed preparation adds
`__preparation_failed` and prevents real-time execution. Existing filenames
are preserved: `X` is appended until a free filename is found.

## Code structure

| Directory | Responsibility |
|---|---|
| `Algorithm/Parameter_Space/` | grouped Boost options, physical parameters, path/name conventions |
| `Algorithm/Mean_Field_Models/` | configuration reader, cluster Hamiltonian, first/second-moment feedback |
| `Algorithm/Contour/` | site-resolved contour types and causal staggered covariance kernels |
| `Algorithm/Functions/` | matrix exponentials, initialization/import, propagation, pCN chains and frozen feedback |
| `Algorithm/Run_Time_Data/` | paired MPI ratio statistics, diagnostics and convergence histories |
| `Algorithm/Storage_Concept/` | standard HDF5 groups, per-site-pair tensors and unique filenames |
| `Algorithm/Time_Measure/` | preparation/real-time iteration durations |
| `Algorithm/Tests/` | exact-physics, options, storage, sampling and MPI regressions |

`main.cpp` contains explicit preparation and real-time `while` loops. Each
loop constructs the mean/covariance, sets up the Gaussian sampler and RNG,
runs pCN burn-in when selected, samples and measures trajectories, reduces
the paired MPI ratios and records convergence. Each self-consistency iteration
uses one `for` loop over its fixed per-core sample count.
There is no `iteration` wrapper or accumulator object. Functions add samples
directly to the next correlation and magnetization containers. Separate functions
record paired statistical blocks in runtime data and reduce/normalize the results. Parameters are fixed after
construction; runtime counters and histories live in `RunTimeData`, separately
for preparation and real time.

## Observable containers

Correlations and magnetizations are independent objects throughout the solver.
There is no combined `State` object.

| Stage | Correlations | Magnetization |
|---|---|---|
| Preparation | `MatsubaraCorrelations`: `[i,j,direction,tau]` | `ClusterMagnetization`: `[site,direction]` |
| Real time | `ContourCorrelations`: `[t,i,j,direction,tau]` | `TimeMagnetization`: `[t,site,direction]` |

The preparation classes have no real-time dimension, including no singleton
time slice. Their public accessors take only spin/site and imaginary-time
indices. The real-time classes add the explicit time index. Each stage uses a
second correlation and magnetization variable for the next iteration, then swaps those buffers after computing the residual. Statistical
errors, closure measurements, sampling diagnostics and their histories live
in `runtime.preparation` and `runtime.real_time`. Storage receives the final
observable variables and reads their statistics from runtime data.

Each class allocates only its symmetry-allowed direction components, using
the shared spinDMFT correlation and magnetization symmetry definitions. A
stores `xx`; B stores `xx,zz`; C stores `xx,xy,yx,zz`; D stores all nine.
Equivalent diagonal components expand on access; forbidden components read
as exact zeros and have no storage. A has no magnetization components; B and
C store only z; D stores x,y,z. Ordered site pairs remain distinct.

The field covariance still uses all spin directions after the local/bath
coupling rotation, so an independently chosen anisotropic bath model is
retained. HDF5 direction order and output shapes are unchanged.

## Options

Run `./executable_DOUBLE.out --help` or `--helpNum` for the full list.
Cluster-specific names follow `CspinDMFT`; contour-specific names follow
`spinDMFT_Keldysh_staged`.

| Area | Options |
|---|---|
| Cluster embedding | `config`, `mfmodel=CorrelationReplica`, `project`, `dstproject` |
| Local/bath models | `spinspinmodel`, `lambda`, `rho`, `spinmfmodel=auto`, `smflambda`, `smfrho` |
| Coupling scales | `rescss`, `rescmf` |
| Fields and local interactions | `unifield`, `stagfield`, `chemshift`, `extraint`, `intstrength` |
| Static noise | `Cname=none/z`, `C` |
| Temperature and symmetry | `beta`, `cstype=A/B/C/D`; cluster alias `stype` |
| Contour grid | `numImagTimeSteps`, `numRealTimeSteps`, `Tmax`, `realTimeSubsteps` |
| Sampling | `numSamplesPerCore`, `preparationSamplesPerCore`, `numSamplesPerSet`, `numBlocks`, `seed` |
| Joint pCN | `samplingStrategy=pcn`, `jointMatsubaraPCN`, `mhStepSize`, `mhBurnIn` |
| Iteration | `reliterror`, `iterlimit`, `preparationIterlimit`, `itermode`, `absiterror` |
| Contour diagnostics | `covarianceTolerance`, `takagiTolerance`, `minimumPhaseMagnitude` |
| Feedback constraints | `constantMagnetization`, `noselfcons` |
| Initialization | `initmode`, `initdcorr`, `initndcorr`, `corrperiods`, `loadinit`, `initcorrfile`, `impcorrfile`, `extrapolate` |
| Prescribed harmonic bath | `bath=harmonic`, `bathOmega`, `bathCoupling`, `bathComponent` |
| Naming and control | `info`, `fileext`, `numPrintDigits`, `runID` |

`spinmodel` is an alias for `spinspinmodel`. Supported coupling matrices are
the shared Physics models `ISO`, `Ising`, `XXZ`, `DRF`, `XY`, and `het`.
The bath model defaults to the local model. Fields enter as `+h_i S_i^z`.
`unifield=h` applies the same z field to every site. `stagfield=h` assigns
`+h` to site 1 and alternates signs using BFS two-coloring of the configuration's
spin-spin coupling graph (edges have `abs(J_ij) > 1e-12`). The graph must be
connected and bipartite; the assignment uses the original couplings before
`rescss` scaling. `chemshift` accepts one comma-separated z field per site.
Choose only one of `unifield`, `stagfield`, or non-`none` `chemshift`, including
when the field amplitude is zero. With none specified, all fields are zero.
The former `Bname`, `Babs`, `Bphi`, and `Btheta` options have been removed.
Generated filenames identify `unifield` or `stagfield`; HDF5 parameters save
that amplitude and the resolved per-site fields in `chemical_shift`.
`extraint=quadrupolar` adds `3*intstrength*sum_i (S_i^z)^2` to the local
Hamiltonian. `intstrength` has energy units and defaults to zero; it is unused
when `extraint=none`. For spin 1/2 this term is a constant energy shift and
does not affect normalized observables.

`rescss` scales local `J`; `rescmf` scales first-moment weights linearly and
covariance weights quadratically. Fields are independent of those coupling
scales. `Cname=z` adds independent static z noise to each site.

The equilibrium solver requires `timeTranslationInvariant=true` and
`realTimeEvolution=staggered-strang`. Pulses and full two-time feedback are
outside this solver's scope. Real covariance eigenvalue clipping options
from the classical cluster solver do not apply to complex Keldysh
pseudo-covariances; these use Takagi factorization. Cluster propagation uses
general complex matrix exponentials, with the staged solver's Padé(13)
scaling-and-squaring scheme.

The immutable run parameters contain the preparation and real-time contour
layouts, constructed once from the configuration's cluster size and time grids.
Both sampling and pCN reuse these layouts.

`self_consistent_equations` takes the cluster correlations and magnetizations
and returns the joint contour covariance, including embedding couplings and
static noise. The sampler factorizes that covariance to draw the fields.

`numBlocks` is a requested count: an equal-sized divisor of the sample count
is selected automatically, with at least two blocks per rank. For pCN the
count is a power of two and errors take the largest paired-ratio jackknife
error under within-chain block coarsening. Use even sample counts.
Independent blocks are pooled over MPI ranks. Sample counts are fixed by
`numSamplesPerCore` and `preparationSamplesPerCore`; there is no adaptive
retry loop. The `adaptive` and `staterrtolerance` options have been removed.

The default seed is `random`; a numeric seed is reproducible for fixed rank
count and options. The generated seed is broadcast to all ranks and saved
exactly as an integer. Create `terminate_<runID>` in the working directory
to stop after an iteration. Keep the terminal's generated run ID when using
this mechanism.

Terminal output follows `spinDMFT_Keldysh_staged`: aligned startup parameters,
Matsubara and mixed/real-time iteration headings, phase/closure diagnostics,
per-step timings, stage totals and the total duration. Each sampling loop uses
the shared duration estimator and progress bar; burn-in is excluded from the
production duration estimate. As in
`CspinDMFT_finite_T`, startup also lists representative site pairs (one-based),
configuration, cluster size, spin lengths, local couplings, coupling models,
site fields and rescaling factors. Clusters below five sites print the whole
local coupling matrix; larger clusters print `J(0,1)`. Output is printed once
on rank zero. The generated run ID is printed for manual termination.

Startup uses `delta_imag_t` for the imaginary-time grid spacing and omits the
fixed propagator name. `delta_real_t = Tmax/numRealTimeSteps` is the measurement
spacing. `delta_real_propagation_t = delta_real_t/realTimeSubsteps` is the
integration substep; refining it also refines the sampled real-time field grid.

The cluster solver reports `complex-weight N_eff` for both sampling modes;
this quantity does not measure pCN autocorrelation. It is not relabeled as
the staged solver's `minimum autocorrelation N_eff`.

## Convergence

The default `itermode=relative` uses the staged solver's rules:

- Preparation compares the largest absolute update to
  `reliterror * max_i[S_i(S_i+1)/3] / sqrt(3*N_preparation)`.
- Real time compares the largest real/imaginary update in current-standard-error
  units to `reliterror`. It also requires the complex absolute update to stay
  below `max_i(max(S_i,S_i²))`.
- Phase magnitude, covariance/factor reconstruction, and imaginary origin
  magnetization diagnostics must pass.

The prepared origin is excluded from real-time residuals. With
`constantMagnetization`, measured magnetizations are saved, but feedback
uses the prepared origin magnetization and excludes real-time magnetization
updates from the residual.

`--itermode=absolute --absiterror=0.005` selects an explicit absolute stopping
threshold. `--absiterror=conservative` uses the preparation-scale expression
above. These are iteration stopping criteria, not bounds on the physical
error or time-discretization error.

Exit codes follow the staged workflow: `0` for normal termination, `2` for
failed/interrupted preparation, `3` for a real-time iteration limit, and
a nonzero process status for unhandled input, numerical or I/O failures.

## Algorithm and statistics

Internally the mixed primitive is

\[
X_{ij}^{ab}(t,\tau)=\langle S_j^b(-i\tau) S_i^a(t)\rangle.
\]

Its `tau=beta` endpoint is the greater real-time function
`<S_i^a(t) S_j^b(0)>`; `tau=0` is the lesser function.
Only site pairs listed in the configuration's `correlation_categories` are
measured, stored and used for convergence. Their order is preserved in memory
and in the category metadata. Symmetry-equivalent source pairs are already
folded into these representatives by the configuration weights. The field
covariance still covers every cluster site pair. Configurations without a
category list default to all ordered pairs.

The mean and connected covariance feedback use configuration weights `A` and
`W`:

\[
\bar V_i^a(t)=\sum_{kc} A_{ik}D_{ac}m_k^c(t),
\]
\[
K_{ij}^{ab}(t,\tau)=\sum_{(k,l)\in\mathcal C}\sum_{cd}s_{kl}W_{ij,kl}D_{ac}D_{bd}
\left[X_{kl}^{cd}(t,\tau)-m_k^c(t)m_l^d(0)\right].
\]

Products are unconjugated. Weights are stored per source-pair orientation:
`s_kl=2` for an off-diagonal representative whose reverse is absent from
`correlation_categories`, and `s_kl=1` for diagonal pairs or when both ordered
pairs are measured explicitly. No manual factor of two is needed. A symmetric
mirror entry for an unlisted reverse pair is allowed and is not counted again.
Other lattice-symmetry multiplicities remain folded into representative weights.
See [the configuration convention](../CspinDMFT_finite_T/Configuration/README.md)
for migration of external legacy files. The joint field
basis is `[V_M(tau_k), eta(t_(n+1/2)), nu(t_n)]`, with `3*N` components per
node. Causal support is enforced at the actual staggered nodes; smooth
kernel values use up to four interpolation knots. Sampling weights `1,2,1/2`
preserve the branch-equivalent norm. Distinct Matsubara endpoints are kept.

Imaginary propagation uses two-Gauss-node CF4 with local cubic field
interpolation. Each real microstep uses

\[
R_n=e^{-ih\nu_n\cdot S/4},\quad
P_n=e^{-ih(H_{\rm cl}+(\bar V+\eta_n)\cdot S)},\qquad
U_n=R_{n+1}P_nR_n,\quad B_n=R_nP_n^{-1}R_{n+1}.
\]

`h=Tmax/(numRealTimeSteps*realTimeSubsteps)`. Adjacent microsteps share
response fields. The complex inverse uses `exp(+ihH)`. The mean is linearly
interpolated to microstep midpoints. Measurements stay on coarse edges.

Independent sampling accumulates raw numerators and divides by complex
`sum Z_M`. The joint pCN sampler updates only the current joint Takagi law's
Matsubara row-space coordinates, targets `p0*abs(Z_M)`, and draws fresh
null-space real coordinates each production step. Both numerator and
partition contributions are divided by `abs(Z_M)` before MPI reduction.
Burn-in is repeated after each field update. No adjoint or absolute partition
normalization replaces the complex denominator.

`bath=harmonic` supplies independent identical oscillator baths at every
site, coupled through `bathComponent`, and performs one measurement iteration
per stage. `noselfcons` performs one iteration per stage without a fixed-point
convergence claim. Its real stage is seeded from measured preparation as
usual; it is intended for diagnostics rather than a converged thermal bath.

## Files and HDF5 layout

The standard top-level groups are `parameters`, `runtimedata`, `results`,
`preparation`, and `timedata`. As in the existing solvers, scalars, lists and
strings are attributes; correlation tensors are datasets.

| Entry | Shape/convention |
|---|---|
| `results/Re_correlation/i-j`, `Im_correlation/i-j` | `[real_time,direction,tau]`, one-based site labels |
| `results/Re_magnetization`, `Im_magnetization` | `[real_time,site,xyz]` |
| `results/correlation_direction_labels` | selected directions in standard tensor order |
| `runtimedata/Re_correlation_sample_stds/i-j`, `Im_correlation_sample_stds/i-j` | standard errors, same shape as correlations |
| `runtimedata/Re_magnetization_sample_stds`, `Im_magnetization_sample_stds` | `[real_time,site,xyz]` |
| `preparation/Re_correlation/i-j`, `Im_correlation/i-j` | `[direction,tau]` |
| `preparation/Re_magnetization`, `Im_magnetization` | `[1,site,xyz]` |
| `timedata` attributes | seconds per iteration, total duration and start time |

Standard direction order is `xx,xy,xz,yx,yy,yz,zx,zy,zz` for D, `xx,xy,yx,zz`
for C, `xx,zz` for B, and `xx` for A. The containers enforce that symmetry during storage and access. Choose a symmetry compatible
with the Hamiltonian; D is the default.

`parameters` stores the configuration, its measured site-pair categories,
physical options and numerical settings. The redundant
`configuration_correlation_categories` dataset is no longer written.
Parameter names have no `(original)` prefix. `filename` and
`configuration_src_file` contain only the respective file names; projects
are stored separately. `correlation_symmetry_type` is the string `A`, `B`,
`C` or `D`. The redundant `algorithm`, `normalization_strategy` and
`real_time_evolution` attributes are omitted.
`runtimedata` stores iteration histories, termination reason, actual block
sizes, errors and closure diagnostics. `preparation` stores the frozen
state and its own convergence histories. Final
measured real-time observables are saved unchanged by origin freezing and
constant-magnetization feedback.

Coordinates are reconstructed from the grid parameters; they are not repeated
in `results`. Preparation omits `num_samples_per_core_history` and
`real_time_started`; its fixed sample count is stored as `samples_per_core`.

For example:

```python
import numpy as np

with h5py.File(filename) as f:
    p = f['parameters'].attrs
    t = np.arange(int(p['num_RealTimePoints'])) * p['delta_real_t']
    tau = np.arange(int(p['num_ImagTimePoints'])) * p['delta_imag_t']
    xx_greater_12 = f['results/Re_correlation/1-2'][:, 0, -1]
```

`corrperiods` sets the number of oscillations of generated `cos`/`sin`
initial guesses on the imaginary-time interval `[0,beta]`: the period is
`beta/corrperiods`. It defaults to one and affects only `initdcorr=cos` and
`initndcorr=sin`. Other initializers and imported data do not use it.

`--loadinit --initcorrfile=Data/.../file.hdf5` imports the saved preparation,
including the configured representative site pairs and magnetizations. Both initialization
generators handle import mode directly, so `main` initializes each observable
through the same function for generated and imported data. `--extrapolate` permits
linear interpolation onto a different imaginary grid within the saved beta
interval. It seeds a new preparation run; it does not skip thermal
self-consistency or resume a Markov chain. With no path, the auto-generated
destination filename is used. `--initmode=import --impcorrfile=...` is the
cluster spelling of the same operation.

Generated initial correlations use the shared Physics functions:
`imagtime`, `const`, `exponential`, `Gaussian`, or `cos`; `zero` is also
available. `initndcorr=zero/sin` controls off-diagonal guesses. The optional
`initdcorr=isolated --initndcorr=isolated` uses exact isolated-cluster Gibbs
correlations instead.

## Migration from the initial implementation

| Previous option | Current convention |
|---|---|
| `--config=path.hdf5` | still accepted; normally use `--project=... --config=basename` |
| `--output=...` | automatic `Data/<dstproject-or-project>/<physics>.hdf5`; use `fileext` |
| `--numTimeSteps` | `--numImagTimeSteps` |
| `--dt` | `--Tmax` and `--numRealTimeSteps` |
| `--h_z=0.4` | `--unifield=0.4` or site-specific `chemshift` |
| `--blocks` | `--numBlocks` |
| `--preparationMaxIterations` | `--preparationIterlimit` |
| `--realTimeMaxIterations` | `--iterlimit` |
| `--tolerance=0.005` | `--itermode=absolute --absiterror=0.005`; default is relative convergence |
| `--rescale_meanfield`, `--rescale_spinspin` | `--rescmf`, `--rescss` |

The old per-iteration-file schema is not written by this version. Existing
files and copied configurations are preserved.

## Validation and limits

The tests cover exact interacting two-spin thermal/mixed dynamics, arbitrary
complex cluster fields, non-Hermitian inverses, response endpoint factors,
mixed spin lengths, all one-site kernel entries against the staged solver,
site exchange and connected feedback. End-to-end tests cover local/bath model
options, per-pair output, filename collisions, pCN, imports/interpolation,
reduced symmetry, harmonic baths, fixed sample counts and termination.

The Python test needs NumPy and h5py; select an interpreter with
`-DPython3_EXECUTABLE=/path/to/python3` if needed. To run the same tests using
MPI from this folder:

```sh
python3 Algorithm/Tests/test_solver.py ./executable_DOUBLE.out mpirun -n 2
```

The solver assumes equilibrium bath time-translation invariance. Dense field
dimension is `3*N*(N_tau+1+2*N_t*q+1)`; matrix storage and factorization grow
quadratically/cubically. The interacting cluster Hilbert space grows as
`product_i(2*S_i+1)`, capped at 4096 states with substantially lower practical
limits. Standard errors describe the current Monte Carlo iteration; they do
not propagate feedback uncertainty. pCN blocking diagnostics require enough
samples to resolve autocorrelation. Finite-coupling continuum accuracy and
long-time stability remain to be established by refinement studies.
