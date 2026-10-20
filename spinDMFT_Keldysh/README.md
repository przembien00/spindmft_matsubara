# spinDMFT_Keldysh

This solver implements the three-branch finite-temperature real-time design in
`Notes/forward_backward_contour_implementation.tex`. It was derived from
`spinDMFT_real_time`; that source tree is unchanged.

## Contour and stochastic field

The contour convention is

```text
0 --(+)-> T --(-)-> 0 --(M)-> -i beta
```

Each sample contains the joint vector `[V_M(tau_k), eta_midpoint, nu_boundary]`.
The imaginary field retains distinct one-sided variables at `tau=0+` and
`tau=beta-`. The real-time fields are sampled on the staggered grid described
below, jointly with both imaginary-time endpoints.

The covariance is the complex symmetric, unconjugated pseudo-covariance
`Gamma = E[V V^T]`. `Contour/Contour_Kernel.cpp` reconstructs all `MM`, mixed,
`++`, `--`, `+-`, and `-+` blocks from the edge-grid mixed primitive before the
upper triangle is mirrored. It is not treated as Hermitian PSD and no
negative-eigenvalue clipping is performed.

Coincident points on a time-ordered or anti-time-ordered branch use
`theta(0)=1/2`, hence `C(t,t)=(G greater(t,t)+G lesser(t,t))/2`. This is fixed
by the contour definition rather than exposed as a selectable prescription.

With `eta=(V_+ + V_-)/2` and `nu=V_+ - V_-`, the branch algebra gives a
symmetrized `eta-eta` block, a causal `eta-nu` response block, and a zero
`nu-nu` pseudo-covariance. A zero `E[nu nu]` does not make sampled `nu` vanish.

### Staggered Strang evolution

Staggered Strang is the only real-time evolution method and is the default.
The optional explicit spelling is `--realTimeEvolution=staggered-strang`. It
implements Solution II of `Notes/magnetization_drift_causal_discretization.tex`.
Other evolution names, including `cfet`, are rejected.
For example, from this directory:

```bash
mpirun -n 1 ./executable_DOUBLE.out \
  --realTimeEvolution=staggered-strang --realTimeSubsteps=1 \
  --numImagTimeSteps=40 --numRealTimeSteps=100 --beta=2 --Tmax=5 \
  --numSamplesPerCore=2500 --cstype=D
```

The solver samples one joint vector `[V_M, eta_midpoint, nu_boundary]`, with
`N = numRealTimeSteps * realTimeSubsteps` real intervals and dimension
`3 * (numImagTimeSteps + 1 + 2*N + 1)`. Its covariance retains the Matsubara
block and midpoint mixed thermal correlations, uses `(Cgreater+Clesser)/2`
for the eta block, and uses `Cgreater-Clesser` for eta--nu only when the
response boundary precedes the midpoint. Future response contractions,
nu--nu, and nu--M are exactly zero before factorization. Smooth primitive
correlations are interpolated with up to four neighboring measurement knots
(linear/quadratic for shorter real grids), then the causal support is applied.
Sampled fields are used directly without interpolation.

For each microstep of width `h`, the maps are

```text
Q_n = exp(-i h nu_n dot S / 4)
P_n = exp(-i h H_eta(t_(n+1/2)))
U_n = Q_(n+1) P_n Q_n
B_n = Q_n P_n^(-1) Q_(n+1)
```

`P_n^(-1)` is the actual complex inverse, with the deterministic magnetic,
mean-field, and extra-interaction terms included only in `H_eta`. Adjacent
microsteps share the same response field and half kick at their boundary.
The existing Matsubara Gauss-CF4 preparation, density and spin-insertion
machinery are retained. Measurements and output correlations stay on the
original edge grid. Both prefix and closed-contour insertions and both
existing correlation normalizations are supported.

Staggered Strang defaults to `--gaussianFactorization=dense`; explicit
`weighted-dense` is also supported with independent sampling. Its existing
M, eta, kappa weights apply directly to M, midpoint eta and boundary
`kappa=nu/2`, so the effective nu weight is `gaussianWeightKappa/4`.
Explicit FFT factorization, FFT cutoff options, and `--realTimeSubsteps=0`
are rejected. The obsolete `--contactCounterterm` option has been removed. At least three
imaginary-time intervals are needed for Gauss-CF4 preparation; one or more
real-time intervals are allowed.

Increasing `realTimeSubsteps` here refines the jointly sampled stochastic
grid and increases factorization cost. Output filenames
include `__prop=staggered-strang` and, when greater than one, `__substeps=q`;
HDF5 metadata records the layout, covariance and maps.

Tests cover the full joint covariance, complex mean subtraction and coupling
rotation, dense and weighted reconstruction, shared half kicks, general spin,
second-order prescribed-field convergence, joint Gaussian ensemble closure,
inserted observables, pCN caching, CLI validation and solver output. The
second-order statement concerns smooth prescribed fields. Weak convergence
and equilibrium magnetization drift of the complete self-consistent
thermal-contour problem still require matched numerical scans.

The standard solver supports dense joint sampling (default) and
`weighted-dense` with independent sampling. FFT interpolation is incompatible
with the causal staggered grid and is rejected.

## Weighted dense noise optimization

For example, append the following options to a run:

```sh
--gaussianFactorization=weighted-dense --samplingStrategy=independent \
--gaussianWeightM=1 --gaussianWeightEta=2 --gaussianWeightKappa=8
```

These are the default weights for `weighted-dense`; the default factorization
is `dense`. All weights must be finite and strictly positive. Explicit
weight options on another factorization are rejected. The weights apply per
native grid point and spin component, not as time-integration weights.

The sampler applies diagonal weights to the directly sampled joint vector
`[V_M, eta_midpoint, nu_boundary]`, factors the weighted complex symmetric
covariance, and scales the draw back to physical fields. The nu weight is
`gaussianWeightKappa/4` because `kappa=nu/2`. The complete unconjugated
covariance, including mixed thermal correlations, is reconstructed after
inverse scaling. The weighted sampler requires independent sampling.

Output names receive `__noiseW=M,eta,kappa`. HDF5 records the attributes
`gaussian_weight_M`, `gaussian_weight_eta`, `gaussian_weight_kappa` on `parameters`,
the weight basis/objective, and the physical reconstruction basis. The existing
per-iteration reconstruction-error and numerical-rank attributes remain in use.

Validation from the repository root:

```sh
cmake --build spinDMFT_Keldysh/Algorithm/build -j 4
ctest --test-dir spinDMFT_Keldysh/Algorithm/build --output-on-failure
/Users/przembien/Projects/python_venv/bin/python \
  spinDMFT_Keldysh/Algorithm/Tests/test_staggered_strang_solver.py \
  spinDMFT_Keldysh/executable_DOUBLE.out
```

The executable/HDF5 test requires NumPy, h5py and working MPI; it isolates all
outputs in temporary directories. C++ tests independently check the weighted
nuclear-norm optimum, full covariance, canonical-weight distribution, batching,
and nonlinear propagation against analytic Gaussian averages using deterministic
quadrature. These checks establish implementation correctness, not a measured
long-time variance improvement.

## Propagation and estimators

Real-time evolution uses the staggered Strang maps above. Imaginary-time
preparation uses two Gauss nodes, cubic interpolation of the Matsubara edge
field and a two-exponential Gauss-CF4 step. It requires at least three
imaginary-time intervals. `realTimeSubsteps=q` must be positive and affects
only the real-time stochastic grid and propagation.

For spin `1/2`, the matrix exponentials use the exact complex Pauli identity;
higher spins and non-scalar local interactions use a general matrix exponential.

Correlations and magnetization use the closed-contour spin insertion

```text
S_closed(t) = B_-(T,0) U_+(t,T) S U_+(t,0).
```

Select this with `--spinInsertionStrategy=closed-contour`; the default is `prefix`.

Here `U_+(t,T)` denotes the forward continuation from `t` to `T`, namely
`U_N ... U_{t+1}`, not an inverse propagator. Every forward and backward step
appears exactly once at every measurement time, including `t=0`. At `t=T`
the continuation is the identity. Equal branch fields recover the prefix
insertion `B_-(t,0) S U_+(t,0)` by cancellation; distinct fields do not.

The measurement builds `B_-(T,0)` once and a right-to-left array of combined
suffixes `B_-(T,0) U_+(t,T)`, then sweeps the forward prefix `U_+(t,0)`.
This requires linear work and linear matrix storage in the real-time grid.
Only the spin directions needed for correlations or magnetization are formed.
They are contracted with fixed imaginary-time insertions (or `rho_M` for
magnetization) through `Tr(A B)=sum_ij A_ij B_ji` without a temporary matrix
product. No inverse or adjoint relation is assumed.

For the mixed correlation currently accumulated by the implementation, the
two equal-time spin insertions are not adjacent on the closed contour. At
`t=0`, `tau=0`, its numerator is
`Tr[rho_M S_b B_-(T,0) U_+(T,0) S_a]`. Even for `a=b` and spin `1/2`, this
does not reduce to `Tr[rho_M B_-(T,0) U_+(T,0)]/4`, because the complete
real-time contour product lies between the two spin operators. Consequently,
closed-contour normalization alone does not enforce `G^{aa}(0,0)=1/4`.

Separately, `R(t+dt)=U_step R(t) B_step` computes
`D(t)=Tr[U_+(t,0) rho_M B_-(t,0)]`. By default this remains only the prefix
closure diagnostic. Select

```text
--correlationNormalization=closed-contour
```

to normalize correlations and magnetization at every real time by the single
fixed endpoint value `D(T)` instead of `Z_M`. The denominator therefore always
contains the complete product `B_-(T,0) U_+(T,0)` and is independent of the
measurement time. The default is `partition-function`. New HDF5 files record
the normalization in `parameters/correlation_normalization` and
`parameters/magnetization_normalization`, and the insertion choice in
`parameters/spin_insertion_strategy`.
Closed-contour-normalized filenames receive `__corrnorm=D`.

The authoritative correlation is represented on the imaginary edge grid:

```text
X_edge^{ab}(t,tau_k)       = <S_b(-i tau_k) S_a(t)>
```

The edge endpoints independently measure greater and lesser functions, and the
edge values supply the mixed field covariances. Independent samples accumulate
raw complex `N`, `Z_M`, the closed-contour trace `D(t)`, and `m_a(t)`. The
default physical bare-prior estimator is

```text
(sum N) / (sum Z_M).
```

There is no trajectorywise division and no absolute-`Z` normalization. With
`--correlationNormalization=closed-contour`, the correlation and magnetization
estimators are instead `(sum N(t))/(sum D(T))` and
`(sum M(t))/(sum D(T))`. They remain ratios of ensemble sums rather than
averages of trajectorywise ratios. In pCN mode this setting also changes the
importance target to `p0(r) Re D_r(T)`. The estimator is formed as
`sum[A/Re D(T)]/sum[D(T)/Re D(T)]`, with `A=N` or `M`. With
partition-function normalization, pCN instead retains the
`p0(r) Re Z_M(r)` target and forms the self-normalized complex ratio
`sum[N/Re Z_M]/sum[Z_M/Re Z_M]`. Keeping the reweighted denominator is
important at finite chain length: it preserves exact samplewise identities
such as `g^{aa}(0)=1/4` instead of leaving them with residual phase noise.
Uncertainty is computed with a delete-one-block jackknife of paired complex
numerator and denominator sums across MPI ranks. `sum |Z_M|` is retained only
for phase and effective-sample-size diagnostics. Iteration totals are packed
into one `MPI_Allreduce`; each rank then evaluates its local delete-one-block
replicates and a second packed `MPI_Allreduce` combines their centered moments.
No rank gathers or replicates the block-resolved correlation tensors.

## Execution and statistical accounting

Independent sampling processes up to 32 trajectories per batch, including a
short final batch when needed. Dense, weighted-dense, and FFT samplers use
BLAS matrix-matrix products for independent block-factor draws. FFT factors
store alternating real and imaginary rows of the same complex Takagi factor:
GEMM applies them to independent batches, and GEMV applies them to single pCN
proposals. Identical factors share this packed storage, while their latent
coordinates remain independent. The original complex factor storage is released
after packing. FFTs then evaluate each sample on its native and Gauss-node grids.
Random coordinates are generated in sample-major order, so batch boundaries do
not change the random sequence. Multiplication order can change floating-point
roundoff; covariance, cutoff, and latent rank are unchanged. Exactly
`numSamplesPerCore` observations are accumulated on every iteration. There is
no adaptive sample-count or convergence change.

Observable measurement uses contiguous cached numerators, fixed-size spin-1/2
matrix contractions, and reusable trajectory and measurement workspaces.
Higher spins retain the general matrix path. Prefix insertion avoids building
closed-contour suffixes, and closed-contour insertion avoids maintaining an
unused backward prefix. Both preserve the independently propagated branches
and all imaginary- and real-time measurement points.

pCN refreshes its cached observables only after acceptance or before the first
production measurement. Every rejected state is still accumulated with its
original weight, sample square, and position in the contiguous statistical
blocks. With partition-function normalization, proposals first propagate the
Matsubara branch; real-time propagation and spin insertions are completed only
after acceptance. Closed-contour normalization still evaluates the complete
trajectory to obtain the proposal likelihood `Re D(T)`.

Covariance construction skips exactly zero rotation coefficients and diagnoses
the raw transpose residual while filling its canonical triangle. FFT setup
uses a read-only covariance source instead of allocating a dense physical-grid
matrix. Existing cutoff, reconstruction, and raw-kernel diagnostics remain in
the output. Gauss-node Fourier phases are precomputed for each FFT grid.

`--antitheticPairs` has been removed. New files contain ordinary trajectory or
pCN-state counting metadata. Existing historical output files remain readable
by analysis scripts. The scan entry points are `slurm_beta_sample_scan.sh` and
`slurm_beta_discretization_scan.sh`.

Symmetry splitting changes the latent basis. Identical seeds therefore need
not reproduce fields from older binaries, even though the factorized
pseudo-covariance and Hermitian covariance are preserved to the numerical rank
tolerance. This implementation is recorded in the HDF5 execution metadata.

Local timings, numerical validation, and instructions for reproducing the
fixed-input benchmark are in [the performance report](Analysis/performance/README.md).

## Prescribed harmonic-bath validation

`--bath=harmonic` replaces the spinDMFT field closure by the exact thermal
contour correlation of one oscillator

```text
H_B = omega_0 a^dagger a,       X = g (a + a^dagger).
```

Set `omega_0`, `g`, and the coupled spin component with `--bathOmega`,
`--bathCoupling`, and `--bathComponent=x|y|z`. Because this is a genuinely
single-component covariance, the mode requires `--cstype=D`. On the stored
mixed grid it uses

```text
X(t,tau) = g^2 (n_B+1) [exp(-omega_0 tau + i omega_0 t)
                       +exp(-omega_0 (beta-tau) - i omega_0 t)].
```

Thus `tau=0` is the analytic lesser function, `tau=beta` is the greater
function, and the real branch blocks contain a nonzero greater-minus-lesser
response. The full `MM`, mixed, `++`, `--`, `+-`, and `-+` covariance is still
built by the ordinary contour kernel and sampled by the selected Gaussian
factorization.

Harmonic-bath mode performs one fixed-bath Monte Carlo measurement and stores
the raw measured estimators without mixing them with the spinDMFT initial
guess. Measured correlations and magnetization, `JQ`, `JL`, and static noise do
not enter the field distribution. The external magnetic field is not rescaled
by `JQ`; it and local on-site interactions remain active. Output is written
under `Data/HARMONIC_BATH/` with the bath parameters in the filename.

For example:

```bash
mpirun -n 1 ./executable_DOUBLE.out \
  --bath=harmonic --bathOmega=1 --bathCoupling=0.2 --bathComponent=x \
  --cstype=D --beta=2 --Bname=z --Babs=0.5 \
  --numImagTimeSteps=50 --Tmax=10 --numRealTimeSteps=200 \
  --numSamplesPerCore=40000 --numBlocks=40
```

## Self-consistency

The field closure is

```text
mean_M = JL D Re m(0)
mean_+(t) = mean_-(t) = JL D Re m(t)
X_conn^{ab}(t,tau) = X^{ab}(t,tau) - m_a(t) m_b(0)
Gamma = JQ^2 D X_conn,contour D^T + N_static.
```

The product is complex and unconjugated. The connected mixed primitive is
formed before the `MM`, mixed, `++`, `--`, `+-`, and `-+` blocks are
reconstructed. The imaginary-branch one-point function is represented by its
equilibrium value at the contour origin. Thus a finite-sample real-time drift
does not re-enter the covariance as a spurious disconnected contribution; in
the stationary limit this reduces exactly to subtraction of `m m^T`.

The real-time mean is evaluated on the stored real-time grid and enters the
same Strang midpoint mean as the fluctuating forward and backward
fields. The same mean trajectory is used on both real branches. The Matsubara preparation
uses the physical real part of `m(0)`; imaginary magnetization remains a
diagnostic and does not enter the Hamiltonian. The external field is common to
all three branches.
The optional `--mixingAlpha` consistently mixes the full complex
magnetization trajectory and the edge-grid correlation primitive, then the
connected primitive is recomputed. The reported fixed-point residual is formed
from the raw, pre-mixing update and includes the full complex magnetization
trajectory as well as the real and imaginary correlation parts.
Convergence requires the pointwise raw residual to satisfy
`|F(x)-x| < q*s`, where `q=--reliterror` and `s` is the current paired-ratio
jackknife standard error. Real and imaginary correlation components are tested
separately, and the real and imaginary magnetization components are likewise
tested against their respective errors. A nonzero residual with zero error, or any
nonfinite residual statistic, cannot converge.

`--constantMagnetization` applies an equilibrium projection to the
self-consistency state after every Monte-Carlo update:

```text
m_a(t) -> m_a(0)
mean_M = mean_+(t) = mean_-(t) = JL D Re m(0)
X_conn^{ab}(t,tau) = X^{ab}(t,tau) - m_a(0)m_b(0).
```

The raw measured `m_a(t)` is still stored. Projected runs receive the filename
suffix `__mag=constant`.

Convergence also checks the configured covariance, branch identity, Takagi,
phase, closed-contour denominator, and imaginary-magnetization tolerances.
KMS endpoint differences and real-time magnetization stationarity are
postprocessing observables rather than convergence gates.

`--loadinit` deliberately imports only the `t=0` Matsubara seed from a file in
this solver's `Data/` directory, including its real and imaginary parts. The
full real-time primitive is then initialized with the documented Gaussian
envelope. This is the explicitly supported Matsubara-seed mode; it is not a
full real-time checkpoint restart.

## Build and run

```bash
cd spinDMFT_Keldysh/Algorithm
cp CMakeLists.txt_ CMakeLists.txt
cmake -S . -B build_keldysh -DCMAKE_BUILD_TYPE=Release -DUSE_EIGEN=OFF
cmake --build build_keldysh -j
ctest --test-dir build_keldysh --output-on-failure
cd ..
mpirun -n 1 ./executable_DOUBLE.out \
  --beta=1 --numImagTimeSteps=20 \
  --Tmax=2 --numRealTimeSteps=40 \
  --numSamplesPerCore=1000 --numBlocks=20
```

Dense joint sampling is the default. `--gaussianFactorization=weighted-dense`
is available with independent sampling. The selected factorization is recorded
in HDF5; existing files are protected by filename collision handling. Dense
factorization scales cubically in the complete staggered field dimension.

## HDF5 output

`results` contains:

- `Re/Im_correlation`, shaped `[real_time, direction_pair, tau_edge]`;
- `Re/Im_magnetization`, shaped `[real_time, spin_component]`;
- explicit real-time and imaginary-edge grid attributes.

`runtimedata` contains matching complex-ratio jackknife errors plus paired
derived jackknife errors for the closed-contour ratio `D(t)-1`,
`gaussian_factor_reconstruction_errors`,
`gaussian_covariance_approximation_errors`,
`gaussian_factor_latent_dimensions`,
`gaussian_largest_factorization_dimensions`, branch identities, phase and effective
sample size, complex partition sums, denominator constancy, and the raw complex
fixed-point residual.
The `standardized_iteration_errors` history stores the largest pointwise
`|F(x)-x|/s`; `parameters/iteration_error_sigma_threshold` stores `q`.

The removed final diagnostics are recovered directly from `results` by forming

```text
X(t,tau) = Re_correlation + i Im_correlation
KMS(t,a,b) = conj(X(t,a,b,beta)) - X(t,a,b,0)
m(t) = Re_magnetization + i Im_magnetization
stationarity(t,a) = m(t,a) - m(0,a).
```

Their values are exact functions of the stored estimates. The stored marginal
errors of `X` and `m` do not retain the covariance needed to reconstruct the
old paired-jackknife error of either difference.

## Validation supplied

```sh
ctest --test-dir Algorithm/build --output-on-failure
python3 Algorithm/Tests/test_staggered_strang_solver.py ./executable_DOUBLE.out
```

The C++ tests cover joint covariance reconstruction, Strang convergence for
smooth prescribed fields, shared response kicks, general spin, contour
measurements and pCN caching. The Python suite checks the default method,
substeps, weighted sampling, pCN, harmonic baths, HDF5 metadata and rejection
of removed evolution options. It requires NumPy and h5py.
