# spinDMFT_Keldysh

This solver implements the three-branch finite-temperature real-time design in
`Notes/forward_backward_contour_implementation.tex`. It was derived from
`spinDMFT_real_time`; that source tree is unchanged.

## Contour and stochastic field

The contour convention is

```text
0 --(+)-> T --(-)-> 0 --(M)-> -i beta
```

Each bare-Gaussian sample contains jointly correlated physical field values

```text
[ V_M(tau_k), V_+(t_n), V_-(t_n) ], k=0,...,N_tau.
```

The imaginary field uses distinct one-sided variables at `tau=0+` and
`tau=beta-`; they are jointly sampled with their separate lesser and greater
mixed covariances and are never identified samplewise. The two real fields are
also distinct at the contour joins. The dense field dimension is

```text
3 * (numImagTimeSteps + 1 + 2 * (numRealTimeSteps + 1)).
```

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

Three sampling algorithms are available:

- `--gaussianFactorization=dense` applies Autonne--Takagi
  factorizations to exactly disconnected blocks of the joint-contour covariance;
- `--gaussianFactorization=weighted-dense` optimizes a weighted noise power
  in the full `(V_M,eta,kappa)` basis, with `kappa=(V_+-V_-)/2`, then returns
  physical branch fields with the original unconjugated covariance. It
  requires `--samplingStrategy=independent`; see the weight options below;
- `--gaussianFactorization=fft` (default) is the frequency-space sampler with optional truncation. It
  mirrors the real-time covariance into a doubled block-circulant grid,
  Fourier-transforms the first `N_tau` Matsubara values and the doubled-real
  axis with FFTW while retaining the distinct `beta-` field as an untransformed
  three-component boundary block. It retains all Matsubara modes and the
  boundary block together with all real modes by default. When a nonnegative
  cutoff is explicitly supplied, it retains real modes satisfying
  `|omega| <= fftCrossFrequencyCutoff` in one dense Takagi block.
  Matsubara coupling to higher real frequencies is set to zero and each
  remaining high-frequency `{omega,-omega}` real block is factorized
  independently. Covariance construction writes directly into these retained
  blocks directly from a read-only contour-covariance source: the real-real blocks come from one FFT of the `6x6` circulant lag
  kernel, and omitted mixed-frequency blocks are never allocated. Each block
  is also sampled independently into persistent FFT buffers. Exact spin-component
  blocks are split within each retained frequency block. No physical-grid dense
  covariance is assembled. Disabling the cutoff keeps the full joint frequency
  covariance within each connected component block.
  The complete frequency field is inverse-transformed and the physical
  real-time interval is retained. The default cutoff is `-1` (disabled).
  Supply, for example, `--fftCrossFrequencyCutoff=3` to enable truncation.

The dense sampler draws independent real latent coordinates for every block. Exactly identical
blocks reuse their factorization while retaining independent random coordinates.
The rank cutoff is evaluated at the original whole-matrix dimension and spectral
scale. No nonzero coupling is discarded to create a block. The
`fft` path preserves the Matsubara--Matsubara and real--real
pseudo-covariances. Only when a nonnegative cutoff is supplied does it remove
the high-frequency tail of the Matsubara--real pseudo-covariance. The discarded relative Frobenius norm and
largest dense block dimension are printed and stored each iteration. Its
doubled-grid canonical Hermitian covariance need not equal the canonical
Hermitian covariance of a direct physical-grid factorization, so it remains a
separately selectable comparison algorithm.

## Weighted dense noise optimization

For example, append the following options to a run:

```sh
--gaussianFactorization=weighted-dense --samplingStrategy=independent \
--gaussianWeightM=1 --gaussianWeightEta=2 --gaussianWeightKappa=8
```

These are the default weights for `weighted-dense`; the default factorization
is `fft`. All weights must be finite and strictly positive. Explicit
weight options on another factorization are rejected. The weights apply per
native grid point and spin component, not as time-integration weights.

With `z=T V=(V_M,eta,kappa)`, the sampler factors
`P=A (T Gamma T^T) A^T=L L^T` and draws `V=T^-1 A^-1 L r` for real standard
normal `r`. Here `A=sqrt(W)` and a harmless common scale sets the largest
weight to 2. Exactly disconnected covariance blocks and batched draws use
the existing dense machinery. Every Matsubara, mixed and real-time block is
retained. Reconstruction errors are measured after transforming back to the
original physical covariance. Excessive weight contrast that causes loss of
physical covariance through numerical rank truncation raises an error.

Weights `(1,2,2)` reproduce the canonical dense Gaussian distribution, though
individual trajectories need not match for the same seed. Weights `(1,2,8)`
match the relative real-time penalties in Schmitz--Stockburger, since their
`nu=2*kappa`. A smaller weighted noise power does not establish smaller
observable errors. See [the derivation](Notes/optimized_gaussian_sampling.tex)
and [the existing frozen-bath analysis](Analysis/noise_covariance/README.md).

The [native-grid sampling audit](Analysis/weighted_sampling_audit/README.md)
finds that `(1,2,8)` worsens late-time magnetization errors at `beta=2.5`,
`h_z=0`, `Tmax=15`: on the same symmetric bath, 64 independent repetitions
give about 51% larger reported endpoint errors than `(1,2,2)`. The physical
covariance and weighted objective pass independent checks. Treat these weights
as experimental; validate the variance of the observable of interest before
using them for production. Final-iteration error bars also do not include all
uncertainty propagated through preceding self-consistency iterations.

A subsequent [observable-variance weight scan](Analysis/weighted_sampling_audit/weight_selection.md)
finds `(1,2,1)` gives about 3% lower combined late-time magnetization variance
than `(1,2,2)` on two frozen zero-field baths at these same parameters. The
endpoint standard-error improvement is only about 1–2%, and closure variance
increases. This is a modest experimental alternative, not a validated change
to the default or a solution to systematic drift.

The new mode supports endpoint propagation (`q=0`) and CF4 (`q>=1`) using
cubic interpolation of the sampled edge fields, and either existing insertion
and normalization choice. It currently rejects pCN: positivity of its selected
real sampling weight has not been established for the changed Gaussian
ensemble. Propagation and ensemble-ratio normalization are unchanged.

Output names receive `__noiseW=M,eta,kappa`. HDF5 records the attributes
`gaussian_weight_M`, `gaussian_weight_eta`, `gaussian_weight_kappa` on `parameters`,
the weight basis/objective, and the physical reconstruction basis. The existing
per-iteration reconstruction-error and numerical-rank attributes remain in use.

Validation from the repository root:

```sh
cmake --build spinDMFT_Keldysh/Algorithm/build -j 4
ctest --test-dir spinDMFT_Keldysh/Algorithm/build --output-on-failure
/Users/przembien/Projects/python_venv/bin/python \
  spinDMFT_Keldysh/Algorithm/Tests/test_weighted_solver.py \
  spinDMFT_Keldysh/executable_DOUBLE.out
```

The executable/HDF5 test requires NumPy, h5py and working MPI; it isolates all
outputs in temporary directories. C++ tests independently check the weighted
nuclear-norm optimum, full covariance, canonical-weight distribution, batching,
and nonlinear propagation against analytic Gaussian averages using deterministic
quadrature. These checks establish implementation correctness, not a measured
long-time variance improvement.

## Propagation and estimators

By default, the imaginary branch uses edge fields and the three-exponential
endpoint CFET composition from `spinDMFT`, with insertions on the same edge
grid. Despite its historical `CFET4-opt` name, using only endpoint Hamiltonians
makes this path globally second order for a general time-dependent field.
Writing that endpoint composition as `C4(H_new,H_old;z)`, the real branches are
propagated independently:

```text
U_+(t+dt) = C4(H_+(t+dt),H_+(t);-i dt) U_+(t)
B_-(t+dt) = B_-(t) C4(H_-(t),H_-(t+dt);+i dt).
```

The reversed endpoint order on the backward branch reverses the CFET4
exponential composition, so equal forward/backward fields still close the
real contour algebraically.

Both branch propagators support general complex non-Hermitian matrices. `B_-` is not
formed from `U_+^{-1}` or `U_+^dagger`.

CF4 is the default propagator. `--realTimeSubsteps=q` selects the method:

- `q=0`: the original endpoint propagator, one full step per measurement
  interval, with no internal field nodes. This also selects the original
  endpoint treatment on the imaginary branch.
- `q=1` (default): one fourth-order CF4 step per measurement interval, using
  two Gauss--Legendre nodes and two exponentials.
- `q>1`: `q` CF4 substeps per measurement interval, with dense, weighted-dense,
  or FFT sampling.

CF4 needs at least three imaginary- and real-time intervals. For shorter grids,
use `--realTimeSubsteps=0`. Dense and weighted-dense sampling support arbitrary
nonnegative `q`; their CF4 nodes use local
four-point cubic interpolation of the sampled fields. There is no separate
CF4 option.

For example:

```text
--numRealTimeSteps=40 --Tmax=2 --realTimeSubsteps=4
```

This measures at 41 points spaced by `0.05` and propagates with steps of
`0.0125`. Each realization samples its frequency amplitudes once. Phase-shifted
inverse FFTs evaluate that same realization at both Gauss nodes per substep.
The signed-frequency convention assigns the Nyquist mode to the positive
frequency. No new random variables or larger covariance factorization are
introduced; this refines propagation of the finite Fourier field and does not
restore frequencies absent from the measurement grid. The deterministic mean
field and one-sided Matsubara field use the original four-point cubic stencil
for CF4. Increasing a positive `q` leaves imaginary-time preparation unchanged.

To use dense sampling with the same substeps, add `--gaussianFactorization=dense`
to the example above. Each substep evaluates the same four-point cubic
interpolant of the sampled native edge fields at its two Gauss nodes. The
covariance matrix, factorization, random draws, and measurement grid stay fixed
as `q` increases. This refines propagation of the interpolated field; it does
not improve the underlying field discretization. Weighted-dense uses the same
interpolation and continues to require independent sampling.

The backward CF4 step exchanges the early and late Gauss nodes and changes the
contour-step sign. Forward microsteps compose from the left and backward
microsteps from the right. Only their full measurement-interval products are
stored, so correlations, magnetization, and self-consistency retain their
original grids. Equal forward/backward fields continue to close algebraically.
Both independent sampling and pCN use the selected trajectories, including the
final trace when closed-contour normalization is selected.

HDF5 records `real_time_substeps`, `delta_real_propagation_t`, and
`real_time_field_interpolation`; `delta_real_t` continues to mean the measurement
spacing. For `q=0`, the propagation spacing equals the measurement spacing and
`propagator=endpoint-cfet4`. Positive `q` uses `propagator=gauss-cf4` and the
`__prop=cf4` filename suffix. For `q>1`, filenames also include `__substeps=q`.

For spin `1/2`, each CFET exponential is evaluated directly from its weighted
complex field with the Pauli identity
`exp(a I + w.sigma)=exp(a)[cosh(q) I + sinh(q)/q w.sigma]`, where
`q^2=w.w` is an unconjugated complex dot product.  This path also includes the
quadrupolar interaction, which is proportional to the identity for spin `1/2`.
Higher-spin or non-scalar local-interaction cases retain the general matrix
exponential.

Correlations and magnetization use the closed-contour spin insertion

```text
S_closed(t) = B_-(T,0) U_+(t,T) S U_+(t,0).
```

This is selected by the default option
`--spinInsertionStrategy=closed-contour`; the alternative remains `prefix`.

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
same CFET4 endpoint Hamiltonians as the fluctuating forward and backward
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

FFT sampling is the default, with no frequency truncation. Use
`--gaussianFactorization=dense` or `weighted-dense` to choose another
factorization. The selected factorization and cutoff are recorded in HDF5;
existing files are protected by the filename collision handling. The dense
factorization scales cubically in the complete field dimension.
With an explicit nonnegative FFT cutoff, only the low-frequency joint block
is dense and high real-frequency pair blocks are factorized independently.
Without truncation, the joint frequency blocks can be substantially larger.

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

CTest covers:

- complex and rank-deficient dense Autonne--Takagi reconstruction and empirical
  `E[V V^T]`, including a nonzero sampled response field with `E[nu nu]=0`;
- full and frequency-truncated blockwise FFT sampling, physical-grid
  round-trip pseudo-covariance, reduced dense-block dimension, and zero-rank
  block handling;
- the complete branch table, transpose identities, Keldysh transform, and
  checked flat layout;
- exact complex-ratio and paired jackknife algebra;
- general complex matrix exponentials, independent real branches, contour
  closure for equal fields, and the analytic `JQ=0` finite-field spin result:
  `Z=2 cosh(beta h_z/2)`, `m_z=-tanh(beta h_z/2)/2`, `G_zz=1/4`, and the complex
  transverse greater correlation.

The `real_time_substeps` CTest checks signed-frequency interpolation (including
Nyquist), unchanged native fields and latent rank, both branch multiplication
orders, equal-field closure, and propagation convergence for FFT nodes and dense
cubic interpolation. To check CLI defaults,
invalid inputs, pCN/independent execution, and HDF5 grids and metadata, run:

```bash
python3 Algorithm/Tests/test_substeps_solver.py ./executable_DOUBLE.out
```

FFT BLAS regression coverage includes single/batched field equivalence, zero-rank
blocks, empty and partial batches, alternating batch sizes, RNG order, and
Gauss nodes for `q=0,1,3`, with and without frequency truncation. A saved solver
executable can also be compared against the current one using:

```bash
python3 Algorithm/Tests/test_fft_blas_solver.py ./executable_DOUBLE.out /path/to/reference.out
```

This comparison fixes the input covariance by using one self-consistency
iteration. In later iterations, roundoff differences can rotate degenerate
factorization bases and therefore change individual same-seed trajectories;
bitwise equality across self-consistency iterations is not an invariant.
