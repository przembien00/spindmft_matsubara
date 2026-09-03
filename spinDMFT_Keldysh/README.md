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

- `--gaussianFactorization=dense` (default) applies one Autonne--Takagi
  factorization to the complete joint-contour covariance;
- `--gaussianFactorization=svd` constructs the same canonical Takagi ensemble
  from a full complex SVD of the joint covariance;
- `--gaussianFactorization=fft` is the frequency-truncated sampler. It
  mirrors the real-time covariance into a doubled block-circulant grid,
  Fourier-transforms the first `N_tau` Matsubara values and the doubled-real
  axis with FFTW while retaining the distinct `beta-` field as an untransformed
  three-component boundary block. It retains all Matsubara modes and the
  boundary block together with real modes satisfying
  `|omega| <= fftCrossFrequencyCutoff` in one dense complex-SVD Takagi block.
  Matsubara coupling to higher real frequencies is set to zero and each
  remaining high-frequency `{omega,-omega}` real block is factorized
  independently. Covariance construction writes directly into these retained
  blocks: the real-real blocks come from one FFT of the `6x6` circulant lag
  kernel, and omitted mixed-frequency blocks are never allocated. Each block
  is also sampled independently into persistent FFT buffers; neither a global
  doubled covariance nor a global dense factor is assembled.
  The complete frequency field is inverse-transformed and the physical
  real-time interval is retained. The default cutoff is `3`; a
  negative cutoff restores the untruncated full-frequency dense test.

For `Gamma=U Sigma V^dagger`, the SVD algorithm corrects the singular-vector
phases (and any degenerate singular-value subspaces) to obtain `L` with

```text
L L^T      = Gamma,
L L^dagger = U Sigma U^dagger = sqrt(Gamma Gamma^dagger).
```

Consequently `dense` and `svd` have identical pseudo-covariance, Hermitian
covariance, and complete real Gaussian distribution. Both draw one independent
real latent vector and perform one factor multiplication per trajectory. The
`fft` path preserves the Matsubara--Matsubara and real--real
pseudo-covariances, but deliberately removes the high-frequency tail of the
Matsubara--real pseudo-covariance. The discarded relative Frobenius norm and
largest dense block dimension are printed and stored each iteration. Its
doubled-grid canonical Hermitian covariance need not equal the canonical
Hermitian covariance of a direct physical-grid factorization, so it remains a
separately selectable comparison algorithm.

## Propagation and estimators

The imaginary branch uses edge fields and the three-exponential CFET4-opt
composition from `spinDMFT`, with insertions on the same edge grid. Writing
that endpoint composition as `C4(H_new,H_old;z)`, the real branches are
propagated independently:

```text
U_+(t+dt) = C4(H_+(t+dt),H_+(t);-i dt) U_+(t)
B_-(t+dt) = B_-(t) C4(H_-(t),H_-(t+dt);+i dt).
```

The reversed endpoint order on the backward branch reverses the CFET4
exponential composition, so equal forward/backward fields still close the
real contour algebraically.

Both exponentials support general complex non-Hermitian matrices. `B_-` is not
formed from `U_+^{-1}` or `U_+^dagger`.

For spin `1/2`, each CFET exponential is evaluated directly from its weighted
complex field with the Pauli identity
`exp(a I + w.sigma)=exp(a)[cosh(q) I + sinh(q)/q w.sigma]`, where
`q^2=w.w` is an unconjugated complex dot product.  This path also includes the
quadrupolar interaction, which is proportional to the identity for spin `1/2`.
Higher-spin or non-scalar local-interaction cases retain the general matrix
exponential.

Correlations and magnetization use the full-contour spin insertion

```text
S_full(t) = B_-(T,0) U_+(t,T) S U_+(t,0).
```

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

Separately, `R(t+dt)=U_step R(t) B_step` is retained only for the existing
prefix closure diagnostic `D(t)=Tr[U_+(t,0) rho_M B_-(t,0)]`. It is not used
as the observable normalization. New HDF5 files record the insertion rule in
`parameters/spin_insertion_strategy` (an attribute).

The authoritative correlation is represented on the imaginary edge grid:

```text
X_edge^{ab}(t,tau_k)       = <S_b(-i tau_k) S_a(t)>
```

The edge endpoints independently measure greater and lesser functions, and the
edge values supply the mixed field covariances. Independent samples accumulate
raw complex `N`, `Z_M`, the closed-contour trace `D(t)`, and `m_a(t)`. The
physical bare-prior estimator is

```text
(sum N) / (sum Z_M).
```

There is no trajectorywise division and no absolute-`Z` normalization.
Uncertainty is computed with a delete-one-block jackknife of paired complex
numerator and denominator sums across MPI ranks. `sum |Z_M|` is retained only
for phase and effective-sample-size diagnostics. Iteration totals are packed
into one `MPI_Allreduce`; each rank then evaluates its local delete-one-block
replicates and a second packed `MPI_Allreduce` combines their centered moments.
No rank gathers or replicates the block-resolved correlation tensors.

Both sampling strategies optionally support antithetic Gaussian pairs:

```text
--samplingStrategy=independent --antitheticPairs
```

For every independent real latent vector `r`, it evaluates the two joint
contour-field fluctuations `L r` and `-L r`. The deterministic mean field is
unchanged, and all Matsubara, forward, and backward fluctuations are negated
together. `numSamplesPerCore` remains the total number of evaluated
trajectories and must be even, so the number of independent latent draws is
half that value. Numerators and `Z_M` from both members remain in the global
complex ratio; no memberwise ratios are formed. Jackknife blocks are adjusted
downward when necessary so that a block never splits an `r,-r` pair. Output
filenames receive the suffix `__antithetic`.

With pCN, the same flag selects a sign-symmetrized Markov chain:

```text
--samplingStrategy=pcn --antitheticPairs
```

One pCN state consists of both `r` and `-r`. Writing
`z(r)=Re Z_M(r)`, the chain targets

```text
pi_pair(r) proportional to p0(r) [z(r) + z(-r)]
```

and accepts the usual pCN proposal with the likelihood ratio
`[z(r')+z(-r')]/[z(r)+z(-r)]`. The measured pair observable is

```text
[N(r) + N(-r)] / [z(r) + z(-r)].
```

This preserves the ordinary pCN expectation while integrating out the binary
sign of the latent state. Both real partition functions must be finite and
positive, consistently with the real-positive pCN target assumption.
`numSamplesPerCore` counts production pair states (not individual contour
trajectories), and pCN blocking lengths and autocorrelation times are expressed
in pair steps. Every proposed state costs two contour-trajectory evaluations;
rejected states reuse the already stored pair. Output filenames retain the
`__antithetic` suffix and store the pair-state counting convention in HDF5.

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

Add `--gaussianFactorization=svd` or `--gaussianFactorization=fft` to select an
alternative. Their files receive the corresponding `__factor=svd` or
`__factor=fft_wcut=<omega>` suffix, while the default dense filenames are unchanged, so
identical runs can be compared without overwriting one another. The dense
factorization scales cubically in the complete field dimension; `svd` uses a
complex matrix of the physical dimension instead of the `2N x 2N` real lift.
For `fft`, only the low-frequency joint block is dense; high real-frequency
pair blocks are factorized and sampled independently.

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

- complex and rank-deficient dense Autonne--Takagi reconstruction, canonical
  SVD reconstruction, dense/SVD Hermitian-covariance equivalence, and empirical
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
