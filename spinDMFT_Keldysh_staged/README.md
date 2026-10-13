# Staged thermal-contour spinDMFT

This is an isolated source copy of `spinDMFT_Keldysh`, with a new two-stage
self-consistency driver. It retains the Matsubara–real preparation correlations.

1. Solve the imaginary-time spinDMFT problem and its initial magnetization.
2. Freeze the converged Matsubara pseudo-covariance and
   the initial mean/magnetization.
3. Iterate the mixed and real-time spin correlations with whole-contour
   Takagi fields. Independent draws are the default; `--jointMatsubaraPCN`
   enables Matsubara-only pCN using the same bare joint Gaussian law.
   The physical Matsubara pseudo-covariance stays frozen; its auxiliary
   Hermitian sampling covariance can change between real-time iterations.

Read [the LaTeX loop note](staged_loop.tex) for the equations, joint sampling,
convergence rules, output schema, and limitations. A [Markdown version](staged_loop.md)
is also provided.

## Build

From the repository root:

```sh
cmake -S spinDMFT_Keldysh_staged/Algorithm \
      -B spinDMFT_Keldysh_staged/Algorithm/build -DCMAKE_BUILD_TYPE=Release
cmake --build spinDMFT_Keldysh_staged/Algorithm/build -j4
ctest --test-dir spinDMFT_Keldysh_staged/Algorithm/build --output-on-failure
```

`Algorithm/CMakeLists.txt_` is also supplied as the build template. The executable
is `spinDMFT_Keldysh_staged/executable_DOUBLE.out`.

## Run

Run from this solver directory so its `Data/` receives the output:

```sh
cd spinDMFT_Keldysh_staged
mpirun -n 1 ./executable_DOUBLE.out \
  --spinmodel=ISO --cstype=C --Bname=z --Babs=0.5 --beta=1 \
  --numImagTimeSteps=20 --numRealTimeSteps=40 --Tmax=2 \
  --numSamplesPerCore=10000 --preparationSamplesPerCore=20000 \
  --numBlocks=32 --iterlimit=30 --preparationIterlimit=30 \
  --mixingAlpha=0.5 --seed=17
```

These are example settings, not validated production convergence parameters.
Defaults select independent sampling, dense whole-contour Takagi factorization, prefix
insertions, and partition-function normalization. The executable also supports Matsubara-only pCN with `--jointMatsubaraPCN`, as described below. It rejects FFT, weighted-dense,
closed-contour insertion/normalization, and a nonnegative FFT cutoff. The inherited machinery remains available to the copied unit tests.

`--preparationSamplesPerCore=0` uses `--numSamplesPerCore`; similarly,
`--preparationIterlimit=0` uses `--iterlimit`.
Preparation uses the `spinDMFT_cplx_field` stopping rule: the largest component
mean absolute change in Re/Im imaginary-time correlations or real origin
magnetization must be strictly below
`reliterror * S*(S+1) / (3*sqrt(3*N))`, where `N` is the total preparation sample
count across MPI ranks. Changes are measured before mixing. Sampling errors are
reported but do not relax this threshold. The existing diagnostic gates still
apply. The real-time stage requires its standardized residual threshold and
an absolute update smaller than `max(S,S*S)`. This additional finite-spin
scale check prevents exploding error bars from turning a divergent update
into apparent convergence; it is a divergence guard, not a precision guarantee.

## Joint Takagi Matsubara pCN

Enable the new sampler with:

```sh
--jointMatsubaraPCN --mhStepSize=0.3 --mhBurnIn=1000
```

`--samplingStrategy=pcn` is an alias. Explicitly selecting independent
sampling together with the flag is rejected. Omit the flag for independent
sampling; `--jointMatsubaraPCN=false` also selects the default when no pCN
strategy is supplied. `--conditionalSampling` and `--matsubaraGaugeScale`
have been removed and are rejected as unknown options. Historical HDF5 data
and diagnostic records remain unchanged.

At each real-time iteration, the code constructs the same weighted joint
Takagi factor L as the independent sampler. A thin real SVD of
`[Re L_M; Im L_M]` identifies an orthonormal basis U for the latent directions
that affect Matsubara fields. The sampler uses

```text
V_M = L_M U a
x   = U a + (I - U U^T) xi
V   = L x
```

pCN updates only a, targeting `p0(a) |Z_M(L_M U a)|`. The real normal vector
xi is freshly drawn after EVERY production step, including rejection.
Its projection does not change the Matsubara field, while preserving the
mixed coupling in the shared contribution. This retains both `L L^T`
and `L L^dagger` of the bare joint Gaussian, up to numerical rank tolerance.
The projector is applied without storing a full null-space basis.

Acceptance evaluates only the imaginary-time trace, using
`min(1, |Z_M(new)| / |Z_M(old)|)`. Real-time propagation never enters the
acceptance step. Negative-real and complex Z_M values remain in support;
zero modulus proposals are rejected and nonfinite traces abort.
All observable numerators and the denominator are divided by |Z_M|:
`sum[N/|Z_M|] / sum[Z_M/|Z_M|]`.

Stage 1 uses the imaginary-only dense Takagi factor with the same pCN
acceptance and estimator. Stage 2 rebuilds the row space and initializes and
burns in a new chain at every bath iteration, independently on each MPI rank.
No full Matsubara probability distribution is frozen across iterations.
`numSamplesPerCore` counts production steps, not accepted proposals.
Preparation has its own sample budget and uses the same pCN controls.

Errors use contiguous multi-scale batch means and paired complex ratio
influences, including denominator fluctuations. Inspect blocking curves and
ESS; preparation stopping still uses the nominal sample count. Real-time
errors are conditional on the frozen preparation. pCN importance-samples
only the thermal weight, so long-time real observable variance can still grow.
No general stability or drift improvement is claimed.

New filenames begin `staged__jointMpcn_absZ__`, distinguishing this method
from historical samplers. `/parameters` records `joint_matsubara_pcn`, the
updated coordinates, null-space refresh, and SVD rank rule. Both stages
store acceptance, phase, autocorrelation and blocking diagnostics. The
saved preparation factor is an imaginary-only reference, not a factor reused
for real-time sampling.

Regression commands:

```sh
python3 -B Algorithm/Tests/test_staged_pcn_solver.py ./executable_DOUBLE.out
python3 -B Algorithm/Tests/test_staged_pcn_solver.py ./executable_DOUBLE.out --mpi-ranks=2
```

Tests cover both Gaussian covariances (including weighted factors), weak and
singular Matsubara blocks, zero bath, mixed coupling with zero MM block,
actual pCN rejections with fresh real noise, exact thermal spins, harmonic
bath comparison with independent sampling, and a short self-consistent run
at `mixingAlpha=1`. These are implementation checks, not a long-time
production stability scan.

Validation on 2026-10-12: the Release build, all 18 registered CTests, and
the two-rank MPI pCN executable/HDF5 suite passed. The updated LaTeX loop
note also compiled successfully.

## Staggered Strang option

Add `--realTimeEvolution=staggered-strang` to use the point-sampled causal
scheme from Solution II of
`../spinDMFT_Keldysh/Notes/magnetization_drift_causal_discretization.tex`.
The existing `--realTimeEvolution=cfet` remains the default. Both independent sampling and `--jointMatsubaraPCN` support this evolution;
the dense/prefix/partition-function requirements still apply.

```sh
mpirun -n 1 ./executable_DOUBLE.out \
  --realTimeEvolution=staggered-strang --realTimeSubsteps=2 \
  --spinmodel=ISO --cstype=C --Bname=z --Babs=0.5 --beta=1 \
  --numImagTimeSteps=20 --numRealTimeSteps=40 --Tmax=2 \
  --numSamplesPerCore=10000 --preparationSamplesPerCore=20000 \
  --numBlocks=32 --iterlimit=30 --preparationIterlimit=30
```

The Matsubara-only preparation is unchanged. The real-time stage uses one
joint draw `[V_M, eta_midpoint, nu_boundary]`, retaining the frozen MM block
and all midpoint eta--M thermal correlations. The eta covariance is
`(Cgreater+Clesser)/2`; the eta--nu response is `Cgreater-Clesser` only when
the response boundary precedes its midpoint. Future response contractions,
nu--nu and nu--M are exactly zero in the target covariance. Smooth
correlation kernels are interpolated before applying this causal support;
the sampled stochastic fields are used directly.

Factorization uses the branch-equivalent basis
`[V_M, sqrt(2)*eta, nu/sqrt(2)]`, then scales each draw back to `[V_M,eta,nu]`.
This preserves the original branch noise objective because
`|V_+|^2+|V_-|^2 = 2|eta|^2+|nu|^2/2`. Using unweighted Takagi directly on
`[V_M,eta,nu]` preserves the unconjugated covariance but changes the
conjugated covariance and can strongly increase complex-trajectory variance.
Joint Matsubara pCN splits the already inverse-scaled physical factor,
so its bare Gaussian distribution retains this metric.

Each microstep of width `h` uses

```text
Q_n = exp(-i h nu_n dot S / 4)
P_n = exp(-i h H_eta(t_(n+1/2)))
U_n = Q_(n+1) P_n Q_n
B_n = Q_n P_n^(-1) Q_(n+1)
```

These maps evolve both the density and thermal spin insertions. Adjacent
steps share the same response boundary and half kick. The precession inverse
is an actual complex inverse. The deterministic magnetic field, mean field
and extra interaction enter the precession Hamiltonian.

`--realTimeSubsteps=q` must be positive. It divides each original measurement
interval into `q` Strang steps with `h=delta_real_t/q`, refining both the
evolution and the jointly sampled field grid. Measurements stay on the
original edge grid. With `N=numRealTimeSteps*q`, the full Gaussian dimension
is `3*(numImagTimeSteps+1+2*N+1)`, so increasing q increases factorization
cost. At least three imaginary-time intervals are needed for the unchanged
Gauss-CF4 preparation; one real-time interval is allowed.

Output filenames include `__prop=staggered-strang` and, for q greater than
one, `__substeps=q`. HDF5 records the evolution, field layout, covariance and
maps. The frozen MM check continues to run each real-time iteration.

Tests cover exact causal support, joint reconstruction,
the branch-equivalent conjugated covariance, shared kicks, general spin,
second-order convergence for smooth prescribed fields, Gaussian ensemble
trace closure, exact zero-bath observables, substeps, frozen-MM feedback,
failed preparation and unchanged preparation output. Second-order weak
convergence and reduced equilibrium magnetization drift still require
matched self-consistent scans.

Output filenames begin with `staged__`. `/results` contains the final raw spin
measurements; `/preparation` records the frozen input, last raw preparation
measurements, fixed physical covariance, reference imaginary factor,
preparation statistics, and invariant checks.
A failed preparation writes an explicitly marked file in
`Data/.../NOT_CONVERGED/` and does not launch the real-time stage.

Exit codes: 0 = both stages converged (or prescribed bath completed),
2 = preparation iteration limit, 3 = real-time iteration limit,
1 = invalid input or numerical failure (MPI abort).

For the executable/HDF5 regression suite, use a Python installation with numpy
and h5py:

```sh
python3 Algorithm/Tests/test_staged_solver.py ./executable_DOUBLE.out
python3 Algorithm/Tests/test_staged_solver.py ./executable_DOUBLE.out --evolution staggered-strang
```

The `--bath=harmonic` path runs prescribed-bath preparation and real-time sampling
without spinDMFT feedback. It provides a controlled test of the staged sampling
construction, not a converged spinDMFT solution.
