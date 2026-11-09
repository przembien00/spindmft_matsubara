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

The loop variables own the observable estimates. At the start of each sampling
sweep, `runtime.begin_iteration_accumulation(correlations, magnetization)` zeros
those variables and binds them as the accumulation destinations. Each sample
adds its unnormalized numerator directly to them.
`runtime.mpi_reduce_and_normalize()` reduces those same sums across MPI ranks,
computes statistical errors, and normalizes them by the complex partition sum.
Preparation uses the analogous preparation variables and passes `true` to the
begin call. Magnetization is one `Magnetization` object with `.Re` and `.Im`
members. Runtime retains partition sums, statistical blocks, errors, and the two
loops' diagnostic histories.

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
  --seed=17
```

These are example settings, not validated production convergence parameters.
Defaults select independent sampling and dense whole-contour Takagi factorization.
Prefix insertions and partition-function normalization are fixed; the
`spinInsertionStrategy` and `correlationNormalization` options have been removed. The executable also supports Matsubara-only pCN with `--jointMatsubaraPCN`, as described below. It always uses dense factorization. The `gaussianFactorization`, `gaussianWeightM`, `gaussianWeightEta`, and `gaussianWeightKappa` options have been removed.

`--dt` and `--partitionImagTolerance` have been removed because they did not
affect this solver. Imaginary-time spacing is `beta/numImagTimeSteps`;
real-time measurement spacing is `Tmax/numRealTimeSteps`.

Covariance construction always produces a dense matrix: the Matsubara matrix
for preparation and the full joint matrix for real-time sampling.

## Full two-time correlations

The default `--timeTranslationInvariant=true` reconstructs the real-time bath
from the endpoints of the mixed primitive at the lag `t1-t2`. For a driven
calculation, select

```sh
--timeTranslationInvariant=false --cstype=D
```

This measures both `Cgreater_ab(t1,t2)=<S_a(t1) S_b(t2)>` and
`Clesser_ab(t1,t2)=<S_b(t2) S_a(t1)>` on the full real-time measurement grid.
The self-consistent eta--eta noise and causal eta--nu response use these
two-time tensors, with interpolation in each time argument before imposing
causal support. Connected correlations subtract `m_a(t1)*m_b(t2)`.
Both tensors participate in MPI reduction, complex partition-sum normalization,
statistical errors and convergence checks, for independent sampling and pCN.
The measured magnetization remains time dependent; combining this mode with
`--constantMagnetization` is rejected. Spin-component symmetry is still
controlled separately by `--cstype`; use `D` for a generic drive that breaks
the reduced spin symmetries.

The initial thermal preparation remains an equilibrium Matsubara calculation.
Only `X(0,tau)` and the origin magnetization are frozen in feedback; real-time
pairs involving the origin remain live. The first real-time guess is seeded
from the stationary preparation. The prescribed harmonic bath retains its
analytic stationary covariance, while spin observables are measured in full.

Filenames in this mode include `two_time__`. The mixed datasets retain their
existing names and axes. Additional `/results/Re_correlation_greater`,
`Im_correlation_greater`, `Re_correlation_lesser` and `Im_correlation_lesser`
have axes `[t1,direction_pair(a,b),t2]`; both time axes use
`/results`' `real_times` attribute. `/runtimedata` contains matching
`*_sample_stds` datasets and, for pCN, `*_tau_int` datasets. The parameter
attribute `time_translation_invariant` records `yes` or `no`.
Observable storage, statistical-block storage and pair measurement now grow
quadratically with the number of real-time points; the Gaussian field grid
and its dense factorization dimension are unchanged.

This option supplies the two-time measurement and feedback needed for drives.
Pulse files automatically select it, as described below. Continuous
time-dependent external-field input is not yet provided.
Regression tests include prescribed pulse maps and nonunitary branch maps,
exact thermal-spin two-time correlations, both samplers, and frozen-MM feedback:

```sh
python3 -B Algorithm/Tests/test_two_time_solver.py ./executable_DOUBLE.out
python3 -B Algorithm/Tests/test_two_time_solver.py ./executable_DOUBLE.out --mpi-ranks=2
```

## Pulse sequences

Supply `--pulseSequence=sequence.hdf5`. The file contains three numeric datasets:

| Dataset | Shape | Meaning |
|---|---|---|
| `/pulses/times` | `[N]` | Pulse times in the same units as `Tmax` |
| `/pulses/axes` | `[N,3]` | Rotation-axis vectors `(x,y,z)`; normalized on input |
| `/pulses/angles` | `[N]` | Rotation angles in radians |

Times must be finite, nondecreasing, and in `[0,Tmax]`. Coincident off-grid
or origin pulses retain row order; coincident measurement-grid pulses add
as Hamiltonian field vectors. Axis vectors must be
finite and nonzero, and angles must be finite. An empty sequence is allowed.
Nonempty sequences automatically select `timeTranslationInvariant=false`;
explicit `true`, `--constantMagnetization`, and `--cstype` values other than
`D` are rejected.

With `hbar=1`, the positive angle convention is `+i angle axis.S`.
For a pulse at a nonzero measurement time, the last propagation substep
ending at that time uses

```text
h = delta_real_t / realTimeSubsteps
H_eff = H - (angle/h) axis.S
P = exp(-i h H_eff)
```

The pulse field is added to the magnetic, mean and sampled precession fields
in the same exponential. It is included on both forward and inverse-precession
branches, with the original response kicks retained. Multiple pulses at the
same measurement time contribute the sum of their `-(angle/h) axis` vectors.
They do not create separate rotation matrices. For noncommuting background
and pulse terms this is a finite-width field pulse over `h`; its finite-step
result differs from a separate instantaneous rotation. Increasing
`realTimeSubsteps` narrows that width while retaining its integrated angle.

Off-grid pulses and pulses at `t=0` use the instantaneous forward rotation
`R=exp(+i angle axis.S)` and backward rotation `R^dagger`.
The spin matrices already contain their spin normalization: for spin 1/2,
`R=cos(angle/2) I + i sin(angle/2) axis.sigma`. No additional factor of two
is applied to the angle. This positive-exponent convention corresponds to
the negative of the angle in the usual active `exp(-i angle axis.S)` convention.

A pulse off the measurement grid and strictly within a Strang microstep splits that step's precession
into the elapsed and remaining durations, retaining its sampled midpoint
field and original response-boundary kicks. A pulse exactly on an internal
substep edge that is not a measurement time acts after the complete preceding
step, including its response kick. Decimal times that differ from an edge
only by floating-point roundoff are snapped to it. Multiple pulses within
one step keep their individual times and order. `--realTimeSubsteps` refines
the propagation grid without changing the listed physical pulse times.

A pulse at `t=0` follows thermal preparation and precedes every real-time
measurement. The real-time origin is therefore `0+`; the imaginary branch
and its mean remain anchored at the pre-pulse state `0-`. Mixed correlations
retain their imaginary-branch spin in that pre-pulse preparation, while both
arguments of real-time greater/lesser correlations refer to the driven state.
Consequently, mixed endpoints need not equal real-time pairs with `t2=0`
when a zero-time pulse is present. A pulse at `Tmax` contributes to the last
substep's Hamiltonian and affects the final measurement.
The Matsubara pCN acceptance weight is unaffected by pulses.

The supplied helper creates a sequence (and refuses to overwrite an existing file):

```sh
python3 write_pulse_sequence.py echo.hdf5 \
  --pulse 0.0 y 1.5707963267948966 \
  --pulse 0.5 x 3.141592653589793
mpirun -n 1 ./executable_DOUBLE.out --pulseSequence=echo.hdf5 --cstype=D \
  --numImagTimeSteps=20 --numRealTimeSteps=40 --Tmax=1 \
  --numSamplesPerCore=10000 --preparationSamplesPerCore=20000
```

The axis argument can also be a vector such as `1,1,0`. To construct files
directly in Python, write the three datasets with h5py; `axes` must have
shape `(N,3)` even for an empty sequence. Output filenames include `pulsed__`.
The result file records the source path, rotation/measurement convention,
and a `/pulses` group containing the normalized, effective sequence; the
result file itself can be passed back as `--pulseSequence` input.

As elsewhere in the solver, bath correlations and the deterministic mean are
interpolated from the measurement grid. Resolve pulse neighborhoods with an
appropriate measurement grid when checking self-consistent discretization errors.
Driven-spin tests compare against the combined Hamiltonian exponentials at
measurement times and timed rotations elsewhere, including noncommuting axes,
spin 1/2 and spin 1, both samplers, origin/final pulses, and substeps. Separate
nonunitary-map and finite-bath tests check branch ordering and frozen preparation:

```sh
python3 -B Algorithm/Tests/test_pulse_solver.py ./executable_DOUBLE.out
python3 -B Algorithm/Tests/test_pulse_solver.py ./executable_DOUBLE.out --mpi-ranks=2
```

`--preparationSamplesPerCore=0` uses `--numSamplesPerCore`; similarly,
`--preparationIterlimit=0` uses `--iterlimit`.
Preparation uses the `spinDMFT_cplx_field` stopping rule: the largest component
mean absolute change in Re/Im imaginary-time correlations or real origin
magnetization must be strictly below
`reliterror * S*(S+1) / (3*sqrt(3*N))`, where `N` is the total preparation sample
count across MPI ranks. Changes are measured between successive estimates. Sampling errors are
reported but do not relax this threshold. The existing diagnostic gates still
apply. `--reliterror` is also the shared standard-error multiplier for requiring
the measured imaginary origin magnetization to be consistent with zero; there
is no separate `imagMagSigma` option. The real-time stage requires its standardized residual threshold and
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
preparation factor is not retained or saved; it can be reconstructed from the
preparation correlations and magnetization. Real-time sampling uses a fresh
full-contour factor.

Regression commands:

```sh
python3 -B Algorithm/Tests/test_staged_pcn_solver.py ./executable_DOUBLE.out
python3 -B Algorithm/Tests/test_staged_pcn_solver.py ./executable_DOUBLE.out --mpi-ranks=2
```

Tests cover both Gaussian covariances (including weighted factors), weak and
singular Matsubara blocks, zero bath, mixed coupling with zero MM block,
actual pCN rejections with fresh real noise, exact thermal spins, harmonic
bath comparison with independent sampling, and a short self-consistent run
without iteration mixing. These are implementation checks, not a long-time
production stability scan.

Validation on 2026-10-12: the Release build, all 18 registered CTests, and
the two-rank MPI pCN executable/HDF5 suite passed. The updated LaTeX loop
note also compiled successfully.

## Staggered Strang evolution

Staggered Strang is the only real-time evolution method and is the default.
The optional explicit spelling is `--realTimeEvolution=staggered-strang`.
It uses the point-sampled causal scheme from Solution II of
`../spinDMFT_Keldysh/Notes/magnetization_drift_causal_discretization.tex`.
Other evolution names, including `cfet`, are rejected. Both independent
sampling and `--jointMatsubaraPCN` support Strang splitting;
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

Output filenames include `__substeps=q` for q greater than one. HDF5 records
the evolution, field layout, covariance and maps. The frozen MM check continues to run each real-time iteration.

Tests cover exact causal support, joint reconstruction,
the branch-equivalent conjugated covariance, shared kicks, general spin,
second-order convergence for smooth prescribed fields, Gaussian ensemble
trace closure, exact zero-bath observables, substeps, frozen-MM feedback,
failed preparation and unchanged preparation output. Second-order weak
convergence and reduced equilibrium magnetization drift still require
matched self-consistent scans.

Output filenames begin with `staged__`. `/results` contains the final raw spin
measurements; `/preparation` records the frozen input, preparation statistics,
and invariant checks. Factor and covariance matrices are not saved. Each stage stores correlations and real/imaginary magnetization separately.
Each iteration replaces the previous estimates directly; `mixingAlpha` is no longer
supported. Preparation `*_sample_stds` describe the saved final estimates.
Preparation correlations, their standard errors, and pCN autocorrelation diagnostics
have axes `[direction, tau]`; real-time output retains `[t, direction, tau]`.
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
