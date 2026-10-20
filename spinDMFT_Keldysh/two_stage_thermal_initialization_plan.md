# Two-stage thermal preparation and real-time spinDMFT: design plan

Status: proposal only. No solver implementation or numerical validation is implied.

## Physical definition

Introduce an opt-in `--thermalPreparation=decoupled` mode; keep the existing
three-branch contour as the default reference. The new mode first converges the
imaginary-time spinDMFT fixed point, then freezes its equilibrium local state
and mean field while solving only real-line correlations on the forward and
backward branches. The real-real bath covariance remains self-consistent. In
the Gaussian contour field, this sets

    Gamma = [ Gamma_MM       0      ]
            [    0      Gamma_RR   ],

where `Gamma = E[delta V delta V^T]` is the **unconjugated**, complex symmetric
pseudo-covariance and `R=(+,-)`. The `++`, `+-`, `-+`, and `--` blocks remain;
the approximation discards only Matsubara-real blocks. It therefore retains
the real-time response field `nu=V_+-V_-` and its correlation with
`eta=(V_++V_-)/2`. Sampling one identical field on both real branches would
define a different approximation and lose that response channel.

For each imaginary-field trajectory let `M[V_M]` be the unnormalized
imaginary-time propagator and `Z_M=Tr M`. Prepare

    rho_eq = E_M[M] / E_M[Z_M],

with a ratio of ensemble means, **not** `E_M[M/Z_M]`. For independent
Matsubara and real fields, linearity makes replacing every trajectory's `M`
by `rho_eq` exactly equivalent within the block-diagonal approximation. The
real-time stage then starts from one fixed density matrix, with `Tr rho_eq=1`.

The equilibrium bath mean is `h_mf = J_L D m_eq`. Its connected real-time
noise primitive uses `G^{>/<}_{ab}(t)-m_eq,a m_eq,b`, the same coupling rotation
and static-noise term as the full solver. Measure the physical `m_a(t)` but
keep `m_eq` in the equilibrium bath closure. Drift of measured `m_a(t)` is an
approximation error, not a quantity to hide by resetting the output.

For every real-time sample, use the existing forward/backward propagation and
construct the two initial operator insertions `S_b rho_eq` and `rho_eq S_b`.
They replace the `tau=beta` and `tau=0` boundary insertions of the full
contour, respectively. Accumulate both greater and lesser `G_ab(t)`, their
complex block errors, measured `m_a(t)`, and `D(t)=Tr[U_+(t) rho_eq B_-(t)]`.
Use a denominator of one for the independent-sampling estimator because
`rho_eq` is normalized; retain the `D(t)` diagnostic. Do not apply an extra
samplewise trace normalization.

## Physics expectation

This is a factorized impurity/bath initialization at the switch to real time.
The local thermal state may be exact for the converged imaginary-time
approximation while its correlations with the initial bath are absent. It is
therefore plausible in weak-coupling/high-temperature regimes or after a
short initial transient, but it has no general guarantee of stationarity,
KMS/detailed balance, or the finite-field long-time plateau. In particular,
correct `rho_eq` and real-real bath covariance alone do not restore the
missing imaginary-real preparation correlations. Treat any speed gain as
conditional on these accuracy checks.

## Implementation sequence (for a later coding task)

1. **Separate data contracts.** Add explicit `MatsubaraCorrelation` and
   `RealCorrelation` (`greater`, `lesser`, magnetization) types. Do not store
   fake `tau` slices in the mixed `CorrelationSet` or write zero-filled mixed
   correlations. Keep the full-contour type and path intact.
2. **Imaginary fixed point.** Add an imaginary-only loop in
   `Algorithm/main.cpp` using the `MM` covariance and current imaginary
   propagator/insertion algebra in `Functions/Functions.cpp`. Preserve the
   Keldysh solver's complex symmetric `MM` covariance, one-sided edge
   contacts, and compatible complex ratio/block statistics. The separate
   `spinDMFT` imaginary-time solver feeds back only the real correlation and
   cannot be plugged in unchanged. Accumulate `sum M` and `sum Z_M` per MPI
   rank and block, reduce, and form `rho_eq` only after a converged imaginary fixed
   point. Check trace, Hermiticity, positivity within Monte Carlo error, and
   consistency of `Tr(rho_eq S_a)` with `m_eq`. For spin 1/2, `rho_eq` is
   determined by `m_eq`, but store the measured matrix as a cross-check;
   higher spins require the full matrix.
3. **Real-only covariance.** In `Contour/Contour_Kernel.*` or a new real-only
   kernel, build only the `+/-` block from greater and lesser values with the
   existing equal-time `theta(0)=1/2` prescription. Apply the current
   connected subtraction and `J_Q^2 D...D` rotation. Use a sampler whose
   dimension is `6*(N_R+1)`, with dense factorization first; adapt FFT and
   weighted-dense only after their real-only covariance reconstruction is
   checked. Preserve `U_+ rho_eq B_-`, never replace `B_-` by `U_+^{-1}`.
4. **Real-only measurement and iteration.** Make a real-only trajectory and
   accumulator so each sample costs O(N_R) measurements rather than O(N_R N_M)
   mixed insertions. Seed `G^{>/<}(0)` from the imaginary-stage contact
   values and initialize the rest of the real-time guess explicitly. Mix and
   converge the real greater/lesser functions; monitor measured magnetization
   without feeding its drift back into the frozen thermal mean. Reuse the
   current real-time propagators, including the chosen substep option.
5. **CLI and storage.** Parse `--thermalPreparation=full|decoupled` in
   `Parameter_Space`; default `full`. Initially support independent sampling
   and dense factorization in the decoupled mode. Record stage-specific
   sample counts, seeds, convergence residuals, `rho_eq`, and timing. Write
   `Re/Im_matsubara_correlation` and `Re/Im_greater`, `Re/Im_lesser` with
   separate standard errors; label the approximation and absent mixed block
   in HDF5 metadata. Keep legacy full-contour keys unchanged for full runs.
   Add a separate decoupled output filename tag to prevent accidental
   comparison or overwrite. Update plotting loaders only after the schema is
   fixed.

`--loadinit` currently imports only an imaginary-time seed from the local
Keldysh HDF5 layout; it does not import an equilibrated density matrix or a
whole real-time fixed point. An external Matsubara checkpoint import can be a
later extension after its normalization, model/field/beta/grid, symmetry,
spin, and statistical-error contracts are verified. It is not part of the
minimal two-stage mode.

## Validation and decision gates

1. **Algebra and zero-bath limits.** With `J_Q=0`, the mode must reproduce
   direct thermal spin evolution, `D(t)=1`, and spin-1/2 equal-time identities
   (for example `G_aa(0)=1/4`). Check nonzero field and both operator
   orderings. Verify dense `L L^T` reconstructs `Gamma_RR` and branch
   transpose identities before self-consistency.
2. **Fixed-bath ablation.** At the same prescribed thermal bath and grids,
   compare (a) the full three-branch calculation and (b) the block-diagonal
   approximation with the same equilibrated `rho_eq`. Compare
   `m_a(t)`, `D(t)`, greater/lesser correlations, and complex-valued errors.
   This isolates the dropped mixed block before nonlinear feedback can mask
   it. The harmonic-bath option is a natural controlled case.
3. **Equilibrium checks.** At zero and nonzero `h_z`, compare real-time
   `G^{>/<}(0)` with the Matsubara boundary, check `m_a(t)-m_eq` and
   `D(t)-1` against block/jackknife uncertainties, and test frequency-domain
   detailed balance `G^<(omega)=exp(-beta*omega) G^>(omega)` where spectral
   weight is resolved. Separate elastic zero-frequency weight before an FDT
   check. For a clustering finite-field case, test `G_zz(t)->m_z^2`.
   KMS failure is an expected possible *physical* bias here; passing
   equal-time checks alone is insufficient.
4. **Self-consistent comparison and cost.** Run the same beta, field, spin
   model, symmetry, real/imaginary grids, sampling quality, and propagation
   settings in both modes. Compare error bars and grid refinement, not just
   raw curves. Profile imaginary preparation, factor construction, real
   propagation, and measurement separately. The setup gain is roughly from
   reducing a dense factor of dimension `3*(N_M+2*N_R)` to `6*N_R`, and the
   per-sample gain comes mainly from omitting repeated imaginary propagation
   and the N_M-fold mixed insertion sweep; actual runtime and statistical
   efficiency must be measured.

The option is useful only in parameter regimes where its deviation from the
full contour is smaller than the user's required accuracy after accounting
for Monte Carlo and grid errors. Otherwise, retain the full contour result.

## Current code touchpoints

- `Algorithm/main.cpp`: one current mixed fixed-point loop and per-sample
  imaginary-plus-real trajectory construction.
- `Algorithm/Functions/Functions.cpp`: `connected_contour_primitive`,
  `covariance_from_primitive`, `build_contour_trajectory`, and
  `measure_contour_observables`.
- `Algorithm/Contour/Contour_Kernel.cpp`: reconstructs `MM`, mixed, and real
  branch blocks from one `(t,tau)` primitive.
- `Algorithm/Run_Time_Data/Run_Time_Data.cpp`: current mixed-tensor accumulation,
  complex-ratio normalization, block errors, and termination diagnostics.
- `Algorithm/Parameter_Space/Parameter_Space.cpp` and
  `Algorithm/Storage_Concept/Storage_Concept.cpp`: mode validation, metadata,
  HDF5 schema, and distinct output name.

Background: G. McCaul, C. Lorenz and L. Kantorovich, *Phys. Rev. B* 95,
125124 (2017), https://arxiv.org/abs/1612.05466; M. Garny and M. M. Mueller,
*Phys. Rev. D* 80, 085011 (2009), https://arxiv.org/abs/0904.3600.
