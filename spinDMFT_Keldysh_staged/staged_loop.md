# Matsubara preparation followed by selectable Keldysh Gaussian sampling

## Purpose and physical scope

The new solver first converges the imaginary-time problem, fixes that thermal
preparation, then iterates the mixed and real-time parts of the same three-branch
contour. The branches are Matsubara M, forward +, and backward −. Unlike the
older block-diagonal preparation proposal, **no mixed covariance block is set to
zero**. The initial spin propagator and the real-time fields remain correlated.

Both sampling modes factorize the whole contour with Takagi. The option
`--jointMatsubaraPCN` updates only Matsubara-active real latent coordinates
with pCN and freshly draws their orthogonal complement. Independent draws
remain the default. Both keep the physical MM pseudo-covariance and mean
fixed in Stage 2; the auxiliary Matsubara marginal may change.

Within the spinDMFT Gaussian-bath approximation, an accurately converged
Matsubara solution and a converged compatible continuation should describe the
same contour fixed point as a simultaneous iteration. This statement assumes
consistent contour identities, the same connected closure, exact Gaussian
sampling and ensemble normalization, and controlled propagation/sampling error.
It is not a theorem about convergence of this discrete Monte Carlo map.

The preparation describes the initial interacting equilibrium state. Its
Matsubara component does not depend on subsequent real-time evolution. This is
the same causal ordering underlying equilibrium initialization of
Kadanoff–Baym equations; it does not by itself prove every approximate bath guess
admits a physical equilibrium continuation.

References: Aoki et al., *Rev. Mod. Phys.* **86**, 779 (2014),
[Eqs. (41)–(42)](https://doi.org/10.1103/RevModPhys.86.779);
McCaul, Lorenz, and Kantorovich, *Phys. Rev. B* **95**, 125124 (2017),
[Sections IV–V](https://arxiv.org/abs/1612.05466).

## Correlations, covariance, and means

The common primitive is

\[
X^{ab}(t,\tau)=\langle S_b(-i\tau)S_a(t)\rangle.
\]

It is stored on real-time and imaginary-time vertices, including distinct
one-sided imaginary endpoints at 0 and beta. At t=0 it contains the Matsubara
spin correlation. The endpoints X(t,0) and X(t,beta) supply lesser and greater
real-time orderings. The covariance builder reconstructs MM, M+/M−, ++, −−,
+−, and −+ from this common primitive, respecting contour order and transpose
symmetry. M+ and M− are related by the contour identities; they are not two
unconstrained fitted functions. Real-real correlations are reconstructed as
functions of the time difference, so this remains the source solver's equilibrium
ansatz rather than a general two-real-time nonequilibrium solver.

For spin components a,b the feedback covariance is

\[
\Gamma^{ab}_{AB}(s,s')=J_Q^2\sum_{cd}D_{ac}D_{bd}
 \langle T_{\mathcal C}\delta S_c(s_A)\delta S_d(s'_B)\rangle
 +\Gamma^{ab}_{\rm static}.
\]

The disconnected subtraction in the primitive is the **unconjugated** product
m_a(t)m_b(0), including both measured real and imaginary parts. The deterministic
mean is J_L D Re m(t) on the real branches, and J_L D Re m(0) on M. The external
field and local extra interaction remain in the propagation Hamiltonian.
Coincident same-branch contacts use the inherited greater/lesser half-sum.

All field covariances mean E[delta V delta V^T], not E[delta V delta V^dagger].
The kernel is complex symmetric and need not be Hermitian or positive
semidefinite. No eigenvalue clipping or Hermitian projection is applied.

Measuring products of sampled fields is **not** the self-consistency update.
Those products reconstruct the input kernel. The update is obtained from
propagated spin observables and the coupling equation above.

## Stage 1: the imaginary-time fixed point

The driver makes a parameter copy with one real-time origin and zero real-time
intervals. The imaginary grid, physical parameters, spin symmetry, propagation
scheme, sample blocking, and ratio estimator are preserved. The preparation
sample count and iteration limit can be set independently.

At preparation iteration k:

1. Construct the covariance A_k=Gamma_MM from X_k(0,tau), m_k(0), the coupling
   rotation, and prescribed static noise. Only this block is sampled.
2. Factor A_k with the dense Autonne–Takagi sampler. Draw independent real
   Gaussian coordinates, or use Matsubara pCN as described below, and build the
   complex Matsubara fields. No forward/backward fields are propagated. The two unused origin
   field slots are only placeholders required by the common contour layout.
3. Propagate the spin in imaginary time and retain its unnormalized operator
   M[V_M] and every imaginary-time insertion. Let Z_M=Tr M[V_M].
4. Accumulate spin correlation and magnetization numerators together with Z_M.
   Reduce across MPI ranks and form ratios of ensemble sums. Calculate paired
   real/imaginary ratio jackknife errors using the inherited blocking code.
5. Evaluate the raw, pre-mixing absolute residual as in `spinDMFT_cplx_field`:
   take the largest mean absolute change over tau for each Re/Im correlation
   component and the largest real origin-magnetization change. Require this to
   be strictly below `reliterror * S*(S+1)/(3*sqrt(3*N))`, using the total
   preparation sample count across MPI ranks. Standardized residuals remain
   diagnostics and do not set the preparation stopping threshold. Update both
   the imaginary correlations and initial magnetization with `mixingAlpha`.
6. Apply the inherited covariance/factorization/phase and imaginary-origin-
   magnetization gates, with the preparation iteration limit. The branch
   transpose diagnostic excludes the unused real placeholders at this stage.

The next preparation iteration uses the updated A and mean. Stage 1 factors
only A and Stage 2 factors the whole contour. Their auxiliary Matsubara
marginals can differ even for the same physical MM kernel. In pCN mode,
weight every numerator and denominator by the inverse |Z_M| in both stages.

If preparation reaches its iteration limit, the code writes a failure artifact
and exits with code 2. It does **not** silently proceed with an unconverged bath.
The frozen preparation is the final mixed fixed-point input, not a claim that
the last raw estimate or its standard error equals the exact solution.

After convergence the driver constructs A_star from the final mixed input.
It fixes X_star(0,tau), m_star(0) (both parts of the numerical complex estimate),
and the imaginary deterministic mean. An imaginary-only reference factor is saved for reproducibility; its columns
are not reused in Stage 2.

## The default full-contour Takagi draw

Write the centered fields as M=delta V_M and R=(delta V_+,delta V_−):

\[
\Gamma_k=\begin{pmatrix}A_\star&B_k\\B_k^T&C_k\end{pmatrix}.
\]

Here B_k includes both Matsubara–real branches and C_k includes all four
real-real blocks. A_star is never iterated in Stage 2. The deterministic field
mean is added separately during propagation; the sampled Gaussian has zero mean.

At each real-time iteration the driver constructs and factorizes the complete
complex-symmetric kernel:

\[
\Gamma_k=U_k\Sigma_k U_k^T,\qquad
L_k=U_k\sqrt{\Sigma_k},\qquad
\Sigma_k=\operatorname{diag}(\sigma_{k,1},\ldots,\sigma_{k,n}),\quad
\sigma_{k,j}\geq0.
\]

The nonnegative Takagi values are singular values, not eigenvalues of Gamma.
This construction does not require Gamma to be Hermitian or positive
semidefinite. Negative or complex eigenvalues are not clipped. The inherited
implementation uses an equivalent real-symmetric eigenproblem and omits
numerically zero Takagi modes using its global roundoff/rank threshold.
Reconstruction is checked against the original physical kernel, so a material
loss of its entries is not accepted.

For each sample it draws one independent real normal vector and maps it through
this joint factor:

\[
x_k\sim\mathcal N(0,I),\qquad
\begin{pmatrix}M\\R\end{pmatrix}=L_k x_k.
\]

Partitioning the rows of the factor gives

\[
M=L_{M,k}x_k,\qquad R=L_{R,k}x_k,
\qquad
\begin{aligned}
L_{M,k}L_{M,k}^T&=A_\star,\\
L_{M,k}L_{R,k}^T&=B_k,\\
L_{R,k}L_{R,k}^T&=C_k.
\end{aligned}
\]

Sharing the latent coordinates produces the mixed correlations. Drawing M and
R independently would erase B_k. There is no separate regression, inverse of A,
or real-time Schur-complement factorization in this version.

The dense sampler may split exactly disconnected covariance components and
reuse exactly identical factors with independent latent coordinates. Such
components are identified by exact zeros, not by discarding small mixed
entries. This is mathematically the same joint Takagi sampling construction;
it does not impose a split between Matsubara and real branches.

### What is frozen, and what can change

The fixed quantity is E[MM^T]=A, together with the deterministic Matsubara mean.
For complex Gaussian fields, this unconjugated covariance does not determine
the full probability distribution on their real and imaginary parts. That also
depends on the Hermitian sampling covariance. For canonical Takagi draws it is

\[
H_k=E[v_k v_k^\dagger]=L_kL_k^\dagger
    =(\Gamma_k\Gamma_k^\dagger)^{1/2}.
\]

This equality is exact before numerical rank truncation. Its Matsubara block
can change with B_k and C_k even though A_star stays fixed. The rows L_M,k need
not equal the imaginary-only preparation factor. Neither the complete
Matsubara marginal distribution nor individual Matsubara draws are held fixed.
The pCN option uses this same changing Gaussian reference.

A two-variable example makes the distinction explicit:

\[
\Gamma(b)=\begin{pmatrix}1&b\\b&-1\end{pmatrix},\quad b\in\mathbb R,
\qquad
H(b)=\sqrt{1+b^2}\,I.
\]

The Matsubara pseudo-covariance remains 1, but its Hermitian variance grows
with the mixed coupling b. Similarly, a zero A with nonzero B is representable
by full complex Gaussian sampling: M has zero unconjugated covariance but
need not vanish in every sample. A completely zero joint kernel does produce
exactly zero fields and zero latent dimension.

For holomorphic contour numerators and partition functions, Gaussian Wick
contractions depend on Gamma, with no conjugated fields. When these averages
exist, changing the auxiliary Hermitian covariance at fixed Gamma preserves
their exact Gaussian expectations. Hence imaginary-only observables still use
the frozen A in the ensemble sense, despite the changing auxiliary marginal.
This argument applies to ratios of ensemble averages. It does not justify
averaging separately normalized trajectories, conjugated-field observables,
or claiming equality of finite-sample estimates.

### Why use this representation

Among factors F satisfying FF^T=Gamma with independent real normal latents,
canonical Takagi minimizes the total field second moment:

\[
E[\|v\|^2]=\operatorname{tr}(FF^\dagger)
\ \geq\ \|\Gamma\|_* =\sum_j\sigma_j,
\qquad
\operatorname{tr}(L L^\dagger)=\sum_j\sigma_j.
\]

The inequality follows from the nuclear-norm bound on FF^T by the squared
Frobenius norm of F. The weights here are one per native grid component, not
quadrature weights. This does not prove minimal spin-observable variance, a better phase average, or faster
fixed-point convergence. A particular observable may favor a different gauge.

The cost is a new joint factorization each iteration, with dimension up to the
whole contour rather than only the real residual. Exact disconnected symmetry
blocks can reduce that dimension. Imaginary propagation still occurs for each
sample; preparation freezes the feedback kernel, not propagators or samples.

## Joint Takagi Matsubara pCN

Enable with `--jointMatsubaraPCN`; `--samplingStrategy=pcn` is an alias.
Let the current joint factor have Matsubara rows L_M and real rows L_R.
Use a thin real SVD of

\[
T_M=\begin{pmatrix}\operatorname{Re}L_M\\\operatorname{Im}L_M\end{pmatrix}
\]

to obtain a real orthonormal basis U for its row space in latent coordinates.
The numerical rank keeps singular values above
`epsilon * max(2*nM,q) * largest_singular_value`. The matrix is scaled before
SVD and the relative reconstruction is checked after projection. With
independent standard real Gaussian a and xi, construct

\[
x=Ua+(I-UU^T)\xi,\qquad M=L_MUa,\qquad R=L_Rx.
\]

The projector is applied through U, without allocating its full square
matrix. Since UU^T is an orthogonal projector, E[xx^T]=I. Thus the full
physical pseudo-covariance and Hermitian covariance of the default joint
sampler are retained, up to numerical rank tolerance. Both real and
imaginary Matsubara rows are necessary to identify the correct null space.

pCN updates only a:

\[
a'=\sqrt{1-\epsilon^2}a+\epsilon\zeta,\qquad
p_{\rm accept}=\min\{1,|Z_M(L_MUa')|/|Z_M(L_MUa)|\}.
\]

Acceptance propagates imaginary time only. After every production step,
including rejection, draw fresh xi and construct the joint contour fields.
The sampled measure is proportional to p0(a)p0(xi)|Z_M|. The estimator is

\[
\langle O\rangle=\frac{\sum_s N_{O,s}/|Z_{M,s}|}
                         {\sum_s Z_{M,s}/|Z_{M,s}|}.
\]

Exactly zero modulus proposals are rejected; nonfinite traces abort. Complex
and negative-real partition functions remain in the support. The equations
for Gaussian covariances describe the bare reference, not tilted pCN moments.

Stage 1 applies this pCN acceptance directly to its imaginary-only Takagi
latents. Each rank rebuilds the chain and burns in after every bath update.
Stage 2 also rebuilds U at each iteration. `mhStepSize` lies in (0,1], and
`mhBurnIn` counts discarded proposals. `numSamplesPerCore` counts production
steps, including rejection. All real fields retain their mixed dependence
on the active chain state; only the null-space contribution is independent.

The physical prepared A and mean stay fixed. No additional condition freezes
H_MM, and no Gaussian gauge completion is introduced. This removes the old
sampler's extra covariance choice but does not guarantee small real-time
observable variance: the importance weight contains only the thermal trace.

## Stage 2: mixed and real-time fixed-point iteration

The initial real-time guess is generated with the source solver's Gaussian
real-time envelope applied to the prepared imaginary primitive, including its
imaginary/transverse components. Initial magnetization is constant at the
prepared value. The exact prepared origin slice is copied back after seeding.
This is a seed, not an assumption that the measured correlations must keep that
shape.

At iteration k:

1. Build the full contour kernel from the current feedback primitive and
   magnetization. Check that its reconstructed MM block equals A_star; record
   the mismatch and abort if it exceeds `covarianceTolerance`.
2. Overwrite the MM block with an exact copy of A_star and build a fresh
   Takagi factor of the complete kernel. For pCN, extract the Matsubara-active
   row space of that same factor.
   Check full joint reconstruction against `takagiTolerance` before drawing.
   Sample correlated M,+,− fields without a mixed-block cutoff.
3. Recompute each sample's imaginary propagator and insertion operators. Evolve
   the forward and backward branches separately with the inherited complex
   propagators; backward evolution is not replaced by an adjoint/inverse.
4. Measure the entire X(t,tau) tensor, m(t), Z_M, and the closed-contour trace.
   Retain the mixed interior and both real-time endpoints. All branches feed
   into the next kernel through the spin measurements.
5. Reduce ensemble sums and compute raw correlations, magnetization, and errors.
   For independent sampling the estimator is sum N / sum Z_M; pCN uses
   the inverse-modulus weights above. It is not an average of N/Z_M and is
   not divided by a trajectory's real-time trace.
6. For the **feedback copy only**, restore X_star(0,tau) and m_star(0). Exclude
   these frozen coordinates from the fixed-point residual. Update all t>0
   mixed correlations/endpoints and real-time magnetization with `mixingAlpha`.
7. Restore the origin exactly after mixing, preventing roundoff drift. Assess
   convergence using the active raw residual and inherited diagnostics. The
   imaginary-origin-magnetization diagnostic still uses the raw measured
   quantity; freezing must not hide a changed ensemble.

The optional inherited `constantMagnetization` projection sets the entire
feedback magnetization to the **prepared** value. Raw m(t) still appears in the
output. This optional projection is not enabled by default.

The spin-correlation endpoint values and contacts are reconstructed together
from one primitive; blocks are not independently optimized. The current code
still assumes real-real dependence on time difference. A general driven/quench
two-real-time solver would require a different data contract and closure.

## Convergence, statistics, and expected numerical behavior

Preparation and real time have separate iteration counters, histories, sampling
budgets, and termination states. RNG streams use the stage tag as well as seed,
rank, and iteration, so Stage 2 does not reuse the preparation latent draws.
Preparation uses the absolute threshold defined above. The real-time
standardized convergence statistic uses current raw output errors; it does not
include the uncertainty of the previous iterate or the frozen preparation.

Freezing the preparation removes one source of noisy feedback, but cannot
stabilize an unstable mixed/real-time fixed-point map by itself. Whole-contour
Takagi keeps the canonical dense sampling choice used by the source solver;
the staged feedback map is separate from the choice of sampler. Its auxiliary
Matsubara marginal can vary across Stage 2 and can differ from Stage 1. Reconstruction gates catch numerical
factorization failures, not large observable variance or cancellation in the
ensemble denominator.

The real-time standard errors are **conditional on the saved preparation**.
For uncertainty including preparation, repeat the complete two-stage calculation
with independent preparation runs, or introduce a properly coupled resampling
analysis. Last-iteration preparation errors are stored for inspection; they are
not automatically propagated to Stage 2, nor are they the uncertainty of the
mixed frozen input.

Mixed measurements and per-sample imaginary propagation are retained. The
whole contour is factorized each real-time iteration. No saved
preparation trajectories are carried into the real-time loop.

A converged residual alone does not establish equilibrium physics. Check
m(t)−m_star, the raw origin estimates against the frozen preparation,
closed-contour trace closure, KMS endpoint conjugacy/detailed balance, timestep
refinement, and agreement with the source solver or a controlled bath. The
closed-contour trace remains a diagnostic and is not an added stopping gate.

## Output and failure behavior

Files use the existing `/results` and `/runtimedata` schema. Independent
filenames begin `staged__`; joint Matsubara pCN filenames begin
`staged__jointMpcn_absZ__`. `/parameters` records `joint_matsubara_pcn`,
`field_sampling`, the pCN weight, updated coordinates, null-space refresh
rule, and SVD rank rule. The preparation factor is an imaginary-only
reference. The frozen physical covariance, frozen input, raw preparation
measurements, and per-iteration frozen-MM residuals remain in `/preparation`.
Zero-rank preparation factors have no matrix dataset and an explicit zero
latent dimension. `/results` always contains raw final-iteration estimates.

For pCN, both stages store acceptance and zero-modulus rejection rates,
phase magnitude, blocking curves, and autocorrelation/ESS diagnostics.
These use paired complex ratio influences including denominator fluctuations.
The physical MM and mean remain frozen while the auxiliary marginal may change.

Exit code 0 denotes convergence or completion of a prescribed bath; 2 is
preparation iteration limit; 3 is real-time iteration limit; 1 denotes invalid
input or numerical failure. Failed preparation never starts real-time sampling.

## Supported options and verification

The staged solver supports dense joint Takagi, prefix insertion and
partition-function normalization with independent or joint Matsubara pCN
sampling. Propagators are endpoint CFET, Gauss-CF4 and staggered Strang.
Strang uses the coordinate weights (M,eta,nu)=(1,2,1/2) before inverse scaling;
pCN splits the resulting physical factor and preserves that bare law.
FFT, weighted-dense, closed-contour insertion/normalization, nonnegative FFT
cutoffs, and self-consistent `noselfcons` are rejected. Harmonic prescribed
baths bypass fixed-point preparation for controlled tests.

`preparationSamplesPerCore` and `preparationIterlimit` fall back to their
ordinary counterparts when zero. The obsolete `conditionalSampling` and
`matsubaraGaugeScale` options are removed. Historical files under Data remain
provenance for old implementations, not current validation claims.

Implementation: `Functions/Joint_Matsubara_PCN.*` contains the split;
`Staged_Solver.*` contains the Matsubara chain; `main.cpp` refreshes the null
space after every production step. The new sampler tests verify both complete
Gaussian covariances for weighted, complex, weak, singular and zero kernels.
The pCN tests cover actual rejected steps, analytic complex-phase averages,
normalization invariance, exact thermal spins, harmonic-bath agreement with
independent draws, and short self-consistency with mixingAlpha=1.
See README.md for build and test commands. These checks do not establish
long-time production stability or reduced magnetization drift.
