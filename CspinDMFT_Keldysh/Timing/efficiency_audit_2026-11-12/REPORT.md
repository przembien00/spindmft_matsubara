# CspinDMFT_Keldysh efficiency and memory audit — 2026-11-12

Investigation of the current working tree. Solver sources were not changed. The normal executable was rebuilt; all three CTest checks passed (`cluster_core`, `cluster_solver`, `algorithm_independence`). The files in this directory are an instrumented audit executable, its generation/build/run scripts, logs and source fingerprints.

## Scope and measurement limits

Local macOS arm64, Release `-O3`, DOUBLE, Accelerate LAPACK/BLAS, one MPI rank and one library thread. The probes use the current sweep grid: 50 imaginary intervals, 150 real intervals, `Tmax=10`, one real-time substep, beta=5, FM cubic J=-1/sqrt(6); A has zero field and C has uniform field 0.5.

Each probe completes one preparation and one real iteration with `--noselfcons`, four samples per stage, two blocks, seed 613. Within the resulting real-stage **frozen field law**, it separately times 128 Gaussian draws for each batch size 1/8/32/64, after warming the shared batch-factor storage. Each size starts from the same sampler/RNG state. It also measures 16 fixed field trajectories with timers around the main phases. These are diagnostic probes, not converged simulations or production-throughput measurements. A first-iteration law can be unusually low rank, especially in A. Later self-consistency iterations, other fields and other machines can have different costs.

Use the `AUDIT` lines for sample timings: the normal main-loop “mixed real-time sampling” timer includes the extra probe work. Peak RSS is the process high-water mark taken after the extra measurements; it includes factorization and small probe buffers, but does **not** measure a 16-rank production reduction. Values below use decimal GB (10^9 bytes).

The older `Timing/cubic_beta_sweep/measurements.json` uses executable hash `908412be...` and the former 80-imaginary/200-real grid, as confirmed by its logs. Those measurements must not be compared with the present grid as an optimization speedup. Its child-RSS figures are cumulative high-water marks across completed subprocesses, whereas this probe obtains its own process RSS.

## Current measurements

| N / symmetry | Joint setup s | Draw ms/sample, batch 1 | Draw ms/sample, batch 32 | Trajectory ms/sample | Peak RSS GB |
|---|---:|---:|---:|---:|---:|
| 2 / A | 0.327 | 0.0169 | 0.0140 | 1.123 | 0.351 |
| 2 / C | 10.540 | 0.4866 | 0.0636 | 2.085 | 0.672 |
| 4 / A | 2.859 | 0.0418 | 0.0310 | 44.674 | 1.310 |
| 4 / C | 88.711 | 1.8782 | 0.1970 | 88.417 | 2.198 |

The four-spin C trajectory breaks down into 3.77 ms imaginary evolution, 4.54 ms Matsubara insertions, 22.31 ms real propagation, and **57.77 ms observable evaluation**. Thus optimize the measurement/propagation path before expecting large overall gains from field batching.

For N=2 C, batching 32 reduces the draw cost about 7.65x; for N=4 C, about 9.53x. Combining measured draw and trajectory costs predicts only approximately **16.5%** and **1.9%** reductions, respectively, in those two per-sample components. These exclude setup, statistics and I/O, and are not measured whole-run speedups.

Across all batches the maximum fixed-factor field difference was **6.01e-15**. The fresh C probe RSS was 2.20 GB for four spins; an earlier cold-cache probe reached 2.45 GB. Treat this as observed allocation/RSS variation, not a measured optimization. The older large-grid 4.74 GB result is not a current-grid measurement.

| N / symmetry | Field dimension | Retained latent dimension | Largest complex block |
|---|---:|---:|---:|
| 2 / A | 2112 | 306 | 702 |
| 2 / C | 2112 | 1506 | 1404 |
| 4 / A | 4224 | 612 | 1404 |
| 4 / C | 4224 | 3012 | 2808 |

## Highest-value memory changes

### 1. Keep jackknife blocks distributed

`Algorithm/Functions/Sampling.cpp:90-91` allocates `pooled` with all ranks' blocks **on every rank**, then repeats the entire jackknife calculation on every rank. This gives per-rank memory proportional to rank count and aggregate memory proportional to its square.

Let P be the number of stored representative pairs, c the stored direction count, nI/nR the imaginary/real interval counts, N the number of sites, B the actual block count per rank, and R the MPI rank count. For C symmetry:

```
Q = P*c*(nR+1)*(nI+1)                    correlation entries
W = 1 + Q + N*(nR+1) + (nR+1)           partition, correlations, moments, closure
local blocks = 16*B*W bytes
pooled blocks = 16*R*B*W bytes PER RANK
```

For A the magnetization term vanishes. The live cubic configurations have P=2/5/13/14 for N=2/4/7/8. The current 16-rank sweep's sample counts and default requested 32 blocks actually select:

| N / symmetry | Samples per rank | Actual B | Extra pooled GB per rank |
|---|---:|---:|---:|
| 2 / A | 368700 | 30 | 0.119 |
| 2 / C | 116700 | 30 | 0.477 |
| 4 / A | 9900 | 30 | 0.297 |
| 4 / C | 2900 | 29 | **1.149** |

The four-spin C buffer alone is 18.39 GB across the node. Its local block array is only 71.8 MB per rank. These are exact payload calculations, excluding allocator overhead. A one-rank/four-sample probe does not test this production memory cost.

Preserve the existing estimator by reducing total numerator/denominator sums first. Each rank then evaluates its own delete-block ratios `(T_j-B_bj)/(Z-B_bz)`. Allreduce the ratio sums to obtain the global jackknife mean, then Allreduce centered squared deviations in a second pass. Apply the same `(R*B-1)/(R*B)` prefactor. For pCN, coarsen blocks locally within each chain and repeat the existing maximum-over-levels rule. Coordinate invalid-denominator checks across ranks. This removes the pooled array with O(W) additional work buffers and distributes the calculation. Floating-point reduction order can change; validate standard errors, not just means.

Gathering only onto rank zero is a smaller patch but leaves the root memory bottleneck. Reducing `numBlocks` would instead change statistical resolution and is not the preferred fix.

### 2. Reduce factorization workspace and matrix lifetimes

The real field dimension is

```
D = 3*N*(nI + 2*nR*q + 2)
```

where q is `realTimeSubsteps`. One complex DOUBLE D-by-D matrix costs `16*D^2` bytes. `main.cpp:100-102` keeps the original covariance alive for the entire real iteration, including sampling and MPI error calculation. It is not used after sampler construction. A construction scope/factory can release it immediately afterwards, without changing the sampled law. This saves sampling/reduction live memory, **not the factorization peak**.

| Sites | D, current grid | One full covariance GB | All Matsubara insertion matrices GB |
|---|---:|---:|---:|
| 2 | 2112 | 0.071 | 0.000078 |
| 4 | 4224 | **0.285** | 0.00251 |
| 7 | 7392 | 0.874 | 0.281 |
| 8 | 8448 | 1.142 | **1.283** |

The last column assumes spin 1/2 and is `16*(3*N)*(nI+1)*(2^N)^2`. It excludes prefix/suffix/step matrices and observables. At N=8 those three trajectory arrays add another 0.159 GB. N=7/8 were not run in this audit.

`cpp_libs/Complex_Gaussian/Complex_Gaussian.cpp:90-105` constructs a real 2b-by-2b lift and a separate eigenvector matrix for a b-by-b complex block. The installed Blaze `Eigen.h:416-423` additionally copies the lift into `Atmp` before DSYEVD and then copies the eigenvectors out. Blaze's `syevd.h:179` allocates about `2*(2b)^2` doubles of work memory.

Use an appropriate in-place symmetric eigensolver on one real lift, consuming its overwritten eigenvectors directly. Removing the separate lift-copy/output buffers eliminates roughly `64*b^2` bytes of explicit allocation: **505 MB for the measured largest four-spin C block b=2808**. The LAPACK workspace remains. Scope the eigensolver storage so it is destroyed before allocating reconstruction matrices. This is an allocation estimate, not a measured new peak-RSS result; lazy allocation and allocator retention affect actual RSS.

Other avoidable quadratic allocations are the full scaled covariance (`Complex_Gaussian.cpp:390`), full diagnostic differences (`:64`, `:404`), and simultaneously retained unique covariance blocks/reconstructions. Evaluate Frobenius residuals directly or in tiles rather than materializing the difference. Verify the physical, unweighted residual and the existing global singular-value cutoff remain unchanged.

The independent batch sampler also retains a real-packed copy of every unique complex factor (`:293-299`), with the same payload size as the complex factor. A common packed representation or separate independent/pCN storage path could avoid that duplication. The saved amount is `16*sum(unique_block_rows*retained_rank)` bytes; it depends on the actual numerical ranks.

### 3. Build exact blocks directly and share factors between ranks

The implementation already finds disconnected covariance blocks, reuses exactly identical blocks and drops numerically zero singular directions (`Complex_Gaussian.cpp:163-246`). Simply suggesting “exploit symmetry” or “use low rank” would duplicate existing work.

The opportunity is to avoid first constructing/scaling/scanning the **full** D-by-D matrix before extracting those blocks. Add a block-oriented covariance interface, preserving exact allowed couplings, block ordering and the global rank threshold. The final nu endpoint has exactly zero covariance and already becomes singleton zero blocks in the current implementation. Known spin-direction block structure can also eliminate full-matrix storage, but must be derived from the actual bath coupling/model, not assumed from the observable symmetry alone.

All MPI ranks currently build and factor the same covariance after global reduction (`main.cpp:100-102`). Factorizing on one rank and broadcasting factors avoids redundant eigensolves. Node-local shared factor storage can further reduce node memory. This does not imply a 16x wall-time speedup: currently ranks factor concurrently; the likely benefits are less contention and lower aggregate resource use. Root peak memory and broadcast costs still need measurement. Preserve per-rank RNG streams and coherent error handling.

## Highest-value runtime changes

### 4. Remove repeated work in trajectory propagation and observables

`Algorithm/Functions/Functions.cpp:157-160` evaluates four dense matrix exponentials per real substep. The right nu-boundary exponential at one substep is exactly the left exponential at the next. Cache this boundary propagator. With 150 substeps this reduces nu exponentials from 300 to 151, and all real-time exponentials from 600 to 451, with no new approximation. It does not remove 25% of total trajectory time because the rest of propagation and measurement remains.

The nu Hamiltonian contains only one-site field terms (`node_h(...,false)`); terms on distinct sites commute. A later optimization can form/apply its exponential through small local spin matrices and tensor operations. Keep the interacting eta propagator distinct. Complex contour Hamiltonians are generally non-Hermitian: backward propagation must not be replaced by an adjoint.

At each stored time, `Functions.cpp:166` and `:171` repeatedly compute `backward*S_i^a*forward` for different pairs/directions. Compute each requested `(i,a)` once per time and reuse it. For the current N=4 C categories, 20 correlation transformations plus 4 magnetization transformations reduce to 8 unique transformations. This removes 32 of 48 dense matrix products in that part of each time slice. The subsequent tau contractions still remain.

`Functions.cpp:146` builds Matsubara insertions for **all 3N operators**, although A only measures xx. Build only the right-hand operators required by the configured pair/direction list. In the cubic A cases this reduces this insertion set from 3N to N, cutting its storage and construction work by two thirds. This becomes substantial at seven/eight spins.

The scalar `contract()` and correlation accessors are also hot-loop candidates: arrange transposed/flattened insertions once and evaluate groups of contractions using contiguous loops or matrix multiplication; precompute storage offsets instead of repeatedly looking up the direction/category for every tau. Benchmark at the actual Hilbert-space size, since small matrices may not benefit from an extra BLAS call. Reusable matrix/observable workspaces can remove per-sample allocation traffic.

The Padé exponential at `Functions.cpp:37` forms `inv(V-U)*(V+U)`. Solve `(V-U)*X=V+U` directly instead. Paired positive/negative eta exponentials can share their matrix powers/Padé intermediates. Both require numerical comparison on arbitrary complex, non-Hermitian fields and full trajectories.

### 5. Use the existing batching option selectively

`numSamplesPerSet` defaults to one (`Parameter_Space.cpp:75`). The code already uses real DGEMM for Gaussian batches. The measurements above quantify the available option without changing the algorithm. C benefits substantially in the Gaussian-draw kernel; A has low-rank first-iteration factors and much less benefit. Propagation/measurement dominates the measured four-spin workloads, so Gaussian-only speedups must not be reported as whole-solver speedups.

Batching preserves the frozen factor and latent stream; the measured field differences are roundoff-sized. This audit does not establish identical complete self-consistent results after changing batch size: tiny changes can feed into subsequent degenerate factorizations and their coordinate bases. Validate complete runs and statistics before changing the default or sweep commands. Batch working memory is linear in batch size; tune it on the cluster.

### 6. Skip calculations whose results are discarded

`generate_initial_correlations()` always computes isolated-cluster correlations (`Functions.cpp:73`), even though the default `imagtime`/`zero` choices overwrite all those entries at `:83-84`. Allocate/fill the requested initial shape directly unless an isolated component is requested. Keep the separately required isolated initial magnetization. When isolated correlations are requested, cache their repeated time exponentials outside the pair/direction loops. At N=8 C on this grid the current initializer evaluates 5715 matrix exponentials, overwhelmingly recomputations; the defaults discard the resulting correlations. This is a startup optimization, not an iteration speedup.

For pCN, `partition_only` returns only after creating every imaginary step, prefix, suffix and observable container (`Functions.cpp:133-144`). The partition needs just the running forward product. Add a dedicated partition evaluator with O(H^2) trajectory storage and no suffix/observable construction. Keep acceptance ratios and RNG semantics unchanged. pCN itself is already using a Matsubara row-space method without constructing a full null-space basis; that part is not a missing optimization.

## Order of implementation and validation

1. Release covariance after setup, skip discarded initialization and make the partition-only path minimal. These are scoped lifetime/control-flow changes.
2. Replace replicated jackknife storage with distributed paired-ratio statistics. This is the largest clear multi-rank memory opportunity.
3. Cache nu-boundary exponentials and unique operator transforms; construct only needed insertions. The timings put these on the dominant per-sample path.
4. Replace redundant eigensolver buffers and full-matrix diagnostics. Then consider a direct block covariance interface and MPI factor sharing.
5. Tune batching and dense kernels on actual cluster hardware, after the above measurements are repeatable.

For claimed behavior preservation, compare complete same-seed old/new runs and HDF5 preparation/results/errors, not builds alone. Use both independent and pCN modes, A/C/D, non-Hermitian fields, and multiple ranks where relevant. Factorization/RNG-coordinate changes also require checking both pseudo-covariance `E[z z^T]` and Hermitian covariance `E[z z^dagger]`, the real joint Gaussian law and observable variance. A different sampler that merely reconstructs `E[z z^T]` is not automatically an equivalent memory optimization.

Measure peak memory separately during factorization, sampling and statistical reduction, including later self-consistency iterations with production block counts. These peaks generally occur at different times and should not simply be added. Preserve grid, sample counts, convergence criteria and exit status in timing comparisons. No new production memory limit or seven/eight-spin feasibility claim is established by this audit.

The four completed audit runs before and after adding trajectory phase timers have exactly equal datasets under `preparation`, `results` and `runtimedata` (25 datasets for each two-spin case; 49 for each four-spin case). This checks instrumentation only, not any proposed optimization. See `instrumentation_comparison.json`.

Reproduce locally after building the solver:

```sh
sh CspinDMFT_Keldysh/Timing/efficiency_audit_2026-11-12/build_probe.sh
python3 CspinDMFT_Keldysh/Timing/efficiency_audit_2026-11-12/run_probes.py
```

The build script records the local dependency paths; adapt them for another machine. `runs.json` records commands, statuses and the probe executable hash; `manifest.json` records solver source hashes.
