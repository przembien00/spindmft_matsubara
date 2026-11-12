# CspinDMFT_Keldysh efficiency fixes — 2026-11-12

Implemented the prioritized changes in the [assessment](../efficiency_audit_2026-11-12/REPORT.md), preserving the dense Gaussian sampling law, per-rank RNG streams, numerical cutoffs, observable schema, and statistical estimator.

## Changes

- `Algorithm/main.cpp`: pass temporary covariances into the sampler constructors. Covariance storage is released before sampling in both stages.
- `Algorithm/Functions/Sampling.cpp`: retain blocks locally, reduce global numerator/partition totals, then reduce delete-block ratio means and centered real/imaginary variances. pCN coarsens adjacent blocks within each chain in place and retains the maximum error across levels. Invalid local denominators/ratios are coordinated across ranks before subsequent collectives. Additional statistical workspace is proportional to observable width, independent of MPI rank count.
- `Algorithm/Functions/Functions.cpp`: reuse adjacent nu-boundary exponentials, compute each requested transformed spin operator once per measurement time, and construct only requested right-hand Matsubara insertions. The forward/backward non-Hermitian propagators remain distinct.
- Initialization constructs isolated correlations only when an isolated initializer is selected, and caches repeated imaginary-time exponentials when they are needed. Initial magnetization still uses its required isolated calculation.
- Partition-only pCN proposals use a running imaginary-time product, with no prefix/suffix/step arrays or observable payloads. The existing partition-only API returns empty observable containers.
- `cpp_libs/Complex_Gaussian/Complex_Gaussian.cpp`: call the same symmetric LAPACK eigensolver in place on one real lift, release it before reconstruction, scale entries during exact block extraction, and compare repeated blocks through their original rows. Full scaled covariance, unique covariance copies, eigensolver input/output copies, and full diagnostic difference matrices are eliminated. Reconstruct each final unique factor once, after the unchanged global rank cutoff, and release each reconstruction before the next one.
- Physical reconstruction diagnostics include couplings lost through underflow at extreme explicit weights. Tolerated transpose asymmetry retains the former symmetric-adapter projection. Tests cover both cases.

The shared Gaussian changes also apply to `spinDMFT_Keldysh_staged`; its full test suite was rebuilt and run. The separate `spinDMFT_Keldysh` Gaussian implementation was not changed. Existing unrelated workspace changes were preserved.

## Validation

- **4/4 cluster CTests passed**, including a new pooled-reference jackknife test and tests for partition-only complex propagation and physical covariance underflow.
- **17/17 staged-solver CTests passed** against the shared Gaussian implementation.
- The new jackknife test also passed with **two MPI ranks**. Its small test-only allgather reproduces the original estimator for an independent oracle. It checks A/B/C/D, both stage types, complex denominators, pCN weighting/coarsening, real/imaginary errors, the blocking envelope, and a denominator failure confined to one rank. Production code performs no block allgather.
- **24 complete old/new run pairs** passed: A/C/D × independent/pCN × generated/isolated initialization × one/two MPI ranks. All 61 preparation/results/runtime datasets and attributes per case were compared, including standard errors, blocking diagnostics, acceptance rates, iteration histories and termination. Maximum absolute difference was **7.11e-15**. On one rank it was **2.75e-26**. These deliberately small-sample tests include multi-iteration feedback. Two one-rank D/pCN cases reached the real-time iteration limit in both versions; the remaining cases converged in both versions. Failed preparation was not counted as a complete real-time comparison.
- **18 Gaussian-map comparisons were exactly equal** between the saved original implementation and the final implementation: physical real-linear factors, `E[z z^T]`, `E[z z^dagger]`, and 19 fixed-seed draws per case. Fixtures cover arbitrary complex dense matrices, repeated blocks, zero modes, the global cutoff, tolerated asymmetry, and all three weighting interfaces. Both implementations also reject the extreme-weight underflow fixture. Equal real-linear maps establish equality of the real joint Gaussian law for these fixtures; this is stronger than pseudo-covariance reconstruction alone.

See `comparison_r1.json`, `comparison_r2.json`, and `gaussian_comparison.json` for results and executable/source fingerprints. Timing and run IDs are excluded from numerical comparisons. No empirical long-time/continuum accuracy claim follows from these regression checks.

## Matched measurements

One fresh process per version and case on local macOS arm64, Release DOUBLE, one MPI rank and one library thread. Both versions use beta=5, FM cubic `J=-1/sqrt(6)`, 50 imaginary/150 real intervals, `Tmax=10`, one substep, 32 samples per stage, eight blocks, batch size one, and seed 613. A has zero field; C has uniform field 0.5. Each run completes one preparation and one real iteration with `--noselfcons`, exit status zero. These are matched diagnostic runs, not converged production benchmarks.

<!-- BENCHMARK_TABLE -->

| Sites / symmetry | Setup s, old → new | Sampling s, old → new | Total wall s, old → new | Peak RSS GB, old → new |
|---|---:|---:|---:|---:|
| 2 / A | 0.353 → 0.291 | 0.0361 → 0.0308 | 0.494 → 0.424 | 0.350 → 0.166 |
| 2 / C | 10.573 → 10.463 | 0.0852 → 0.0771 | 10.771 → 10.636 | 0.667 → 0.392 |
| 4 / A | 2.903 → 2.735 | 1.4182 → 1.0781 | 4.707 → 4.073 | 1.308 → 0.572 |
| 4 / C | 88.071 → 88.415 | 2.9218 → 2.1934 | 91.509 → 91.031 | 2.266 → 1.401 |

For four-spin C, peak RSS fell from 2.266 to 1.401 GB (38.2%), and mixed/real-time sampling time fell from 2.922 to 2.193 s (24.9%). Total runtime was nearly unchanged because setup still dominates. Four-spin A sampling time fell by 24.0% and its total wall time by 13.5% in this measurement. All four matched benchmark output pairs also passed numerical comparison; the largest stored difference was 3.15e-26 (`benchmark_comparison.json`).

Sampling is the normal mixed/real-time sampling phase, including Gaussian draws, trajectory measurement and progress reporting. Setup is the joint field-distribution phase. Peak RSS is the fresh process high-water mark from macOS `/usr/bin/time -l`, in decimal GB. It includes all phases and is not a phase-isolated measurement. Allocator retention and run-to-run timing variation remain relevant; these single-run measurements do not establish a general end-to-end speedup.

An initial implementation reconstructed each factor twice and regressed setup time. That duplicate product was removed before final validation. `exploratory_benchmarks.json` preserves the exploratory measurements; `benchmarks.json` records the final executable measurements against the saved baseline.

The original four-spin C sweep's 16-rank, 29-block pooled array occupied 1.149 GB **per rank**. That array is removed. Two new complex vectors of observable width add approximately 4.95 MB per rank for that shape, giving approximately **1.144 GB net statistical payload savings per rank**, in addition to other phase-specific changes. This is an allocation calculation, not a measured 16-rank peak RSS. Do not add this estimate to the one-rank peak-RSS reduction: statistical reduction and factorization have different lifetimes. No new seven/eight-spin production feasibility claim is made.

## Scope of follow-up candidates

`numSamplesPerSet` and sweep defaults remain unchanged. The assessment's existing frozen-factor batching measurements support selective tuning; cluster-hardware throughput and later-iteration ranks still need measurement before changing production settings.

The assessment's broader candidates remain separate work: constructing covariance blocks directly from the model, cross-rank/node factor sharing, a common independent/pCN packed-factor representation, tensor-product local nu exponentials, alternative contraction kernels, and paired Padé/linear-solve exponential kernels. This change implements the prioritized lifetime, statistics, reuse, initialization and eigensolver fixes without changing those numerical kernels or introducing a sampler option.

## Reproduction

The original executable was saved as `before.out` before editing, after confirming the baseline build was current. `shared_before.cpp` is the unchanged original shared Gaussian source. Both baseline snapshots are retained in this directory. Build local solvers as usual, then from the repository root:

```sh
ctest --test-dir CspinDMFT_Keldysh/Algorithm/build --output-on-failure
ctest --test-dir spinDMFT_Keldysh_staged/Algorithm/build --output-on-failure
mpirun -n 2 CspinDMFT_Keldysh/Algorithm/build/test_sampling
python3 CspinDMFT_Keldysh/Timing/efficiency_fixes_2026-11-12/compare_runs.py
python3 CspinDMFT_Keldysh/Timing/efficiency_fixes_2026-11-12/compare_runs.py --ranks 2
python3 CspinDMFT_Keldysh/Timing/efficiency_fixes_2026-11-12/compare_gaussian.py
python3 CspinDMFT_Keldysh/Timing/efficiency_fixes_2026-11-12/benchmark.py
```

Python needs NumPy/h5py. The Gaussian comparison compiler paths and benchmark RSS collection describe this local macOS installation; adapt them on other machines. `--after-only` reruns the final benchmark executable while retaining saved baseline measurements. Raw HDF5 outputs and logs are under the per-case directories; `manifest.json` fingerprints the final source and executable.
