# Release covariance blocks before eigensolver workspace allocation

The internal dense sampler now constructs and validates the real lift, releases the copied complex `block` allocation, and then calls the eigensolver. A scoped swap with an empty matrix actually releases the allocation; changing dimensions alone would not establish that memory was freed.

`takagi_lift()` constructs the real matrix. `factor_from_lift()` takes ownership of it and computes the factor. The sampler's later reconstruction diagnostics continue to read the original full covariance through the existing row indices. The public `autonne_takagi(const Matrix&)` interface retains its caller-owned input for its own reconstruction diagnostic and does not consume or mutate that input.

For the largest four-spin C block (`b=2808`), the released complex matrix contains `16*b*b = 126,157,824` bytes. Its lifetime no longer overlaps the approximately 505 MB LAPACK work array. The main explicit covariance/block/lift/work payload at the largest eigensolve falls from approximately 1.169 GB to 1.043 GB. Other allocations and allocator retention affect actual RSS.

## Validation

- All 4 cluster CTests and all 17 staged-solver CTests passed.
- All 24 complete old/new run pairs (A/C/D, independent/pCN, default/isolated initialization, one/two MPI ranks) had exactly equal stored preparation, result and runtime data. Termination and iteration counts were unchanged, including the same two one-rank pCN iteration-limit cases.
- All 18 Gaussian comparisons were exactly equal for the factors, pseudo-covariance, Hermitian covariance and fixed-seed draws. The extreme-weight physical-residual rejection also passed in both versions.

## Memory measurements

<!-- BENCHMARK_TABLE -->

| Sites / symmetry | Setup s, old → new | Sampling s, old → new | Total wall s, old → new | Peak RSS GB, old → new |
|---|---:|---:|---:|---:|
| 2 / A | 0.291 → 0.289 | 0.0308 → 0.0298 | 0.424 → 0.449 | 0.166 → 0.166 |
| 2 / C | 10.463 → 10.634 | 0.0771 → 0.0769 | 10.636 → 10.809 | 0.392 → 0.360 |
| 4 / A | 2.735 → 2.946 | 1.0781 → 1.0562 | 4.073 → 4.268 | 0.572 → 0.572 |
| 4 / C | 88.415 → 87.750 | 2.1934 → 2.1921 | 91.031 → 90.363 | 1.401 → 1.337 |

For four-spin C, measured peak RSS fell from 1.4006 to 1.3368 GB, about 63.8 MB (4.6%). This is smaller than the 126.2 MB reduction in overlapping live payload: process high-water RSS is affected by other phases and allocator retention. These measurements do not isolate the reason for that difference. Four-spin A's whole-process peak was essentially unchanged. All four benchmark output pairs were exactly equal.

These are one-rank fresh-process high-water RSS measurements, with 50 imaginary/150 real intervals, one substep, beta=5, 32 samples per stage, eight blocks, batch size one and seed 613. Each run completes one preparation and one real-time iteration with `--noselfcons`, exit status zero. They do not establish converged production throughput or multi-rank peak memory.

The baseline measurements are reused from the immediately preceding implementation's `efficiency_fixes_2026-11-12/benchmarks.json`; their executable hashes were verified against the saved `before.out` before reuse. New measurements and commands are in `benchmarks.json`. Raw outputs are checked in `benchmark_comparison.json`; source and binary hashes are in `manifest.json`.

The snapshots and comparison scripts in this directory compare this isolated follow-up against the previous implementation. The compiler paths in `compare_gaussian.py` and RSS collection in `benchmark.py` describe the local macOS installation.
