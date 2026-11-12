# On-demand covariance blocks and optional reconstruction — 2026-11-12

The cluster solver no longer constructs a full D-by-D covariance matrix in either stage. It evaluates entries from the bath primitive, discovers exact connected components, and constructs one complex block at a time. After conversion to the real lift, that block is released before the eigensolve. The sampler retains factors and row mappings after construction; the provider and bath primitive are released.

The full matrix was an input-interface requirement of the earlier implementation, not a mathematical requirement. Disconnected blocks can be factorized independently. Blocks must include every nonzero coupling, including Matsubara/eta/nu cross-couplings. A fully connected problem can still have a single D-by-D block. Discovery still scans entries with O(D²) work; this change removes full-matrix storage, not all quadratic work.

The entry provider uses the same kernel, coupling rotations, connected correlations, static noise, interpolation, and canonical symmetric entries as the dense reference helper. It does not infer field blocks from observable symmetry. Existing block order, identical-block sharing, local/global singular-value thresholds and latent RNG ordering are retained. The dense helper remains available for reference tests, while production `main.cpp` uses the provider.

## Reconstruction option

Gaussian reconstruction checking is **off by default** in `CspinDMFT_Keldysh`.

```sh
./executable_DOUBLE.out ... --checkGaussianReconstruction
```

The flag also accepts explicit `true`/`false`. When enabled, the sampler reconstructs final factors after the global cutoff, evaluates the physical unweighted residual, and enforces `takagiTolerance` as before. When disabled, it skips all covariance-factor reconstruction products and residual scans. Eigensolver failures and nonfinite covariance inputs remain errors.

The HDF5 parameter attribute `check_gaussian_reconstruction` records yes/no. Skipped `gaussian_factor_reconstruction_errors` entries are NaN rather than a fabricated zero and are excluded from convergence conditions. Both preparation and real time use the flag. The shared library keeps checking enabled by default for existing callers; staged-solver defaults are unchanged. The standalone checked `autonne_takagi` API is preserved.

## Validation

- All **4 cluster CTests** and **17 staged-solver CTests** passed.
- **48 complete old/new run pairs** passed: A/C/D × independent/pCN × generated/isolated initialization × one/two MPI ranks × check disabled/enabled. Observable values, standard errors, blocking statistics, pCN acceptance, iteration histories, and termination outcomes were exactly equal. Disabled reconstruction histories intentionally changed to NaN; enabled histories were exactly equal to the baseline. Two one-rank D/pCN cases reached the same real-time iteration limit in both versions; the remaining cases converged.
- The Gaussian comparison fixtures produce exactly equal factors, pseudo-covariance `E[z z^T]`, Hermitian covariance `E[z z^dagger]`, and 19 fixed-seed draws in each of 18 cases. Additional comparisons inside those fixtures check the provider and disabled-check paths against the checked dense path. Tests include repeated blocks, the global cutoff, zero modes, tolerated asymmetry, and physical-residual rejection at extreme weights.
- A cluster-core test compares every provider entry against the dense reference with a mixed bath-coupling matrix and static noise. It checks exact physical factor maps for checked/unchecked providers. CLI tests cover flag enablement, persisted yes/no settings, unmeasured NaN residuals, and unchanged convergence. Nonfinite covariance rejection remains tested with checking disabled.

See `comparison_r{1,2}_{default,checked}.json` and `gaussian_comparison.json` for results and fingerprints.

## Matched measurements

<!-- BENCHMARK_TABLE -->

| Sites / symmetry | Setup s, old → new | Sampling s, old → new | Total wall s, old → new | Peak RSS GB, old → new |
|---|---:|---:|---:|---:|
| 2 / A | 0.296 → 0.184 | 0.0306 → 0.0303 | 0.418 → 0.301 | 0.162 → 0.091 |
| 2 / C | 10.517 → 1.525 | 0.0757 → 0.0754 | 10.699 → 1.699 | 0.358 → 0.286 |
| 4 / A | 2.782 → 1.613 | 1.0599 → 1.0699 | 4.107 → 2.949 | 0.569 → 0.284 |
| 4 / C | 88.001 → 14.901 | 2.2326 → 2.1904 | 90.665 → 17.515 | 1.334 → 1.060 |

For four-spin C, peak RSS fell from 1.334 to 1.060 GB (20.6%), and total wall time fell from 90.7 to 17.5 s. The sampling phase was essentially unchanged; the gain is in distribution setup. All four matched large-grid output pairs had exactly equal observables, error estimates and runtime diagnostics, apart from the deliberately unmeasured reconstruction histories.

These are fresh-process runs of the saved previous executable and the new default solver, with beta=5, FM cubic J=-1/sqrt(6), 50 imaginary/150 real intervals, one substep, 32 samples per stage, eight blocks, batch size one and seed 613. A has zero field; C has uniform field 0.5. Both versions complete one preparation and one real-time iteration with `--noselfcons`, exit status zero. The baseline has reconstruction checking enabled, as it was previously mandatory. These results combine optional-check removal and the provider change; they do not isolate the speedup of each change or establish production convergence.

Peak RSS is collected with `RUSAGE_CHILDREN` in a separate fresh Python worker for each solver invocation. Each worker launches exactly one solver child, avoiding cumulative high-water marks across previous benchmark cases. macOS reports bytes. Raw logs/HDF5 are under `benchmarks/`; `benchmarks.json` records commands, statuses and binary hashes. `benchmark_comparison.json` checks numerical output, excluding the intentionally disabled reconstruction residuals.

For four-spin C, the provider's full bath primitive occupies about 17.7 MB. It replaces the 285.5 MB full covariance during construction. The largest real lift and LAPACK workspace remain necessary for the chosen eigensolver. Actual peak RSS includes these arrays, factors, other solver data, libraries, padding and allocator retention. No production multi-rank or seven/eight-spin feasibility claim follows from these measurements.

## Reproduction

The saved baseline executable is `before.out`; shared Gaussian source/header snapshots are included for the comparison probe. Build the solvers normally, then from the repository root:

```sh
ctest --test-dir CspinDMFT_Keldysh/Algorithm/build --output-on-failure
ctest --test-dir spinDMFT_Keldysh_staged/Algorithm/build --output-on-failure
python3 CspinDMFT_Keldysh/Timing/block_covariance_2026-11-12/compare_runs.py
python3 CspinDMFT_Keldysh/Timing/block_covariance_2026-11-12/compare_runs.py --check
python3 CspinDMFT_Keldysh/Timing/block_covariance_2026-11-12/compare_runs.py --ranks 2
python3 CspinDMFT_Keldysh/Timing/block_covariance_2026-11-12/compare_runs.py --ranks 2 --check
python3 CspinDMFT_Keldysh/Timing/block_covariance_2026-11-12/compare_gaussian.py
python3 CspinDMFT_Keldysh/Timing/block_covariance_2026-11-12/benchmark.py
python3 CspinDMFT_Keldysh/Timing/block_covariance_2026-11-12/summarize.py
```

Python requires NumPy/h5py. The Gaussian comparison compiler paths and RSS units describe this local macOS setup. Final source and binary fingerprints are in `manifest.json`.
