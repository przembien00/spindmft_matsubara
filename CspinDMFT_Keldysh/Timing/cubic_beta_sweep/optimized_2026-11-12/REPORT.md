# Optimized cubic sweep calibration, 2026-11-12

The updated plan fits ten real-time iterations plus up to twenty preparation
iterations within a conservative 6.4-hour estimate and an eight-hour Slurm
allocation. Counts are per MPI rank per iteration; production uses 16 ranks.
Both coupling signs use the same conservative size/symmetry budget. C covers
both nonzero-field cases.

| Sites / symmetry | Real samples/rank | Preparation samples/rank | Max measured real sample (ms) | Max setup/iteration (s) | Inflated total (h) |
|---|---:|---:|---:|---:|---:|
| 2 / A | 1,066,900 | 106,690 | 1.050 | 0.208 | 6.400 |
| 2 / C | 454,400 | 45,440 | 2.494 | 1.517 | 6.399 |
| 4 / A | 30,000 | 3,000 | 37.131 | 1.850 | 6.380 |
| 4 / C | 15,900 | 1,590 | 69.591 | 15.657 | 6.394 |

The planner takes the maximum measured cost for each phase over signs, cases,
betas and measured iterations, doubles it, and budgets at most 80% of eight
hours. With real sample count S and preparation count floor(S/10), the model is
`2*initialization + 10*(2*real_setup + 2*real_sample*S)
+ 20*(2*prep_setup + 2*prep_sample*floor(S/10)) <= 23040 seconds`.
The TSV timing columns already contain the factor of two. Sample counts are
rounded down. The remaining 1.6 hours covers unmeasured MPI/I/O and machine
variation; this is a local estimate, not an eight-hour cluster measurement.

## Measurements

`measurements.json` contains 16 completed one-rank, one-thread probes, full
commands and executable hashes. Twelve probes cover both signs, N=2/4 and all
three field cases at beta=5. Four additional FM A/C probes at beta=0.2 each
complete three self-consistent real-time iterations, measuring later field
laws. They use the normal convergence threshold of five; three end at the
requested real-time cap (exit 3), and one converges on iteration three (exit 0).
All preparations succeed. No failed preparation or crashed run enters the plan.

All probes use the production grid (50 imaginary/150 real intervals, Tmax=10),
one real-time substep, independent sampling, batch size one, 32 requested
blocks, time-translation invariance and seed 613. Each uses 256 samples per
stage per iteration. The maximum recorded peak RSS is
1.172 GB, below the 5 GB/rank cap.
RSS uses the largest completed-child high-water mark so far in each benchmark
invocation, conservatively bounding later probes. Production streams samples;
its sample count does not multiply the retained field/sample storage.

The current Release build was checked with `cmake --build Algorithm/build -j 4`.
Executable SHA-256: `2ae372d7f3d5454e603e1f6bcf2dcc73429040de22ed31f0e3bb6c1e11eb1eb7`.

Reproduce these probes from the repository root:

```sh
python3 CspinDMFT_Keldysh/benchmark_cubic_sweep.py --output CspinDMFT_Keldysh/Timing/cubic_beta_sweep/optimized_2026-11-12
python3 CspinDMFT_Keldysh/benchmark_cubic_sweep.py --output CspinDMFT_Keldysh/Timing/cubic_beta_sweep/optimized_2026-11-12 --feedback-iterations 3 --betas 0.2 --sign FM --cases zero hz0p5
```

Use a fresh output directory on a different machine. Measurements from older
executables or grids cannot enter the regenerated plan. Old pre-optimization
measurements remain alongside the default plan as historical records; this
calibration uses only this directory's measurements.

## Launcher verification

- All 132 final dry-run commands use the revised sample plan, unique concise
  names, ten real-time iterations and twenty preparation iterations. All 83
  existing sweep files match the resulting naming convention.
- Three completed solver smoke runs verify actual zero-field, uniform-field
  and staggered-field filenames, with eight samples per stage and temporary
  output; existing data were untouched.
- Submission stubs verify both 66-task arrays, 2G/5G per-rank requests, subset
  submission and Slurm spool paths. Exit-status stubs verify that status 3 is
  accepted while preparation failure status 2 remains an error.
- All twelve plan rows are ready, fit the 6.4-hour inflated budget and use even
  sample counts. The plan is reproducible with `--plan-only`.

No Slurm jobs were submitted. Short timing probes do not establish production
convergence or observable accuracy at the chosen counts.
