# Expanded cubic sweep, 2026-11-12

The production plan now includes N=2,4,7,8, both coupling signs and all three
field cases over eleven betas: 264 jobs, submitted as four 66-job arrays.
Every job requests 16 MPI ranks, one core per rank, 4G per core and eight hours.
N=2/4 rows are exactly unchanged from the prior optimized sample plan.

| Sites / symmetry | Real samples/rank | Preparation samples/rank | Inflated total hours |
|---|---:|---:|---:|
| 2 / A | 1,066,900 | 106,690 | 6.400 |
| 2 / C | 454,400 | 45,440 | 6.399 |
| 4 / A | 30,000 | 3,000 | 6.380 |
| 4 / C | 15,900 | 1,590 | 6.394 |
| 7 / A | 74 | 6 | 6.267 |
| 7 / C | 36 | 4 | 6.228 |
| 8 / A | 8 | 4 | 6.058 |
| 8 / C | 4 | 4 | 7.596 |

Counts are per rank and per iteration. Estimates double the measured costs and
include ten real-time plus up to twenty preparation iterations. Preparation
receives one tenth of the real count rounded down to an even count, with a
minimum of four. All counts admit at least two equal statistical blocks.
N=8 C needs 7.596 hours at its minimum four samples; it uses the additional
20% reserve, and the launcher prints a notice. All other cases fit 6.4 hours.
The memory planning cap is 4,000,000,000 bytes per rank; the largest corrected
measurement is 3,245,772,992 bytes. These are local timing estimates and do not
validate convergence or accuracy at the small N=7/8 counts.

`measurements.json` consolidates the 16 N=2/4 optimized measurements from
`../optimized_2026-11-12/measurements.json` and the four successful N=7/8
feasibility measurements from `../../large_cluster_feasibility_2026-11-12/`.
Each imported large-cluster row records its original command, executable hash
and measurement source. Its production RSS bound includes the exact extra
payload for 30 preparation and real-time statistical blocks. The failed initial
sandbox memory wrapper is excluded. No simulations were rerun for this expansion.

The planner pools signs and field amplitudes within each size and A/C symmetry.
The N=7/8 FM zero/uniform-field measurements therefore supply estimates for
both signs and both nonzero-field cases. Native self consistency remains enabled
in production.

Reproduce the complete sample plan without running simulations:

```sh
python3 CspinDMFT_Keldysh/benchmark_cubic_sweep.py --plan-only --output CspinDMFT_Keldysh/Timing/cubic_beta_sweep/expanded_2026-11-12
```

The default launcher uses the identical copy in
`Timing/cubic_beta_sweep/sample_plan.tsv`. A fresh calibration uses 256 samples
per stage for N=2/4 and four for N=7/8 unless `--samples` explicitly overrides
that choice, avoiding an implicit 256-sample eight-spin probe.

All 264 dry runs were checked for unique concise filenames, correct fields,
counts, grids and timing bounds. Submission stubs verified four disjoint
66-task arrays at 4G/core, explicit subsets and Slurm spool paths. Solver-status
stubs verified real-time cap handling and preservation of preparation errors,
including a runnable N=8 C row. Plan-only regeneration is deterministic. Details
are recorded in `validation.json`. No Slurm jobs were submitted.
