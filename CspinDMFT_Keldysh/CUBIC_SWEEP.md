# Cubic-lattice beta sweep

From the repository root on the cluster, with the usual solver executable built:

```sh
bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --submit
```

This creates the log directory before submission and submits four arrays of
66 jobs, one per cluster size. Each job uses 16 MPI ranks on one node in the
existing `med` partition, one thread per rank, **4 GB per core for every size**,
and an eight-hour wall-time limit. No jobs have been submitted during preparation.

The sweep contains **264 jobs**: N=2,4,7,8, both signs of J=1/sqrt(6), three
field cases, and beta=0.2,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5. The largest measured
peak plus the statistical-block allowance is about 3.25 GB/rank (N=8, C),
within the 4 GB planning cap. See the
[large-cluster feasibility report](Timing/large_cluster_feasibility_2026-11-12/REPORT.md)
for the raw memory and runtime measurements. N=7/8 receive small sample counts
to fit the time budget; their statistical accuracy has not been established.

Positive J is AFM. The selected configuration filenames round the coupling
to 0.40825, but the files contain the unrounded value. They are not modified.

| Case | Symmetry | FM field | AFM field |
|---|---|---|---|
| `zero` | A | uniform 0 | uniform 0 |
| `hz0p5` | C | uniform 0.5 | uniform 0.5 |
| `ordering` | C | uniform 0.005 | staggered 0.005 |

Both local and bath models are ISO. The real grid has 150 intervals up to
Tmax=10 (151 stored points); the imaginary grid has 50 intervals (51 points).
Evolution uses one staggered Strang step per real interval. Sampling is
independent, matching the existing solver default.

Each run permits ten real-time iterations and twenty preparation iterations,
with the default convergence threshold of 5. Preparation uses one tenth of
the production samples per rank, rounded down to an even count, with a minimum
of four. Its convergence is required before
real time starts. The budget is a timing estimate, not a preparation-accuracy
or physical-convergence validation.

Results are separated into `Data/Cubic_beta_sweep/{FM,AFM}/{zero,hz0p5,ordering}`.
Filenames match the existing data: model, configuration, beta, and nonzero
field only, for example
`spinmodel=ISO__config=Cubic_3D_N=2_NN_J=-0.40825__beta=1.hdf5` or
`spinmodel=ISO__config=Cubic_3D_N=4_NN_J=0.40825__beta=1__stagfield=0.005.hdf5`.
Grid and sample counts remain in the HDF5 parameters. Reaching the ten-iteration
real-time cap saves under `NOT_CONVERGED` and is treated as a completed batch
job. Preparation failure remains an error and prevents real-time execution.

## Timings and sample counts

The selected counts are in `Timing/cubic_beta_sweep/sample_plan.tsv`:

| Sites / symmetry | Real samples per rank | Preparation samples per rank |
|---|---:|---:|
| 2 / A | 1,066,900 | 106,690 |
| 2 / C | 454,400 | 45,440 |
| 4 / A | 30,000 | 3,000 |
| 4 / C | 15,900 | 1,590 |
| 7 / A | 74 | 6 |
| 7 / C | 36 | 4 |
| 8 / A | 8 | 4 |
| 8 / C | 4 | 4 |

Both signs share these counts; C covers both nonzero-field cases. The complete
24-row plan and the 20 measurements used to generate it are saved in
`Timing/cubic_beta_sweep/expanded_2026-11-12/`. Regenerate that plan without
running new probes with:

```sh
python3 CspinDMFT_Keldysh/benchmark_cubic_sweep.py --plan-only \
    --output CspinDMFT_Keldysh/Timing/cubic_beta_sweep/expanded_2026-11-12
```

Detailed N=2/4 measurements and validation are in
[`optimized_2026-11-12/REPORT.md`](Timing/cubic_beta_sweep/optimized_2026-11-12/REPORT.md).
Fresh optimized-solver probes, commands, executable hash and raw logs are in
`Timing/cubic_beta_sweep/optimized_2026-11-12/`. These use the production grid
(50 imaginary and 150 real intervals, Tmax=10), beta=5, one substep, 32 blocks,
independent sampling with batch size one and time-translation invariance.
Each probe completes one preparation and one real-time iteration with 256
samples per stage, seed 613, and `--noselfcons`. Both coupling signs and all
three field cases are measured. Four additional FM probes at beta=0.2 (N=2,4,
A/C) retain feedback for up to three real-time iterations, with the normal
convergence threshold and preparation limit. The planner also includes their
maximum per-iteration costs, so later field laws inform the estimate.
Production retains self consistency and native initialization. These probes
estimate cost; they do not establish preparation convergence or accuracy at the
chosen production counts. Four additional N=7/8 probes from the feasibility
report each complete one preparation and one real-time iteration at beta=5,
with four samples per stage and two actual statistical blocks. Their memory
bounds include an allowance for 32 blocks. These FM zero-field/uniform-field
measurements supply the large-cluster A/C timing estimates for both signs
and field cases.

For each size and symmetry, the planner takes the maximum measured sampling,
covariance setup and initialization costs over the available signs, temperatures
and field cases. It doubles those costs and reserves a further 20% of the
eight-hour allocation for MPI, storage and timing variation. It solves the
remaining budget jointly for **ten real-time plus up to twenty preparation
iterations**, including setup on every iteration. Counts are rounded down to
valid even counts, and preparation retains at least four samples.

The inflated estimate for both stages is normally at most 6.4 hours, leaving
1.6 hours of wall-time headroom. **N=8 C is the exception:** its minimum four
samples per rank per stage require about **7.60 hours** in the doubled-cost
model. This case uses part of the extra reserve while remaining within eight
hours, and the launcher prints a notice. Cases whose minimum count exceeds
eight hours remain rejected. The `estimated_seconds` column includes both stages.

**Counts are per rank and are not divided by 16.** Ranks sample in parallel;
16 ranks provide 16 times the total samples per iteration. Covariance setup
is repeated independently on every rank.

Local timing cannot guarantee cluster performance. Slurm enforces the
eight-hour cap. Recalibration on a cluster compute node is available without
third-party Python packages:

```sh
python3 CspinDMFT_Keldysh/benchmark_cubic_sweep.py --output /path/to/fresh_timings
SAMPLE_PLAN=/path/to/fresh_timings/sample_plan.tsv \
    bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --submit
```

The benchmark runs directly as one MPI rank with one thread, saves each
completed probe, and skips identical completed probes when restarted. Use a
fresh output directory when changing machines. Only measurements from the
current executable, grid and sampling settings enter the plan; old pre-optimization
measurements are excluded. `--betas`, `--samples`, `--sites`, `--cases`, and
`--sign` allow additional probes. `--feedback-iterations 3` measures later
self-consistent iterations; exit status 3 is accepted for that iteration cap.
Probe defaults are 256 samples per stage for N=2/4 and four for N=7/8; an
explicit `--samples` overrides these and must be even and at least four.
`--plan-only` regenerates counts from saved measurements. `--hours` defaults to eight; if changed, adjust the Slurm wall-time
request to match. Memory records are the maximum completed-child peak seen so
far, so they conservatively bound the current probe. Small-sample probes add
an explicit statistical-block memory allowance for production. Plans exceeding
4 GB per rank or the full timing allowance are rejected by the launcher.

## Inspection and subsets

```sh
bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --dry-run 0
# All twenty-four field/config cases at beta=0.2:
bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --submit --array=0-23
```

Jobs are ordered by beta, then by the twenty-four rows of the sample plan.
Task 0 uses beta=0.2 and plan row 1; task 24 uses beta=0.5 and plan row 1;
task 263 uses beta=5 and plan row 24. An explicit array selection submits one
array with the same 4 GB-per-core default. Slurm options may be passed after
`--submit`. `KELDYSH_EXECUTABLE` optionally selects another cluster-built
executable; its performance should be checked before reusing this plan.
