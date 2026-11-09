# Cubic-lattice beta sweep

From the repository root on the cluster, with the usual solver executable built:

```sh
bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --submit
```

This creates the log directory before submission and submits two arrays of
66 jobs. Each job uses 16 MPI ranks on one node in the existing `med` partition,
one thread per rank, and a six-hour wall-time limit. N=2 requests 2 GB per rank;
N=4 requests 5 GB per rank. No jobs have been submitted during preparation.

The sweep contains **132 jobs**: N=2,4, both signs of J=1/sqrt(6), three field
cases, and beta=0.2,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5. N=7 and N=8 are excluded
as requested. The seven-spin C factorization alone requires at least 6.9 GB
per rank in simultaneously allocated dense arrays, exceeding the 5 GB cap.

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
with the default convergence threshold of 5. Preparation uses at most one
tenth of the production samples per rank. Its convergence is required before
real time starts. The budget is a timing estimate, not a preparation-accuracy
or physical-convergence validation.

Results are separated into `Data/Cubic_beta_sweep/{FM,AFM}/{zero,hz0p5,ordering}`.
Filenames include the config, beta, and sample counts. Reaching the ten-iteration
real-time cap saves under `NOT_CONVERGED` and is treated as a completed batch
job. Preparation failure remains an error and prevents real-time execution.

## Timings and sample counts

The selected counts are in `Timing/cubic_beta_sweep/sample_plan.tsv`; raw
one-core measurements and logs are alongside it. The existing counts are
retained from the original 80-imaginary/200-real-interval probes, so the
timing estimates refer to that larger grid. New probes use the current
50-imaginary/150-real-interval grid with four samples per stage, at beta=5,
for N=2 and N=4 with A and C symmetry. Grid sizes are part of new probe
identities, preventing reuse of the old probes as new measurements.
Only the probes disable self-consistent feedback.
Production uses the normal self-consistent solver and native initialization.

For each size and symmetry, the planner uses the maximum measured costs,
pooling any available signs, temperatures and C-field amplitudes. The
`timing_source` column identifies the probes used. Sampling, covariance setup,
and initialization costs are inflated by a factor of two. A further 20% of
the six hours is reserved for MPI, storage and timing variation. Of the
remaining budget, 80% is allocated to ten real-time iterations and at most
20% to preparation; the sample-count cap makes preparation substantially
cheaper in the measured cases. Counts are rounded down.

**Counts are per rank and are not divided by 16.** Ranks sample in parallel;
16 ranks provide 16 times the total samples per iteration. Covariance setup
is repeated independently on every rank.

Local timing cannot guarantee cluster performance. Slurm enforces the
six-hour cap. Recalibration on a cluster compute node is available without
third-party Python packages:

```sh
python3 CspinDMFT_Keldysh/benchmark_cubic_sweep.py --output /path/to/fresh_timings
SAMPLE_PLAN=/path/to/fresh_timings/sample_plan.tsv \
    bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --submit
```

The benchmark runs directly as one MPI rank with one thread, saves each
completed probe, and skips identical completed probes when restarted. The
executable hash is part of each probe identity. Use a fresh output directory
when changing machines. `--betas`, `--samples`, `--sites`, `--cases`, and
`--sign` allow additional probes. `--plan-only` regenerates counts from the
saved measurements. Memory records are the maximum completed-child peak
seen so far, so they conservatively bound the current probe. Plans exceeding
5 GB per rank or the timing allowance are rejected by the launcher.

## Inspection and subsets

```sh
bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --dry-run 0
# All twelve field/config cases at beta=0.2:
bash CspinDMFT_Keldysh/slurm_cubic_beta_sweep.sh --submit --array=0-11
```

Jobs are ordered by beta, then by the twelve rows of the sample plan. Task
0 uses beta=0.2 and plan row 1; task 12 uses beta=0.5 and plan row 1; task
131 uses beta=5 and plan row 12. An explicit array selection submits one
array with the 5 GB-per-rank default. Slurm options may be passed after
`--submit`. `KELDYSH_EXECUTABLE` optionally selects another cluster-built
executable; its performance should be checked before reusing this plan.
