# Seven/eight-spin feasibility, 2026-11-12

Full-grid, one-rank measurements of the current optimized solver. This assessment
does not add configurations to the production Slurm sweep.

## Measurements and conservative sample budgets

| N / symmetry | Measured peak GB/rank | Peak plus 32-block payload GB/rank | Real seconds/sample | Real samples/rank | Preparation samples/rank |
|---|---:|---:|---:|---:|---:|
| 7 / A | 0.750 | 0.799 | 14.68 | 74 | 6 |
| 7 / C | 2.282 | 2.476 | 28.37 | 36 | 4 |
| 8 / A | 1.526 | 1.578 | 110.05 | 8 | 4 |
| 8 / C | 3.037 | 3.246 | 211.10 | 0 | — |

All four measured cases fit the 5 GB/rank memory cap, including the statistical
payload allowance. Seven spins are now a reasonable pilot candidate. Eight
spins also fit memory, but their trajectory cost permits very little sampling:
eight-spin C requires about 7.60 inflated hours even for four samples per rank
per stage, including up to twenty preparation iterations. That fits eight hours
only with the additional 20% reserve removed. The existing N=2/4 sample plan is
unchanged; these low-count estimates have not been validated for accuracy.

The sample counts budget ten real-time plus up to twenty preparation iterations.
They double measured costs and use at most 6.4 hours of the eight-hour allocation,
matching the existing sweep's 20% wall-time reserve. Counts are per rank and per
iteration, without division by 16. Preparation uses one tenth of the real
samples, rounded down to an even count, with a minimum of four.
A zero real count means even four samples per stage exceed that budget. The
`summary.json` file also gives counts with the 20% reserve removed, retaining
the factor-of-two cost inflation. These counts are throughput estimates, not
validated recommendations for statistical accuracy.

Each completed probe uses FM J=-1/sqrt(6), beta=5, 50 imaginary/150 real intervals,
Tmax=10, one substep, independent sampling, batch size one, time-translation
invariance and seed 613. A has zero field; C has uniform field 0.5. Both stages
complete four samples with `--noselfcons`. This measures the initial prepared
field law, not converged production or later self-consistent laws. Commands,
termination statuses, hashes, phase timings and actual peak RSS are saved in
`measurements.json`; raw logs and HDF5 output are under each case directory.

## Memory accounting

Memory values are decimal GB. Actual process high-water RSS is measured by
macOS `/usr/bin/time -l` in a fresh process. Four samples permit two statistical
blocks despite the request for 32. The adjusted column adds the exact payload
of 30 more preparation and real-time blocks, using the configured representative
site pairs and selected spin directions. Both stages' block buffers persist
into real time. This is conservative extra statistical payload for the measured
field law, not a measured 32-block production peak. Larger sample counts stream
trajectories rather than retaining them all.

For the current symmetry blocks, covariance storage grows roughly quadratically
with the number of spins. On-demand construction removes the full covariance,
but the real lift and LAPACK workspace still overlap. The eigensolve allocation
for a complex block of size b is approximately `96*b*b` bytes: `32*b*b` for the
real lift and `64*b*b` for the leading LAPACK work array. Other live arrays,
stored factors and allocator retention are additional. For an eight-spin C
block of size 5616, those two arrays alone occupy about 3.03 GB.

Hilbert-space matrices still have dimension 2^N. Their storage grows as 4^N,
and dense matrix products grow roughly as 8^N. The covariance optimization
therefore does not remove the much larger per-trajectory computation at N=7/8.
The seven-spin configuration has six bath-active sites (the center has no
external bath coupling), giving largest measured C blocks of size 4212. All
eight sites are bath-active in the eight-spin configuration, giving size 5616.

## Limits

All production decisions still use the 5 GB/rank cap. Short no-feedback probes
cannot guarantee later-iteration or 16-rank node peaks, preparation convergence,
or observable accuracy. In particular, a handful of samples per rank does not
establish that the larger clusters give useful error bars. No cluster jobs
were submitted, and the existing N=2/4 sample plan remains unchanged.

The first sandboxed N=7 A timing wrapper could not collect memory because
macOS system-statistics access was blocked. That record has no RSS value and
is excluded; the successful fresh-process rerun supplies the table entry.

Reproduce with:

```sh
python3 CspinDMFT_Keldysh/Timing/large_cluster_feasibility_2026-11-12/benchmark.py
/Users/przembien/Projects/python_venv/bin/python CspinDMFT_Keldysh/Timing/large_cluster_feasibility_2026-11-12/summarize.py
```
