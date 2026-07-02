# Example run commands

pCN / Metropolis algorithms, `executable_DOUBLE.out`

## Notes

- Build first (per-algorithm):

  ```bash
  cd <algorithm>/Algorithm
  cp CMakeLists.txt_ CMakeLists.txt      # only if CMakeLists.txt is missing
  mkdir -p build && cd build && cmake .. && make
  ```

  This produces `<algorithm>/executable_DOUBLE.out` (no "MH" suffix anymore;
  the pCN sampler is now the default algorithm).

- **Run from the algorithm's top-level folder** (e.g. `CspinDMFT_finite_T/` or
  `spinDMFT/`), NOT from the repo root. The executable resolves
  `Configuration_Data/` and `Data/` paths relative to the current working
  directory.

- Output HDF5 is written under `Data/<project>/`; non-converged runs go to
  `Data/<project>/NOT_CONVERGED/`. Inspect with `quickreader.py` / `h5py`.

- `-n <N>` sets the number of MPI cores; total samples = `N * numSamplesPerCore`.

## pCN sampler parameters (both algorithms)

| Option | Meaning |
|--------|---------|
| `--mhStepSize=<beta>` | pCN step `beta` in (0,1]. `V' = sqrt(1-beta^2) V + beta*xi`; `beta=1` -> independence sampler; `beta->0` -> tiny moves. |
| `--mhBurnIn=<N>` | burn-in steps per core in the first self-consistency iteration. |
| `--mhWarmBurnIn=<frac>` | burn-in fraction (0..1) reused from iteration 2 onward (warm-started from the previous iteration's final state). |

## Error-bar estimation (both algorithms)

| Option | Meaning |
|--------|---------|
| `--errmethod=<method>` | standard-error estimator for the correlations / spin expectations. `blocking`: batch means pooled over blocks and cores; robust to autocorrelation (default). `ar1`: legacy acceptance-based single-mode factor; under-/over-estimates the error for slow or very-high-acceptance chains. Anything other than `ar1` selects blocking. |
| `--numBlocks=<N>` | number of batch-mean blocks per core for the blocking estimator (default 32). Each block must be much longer than the chain autocorrelation time; the `blocking_curve` diagnostic in the output reports whether that holds. With blocking, the output also stores the integrated autocorrelation time `tau_int` (in pCN steps) per point. |

## CspinDMFT_finite_T (finite-T cluster, Matsubara)

Run from `CspinDMFT_finite_T/`. Needs a cluster config in
`Configuration_Data/<lattice>/` (here `Square_2D_N=2_NN_J=0.5`).

```bash
mpirun -n 4 executable_DOUBLE.out \
  --config=Square_2D_N=2_NN_J=0.5 \
  --spinspinmodel=ISO \
  --stype=A \
  --beta=1 \
  --numTimeSteps=99 \
  --numSamplesPerCore=10000 \
  --mhStepSize=0.3 \
  --mhBurnIn=500 \
  --iterlimit=20 \
  --numBlocks=32 \
  --project=SquareLattice --dstproject=test \
  --fileext=test 


```

Common variants:

```text
--uncoupled_spins        simulate as independent single spins (per-site Z_i)
--chemshift=0.5,-0.5     per-site magnetic field in z (or --unifield / --stagfield)
--stype=B|C|D            correlation symmetry type (A is most constrained, D least)
--errmethod=ar1          legacy error bars instead of blocking (default)
--adaptive --staterrtolerance=0.0025   grow sample size to hit a target error
```

See all options: `./executable_DOUBLE.out --help`

## spinDMFT (single-site spinDMFT)

Run from `spinDMFT/`. Single-site: no lattice `--config` needed.

```bash
mpirun -n 4 executable_DOUBLE.out \
  --spinmodel=ISO \
  --cstype=A \
  --beta=5 \
  --numTimeSteps=99 \
  --numSamplesPerCore=10000 \
  --mhStepSize=0.9 \
  --mhBurnIn=500 \
  --iterlimit=20 \
  --fileext=testing \
  --numBlocks=20 \
  --critneg=100 \
  --loadinit \
  --initcorrfile=abc
  
```

Common variants:

```text
--Bname=<field> --Babs=<h> --Btheta=<deg> --Bphi=<deg>   external field
--spin=1   (spin-1 instead of spin-1/2) --project
--cstype=A|B|C|D                                         correlation symmetry
```

See all options: `./executable_DOUBLE.out --help` (and `--helpNum` for numerics)
