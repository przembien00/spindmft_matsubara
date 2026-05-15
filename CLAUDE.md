# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a C++17/C++20 computational physics framework implementing **Spin Dynamical Mean-Field Theory (spinDMFT)** and variants. It simulates spin dynamics in quantum/statistical physics using self-consistent Monte-Carlo iteration with MPI parallelization and HDF5 data storage.

## Building

Each algorithm lives in its own subdirectory (`spinDMFT/`, `CspinDMFT/`, `CspinDMFT_finite_T/`, `nl-spinDMFT/`, `ExactDiagonalization/`). The build template must be copied before use:

```bash
cd <algorithm>/Algorithm
cp CMakeLists.txt_ CMakeLists.txt   # template copy required
mkdir build && cd build
cmake ..
make
```

Executables are named `executable_DOUBLE.out` or `executable_FLOAT.out` depending on the precision set in CMakeLists.txt.

**Dependencies:** MPI 4.1.6, HDF5 1.10.10, LAPACK 3.12.0, Boost 1.83, Blaze 3.8 (header-only matrix library). Eigen is optional (alternative diagonalization backend).

## Running

```bash
mpirun -n <num_cores> executable_DOUBLE.out \
  --numSamplesPerCore=2500 --numTimeSteps=100 --dt=0.02 \
  --spinmodel=ISO --cstype=A
```

Parameters are parsed at runtime via Boost program options. Output is written as HDF5 files to `Data/`. Use `quickreader.py` to inspect HDF5 output.

## Architecture

### Shared Libraries (`cpp_libs/`)

All algorithms share these libraries — changes here affect every algorithm:

| Library | Role |
|---------|------|
| `Physics/` | Spin models (ISO, DRF), magnetic fields, lattice geometries, coupling constants |
| `Observables/` | Spin correlations, tensors, cluster structure definitions |
| `Multivariate_Gaussian/` | Covariance matrix construction, Gaussian noise sampling, symmetry blocks |
| `Frequency_Multivariate_Gaussian/` | Frequency-domain MVG for Matsubara/frequency formalisms |
| `Matrices/` | LAPACK/Eigen diagonalization wrappers |
| `HDF5/` | HDF5 read/write abstraction |
| `Globals/` | Shared type definitions (`Types.h`, `MPI_Types.h`) |

### Per-Algorithm Components

Each algorithm's `Algorithm/` directory contains:

- **`Parameter_Space/`** — CLI argument parsing, physics parameter storage, simulation configuration
- **`Run_Time_Data/`** — Statistics accumulation, eigenvalue tracking across iterations
- **`Functions/`** — Core computational kernels: time propagators, spin correlation integrals, convergence checks
- **`Storage_Concept/`** — HDF5 I/O for simulation results and metadata
- **`Time_Measure/`** — Per-iteration performance profiling

`CspinDMFT` additionally has `Mean_Field_Models/` for physics-specific mean-field equations parameterized by cluster geometry.

### Self-Consistent Loop (main.cpp)

The core algorithm follows this pattern in every variant:

1. MPI init → parse parameters → initialize operators/correlations/RNG
2. **Outer loop** (self-consistency iterations):
   - Compute mean-field from current correlations
   - Diagonalize covariance matrix (LAPACK/Eigen)
   - Sample Gaussian noise
   - **Inner loop** (Monte-Carlo sampling, per MPI core): time-evolve spins, accumulate correlation functions
   - MPI reduce → check convergence → write checkpoint to HDF5

### Algorithm Variants

| Directory | Description |
|-----------|-------------|
| `spinDMFT/` | Main single-site spinDMFT algorithm |
| `CspinDMFT/` | Cluster spinDMFT (finite clusters, multiple correlation symmetry types A/B/C/D) |
| `CspinDMFT_finite_T/` | Finite-temperature cluster variant (Matsubara formalism) |
| `nl-spinDMFT/` | Non-local spinDMFT (environment + embedded cluster) |
| `ExactDiagonalization/` | Full Hilbert space reference solver for benchmarking |

### Configuration Files

`CspinDMFT/Configuration/` and `CspinDMFT_finite_T/Configuration/` contain Python scripts (`config_generator.py`, `Square2D_NSpin.py`) that generate cluster geometry configs used at runtime.

## Data Format

All output is HDF5. Files are named by physics parameters, e.g.:
`spinmodel=ISO__config=Square_2D_N=2_NN_J=0.5__beta=1.5.hdf5`

## Plotting

Use the matplotlibrc style provided in the directory you're working in. Don't add titles to plots. Use exactly the range of x values on the x axis. When plotting correlations, normalize tau values by beta so they lie in [0,1]. Save plots in dedicated Plots/ folders. Use the name $g^{ab}_{ij}(\tau)$ for correlations, where a and b are the spin indices and i, j the site indices. Use already available scripts for plotting if they are available, otherwise create new ones.