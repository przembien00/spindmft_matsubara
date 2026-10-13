# Staggered Strang sampling-basis correction

The new production files at beta=2.5, JL=-2, h_z=0.5 and Tmax=15 show
exploding correlations despite small covariance reconstruction residuals.
Their reported `by convergence` status is unreliable: the growing sampling
errors can reduce the standardized fixed-point residual even when the
absolute update is enormous. The original files were left intact.

## Cause and correction

The initial Strang implementation factorized the target covariance directly
in `(M,eta,nu)` with unit coordinate weights. Takagi factorization in that
basis preserves `E[z z^T]`, including all causal zeros and thermal mixed
blocks, but changes the auxiliary conjugated covariance compared with the
original branch-field norm. The original norm is

```
|V_+|^2 + |V_-|^2 = 2|eta|^2 + |nu|^2/2.
```

The corrected staged driver factors `(M,sqrt(2)*eta,nu/sqrt(2))`, then
inverse-scales the draw. Both pseudo-covariance and branch-equivalent
conjugated covariance are covered by tests. In conditional mode this
scaling is applied to the real residual factor; the frozen Matsubara factor
and the mixed conditional mean are unchanged. Preparation and Strang maps
are unchanged.

The real-time stopping rule also requires the absolute update to be below
`max(S,S*S)` (0.5 for spin 1/2). This blocks the observed false-convergence
mechanism while retaining the standardized residual test. This ceiling is
a divergence guard, not a precision guarantee.

## Matched executable comparison

Both executables used the same settings and seed, with independent
iteration streams from the normal staged driver:

```
--realTimeEvolution=staggered-strang --spinmodel=ISO --cstype=C
--JL=-2 --beta=2.5 --Bname=z --Babs=0.5 --dontrescaleB
--numImagTimeSteps=30 --numRealTimeSteps=100 --Tmax=15
--numSamplesPerCore=60000 --preparationSamplesPerCore=60000
--numBlocks=24 --iterlimit=12 --preparationIterlimit=20 --seed=734
```

One MPI rank was used. The total sample count matches the 6 x 10000
budget in the reported coarse-grid runs, though the rank streams differ.

| Real iteration | Original absolute update | Corrected absolute update |
|---|---:|---:|
| 1 | 0.7311 | 0.7450 |
| 2 | 0.6376 | 0.5270 |
| 3 | 11.5602 | 3.7000 |
| 4 | 312550.45 | 4.0245 |
| 5 | 8.6771e37 | 0.3608 |
| 6 | non-finite sample; aborted | 0.0811 |
| 10 | — | 0.0649443; converged |

The corrected final output has:

- maximum complex correlation magnitude: 0.25;
- maximum Re/Im correlation standard errors: 0.0189305941 / 0.0196471005;
- final standardized residual: 4.48532269;
- effective complex-weight sample count: 55973.0368 of 60000;
- average phase magnitude: 0.998823169;
- frozen MM residual: exactly zero;
- maximum closed-contour trace residual: 0.0600478503.

The transient remains noisy in both the old CFET solver and corrected
Strang solver. This single seeded comparison establishes the reported
sampling instability and verifies recovery for this case; it does not
establish second-order weak convergence, general long-time stability,
or reduced equilibrium magnetization drift.

Logs are in `unbalanced/run.log` and `balanced/run.log`. The corrected
HDF5 output is under `balanced/Data/`; `validation.json` records exact
final metrics. Earlier production files should be rerun with the rebuilt
executable. New files record `staggered_gaussian_basis` and
`staggered_gaussian_objective`, distinguishing them from the original
unbalanced Strang outputs.

## Regression checks

The default CTest suite includes joint/conditional reconstruction, the
branch-equivalent conjugated covariance, frozen conditional marginals,
unchanged preparation, and rejection of huge absolute updates with small
standardized residuals. The optional long-time regression reproduces the
settings above and takes a few minutes:

```
python3 Algorithm/Tests/test_strang_stability.py ./executable_DOUBLE.out
```

A saved pre-fix executable can be supplied with `--baseline` for a matched
before/after comparison. `--output-dir` preserves logs and outputs.
