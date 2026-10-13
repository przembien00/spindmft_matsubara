# Mixed-stage instability audit

These controls use one MPI rank, seed 17 (with a separate seed-29 mixing control),
beta=2.5, ISO/C, h_z=0.5, JL=-2, 30 imaginary-time steps, and staggered Strang.
The long-time cases use 100 real-time steps over T=15 (dt=0.15),
10,000 real-time samples and 50,000 preparation samples per rank.
The reported user run has six ranks; these tests isolate behavior on the same
grids and physics but do not reproduce its total sampling budget.
Every case has its launch command, console log and HDF5 output in its own folder.

## Observations

- With Matsubara pCN and the default mixing alpha=1, mixed iteration 3 reaches
  max |g|=1069.23 with max quadrature-combined Re/Im component standard error
  1222.05. The sampled phase magnitude remains 0.9643. The frozen-MM residual is
  zero. The first mixed iteration has max |g|=0.25 and maximum error 0.01244.
- Independent conditional sampling also reaches max |g|=2152.48 by iteration 3,
  while its phase magnitude is 0.9635. The effect is not specific to pCN.
- The joint-Takagi independent control also becomes unstable at the same horizon
  and seed: max |g|=122.48 and maximum error 121.62. The conditional gauge is not
  the only possible source of growth.
- Lowering the conditional Matsubara gauge scale to 0.3 improves the sampled
  phase to 0.9994 but worsens mixed-stage growth to 1.07e17. Improving the
  imaginary-branch phase is insufficient to control real-time estimates.
- With mixingAlpha=0.5, the seed-17 run terminates by the existing convergence
  criterion after 10 mixed iterations with max |g|=0.25. The largest late-time
  correlation error is still 0.09632. This is absence of explosive growth in this
  test, not an accurate or seed-independent long-time solution.
- A second seed (29) with mixingAlpha=0.5 also terminates without explosive
  growth after 12 mixed iterations: max |g|=0.26730,
  maximum correlation error 0.09537, and sampled phase
  0.96281. These two seeds support this setting as a
  practical stabilization control, not a production precision guarantee.
- With T=3 and 20 real-time steps (the same dt=0.15), default mixing alpha=1
  terminates after 5 mixed iterations with max |g|=0.25 and maximum error 0.00610.
- mixingAlpha=0.2 slows preparation enough that it does not finish within the
  configured 30 preparation iterations. It therefore supplies no mixed-stage
  stabilization evidence.

## Interpretation and next run

The controls are consistent with real-time sampling fluctuations being amplified
by undamped feedback of the whole mixed-time primitive and magnetization into
the next Gaussian covariance. The healthy phase magnitude argues against division
by a near-zero phase denominator as the primary cause. Errors comparable to the
huge estimates show that the late-time values are unresolved Monte Carlo outputs.
This does not establish a unique microscopic explanation or exclude propagation
bias; no long-time physical accuracy claim is made.

Matsubara-only pCN weights |Z_M|. It does not importance-sample the direct real
residuals or the real-time observable numerators. Consequently good pCN acceptance
and a large average phase do not bound their variance.

For a controlled retry, add --mixingAlpha=0.5 to the user command. For an early-time
validation, use --Tmax=3 --numRealTimeSteps=20. Inspect late-time uncertainties and
repeat seeds before interpreting a termination message as a precise solution.
No sampler or normalization code was changed during this audit.

The machine-readable results are in audit_results.json. Recreate the summary with:

    python3 -B audit_instability.py
