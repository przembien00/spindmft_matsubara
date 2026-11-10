# Normalized cubic configurations

`Configuration_Data/CubicLattice/` includes both signs of
`J=0.40825` for clusters with `N=2,4,7,8`. The filenames round the coupling
to five decimal places; stored values retain full precision at
`|J|=1/sqrt(6)`, so the six nearest neighbors give `sum_i J_i0^2=1`.
For example, from `CspinDMFT_Keldysh/`, select the AFM dimer with:

```sh
./executable_DOUBLE.out --project=CubicLattice --config=Cubic_3D_N=2_NN_J=0.40825
```

The original `J=±0.5` configurations remain available. The normalized copies
preserve their cluster geometry and correlation categories, scale internal
couplings and expectation weights by `2/sqrt(6)`, and scale correlation
weights by `2/3`. The normalization includes both internal and bath neighbors
and holds separately at every cluster site.

To generate the normalized files from the originals, run
`python3 Configuration/make_normalized_cubic.py`. This requires NumPy and h5py
and refuses to overwrite existing target files.

## Correlation weight convention

Correlation weights are stored per source-pair orientation. For an off-diagonal
representative whose reverse is absent from `correlation_categories`, the solver
counts the stored weight twice. Diagonal pairs and explicitly listed reverse
pairs are counted once. Configuration files identify this as
`correlation_weight_convention_version=2`.

External legacy configurations must halve off-diagonal weights when the reverse
pair is absent, then set that version attribute to 2. Do not modify diagonal
weights, explicitly ordered pairs, or expectation weights.
