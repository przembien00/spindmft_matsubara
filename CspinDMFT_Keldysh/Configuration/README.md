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
