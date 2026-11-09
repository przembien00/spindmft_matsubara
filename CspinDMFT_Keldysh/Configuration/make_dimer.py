#!/usr/bin/env python3
"""Write a two-spin cluster with independent on-site replica baths.

This small benchmark embedding is not a square-lattice configuration. For a
lattice embedding, pass a configuration from CspinDMFT_finite_T instead.
"""
import argparse
from pathlib import Path
import h5py
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output', type=Path)
parser.add_argument('--J', type=float, default=0.7, help='intra-cluster coupling')
parser.add_argument('--bath', type=float, default=0.2, help='bath coupling scale (squared in covariance)')
parser.add_argument('--mean', type=float, default=0.0, help='on-site first-moment feedback weight')
p = parser.parse_args()
if not np.all(np.isfinite([p.J,p.bath,p.mean])):
    parser.error('couplings must be finite')
p.output.parent.mkdir(parents=True, exist_ok=True)
with h5py.File(p.output, 'x') as f:
    f.attrs['correlation_weight_convention_version'] = 2
    f.attrs['num_Spins'] = 2
    f.attrs['spin lengths'] = np.array([0.5,0.5])
    f.attrs['num_Categories'] = 3
    f['spin-spin couplings'] = [[0.,p.J],[p.J,0.]]
    f['expectation weights'] = p.mean*np.eye(2)
    weights = np.zeros((2,2,2,2))
    weights[0,0,0,0] = weights[1,1,1,1] = p.bath**2
    f['correlation weights'] = weights
    f['correlation_categories'] = np.array([[0,0],[0,1],[1,1]],dtype=np.uint32)
print(p.output)
