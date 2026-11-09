#!/usr/bin/env python3
"""Add unit-second-moment copies of the supplied cubic-lattice configurations."""

from pathlib import Path
import shutil

import h5py
import numpy as np


def main():
    config_dir = Path(__file__).resolve().parent.parent / "Configuration_Data" / "CubicLattice"
    coupling = 1.0 / np.sqrt(6.0)
    scale = coupling / 0.5
    for sites in (2, 4, 7, 8):
        for sign in (-1, 1):
            source = config_dir / f"Cubic_3D_N={sites}_NN_J={sign * 0.5}.hdf5"
            target = config_dir / f"Cubic_3D_N={sites}_NN_J={sign * coupling:.5f}.hdf5"
            # Preserve the geometry, categories and convention metadata exactly.
            # Exclusive creation prevents overwriting an existing configuration.
            with source.open("rb") as src, target.open("xb") as dst:
                shutil.copyfileobj(src, dst)
            with h5py.File(target, "r+") as config:
                for name in ("spin-spin couplings", "expectation weights"):
                    config[name][...] = config[name][...] * scale
                config["correlation weights"][...] = config["correlation weights"][...] * scale**2
                local = np.square(config["spin-spin couplings"][...]).sum(axis=1)
                bath = np.einsum("iikk->i", config["correlation weights"][...])
                np.testing.assert_allclose(local + bath, 1.0, rtol=0, atol=1e-14)
            print(target)


if __name__ == "__main__":
    main()
