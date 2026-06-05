from __future__ import annotations

"""Generate cubic-lattice nearest-neighbor cluster configurations for 3D.

Two cluster geometries are supported:

- **7-spin star**: one central site at (0,0,0) and its 6 nearest neighbors along
  the ±x, ±y, ±z axes.  The cluster does not tile the cubic lattice periodically,
  so external neighbors are mapped to the nearest cluster site by Euclidean distance
  (ties broken by lowest site index).

- **8-spin cube**: 2×2×2 sites at all corners of a unit cube.  This tiles the
  cubic lattice with period 2 in every direction, so the mapping is periodic
  modulo 2 in each coordinate.
"""

import numpy as np
from typing import Callable

try:
    from config_generator import create_config_file, generate_nn_cluster_weights
except ModuleNotFoundError:
    from .config_generator import create_config_file, generate_nn_cluster_weights


LatticePoint3D = tuple[int, int, int]

NN_DISPLACEMENTS_3D: list[LatticePoint3D] = [
    (1, 0, 0), (-1, 0, 0),
    (0, 1, 0), (0, -1, 0),
    (0, 0, 1), (0, 0, -1),
]


def nearest_neighbor_couplings_3d(
    cluster_positions: list[LatticePoint3D],
    coupling: float,
) -> np.ndarray:
    """Build the NN coupling matrix for a 3D cubic-lattice cluster."""

    num_spins = len(cluster_positions)
    position_to_index = {pos: idx for idx, pos in enumerate(cluster_positions)}
    J = np.zeros((num_spins, num_spins), dtype=np.float64)

    for idx, (x, y, z) in enumerate(cluster_positions):
        for dx, dy, dz in NN_DISPLACEMENTS_3D:
            nb = (x + dx, y + dy, z + dz)
            nb_idx = position_to_index.get(nb)
            if nb_idx is not None and nb_idx > idx:
                J[idx, nb_idx] = coupling
                J[nb_idx, idx] = coupling

    return J


# ---------------------------------------------------------------------------
# 7-spin star cluster
# ---------------------------------------------------------------------------

def star_cluster_positions() -> list[LatticePoint3D]:
    """Return the 7 sites of the star cluster: center + 6 NN arms."""

    return [(0, 0, 0)] + list(NN_DISPLACEMENTS_3D)


def star_cluster_map(
    cluster_positions: list[LatticePoint3D],
) -> Callable[[LatticePoint3D], int]:
    """Map any lattice point to the nearest star-cluster site (squared distance).

    Ties are broken by lowest site index.  For on-axis external sites (e.g.
    (2,0,0)) this maps unambiguously to the same-arm site.  For diagonal sites
    (e.g. (1,1,0)) the nearest arm wins.
    """

    def map_to_cluster_index(point: LatticePoint3D) -> int:
        dists = [sum((a - b) ** 2 for a, b in zip(point, pos)) for pos in cluster_positions]
        min_dist = min(dists)
        return next(i for i, d in enumerate(dists) if d == min_dist)

    return map_to_cluster_index


def write_star_config(coupling: float, label: str) -> None:
    """Write a 7-spin star configuration to HDF5."""

    project_name = "CubicLattice"
    config_file = f"Cubic_3D_N=7_star_NN_{label}"

    cluster_positions = star_cluster_positions()
    J = nearest_neighbor_couplings_3d(cluster_positions, coupling)
    map_fn = star_cluster_map(cluster_positions)

    base_weights = generate_nn_cluster_weights(
        cluster_positions=cluster_positions,
        nn_displacements=NN_DISPLACEMENTS_3D,
        map_to_cluster_index=map_fn,
    )

    eff = 0.5 if coupling == 0.0 else coupling
    mf_weights = eff * base_weights.mf_expectation_weights
    corr_weights = (eff * eff) * base_weights.correlation_weights

    output_path = create_config_file(
        project_name=project_name,
        config_file=config_file,
        J=J,
        mf_expectation_weights=mf_weights,
        correlation_weights=corr_weights,
        correlation_categories=base_weights.canonical_pair_categories,
    )
    print(f"wrote {output_path}")


# ---------------------------------------------------------------------------
# 8-spin cube cluster (2×2×2)
# ---------------------------------------------------------------------------

def cube_cluster_positions() -> list[LatticePoint3D]:
    """Return the 8 corner sites of the 2×2×2 cube cluster."""

    return [(x, y, z) for z in range(2) for y in range(2) for x in range(2)]


def cube_cluster_map() -> Callable[[LatticePoint3D], int]:
    """Map any lattice point into the 2×2×2 cube via periodic images (mod 2)."""

    position_to_index = {pos: idx for idx, pos in enumerate(cube_cluster_positions())}

    def map_to_cluster_index(point: LatticePoint3D) -> int:
        return position_to_index[(point[0] % 2, point[1] % 2, point[2] % 2)]

    return map_to_cluster_index


def write_cube_config(coupling: float, label: str) -> None:
    """Write an 8-spin cube configuration to HDF5."""

    project_name = "CubicLattice"
    config_file = f"Cubic_3D_N=8_cube_NN_{label}"

    cluster_positions = cube_cluster_positions()
    J = nearest_neighbor_couplings_3d(cluster_positions, coupling)
    map_fn = cube_cluster_map()

    base_weights = generate_nn_cluster_weights(
        cluster_positions=cluster_positions,
        nn_displacements=NN_DISPLACEMENTS_3D,
        map_to_cluster_index=map_fn,
    )

    eff = 0.5 if coupling == 0.0 else coupling
    mf_weights = eff * base_weights.mf_expectation_weights
    corr_weights = (eff * eff) * base_weights.correlation_weights

    output_path = create_config_file(
        project_name=project_name,
        config_file=config_file,
        J=J,
        mf_expectation_weights=mf_weights,
        correlation_weights=corr_weights,
        correlation_categories=base_weights.canonical_pair_categories,
    )
    print(f"wrote {output_path}")


def main() -> None:
    """Write FM and AFM configs for both 3D cluster geometries."""

    # 7-spin star: AFM (J=0.5) and FM (J=-0.5)
    write_star_config(coupling=0.5, label="J=0.5")
    write_star_config(coupling=-0.5, label="J=-0.5")

    # 8-spin cube: AFM (J=0.5) and FM (J=-0.5)
    write_cube_config(coupling=0.5, label="J=0.5")
    write_cube_config(coupling=-0.5, label="J=-0.5")


if __name__ == "__main__":
    main()
