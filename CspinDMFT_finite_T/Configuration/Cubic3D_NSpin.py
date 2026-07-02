from __future__ import annotations

"""Generate cubic-lattice nearest-neighbor cluster configurations for 3D.

Three cluster geometries are supported:

- **2-spin bond**: two sites at (0,0,0) and (1,0,0), one NN bond.  External
  neighbors are mapped by the bipartite sublattice parity ``(x+y+z) mod 2``,
  so every external neighbor of a cluster site is represented by the *other*
  cluster site (the only one on the opposite sublattice).  Each cubic site
  has 6 NN; one is the in-cluster bond partner, so each cluster spin sees 5
  out-of-cluster neighbors, all mapped to its bond partner.

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
# 2-spin bond cluster
# ---------------------------------------------------------------------------

def two_site_cluster_positions() -> list[LatticePoint3D]:
    """Return the 2 sites of the bond cluster: (0,0,0) and its +x neighbor."""

    return [(0, 0, 0), (1, 0, 0)]


def bipartite_two_site_map_3d() -> Callable[[LatticePoint3D], int]:
    """Map the infinite cubic lattice onto the two-site A/B sublattices."""

    def map_to_cluster_index(point: LatticePoint3D) -> int:
        return (point[0] + point[1] + point[2]) & 1

    return map_to_cluster_index


def write_two_site_config(coupling: float, label: str) -> None:
    """Write a 2-spin bond configuration to HDF5."""

    project_name = "CubicLattice"
    config_file = f"Cubic_3D_N=2_NN_{label}"

    cluster_positions = two_site_cluster_positions()
    J = nearest_neighbor_couplings_3d(cluster_positions, coupling)
    map_fn = bipartite_two_site_map_3d()

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


def write_two_site_uncoupled_config(mf_coupling: float, label: str) -> None:
    """Write a 2-spin bond configuration with J=0 and all-neighbor mean-field weights.

    Mirrors :func:`write_two_site_config` but with ``include_cluster_neighbors=True``
    so the bond partner itself also contributes to the Weiss field, matching the
    uncoupled-cluster convention used by the square-lattice N=2 config.
    """

    project_name = "CubicLattice"
    config_file = f"Cubic_3D_N=2_NN_{label}"

    cluster_positions = two_site_cluster_positions()
    map_fn = bipartite_two_site_map_3d()

    J = np.zeros((2, 2), dtype=np.float64)
    base_weights = generate_nn_cluster_weights(
        cluster_positions=cluster_positions,
        nn_displacements=NN_DISPLACEMENTS_3D,
        map_to_cluster_index=map_fn,
        include_cluster_neighbors=True,
    )
    mf_weights = mf_coupling * base_weights.mf_expectation_weights
    corr_weights = (mf_coupling * mf_coupling) * base_weights.correlation_weights

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
# 7-spin star cluster
# ---------------------------------------------------------------------------

def star_cluster_positions() -> list[LatticePoint3D]:
    """Return the 7 sites of the star cluster: center + 6 NN arms."""

    return [(0, 0, 0)] + list(NN_DISPLACEMENTS_3D)


def cubic_sublattice_parity(point: LatticePoint3D) -> int:
    """Return the bipartite sublattice ``(x + y + z) mod 2`` of a cubic site."""

    return (point[0] + point[1] + point[2]) & 1


def star_cluster_map(
    cluster_positions: list[LatticePoint3D],
) -> Callable[[LatticePoint3D], int]:
    """Map any lattice point to the nearest *same-sublattice* star-cluster site.

    The cubic lattice is bipartite with sublattice ``(x + y + z) mod 2``, and a
    nearest neighbor always flips that parity.  In the 7-spin star the center
    (0,0,0) is the only A-site and the six arms are all B-sites, so every
    external neighbor of an arm lies on sublattice A.  A plain nearest-distance
    map would send such a neighbor to the arm itself (distance 1) — a B-site —
    mixing the sublattices.  Restricting the nearest-site search to cluster
    sites of the *same* parity keeps the mean-field embedding sublattice-
    consistent: external A-neighbors are represented by the A center, so a B-arm
    only ever feels A contributions in its Weiss field.  Ties are broken by
    lowest site index.
    """

    def map_to_cluster_index(point: LatticePoint3D) -> int:
        target_parity = cubic_sublattice_parity(point)
        candidates = [
            (sum((a - b) ** 2 for a, b in zip(point, pos)), i)
            for i, pos in enumerate(cluster_positions)
            if cubic_sublattice_parity(pos) == target_parity
        ]
        return min(candidates)[1]

    return map_to_cluster_index


def write_star_config(coupling: float, label: str) -> None:
    """Write a 7-spin star configuration to HDF5."""

    project_name = "CubicLattice"
    config_file = f"Cubic_3D_N=7_NN_{label}"

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
    config_file = f"Cubic_3D_N=8_NN_{label}"

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
    """Write FM and AFM configs for all 3D cluster geometries."""

    # 2-spin bond: AFM (J=0.5), FM (J=-0.5), and J=0 variants
    write_two_site_config(coupling=0.5, label="J=0.5")
    write_two_site_config(coupling=-0.5, label="J=-0.5")
    write_two_site_config(coupling=0.0, label="J=0.0")
    write_two_site_uncoupled_config(mf_coupling=0.5, label="J=0.0_uncoupled")

    # 7-spin star: AFM (J=0.5) and FM (J=-0.5)
    write_star_config(coupling=0.5, label="J=0.5")
    write_star_config(coupling=-0.5, label="J=-0.5")

    # 8-spin cube: AFM (J=0.5) and FM (J=-0.5)
    write_cube_config(coupling=0.5, label="J=0.5")
    write_cube_config(coupling=-0.5, label="J=-0.5")


if __name__ == "__main__":
    main()
