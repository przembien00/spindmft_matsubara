from __future__ import annotations

"""Generate square-lattice nearest-neighbor cluster configurations.

The generator is written for rectangular ``Lx x Ly`` clusters embedded into the
infinite square lattice. ``main()`` writes the one-site and 2x2 plaquette cases.
"""

import numpy as np
from typing import Callable

try:
    from config_generator import (
        canonicalize_hypercubic_displacement,
        create_config_file,
        generate_nn_cluster_weights,
    )
except ModuleNotFoundError:
    from .config_generator import (
        canonicalize_hypercubic_displacement,
        create_config_file,
        generate_nn_cluster_weights,
    )


LatticePoint2D = tuple[int, int]


def rectangular_cluster_positions(lx: int, ly: int) -> list[LatticePoint2D]:
    """Return lattice sites of an ``lx x ly`` square-lattice cluster."""

    return [(x, y) for y in range(ly) for x in range(lx)]


def nearest_neighbor_couplings(
    cluster_positions: list[LatticePoint2D],
    coupling: float,
) -> np.ndarray:
    """Build the nearest-neighbor coupling matrix for the given cluster."""

    num_spins = len(cluster_positions)
    position_to_index = {position: index for index, position in enumerate(cluster_positions)}
    J = np.zeros((num_spins, num_spins), dtype=np.float64)

    for index, (x, y) in enumerate(cluster_positions):
        for neighbor in ((x + 1, y), (x, y + 1)):
            neighbor_index = position_to_index.get(neighbor)
            if neighbor_index is None:
                continue
            J[index, neighbor_index] = coupling
            J[neighbor_index, index] = coupling

    return J


def rectangular_cluster_map(lx: int, ly: int):
    """Map any square-lattice site back into the rectangular cluster cell."""

    position_to_index = {
        (x, y): index
        for index, (x, y) in enumerate(rectangular_cluster_positions(lx, ly))
    }

    def map_to_cluster_index(point: LatticePoint2D) -> int:
        return position_to_index[(point[0] % lx, point[1] % ly)]

    return map_to_cluster_index


def reflective_cluster_map(lx: int, ly: int):
    """Map any square-lattice site back into the cluster by reflective folding.

    Unlike :func:`rectangular_cluster_map`, which folds with a periodic
    ``x % lx`` wrap, this folds the external site back through the cluster
    boundaries by reflection.  Reflection about an integer axis is
    ``x' = 2 c - x``, which preserves the parity of each coordinate and hence
    the bipartite sublattice ``(x + y) mod 2``.  This matters for clusters with
    an odd side length: there the periodic wrap (period = odd) flips the
    sublattice, so an A-site would pick up A-site representatives for its
    external mean field instead of the physically correct B-site ones.

    With reflective folding, the missing outside neighbor of an edge spin is
    represented by its mirror image just inside the cluster, which sits on the
    same sublattice as the true external neighbor.
    """

    position_to_index = {
        (x, y): index
        for index, (x, y) in enumerate(rectangular_cluster_positions(lx, ly))
    }

    def reflect(coordinate: int, length: int) -> int:
        if length == 1:
            return 0
        period = 2 * (length - 1)
        coordinate %= period
        return coordinate if coordinate < length else period - coordinate

    def map_to_cluster_index(point: LatticePoint2D) -> int:
        return position_to_index[(reflect(point[0], lx), reflect(point[1], ly))]

    return map_to_cluster_index


def sublattice_preserving_cluster_map(lx: int, ly: int):
    """Return a cluster fold map that preserves the bipartite sublattice.

    For even side lengths the periodic wrap already preserves parity, so the
    cheaper :func:`rectangular_cluster_map` is used.  For any odd side length
    the periodic period is odd and flips the sublattice, so the reflective fold
    is used instead.
    """

    if lx % 2 == 1 or ly % 2 == 1:
        return reflective_cluster_map(lx, ly)
    return rectangular_cluster_map(lx, ly)


def bipartite_two_site_map() -> Callable[[LatticePoint2D], int]:
    """Map the infinite square lattice onto the two-site A/B sublattices."""

    def map_to_cluster_index(point: LatticePoint2D) -> int:
        return (point[0] + point[1]) & 1

    return map_to_cluster_index


def write_config(lx: int, ly: int, coupling: float, label: str) -> None:
    """Write one square-lattice nearest-neighbor configuration to HDF5."""

    project_name = "SquareLattice"
    num_spins = lx * ly
    config_file = f"Square_2D_N={num_spins}_NN_{label}"

    cluster_positions = rectangular_cluster_positions(lx, ly)
    nn_displacements = [(1, 0), (-1, 0), (0, 1), (0, -1)]

    map_to_cluster_index = sublattice_preserving_cluster_map(lx, ly)
    if num_spins == 2 and lx == 2 and ly == 1:
        map_to_cluster_index = bipartite_two_site_map()

    J = nearest_neighbor_couplings(cluster_positions, coupling)
    base_weights = generate_nn_cluster_weights(
        cluster_positions=cluster_positions,
        nn_displacements=nn_displacements,
        map_to_cluster_index=map_to_cluster_index,
    )
    effective_mf_coupling = 0.5 if coupling == 0.0 else coupling
    mf_expectation_weights = effective_mf_coupling * base_weights.mf_expectation_weights
    correlation_weights = (effective_mf_coupling * effective_mf_coupling) * base_weights.correlation_weights

    output_path = create_config_file(
        project_name=project_name,
        config_file=config_file,
        J=J,
        mf_expectation_weights=mf_expectation_weights,
        correlation_weights=correlation_weights,
        correlation_categories=base_weights.canonical_pair_categories,
    )
    print(f"wrote {output_path}")




def write_uncoupled_config(lx: int, ly: int, mf_coupling: float, label: str) -> None:
    """Write a square-lattice configuration with J=0 and all-neighbor mean-field weights.

    The cluster Hamiltonian has no spin-spin coupling (J=0), but the mean-field
    weights are computed by treating every nearest neighbor—including the other
    cluster sites—as part of the environment.  This is achieved by passing
    ``include_cluster_neighbors=True`` to ``generate_nn_cluster_weights``.
    """

    project_name = "SquareLattice"
    num_spins = lx * ly
    config_file = f"Square_2D_N={num_spins}_NN_{label}"

    cluster_positions = rectangular_cluster_positions(lx, ly)
    nn_displacements = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    map_to_cluster_index = sublattice_preserving_cluster_map(lx, ly)

    J = np.zeros((num_spins, num_spins), dtype=np.float64)
    base_weights = generate_nn_cluster_weights(
        cluster_positions=cluster_positions,
        nn_displacements=nn_displacements,
        map_to_cluster_index=map_to_cluster_index,
        include_cluster_neighbors=True,
        canonical_displacement=canonicalize_hypercubic_displacement,
    )
    mf_expectation_weights = mf_coupling * base_weights.mf_expectation_weights
    correlation_weights = (mf_coupling * mf_coupling) * base_weights.correlation_weights

    output_path = create_config_file(
        project_name=project_name,
        config_file=config_file,
        J=J,
        mf_expectation_weights=mf_expectation_weights,
        correlation_weights=correlation_weights,
        correlation_categories=base_weights.canonical_pair_categories,
    )
    print(f"wrote {output_path}")


def write_uncoupled_nofield_config(lx: int, ly: int, mf_coupling: float, label: str) -> None:
    """Write an uncoupled config with the first-moment (Weiss) field switched off.

    Identical to :func:`write_uncoupled_config` (J=0, all-neighbor mean-field
    weights, ``include_cluster_neighbors=True``) except that the first-moment
    expectation weights are set to zero while the second-moment correlation
    weights are kept unchanged.  This removes the static Weiss field that lets
    ``<S>`` order spontaneously, leaving only the (sign-blind) bath fluctuations,
    and is meant to test whether the low-T non-convergence is the first-moment
    runaway: with W=0 no magnetization can grow.
    """

    project_name = "SquareLattice"
    num_spins = lx * ly
    config_file = f"Square_2D_N={num_spins}_NN_{label}"

    cluster_positions = rectangular_cluster_positions(lx, ly)
    nn_displacements = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    map_to_cluster_index = sublattice_preserving_cluster_map(lx, ly)

    J = np.zeros((num_spins, num_spins), dtype=np.float64)
    base_weights = generate_nn_cluster_weights(
        cluster_positions=cluster_positions,
        nn_displacements=nn_displacements,
        map_to_cluster_index=map_to_cluster_index,
        include_cluster_neighbors=True,
        canonical_displacement=canonicalize_hypercubic_displacement,
    )
    mf_expectation_weights = np.zeros_like(base_weights.mf_expectation_weights)
    correlation_weights = (mf_coupling * mf_coupling) * base_weights.correlation_weights

    output_path = create_config_file(
        project_name=project_name,
        config_file=config_file,
        J=J,
        mf_expectation_weights=mf_expectation_weights,
        correlation_weights=correlation_weights,
        correlation_categories=base_weights.canonical_pair_categories,
    )
    print(f"wrote {output_path}")


def main() -> None:
    """Write square-lattice configurations for N=1, 2, 4, 9 clusters."""

    # N=1 single-site
    write_config(lx=1, ly=1, coupling=0.5, label="J=0.5")
    write_config(lx=1, ly=1, coupling=-0.5, label="J=-0.5")

    # N=2 two-site chain (bipartite A/B mapping)
    write_config(lx=2, ly=1, coupling=0.5, label="J=0.5")
    write_config(lx=2, ly=1, coupling=-0.5, label="J=-0.5")
    write_config(lx=2, ly=1, coupling=0.0, label="J=0.0")
    write_uncoupled_config(lx=2, ly=1, mf_coupling=0.5, label="J=0.0_uncoupled")

    # N=4 2x2 plaquette
    write_config(lx=2, ly=2, coupling=0.5, label="J=0.5")
    write_config(lx=2, ly=2, coupling=-0.5, label="J=-0.5")
    write_config(lx=2, ly=2, coupling=0.0, label="J=0.0")
    write_uncoupled_config(lx=2, ly=2, mf_coupling=0.5, label="J=0.0_uncoupled")
    write_uncoupled_nofield_config(lx=2, ly=2, mf_coupling=0.5, label="J=0.0_uncoupled_nofield")

    # N=9 3x3 plaquette
    write_config(lx=3, ly=3, coupling=0.5, label="J=0.5")
    write_config(lx=3, ly=3, coupling=-0.5, label="J=-0.5")
    write_uncoupled_config(lx=3, ly=3, mf_coupling=0.5, label="J=0.0_uncoupled")
    write_uncoupled_nofield_config(lx=3, ly=3, mf_coupling=0.5, label="J=0.0_uncoupled_nofield")



if __name__ == "__main__":
    main()
