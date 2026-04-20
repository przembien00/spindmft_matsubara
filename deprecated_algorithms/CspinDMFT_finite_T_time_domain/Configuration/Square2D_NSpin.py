from __future__ import annotations

"""Generate square-lattice nearest-neighbor cluster configurations.

The generator is written for rectangular ``Lx x Ly`` clusters embedded into the
infinite square lattice. For now, ``main()`` writes the 2x2 plaquette case.
"""

import numpy as np

try:
    from config_generator import create_config_file, generate_nn_cluster_weights
except ModuleNotFoundError:
    from .config_generator import create_config_file, generate_nn_cluster_weights


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


def write_config(lx: int, ly: int, coupling: float, label: str) -> None:
    """Write one square-lattice nearest-neighbor configuration to HDF5."""

    project_name = "SquareLattice"
    num_spins = lx * ly
    config_file = f"Square_2D_N={num_spins}_NN_{label}"

    cluster_positions = rectangular_cluster_positions(lx, ly)
    nn_displacements = [(1, 0), (-1, 0), (0, 1), (0, -1)]

    J = nearest_neighbor_couplings(cluster_positions, coupling)
    base_weights = generate_nn_cluster_weights(
        cluster_positions=cluster_positions,
        nn_displacements=nn_displacements,
        map_to_cluster_index=rectangular_cluster_map(lx, ly),
    )
    mf_expectation_weights = coupling * base_weights.mf_expectation_weights
    correlation_weights = (coupling * coupling) * base_weights.correlation_weights

    output_path = create_config_file(
        project_name=project_name,
        config_file=config_file,
        J=J,
        mf_expectation_weights=mf_expectation_weights,
        correlation_weights=correlation_weights,
    )
    print(f"wrote {output_path}")


def main() -> None:
    """Write the 2x2 square plaquette for antiferromagnetic and ferromagnetic J."""

    write_config(lx=2, ly=2, coupling=0.5, label="J=0.5")
    write_config(lx=2, ly=2, coupling=-0.5, label="J=-0.5")


if __name__ == "__main__":
    main()
