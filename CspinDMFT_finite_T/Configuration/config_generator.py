from __future__ import annotations

"""Helpers for building cluster configuration HDF5 files.

The module provides small lattice-geometry utilities plus the logic that turns
an embedded finite cluster into the weight tensors expected by the C++ code.
"""

from dataclasses import dataclass
import os
from pathlib import Path
from typing import Callable

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as np

LatticePoint = tuple[int, ...]
ClusterIndexMap = Callable[[LatticePoint], int]


@dataclass
class NNClusterWeights:
    """Container for mean-field and correlation weights of a cluster embedding."""

    mf_expectation_weights: np.ndarray
    correlation_weights: np.ndarray


def create_whole_categories_list(num_spins: int) -> np.ndarray:
    """Return all unique spin-pair categories ``(i, j)`` with ``i <= j``."""

    categories = [(i, j) for i in range(num_spins) for j in range(i, num_spins)]
    return np.asarray(categories, dtype=np.uint32)


def add_lattice_points(p1: LatticePoint, p2: LatticePoint) -> LatticePoint:
    """Add two lattice-coordinate tuples component-wise."""

    return tuple(a + b for a, b in zip(p1, p2))


def subtract_lattice_points(p1: LatticePoint, p2: LatticePoint) -> LatticePoint:
    """Subtract two lattice-coordinate tuples component-wise as ``p1 - p2``."""

    return tuple(a - b for a, b in zip(p1, p2))


def generate_nn_cluster_weights(
    cluster_positions: list[LatticePoint],
    nn_displacements: list[LatticePoint],
    map_to_cluster_index: ClusterIndexMap,
) -> NNClusterWeights:
    """Build embedding weights for a nearest-neighbor cluster problem.

    ``cluster_positions`` defines the sites kept explicitly in the finite
    cluster. For every cluster site, the function enumerates nearest neighbors
    in ``nn_displacements`` and separates those that lie outside the cluster.
    Each external site is mapped back onto a representative cluster site with
    ``map_to_cluster_index``.

    The returned tensors encode:
    - ``mf_expectation_weights``: how many external neighbors of site ``i``
      belong to representative site ``k``.
    - ``correlation_weights``: how many pairs of external neighbors reproduce
      the same relative displacement as a representative cluster pair.
    """

    num_spins = len(cluster_positions)
    cluster_position_set = set(cluster_positions)

    mf_expectation_weights = np.zeros((num_spins, num_spins), dtype=np.float64)
    external_neighbors_per_site: list[list[LatticePoint]] = [[] for _ in range(num_spins)]

    for i, cluster_site in enumerate(cluster_positions):
        for displacement in nn_displacements:
            external_neighbor = add_lattice_points(cluster_site, displacement)
            if external_neighbor in cluster_position_set:
                continue

            external_neighbors_per_site[i].append(external_neighbor)
            representative_site = map_to_cluster_index(external_neighbor)
            mf_expectation_weights[i, representative_site] += 1.0

    correlation_weights = np.zeros((num_spins, num_spins, num_spins, num_spins), dtype=np.float64)
    for i in range(num_spins):
        for j in range(i, num_spins):
            for neighbor_i in external_neighbors_per_site[i]:
                for neighbor_j in external_neighbors_per_site[j]:
                    k = map_to_cluster_index(neighbor_i)
                    l = map_to_cluster_index(neighbor_j)
                    external_displacement = subtract_lattice_points(neighbor_j, neighbor_i)
                    representative_displacement = subtract_lattice_points(cluster_positions[l], cluster_positions[k])
                    if external_displacement == representative_displacement:
                        correlation_weights[i, j, k, l] += 1.0
                        correlation_weights[j, i, k, l] = correlation_weights[i, j, k, l]
                    elif external_displacement == tuple(-x for x in representative_displacement):
                        correlation_weights[i, j, l, k] += 1.0
                        correlation_weights[j, i, l, k] = correlation_weights[i, j, l, k]

    return NNClusterWeights(
        mf_expectation_weights=mf_expectation_weights,
        correlation_weights=correlation_weights,
    )


def create_config_file(
    project_name: str,
    config_file: str,
    J: np.ndarray,
    mf_expectation_weights: np.ndarray,
    correlation_weights: np.ndarray,
    correlation_categories: np.ndarray | None = None,
) -> Path:
    """Write a configuration HDF5 file consumed by the DMFT executable.

    Parameters are written into ``Configuration_Data/<project_name>`` using the
    dataset and attribute names expected by the existing C++ reader.
    """

    config_root = Path(__file__).resolve().parent.parent / "Configuration_Data"
    if project_name:
        config_root = config_root / project_name
    config_root.mkdir(parents=True, exist_ok=True)

    output_path = config_root / f"{config_file}.hdf5"
    if correlation_categories is None:
        correlation_categories = create_whole_categories_list(J.shape[0])

    if output_path.exists():
        output_path.unlink()

    with h5py.File(output_path, "w") as h5file:
        h5file.attrs["num_Spins"] = np.int64(J.shape[0])
        h5file.create_dataset("spin-spin couplings", data=np.asarray(J, dtype=np.float64))
        h5file.create_dataset("expectation weights", data=np.asarray(mf_expectation_weights, dtype=np.float64))
        h5file.create_dataset("correlation weights", data=np.asarray(correlation_weights, dtype=np.float64))
        h5file.attrs["num_Categories"] = np.int64(len(correlation_categories))
        h5file.create_dataset("correlation_categories", data=np.asarray(correlation_categories, dtype=np.uint32))

    return output_path
