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
    canonical_pair_categories: np.ndarray  # shape (n_canonical, 2), dtype uint32


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


def canonicalize_hypercubic_displacement(d: LatticePoint) -> LatticePoint:
    """Canonicalize a displacement under the full point group of a hypercubic lattice.

    A square (or cubic) Bravais lattice is invariant under permuting and
    sign-flipping its axes, so two displacements that agree after sorting
    their absolute component values (e.g. ``(1, 0)`` and ``(0, -1)``, or
    ``(1, 1)`` and ``(-1, 1)``) connect physically equivalent site pairs.
    Used to merge correlation categories that are only distinguished by
    lattice orientation, not by physics.
    """

    return tuple(sorted(abs(c) for c in d))


def generate_nn_cluster_weights(
    cluster_positions: list[LatticePoint],
    nn_displacements: list[LatticePoint],
    map_to_cluster_index: ClusterIndexMap,
    include_cluster_neighbors: bool = False,
    canonical_displacement: Callable[[LatticePoint], LatticePoint] | None = None,
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

    ``canonical_displacement``, if given, is applied to every displacement
    before it is used to identify or look up a representative cluster pair.
    This additionally merges cluster pairs whose displacements are related by
    a lattice symmetry beyond plain translation (e.g. two nearest-neighbor
    pairs pointing along ``x`` and along ``y``), on top of the translational
    merging that always happens. Pass
    :func:`canonicalize_hypercubic_displacement` to fully exploit the point
    group of a square/cubic lattice; leave as ``None`` to only merge by exact
    (translational) displacement match, as before.
    """

    def _canon(d: LatticePoint) -> LatticePoint:
        return canonical_displacement(d) if canonical_displacement is not None else d

    num_spins = len(cluster_positions)
    cluster_position_set = set(cluster_positions)

    mf_expectation_weights = np.zeros((num_spins, num_spins), dtype=np.float64)
    external_neighbors_per_site: list[list[LatticePoint]] = [[] for _ in range(num_spins)]

    for i, cluster_site in enumerate(cluster_positions):
        for displacement in nn_displacements:
            external_neighbor = add_lattice_points(cluster_site, displacement)
            if (not include_cluster_neighbors) and external_neighbor in cluster_position_set:
                continue

            external_neighbors_per_site[i].append(external_neighbor)
            representative_site = map_to_cluster_index(external_neighbor)
            mf_expectation_weights[i, representative_site] += 1.0

    # Build a map from displacement to canonical cluster pair.  The canonical
    # pair for a given displacement d is the first (lowest-index) cluster pair
    # whose separation equals d.  Both d and -d map to the same canonical pair
    # so that the direction of traversal does not matter.
    #
    # This ensures that a pair of external neighbors contributes to the cluster
    # pair that is *geometrically equivalent* to them (same separation vector),
    # rather than to the pair obtained by mapping each neighbor independently
    # via translational symmetry (which can assign geometrically mismatched
    # cluster pairs, e.g. NNN external pairs to NN cluster pairs).
    ndim = len(cluster_positions[0])
    displacement_to_canonical_pair: dict[tuple[int, ...], tuple[int, int]] = {}
    for k_ in range(num_spins):
        for l_ in range(k_, num_spins):
            d = subtract_lattice_points(cluster_positions[l_], cluster_positions[k_])
            key = _canon(d)
            if key not in displacement_to_canonical_pair:
                displacement_to_canonical_pair[key] = (k_, l_)
            rev_key = _canon(tuple(-x for x in d))
            if rev_key not in displacement_to_canonical_pair:
                displacement_to_canonical_pair[rev_key] = (k_, l_)

    correlation_weights = np.zeros((num_spins, num_spins, num_spins, num_spins), dtype=np.float64)
    for i in range(num_spins):
        for j in range(i, num_spins):
            for neighbor_i in external_neighbors_per_site[i]:
                for neighbor_j in external_neighbors_per_site[j]:
                    external_displacement = subtract_lattice_points(neighbor_j, neighbor_i)
                    canonical = displacement_to_canonical_pair.get(_canon(external_displacement))
                    if canonical is None:
                        continue
                    k_ord, l_ord = canonical
                    correlation_weights[i, j, k_ord, l_ord] += 1.0
                    correlation_weights[j, i, k_ord, l_ord] = correlation_weights[i, j, k_ord, l_ord]

    canonical_pairs = sorted(set(displacement_to_canonical_pair.values()))
    canonical_pair_categories = np.asarray(canonical_pairs, dtype=np.uint32)

    return NNClusterWeights(
        mf_expectation_weights=mf_expectation_weights,
        correlation_weights=correlation_weights,
        canonical_pair_categories=canonical_pair_categories,
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
