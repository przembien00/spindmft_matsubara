"""Repository-wide generator and migration checks; run from shared test support."""
from pathlib import Path
import importlib.util
import sys
import tempfile
import unittest
import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[2]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


migration = load(REPO / 'CspinDMFT_finite_T/Configuration/migrate_correlation_weights.py',
                 'migration')
migrate = migration.migrate
VERSION_ATTRIBUTE = migration.VERSION_ATTRIBUTE

class CorrelationWeightsTest(unittest.TestCase):
    def test_generators_count_per_orientation(self):
        for algorithm in ('CspinDMFT_finite_T', 'CspinDMFT_finite_T_conditional_pcn'):
            generator = load(REPO / algorithm / 'Configuration/config_generator.py', algorithm)
            positions = [(0, 0), (1, 0), (0, 1), (1, 1)]
            nn = [(1, 0), (-1, 0), (0, 1), (0, -1)]
            # Independent ordered-neighbor enumeration, grouped only by distance class.
            canonical = {}
            for k in range(4):
                for l in range(k, 4):
                    key = tuple(sorted(abs(a-b) for a, b in zip(positions[k], positions[l])))
                    canonical.setdefault(key, (k, l))
            for include_cluster in (False, True):
                result = generator.generate_nn_cluster_weights(
                    positions, nn, lambda p: (p[0] % 2) + 2*(p[1] % 2),
                    include_cluster_neighbors=include_cluster,
                    canonical_displacement=generator.canonicalize_hypercubic_displacement)
                expected = np.zeros((4, 4, 4, 4))
                neighbors = [[(p[0]+d[0], p[1]+d[1]) for d in nn
                              if include_cluster or (p[0]+d[0], p[1]+d[1]) not in positions]
                             for p in positions]
                for i in range(4):
                    for j in range(4):
                        for a in neighbors[i]:
                            for b in neighbors[j]:
                                key = tuple(sorted(abs(x-y) for x, y in zip(a, b)))
                                if key in canonical:
                                    k, l = canonical[key]
                                    expected[i, j, k, l] += 1
                effective = result.correlation_weights.copy()
                for k, l in result.canonical_pair_categories:
                    effective[:, :, k, l] *= 1 if k == l else 2
                np.testing.assert_array_equal(effective, expected)

    def test_periodic_normalization(self):
        folder = REPO / 'CspinDMFT_finite_T_conditional_pcn/Configuration'
        # Supply the same generator import used by the standalone periodic script.
        load(folder / 'config_generator.py', 'config_generator')
        periodic = load(folder / 'Square2D_Uncoupled_Periodic.py', 'periodic')
        for size in (2, 4):
            for coupling in (-0.5, 0.5):
                _, weights, categories = periodic.build_periodic_uncoupled_weights(size, coupling)
                effective = sum((1 if k == l else 2)*weights[:, :, k, l] for k, l in categories)
                np.testing.assert_allclose(effective, (4*coupling)**2)

    def test_migration_and_idempotence(self):
        for categories in ([[0, 0], [0, 1]], [[0, 0], [0, 1], [1, 0]]):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'legacy.hdf5'
                old = np.arange(16, dtype=float).reshape(2, 2, 2, 2)
                with h5py.File(path, 'w') as f:
                    f['correlation weights'] = old
                    f['correlation_categories'] = categories
                    f['expectation weights'] = np.eye(2)
                self.assertTrue(migrate(path))
                self.assertFalse(migrate(path))
                with h5py.File(path) as f:
                    self.assertEqual(f.attrs[VERSION_ATTRIBUTE], 2)
                    np.testing.assert_array_equal(f['expectation weights'][:], np.eye(2))
                    np.testing.assert_array_equal(f['correlation weights'][:, :, 0, 0], old[:, :, 0, 0])
                    factor = 2 if len(categories) == 2 else 1
                    np.testing.assert_array_equal(factor*f['correlation weights'][:, :, 0, 1], old[:, :, 0, 1])
