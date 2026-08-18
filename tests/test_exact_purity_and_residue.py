from __future__ import annotations
from itertools import combinations
import unittest
import numpy as np
from subset_states.core import (
    entropy_from_support,
    fixed_cardinality_average_purity,
    matrix_from_support,
    residue_class_entropy_ceiling,
)

class ExactPurityTests(unittest.TestCase):
    def test_exhaustive_n4_average_purity(self) -> None:
        n = 4
        N = 1 << n
        for m in range(1, N + 1):
            purities = []
            for support in combinations(range(N), m):
                c = matrix_from_support(n, support)
                rho = c.T @ c
                purities.append(float(np.trace(rho @ rho)))
            self.assertAlmostEqual(np.mean(purities), fixed_cardinality_average_purity(n, m), places=12)

class ResidueCeilingTests(unittest.TestCase):
    def test_exhaustive_n4_ceiling(self) -> None:
        n = 4
        N = 1 << n
        for mask in range(1, 1 << N):
            support = np.asarray([x for x in range(N) if mask >> x & 1], dtype=np.int64)
            entropy = entropy_from_support(n, support)
            for t in (1, 2):
                self.assertLessEqual(entropy, residue_class_entropy_ceiling(n, support, t) + 2e-12)

if __name__ == "__main__":
    unittest.main()
