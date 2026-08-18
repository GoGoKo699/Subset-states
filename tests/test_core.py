from __future__ import annotations

from itertools import combinations
import unittest

import numpy as np

from subset_states.core import (
    dense_bulk_approximation,
    entropy_from_state_vector,
    fixed_cardinality_mean_spectrum,
    hypergeometric_occupancy_approximation,
    matrix_from_support,
    qft_state_from_support,
)


class FixedCardinalityMeanTests(unittest.TestCase):
    def test_spectrum_trace_and_endpoints(self) -> None:
        for n in (4, 6, 8):
            d = 1 << (n // 2)
            N = 1 << n
            for m in (1, max(2, N // 7), N):
                uniform, orthogonal = fixed_cardinality_mean_spectrum(n, m)
                self.assertAlmostEqual(uniform + (d - 1) * orthogonal, 1.0, places=13)
                self.assertGreaterEqual(uniform, orthogonal)
                self.assertGreaterEqual(orthogonal, -1e-14)
            uniform, orthogonal = fixed_cardinality_mean_spectrum(n, N)
            self.assertAlmostEqual(uniform, 1.0, places=13)
            self.assertAlmostEqual(orthogonal, 0.0, places=13)

    def test_exact_small_ensemble_mean_matrix(self) -> None:
        n, m = 4, 3
        N = 1 << n
        matrices = []
        for support in combinations(range(N), m):
            c = matrix_from_support(n, support)
            matrices.append(c.T @ c)
        empirical = np.mean(matrices, axis=0)
        uniform, orthogonal = fixed_cardinality_mean_spectrum(n, m)
        d = 1 << (n // 2)
        analytic = orthogonal * np.eye(d) + (uniform - orthogonal) * np.ones((d, d)) / d
        np.testing.assert_allclose(empirical, analytic, rtol=0, atol=2e-14)


class ApproximationTests(unittest.TestCase):
    def test_hypergeometric_diagonal_entropy_against_enumeration(self) -> None:
        n, m = 4, 3
        N = 1 << n
        entropies = []
        for support in combinations(range(N), m):
            c = matrix_from_support(n, support)
            diagonal = np.diag(c.T @ c)
            positive = diagonal[diagonal > 0]
            entropies.append(float(-np.sum(positive * np.log2(positive))))
        exact_average = float(np.mean(entropies))
        formula = hypergeometric_occupancy_approximation(n, m)
        self.assertAlmostEqual(formula, exact_average, places=13)

    def test_dense_approximation_is_finite_and_has_full_support_endpoint(self) -> None:
        for n in (4, 8, 14):
            N = 1 << n
            for m in (1, max(2, N // 20), max(2, N // 3), N):
                self.assertTrue(np.isfinite(dense_bulk_approximation(n, m)))
            self.assertEqual(dense_bulk_approximation(n, N), 0.0)


class QFTConventionTests(unittest.TestCase):
    def test_positive_exponent_qft(self) -> None:
        n = 4
        N = 1 << n
        support = np.asarray([0, 1, 3, 8], dtype=int)
        state = np.zeros(N, dtype=complex)
        state[support] = 1 / np.sqrt(len(support))
        j = np.arange(N)
        direct = np.exp(2j * np.pi * np.outer(j, j) / N) @ state / np.sqrt(N)
        computed = qft_state_from_support(n, support)
        np.testing.assert_allclose(computed, direct, rtol=0, atol=2e-14)

    def test_fourier_sign_does_not_change_entropy_for_real_input(self) -> None:
        n = 4
        N = 1 << n
        support = np.asarray([0, 2, 5, 7, 10], dtype=int)
        state = np.zeros(N, dtype=complex)
        state[support] = 1 / np.sqrt(len(support))
        positive = np.fft.ifft(state, norm="ortho")
        negative = np.fft.fft(state, norm="ortho")
        self.assertAlmostEqual(
            entropy_from_state_vector(n, positive),
            entropy_from_state_vector(n, negative),
            places=13,
        )


if __name__ == "__main__":
    unittest.main()
