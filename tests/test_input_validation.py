from __future__ import annotations

from decimal import Decimal, localcontext
from itertools import permutations
import unittest

import numpy as np

from subset_states.core import (
    fixed_cardinality_mean_spectrum,
    hypergeometric_occupancy_approximation,
    matrix_from_state_vector,
    matrix_from_support,
    qft_state_from_support,
    renyi_entropy_from_spectrum,
    residue_class_entropy_ceiling,
    spectrum_from_matrix,
)


class SupportValidationTests(unittest.TestCase):
    def test_support_apis_reject_malformed_labels_before_casting(self) -> None:
        invalid = (
            [], [[1, 2]], 1, [1, 1], [-1], [16], [1.5], [1.0],
            [True], [1j], [np.nan], [np.inf], ["1"],
            np.asarray([2**64 - 1], dtype=np.uint64),
        )
        functions = (
            lambda support: matrix_from_support(4, support),
            lambda support: qft_state_from_support(4, support),
            lambda support: residue_class_entropy_ceiling(4, support, 1),
        )
        for support in invalid:
            for function in functions:
                with self.subTest(support=support, function=function):
                    with self.assertRaises(ValueError):
                        function(support)

    def test_qft_accepts_odd_qubit_count_and_preserves_norm(self) -> None:
        support = np.asarray([0, 3, 7], dtype=np.uint64)
        state = qft_state_from_support(3, support)
        np.testing.assert_allclose(np.vdot(state, state), 1.0, atol=1e-14)

    def test_integer_parameters_are_not_silently_truncated(self) -> None:
        for n in (0, -2, 4.0, True):
            with self.subTest(n=n), self.assertRaises(ValueError):
                qft_state_from_support(n, [0])
        for m in (1.5, True, -1):
            for function in (fixed_cardinality_mean_spectrum, hypergeometric_occupancy_approximation):
                with self.subTest(m=m, function=function), self.assertRaises(ValueError):
                    function(4, m)
        for t in (1.5, True, -1, 3):
            with self.subTest(t=t), self.assertRaises(ValueError):
                residue_class_entropy_ceiling(4, [0, 1], t)

    def test_normalized_integer_dtype_cannot_silently_erase_amplitudes(self) -> None:
        with self.assertRaises(ValueError):
            matrix_from_support(4, [0, 1], dtype=int)
        binary = matrix_from_support(4, [0, 1], normalize=False, dtype=int)
        self.assertEqual(binary.sum(), 2)


class CoefficientMatrixTests(unittest.TestCase):
    def test_vector_and_support_use_identical_basis_order(self) -> None:
        support = [1, 4, 9, 15]
        state = np.zeros(16)
        state[support] = 1 / np.sqrt(len(support))
        for right in permutations(range(4), 2):
            with self.subTest(right=right):
                np.testing.assert_array_equal(
                    matrix_from_state_vector(4, state, right),
                    matrix_from_support(4, support, right),
                )

    def test_complex_coefficients_and_partial_trace_in_computational_basis(self) -> None:
        state = np.arange(16) + 1j * np.arange(16)[::-1] ** 2
        state = state / np.linalg.norm(state)
        expected_matrix = state.reshape(4, 4)
        matrix = matrix_from_state_vector(4, state)
        np.testing.assert_array_equal(matrix, expected_matrix)
        expected_rho = np.zeros((4, 4), dtype=complex)
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    expected_rho[b, c] += state[4 * a + b] * state[4 * a + c].conjugate()
        np.testing.assert_allclose(matrix.T @ matrix.conj(), expected_rho, atol=1e-15)
        self.assertGreater(np.max(abs(matrix.conj().T @ matrix - expected_rho)), 0.01)
        expected_spectrum = np.linalg.eigvalsh(expected_rho)
        expected_spectrum = expected_spectrum[expected_spectrum > 1e-14]
        np.testing.assert_allclose(spectrum_from_matrix(matrix), expected_spectrum, atol=1e-14)

    def test_partitions_and_state_values_are_validated(self) -> None:
        for right in ((0, 0), (0, 4), (0, 1.5), (0, True)):
            for function, data in ((matrix_from_state_vector, np.ones(16)), (matrix_from_support, [0])):
                with self.subTest(right=right, function=function), self.assertRaises(ValueError):
                    function(4, data, right)
        for state in (np.ones((4, 4)), np.full(16, np.nan), np.full(16, np.inf)):
            with self.subTest(state=state), self.assertRaises(ValueError):
                matrix_from_state_vector(4, state)


class SpectrumValidationTests(unittest.TestCase):
    def test_invalid_spectra_and_orders_are_rejected(self) -> None:
        for spectrum in ([], [0, 0], [-0.1, 1.1], [np.nan, 1], [np.inf, 1], [[0.5, 0.5]], [1j]):
            with self.subTest(spectrum=spectrum), self.assertRaises(ValueError):
                renyi_entropy_from_spectrum(spectrum)
        for order in (0, -1, -np.inf, np.nan, 1j):
            with self.subTest(order=order), self.assertRaises(ValueError):
                renyi_entropy_from_spectrum([0.25, 0.75], order)

    def test_tolerated_roundoff_and_positive_infinite_order(self) -> None:
        self.assertAlmostEqual(renyi_entropy_from_spectrum([-1e-16, 0.5, 0.5]), 1)
        self.assertAlmostEqual(renyi_entropy_from_spectrum([0.25, 0.75], np.inf), -np.log2(0.75))

    def test_finite_orders_do_not_underflow_or_collapse_to_order_one(self) -> None:
        self.assertAlmostEqual(renyi_entropy_from_spectrum([0.5, 0.5], 2000), 1)
        p = np.asarray([0.2, 0.8])
        order = 1.0 + 1e-6
        with localcontext() as context:
            context.prec = 60
            alpha = Decimal.from_float(order)
            expected = (Decimal("0.2") ** alpha + Decimal("0.8") ** alpha).ln() / ((1 - alpha) * Decimal(2).ln())
        computed = renyi_entropy_from_spectrum(p, order)
        self.assertAlmostEqual(computed, float(expected), places=11)
        self.assertGreater(abs(computed - renyi_entropy_from_spectrum(p, 1)), 1e-8)
        self.assertAlmostEqual(renyi_entropy_from_spectrum(p, 1 + 1e-12), renyi_entropy_from_spectrum(p), places=11)

    def test_invalid_matrices_and_cutoffs_are_rejected(self) -> None:
        for matrix in ([], [1, 2], np.zeros((2, 2)), [[np.nan]], [[np.inf]]):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                spectrum_from_matrix(matrix)
        for atol in (-1, np.nan, np.inf):
            for function, argument in ((spectrum_from_matrix, [[1]]), (renyi_entropy_from_spectrum, [1])):
                with self.subTest(atol=atol, function=function), self.assertRaises(ValueError):
                    function(argument, atol=atol)


if __name__ == "__main__":
    unittest.main()
