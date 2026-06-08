"""Numerical tools for subset-state entanglement experiments."""

from .peak_fitting import local_m_grid, quadratic_peak_fit
from .tables import PEAK_M, PEAK_S, TABLE_I_PEAKS, read_table_i, table_i_array

from .core import (
    balanced_bipartitions,
    dnm_approximation,
    entropy_from_state_vector,
    entropy_from_support,
    exact_page_entropy_bits,
    matrix_from_state_vector,
    matrix_from_support,
    natural_right_bits,
    omega_prime_factors_sieve,
    qft_state_from_support,
    random_balanced_partition,
    random_subset,
    renyi_entropy_from_spectrum,
    spectrum_from_matrix,
    summary_stats,
    tnm_approximation,
)

__all__ = [
    "PEAK_M",
    "PEAK_S",
    "TABLE_I_PEAKS",
    "local_m_grid",
    "quadratic_peak_fit",
    "read_table_i",
    "table_i_array",
    "balanced_bipartitions",
    "dnm_approximation",
    "entropy_from_state_vector",
    "entropy_from_support",
    "exact_page_entropy_bits",
    "matrix_from_state_vector",
    "matrix_from_support",
    "natural_right_bits",
    "omega_prime_factors_sieve",
    "qft_state_from_support",
    "random_balanced_partition",
    "random_subset",
    "renyi_entropy_from_spectrum",
    "spectrum_from_matrix",
    "summary_stats",
    "tnm_approximation",
]
