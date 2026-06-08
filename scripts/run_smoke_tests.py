#!/usr/bin/env python3
"""Small deterministic checks for the subset-state utilities."""
from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from subset_states.peak_fitting import local_m_grid, quadratic_peak_fit
from subset_states.tables import PEAK_M, read_table_i

from subset_states.core import (
    balanced_bipartitions,
    dnm_approximation,
    entropy_from_state_vector,
    entropy_from_support,
    exact_page_entropy_bits,
    omega_prime_factors_sieve,
    qft_state_from_support,
    random_subset,
    tnm_approximation,
)


def main() -> None:
    rng = np.random.default_rng(123)
    n = 4
    N = 1 << n

    # M=1 and M=N are product states for the natural bipartition.
    assert abs(entropy_from_support(n, np.array([3]))) < 1e-12
    assert abs(entropy_from_support(n, np.arange(N))) < 1e-12

    # A two-Bell-pair/rainbow support has entropy n/2 for n=4.
    rainbow = np.array([0b0011, 0b0110, 0b1001, 0b1100])
    assert abs(entropy_from_support(n, rainbow) - 2.0) < 1e-12

    support = random_subset(n, 5, rng)
    qft_state = qft_state_from_support(n, support)
    assert abs(np.linalg.norm(qft_state) - 1.0) < 1e-12
    _ = entropy_from_state_vector(n, qft_state)

    partitions = list(balanced_bipartitions(10))
    assert len(partitions) == 126

    omega = omega_prime_factors_sieve(16)
    assert omega[0] == 0 and omega[1] == 0 and omega[8] == 3 and omega[12] == 3

    assert abs(dnm_approximation(n, 1)) < 1e-12
    assert tnm_approximation(n, 0) > 0
    assert abs(tnm_approximation(n, N)) < 1e-12
    assert exact_page_entropy_bits(10) > 0

    table = read_table_i(ROOT / "data" / "table_i_peaks.csv")
    assert table[0]["n"] == 10 and PEAK_M[14] == 716

    grid = local_m_grid(10, PEAK_M[10], points=5, span=0.1)
    assert PEAK_M[10] in grid and grid.size >= 5
    fit = quadratic_peak_fit(np.array([80, 107, 140]), np.array([4.0, 4.1, 4.0]), np.array([0.01, 0.01, 0.01]))
    assert "fit_status" in fit and fit["M_peak"] > 0

    print("smoke tests passed")


if __name__ == "__main__":
    main()
