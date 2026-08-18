#!/usr/bin/env python3
"""Independent checks for the final subset-state manuscript candidate."""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.special import digamma

ROOT = Path(__file__).resolve().parents[1]


def entropy_bits(rho: np.ndarray, tol: float = 1e-14) -> float:
    vals = np.linalg.eigvalsh(rho).real
    vals = vals[vals > tol]
    vals = vals / vals.sum()
    return float(-np.sum(vals * np.log2(vals)))


def h_shannon(probs: np.ndarray) -> float:
    p = probs[probs > 0]
    return float(-np.sum(p * np.log2(p))) if p.size else 0.0


def exact_mean_rho(n: int, M: int) -> np.ndarray:
    d = 1 << (n // 2)
    N = 1 << n
    beta = (M - 1) / (d * (N - 1))
    return (1 / d - beta) * np.eye(d) + beta * np.ones((d, d))


def exact_avg_purity(n: int, M: int) -> float:
    d = 1 << (n // 2)
    N = 1 << n
    return (
        1 / M
        + 2 * (d - 1) * (M - 1) / (M * (N - 1))
        + (d - 1) ** 2 * (M - 1) * (M - 2) * (M - 3)
        / (M * (N - 1) * (N - 2) * (N - 3))
    )


def exhaustive_n4_checks() -> dict[str, float]:
    n = 4
    N = 1 << n
    d = 1 << (n // 2)
    m = n // 2
    rho_sum = np.zeros((N + 1, d, d), dtype=float)
    purity_sum = np.zeros(N + 1, dtype=float)
    counts = np.zeros(N + 1, dtype=np.int64)
    max_residue_violation = {1: -math.inf, 2: -math.inf}
    max_entropy_routine_diff = 0.0

    labels = np.arange(N, dtype=np.uint32)
    batch_size = 4096
    for lo in range(1, 1 << N, batch_size):
        hi = min(1 << N, lo + batch_size)
        masks = np.arange(lo, hi, dtype=np.uint32)
        bits = ((masks[:, None] >> labels[None, :]) & 1).astype(float)
        Ms = bits.sum(axis=1).astype(np.int64)
        X = bits.reshape(-1, d, d)
        rho = np.einsum('bai,baj->bij', X, X) / Ms[:, None, None]
        evals = np.linalg.eigvalsh(rho).real
        positive = np.where(evals > 1e-14, evals, 1.0)
        entropy = -np.sum(np.where(evals > 1e-14, evals * np.log2(positive), 0.0), axis=1)
        if lo == 1:
            entropy_direct = np.asarray([entropy_bits(r) for r in rho[:128]])
            max_entropy_routine_diff = max(
                max_entropy_routine_diff,
                float(np.max(np.abs(entropy[:128] - entropy_direct))),
            )
        purity = np.sum(evals * evals, axis=1)

        np.add.at(rho_sum, Ms, rho)
        purity_sum += np.bincount(Ms, weights=purity, minlength=N + 1)
        counts += np.bincount(Ms, minlength=N + 1)

        for t in (1, 2):
            modulus = 1 << t
            residue_counts = np.stack(
                [bits[:, (labels % modulus) == r].sum(axis=1) for r in range(modulus)],
                axis=1,
            )
            probs = residue_counts / Ms[:, None]
            p_safe = np.where(probs > 0, probs, 1.0)
            H = -np.sum(np.where(probs > 0, probs * np.log2(p_safe), 0.0), axis=1)
            bound = m - t + H
            max_residue_violation[t] = max(
                max_residue_violation[t],
                float(np.max(entropy - bound)),
            )

    max_mean_rho_error = 0.0
    max_purity_error = 0.0
    for M in range(1, N + 1):
        empirical_rho = rho_sum[M] / counts[M]
        empirical_purity = purity_sum[M] / counts[M]
        max_mean_rho_error = max(
            max_mean_rho_error,
            float(np.max(np.abs(empirical_rho - exact_mean_rho(n, M)))),
        )
        max_purity_error = max(
            max_purity_error,
            abs(empirical_purity - exact_avg_purity(n, M)),
        )

    return {
        "supports_enumerated": float((1 << N) - 1),
        "max_mean_rho_error": max_mean_rho_error,
        "max_average_purity_error": max_purity_error,
        "max_residue_bound_violation_t1": max_residue_violation[1],
        "max_residue_bound_violation_t2": max_residue_violation[2],
        "max_entropy_routine_difference": max_entropy_routine_diff,
    }

def page_entropy(n: int) -> float:
    d = 1 << (n // 2)
    return float((digamma(d * d + 1) - digamma(d + 1) - (d - 1) / (2 * d)) / math.log(2))


def table_fit_checks() -> dict[str, float]:
    rows = np.asarray(
        [
            (10, 107, 4.072),
            (12, 276, 5.108),
            (14, 716, 6.143),
            (16, 1873, 7.176),
            (18, 4934, 8.196),
            (20, 13091, 9.215),
            (22, 34771, 10.231),
            (24, 93018, 11.242),
            (26, 250660, 12.251),
            (28, 672556, 13.258),
            (30, 1836685, 14.263),
        ],
        dtype=float,
    )
    n, M, S = rows.T
    slope_S, intercept_S = np.polyfit(n, S, 1)
    slope_M, intercept_M = np.polyfit(n, np.log2(M), 1)
    page = np.asarray([page_entropy(int(x)) for x in n])
    gaps = page - S
    fractions = M / np.power(2.0, n)
    assert np.all(np.diff(gaps) < 0), "Page gaps are not strictly decreasing"
    assert np.all(np.diff(fractions) < 0), "support fractions are not strictly decreasing"
    return {
        "S_slope": float(slope_S),
        "S_intercept": float(intercept_S),
        "log2M_slope": float(slope_M),
        "log2M_intercept": float(intercept_M),
        "page_gap_n10": float(gaps[0]),
        "page_gap_n30": float(gaps[-1]),
        "support_fraction_n10": float(fractions[0]),
        "support_fraction_n30": float(fractions[-1]),
    }


def dense_ansatz_equivalence() -> float:
    n = 14
    N = 1 << n
    d = 1 << (n // 2)
    max_diff = 0.0
    for M in range(1, N):
        lam = 1 / d + (d - 1) * (M - 1) / (d * (N - 1))
        old = (
            (1 - lam) * (math.log(d) - math.log(1 - lam) - 0.5)
            - lam * math.log(lam)
        ) / math.log(2)
        h2 = -lam * math.log2(lam) - (1 - lam) * math.log2(1 - lam)
        new = h2 + (1 - lam) * (math.log2(d) - 1 / (2 * math.log(2)))
        max_diff = max(max_diff, abs(old - new))
    return max_diff


def residue_checks() -> dict[str, float]:
    path = ROOT / "data" / "residue_matched_summary.csv"
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))

    by_key = {(int(r["k"]), int(r["matched_low_bits"])): r for r in rows}
    strongest = [(1, 1), (2, 2), (3, 3)]
    expected = {
        (1, 1): (0.124881, 0.114442, 0.885607, 0.894136),
        (2, 2): (0.230201, 0.212232, 0.841422, 0.851975),
        (3, 3): (0.111500, 0.103991, 0.910577, 0.916102),
    }
    max_claim_error = 0.0
    for key in strongest:
        r = by_key[key]
        vals = (
            float(r["position_residual_deficit"]),
            float(r["fourier_residual_deficit"]),
            float(r["cardinality_deficit_reduction_position"]),
            float(r["cardinality_deficit_reduction_fourier"]),
        )
        max_claim_error = max(max_claim_error, max(abs(a - b) for a, b in zip(vals, expected[key])))
        assert int(r["position_null_count_at_or_below_structured"]) == 0
        assert int(r["fourier_null_count_at_or_below_structured"]) == 0

    prime = by_key[(1, 1)]
    forced_gap = float(prime["forced_gap_from_balanced_max_bits"])
    ceiling = float(prime["entropy_ceiling_bits"])
    return {
        "max_residue_claim_rounding_error": max_claim_error,
        "prime_parity_ceiling_bits": ceiling,
        "prime_forced_gap_bits": forced_gap,
        "all_six_strongest_controls_below_all_1000_nulls": 1.0,
    }


def main() -> None:
    results = {
        "exhaustive_n4": exhaustive_n4_checks(),
        "table_fits": table_fit_checks(),
        "dense_ansatz": {"max_old_new_difference": dense_ansatz_equivalence()},
        "residue_controls": residue_checks(),
    }

    lines = ["FINAL SCIENTIFIC VALIDATION", "=" * 29, ""]
    for section, values in results.items():
        lines.append(section)
        lines.append("-" * len(section))
        for key, value in values.items():
            lines.append(f"{key}: {value:.16g}")
        lines.append("")

    # Hard acceptance thresholds.
    ex = results["exhaustive_n4"]
    assert ex["max_mean_rho_error"] < 5e-13
    assert ex["max_average_purity_error"] < 5e-13
    assert ex["max_residue_bound_violation_t1"] < 5e-13
    assert ex["max_residue_bound_violation_t2"] < 5e-13
    assert results["dense_ansatz"]["max_old_new_difference"] < 5e-13
    tf = results["table_fits"]
    assert abs(tf["S_slope"] - 0.5093) < 1e-12
    assert abs(tf["S_intercept"] + 0.9900909090909096) < 1e-12
    assert abs(tf["log2M_slope"] - 0.7035409513877945) < 1e-12
    assert abs(tf["log2M_intercept"] + 0.35773436814817394) < 1e-12

    lines.append("STATUS: PASS")
    out = ROOT / "validation" / "FINAL_SCIENTIFIC_VALIDATION.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    print(out.read_text())


if __name__ == "__main__":
    main()
