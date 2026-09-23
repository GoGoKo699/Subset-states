#!/usr/bin/env python3
"""Validate exact small-system mathematics and consistency of retained CSV evidence.

Prints Markdown without changing the repository. This check is not a peak-search
rerun, an asymptotic theorem, a novelty assessment, or a publication certificate.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states import core, tables
from scripts.fig2_peak_scaling import fit_line
from scripts.validate_residue_results import (
    close, finite, markdown_report, read_rows, require, residue_checks, write_report,
)


def exhaustive_checks(n_values: tuple[int, ...] = (2, 4)) -> dict[str, float | int]:
    """Independently enumerate supports and compare the production implementations."""
    errors = {name: 0.0 for name in (
        "coefficient_matrix", "entropy", "residue_ceiling", "mean_spectrum", "average_purity",
        "diagonal_entropy", "entropy_lower_bound", "qft_amplitude", "qft_entropy",
    )}
    supports_enumerated = 0
    max_bound_violation = 0.0
    for n in n_values:
        require(n in (2, 4), "exhaustive checks are deliberately limited to n=2 and n=4")
        N, d = 1 << n, 1 << (n // 2)
        rho_sum = np.zeros((N + 1, d, d))
        purity_sum = np.zeros(N + 1)
        diagonal_sum = np.zeros(N + 1)
        entropy_sum = np.zeros(N + 1)
        counts = np.zeros(N + 1, dtype=int)
        labels = np.arange(N, dtype=np.uint32)
        # Batch eigensolvers provide an independent route from the production SVD.
        for lo in range(1, 1 << N, 4096):
            masks = np.arange(lo, min(1 << N, lo + 4096), dtype=np.uint32)
            bits = ((masks[:, None] >> labels[None, :]) & 1).astype(float)
            sizes = bits.sum(axis=1).astype(int)
            matrices = bits.reshape(-1, d, d) / np.sqrt(sizes[:, None, None])
            rho = np.einsum("bai,baj->bij", matrices, matrices)
            eigenvalues = np.linalg.eigvalsh(rho)
            safe = np.where(eigenvalues > 1e-14, eigenvalues, 1.0)
            entropies = -np.sum(np.where(eigenvalues > 1e-14, eigenvalues * np.log2(safe), 0), axis=1)
            diagonal = np.diagonal(rho, axis1=1, axis2=2)
            safe_diag = np.where(diagonal > 0, diagonal, 1.0)
            diagonal_entropies = -np.sum(diagonal * np.log2(safe_diag), axis=1)
            np.add.at(rho_sum, sizes, rho)
            purity_sum += np.bincount(sizes, weights=np.sum(rho * rho, axis=(1, 2)), minlength=N + 1)
            diagonal_sum += np.bincount(sizes, weights=diagonal_entropies, minlength=N + 1)
            entropy_sum += np.bincount(sizes, weights=entropies, minlength=N + 1)
            counts += np.bincount(sizes, minlength=N + 1)
            for index, size in enumerate(sizes):
                support = np.flatnonzero(bits[index])
                production_matrix = core.matrix_from_support(n, support)
                matrix_error = float(np.max(np.abs(production_matrix - matrices[index])))
                require(math.isfinite(matrix_error), "non-finite production coefficient matrix")
                errors["coefficient_matrix"] = max(errors["coefficient_matrix"], matrix_error)
                errors["entropy"] = max(errors["entropy"], close(
                    core.entropy_from_support(n, support), float(entropies[index]), "production entropy"))
                for t in range(n // 2 + 1):
                    population = np.bincount(support % (1 << t), minlength=1 << t)
                    probs = population[population > 0] / size
                    bound = n / 2 - t - float(np.sum(probs * np.log2(probs)))
                    errors["residue_ceiling"] = max(errors["residue_ceiling"], close(
                        core.residue_class_entropy_ceiling(n, support, t), bound, "production residue ceiling"))
                    max_bound_violation = max(max_bound_violation, float(entropies[index]) - bound)
                if n == 2:
                    # Explicit DFT checks sign and normalization, independently of FFT.
                    phases = np.exp(2j * np.pi * np.outer(np.arange(N), support) / N)
                    direct_state = phases.sum(axis=1) / math.sqrt(N * size)
                    actual_state = core.qft_state_from_support(n, support)
                    amplitude_error = float(np.max(np.abs(actual_state - direct_state)))
                    require(math.isfinite(amplitude_error), "non-finite production QFT")
                    errors["qft_amplitude"] = max(errors["qft_amplitude"], amplitude_error)
                    c = direct_state.reshape(d, d)
                    eig = np.linalg.eigvalsh(c.conj().T @ c)
                    positive = eig[eig > 1e-14]
                    direct_entropy = float(-np.sum(positive * np.log2(positive)))
                    errors["qft_entropy"] = max(errors["qft_entropy"], close(
                        core.entropy_from_state_vector(n, actual_state), direct_entropy, "production QFT entropy"))
            supports_enumerated += len(masks)
        for size in range(1, N + 1):
            require(counts[size] == math.comb(N, size), f"n={n}, M={size}: enumeration count")
            mean_rho = rho_sum[size] / counts[size]
            uniform, orthogonal = core.fixed_cardinality_mean_spectrum(n, size)
            predicted_rho = orthogonal * np.eye(d) + (uniform - orthogonal) * np.ones((d, d)) / d
            error = float(np.max(np.abs(predicted_rho - mean_rho)))
            require(math.isfinite(error), "non-finite production mean spectrum")
            errors["mean_spectrum"] = max(errors["mean_spectrum"], error)
            close(core.mean_matrix_uniform_eigenvalue(n, size), uniform, "production uniform eigenvalue")
            mean_purity = float(purity_sum[size] / counts[size])
            errors["average_purity"] = max(errors["average_purity"], close(
                core.fixed_cardinality_average_purity(n, size), mean_purity, "production average purity"))
            errors["diagonal_entropy"] = max(errors["diagonal_entropy"], close(
                core.hypergeometric_occupancy_approximation(n, size), float(diagonal_sum[size] / counts[size]),
                "production diagonal entropy"))
            lower = core.average_entropy_lower_bound_from_purity(n, size)
            errors["entropy_lower_bound"] = max(errors["entropy_lower_bound"], close(
                lower, -math.log2(mean_purity), "production purity entropy bound"))
            max_bound_violation = max(max_bound_violation, lower - float(entropy_sum[size] / counts[size]))
    for name, error in errors.items():
        require(error < 5e-12, f"{name} maximum error {error} exceeds tolerance")
    require(max_bound_violation < 5e-12, f"entropy bound violated by {max_bound_violation}")
    return {"supports_enumerated": supports_enumerated,
            **{f"max_{name}_error": value for name, value in errors.items()},
            "max_entropy_bound_violation": max_bound_violation}


def independent_page_entropy(n: int) -> float:
    """Page expression using harmonic sums/expansion, independently of digamma."""
    d = 1 << (n // 2)
    def harmonic(number: int) -> float:
        if number < 32:
            return math.fsum(1 / k for k in range(1, number + 1))
        x = 1.0 / number
        return (math.log(number) + 0.5772156649015328606 + x / 2 - x**2 / 12
                + x**4 / 120 - x**6 / 252 + x**8 / 240 - x**10 / 132)
    return (harmonic(d * d) - harmonic(d) - (d - 1) / (2 * d)) / math.log(2)


def table_fit_checks(root: Path = ROOT) -> dict[str, float | int]:
    rows = read_rows(root / "data/table_i_peaks.csv", {"n", "M_n", "S_n"})
    require([int(row["n"]) for row in rows] == list(range(10, 31, 2)), "Table I n grid is incomplete or duplicated")
    table = np.asarray([(int(row["n"]), int(row["M_n"]), finite(row, "S_n")) for row in rows], dtype=float)
    require(np.array_equal(table, np.asarray([[r["n"], r["M_n"], r["S_n"]]
            for r in tables.read_table_i(root / "data/table_i_peaks.csv")])), "production Table I parser disagrees")
    require(np.array_equal(table, tables.table_i_array()), "Table I CSV and production defaults disagree")
    n, sizes, entropies = table.T
    require(np.all((sizes >= 1) & (sizes <= 2**n)), "Table I support sizes out of range")
    require(np.all((entropies >= 0) & (entropies <= n / 2)), "Table I entropy out of range")
    page = np.asarray([independent_page_entropy(int(value)) for value in n])
    for ni, value in zip(n, page):
        close(core.exact_page_entropy_bits(int(ni)), float(value), "production Page entropy")
    derived = read_rows(root / "outputs/fig2/fig2_table_with_page.csv", {
        "n", "M_n", "log2_M_n", "S_n", "Page_exact_bits", "Page_minus_S_n"})
    require(len(derived) == len(rows), "derived Table I row count mismatch")
    max_error = 0.0
    for row, ni, size, entropy, page_value in zip(derived, n, sizes, entropies, page):
        for field, value in {"n": ni, "M_n": size, "log2_M_n": math.log2(size), "S_n": entropy,
                             "Page_exact_bits": page_value, "Page_minus_S_n": page_value - entropy}.items():
            max_error = max(max_error, close(finite(row, field), float(value), f"derived Table I {field}"))
    fit_fields = {"slope", "intercept", "r_value", "r_squared", "p_value", "slope_stderr", "intercept_stderr", "residual_std"}
    fits = read_rows(root / "outputs/fig2/fig2_linear_fit_summary.csv", {"quantity"} | fit_fields)
    by_quantity = {row["quantity"]: row for row in fits}
    require(len(fits) == 2 and set(by_quantity) == {"log2_M_n", "S_n"}, "fit quantities missing or duplicated")
    results: dict[str, float | int] = {"table_rows_checked": len(rows)}
    for quantity, values in (("log2_M_n", np.log2(sizes)), ("S_n", entropies)):
        slope, intercept = np.polyfit(n, values, 1)
        production = fit_line(n, values)
        close(production["slope"], float(slope), f"independent {quantity} slope")
        close(production["intercept"], float(intercept), f"independent {quantity} intercept")
        for field in fit_fields:
            # Relative tolerance is needed for tiny p-values; do not accept zero in their place.
            expected = production[field]
            max_error = max(max_error, close(finite(by_quantity[quantity], field), expected,
                f"stored {quantity} fit {field}", atol=max(abs(expected) * 1e-8, 1e-30)))
        results[f"{quantity}_slope"] = float(slope)
        results[f"{quantity}_intercept"] = float(intercept)
    gaps = page - entropies
    require(np.all(np.diff(gaps) < 0), "recorded Page gaps are not decreasing")
    require(np.all(np.diff(sizes / 2**n) < 0), "recorded support fractions are not decreasing")
    results["max_derived_table_or_fit_error"] = max_error
    return results


def dense_ansatz_checks() -> dict[str, float | int]:
    max_error, tested = 0.0, 0
    for n in (2, 4, 8, 14):
        N, d = 1 << n, 1 << (n // 2)
        for size in range(1, N + 1):
            lam = (N - size + d * (size - 1)) / (d * (N - 1))
            if size == N:
                independent = 0.0
            else:
                h2 = -lam * math.log2(lam) - (1 - lam) * math.log2(1 - lam)
                independent = h2 + (1 - lam) * (math.log2(d) - 1 / (2 * math.log(2)))
            max_error = max(max_error, close(core.dense_bulk_approximation(n, size), independent,
                                             "production dense-bulk ansatz", atol=5e-13))
            tested += 1
    return {"parameter_pairs_checked": tested, "max_formula_error": max_error}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="save the report to a new .md file")
    args = parser.parse_args()
    try:
        # Check files first so malformed or stale evidence fails promptly.
        results = {"Retained Table I consistency": table_fit_checks(),
                   "Exhaustive small systems": exhaustive_checks(),
                   "Dense-bulk formula implementation": dense_ansatz_checks(),
                   "Residue-control evidence": residue_checks()}
        report = markdown_report("Scientific sanity checks", results,
            "Exhaustive checks cover every nonempty support at $`n=2`$ and $`n=4`$ against production routines. "
            "Table I checks read the actual retained CSV and its derived outputs; they do not rerun the "
            "peak searches or validate an asymptotic law. The residue checks reconstruct all 12,000-row "
            "statistics and replay three seeded supports per group (36 total). Formula agreement does "
            "not establish the dense-bulk ansatz's approximation accuracy. Novelty remains a separate question.")
        write_report(report, args.output)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"# Scientific sanity checks\n\n**Status: FAIL**\n\n{exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
