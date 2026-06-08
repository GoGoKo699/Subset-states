#!/usr/bin/env python3
"""Figure 2: finite-size scaling of the Table-I peak values.

This script does not rediscover the peaks. It takes the production Table-I values
as input, computes log2(M_n), the Page benchmark, and the two linear regressions
reported in the manuscript.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from subset_states.core import exact_page_entropy_bits
from subset_states.experiments import write_rows
from subset_states.plotting import apply_journal_style, save_figure
from subset_states.tables import read_table_i


def fit_line(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    result = stats.linregress(x, y)
    residuals = y - (result.slope * x + result.intercept)
    return {
        "slope": float(result.slope),
        "intercept": float(result.intercept),
        "r_value": float(result.rvalue),
        "r_squared": float(result.rvalue**2),
        "p_value": float(result.pvalue),
        "slope_stderr": float(result.stderr),
        "intercept_stderr": float(result.intercept_stderr),
        "residual_std": float(residuals.std(ddof=2)) if len(x) > 2 else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=Path, default=ROOT / "data" / "table_i_peaks.csv", help="CSV with columns n,M_n,S_n")
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig2")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows_table = read_table_i(args.table)
    n = np.asarray([row["n"] for row in rows_table], dtype=float)
    M_n = np.asarray([row["M_n"] for row in rows_table], dtype=float)
    S_n = np.asarray([row["S_n"] for row in rows_table], dtype=float)
    log2_M = np.log2(M_n)
    page = np.asarray([exact_page_entropy_bits(int(x)) for x in n], dtype=float)

    fit_logM = fit_line(n, log2_M)
    fit_S = fit_line(n, S_n)

    rows = []
    for ni, mi, si, logmi, pagei in zip(n, M_n, S_n, log2_M, page):
        rows.append(
            {
                "n": int(ni),
                "M_n": int(mi),
                "log2_M_n": float(logmi),
                "S_n": float(si),
                "Page_exact_bits": float(pagei),
                "Page_minus_S_n": float(pagei - si),
            }
        )
    write_rows(args.outdir / "fig2_table_with_page.csv", rows, rows[0].keys())

    fit_rows = [
        {"quantity": "log2_M_n", **fit_logM},
        {"quantity": "S_n", **fit_S},
    ]
    write_rows(args.outdir / "fig2_linear_fit_summary.csv", fit_rows, fit_rows[0].keys())

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    x = np.linspace(n.min() - 1, n.max() + 1, 200)
    ax.plot(n, log2_M, "o", label=r"$\log_2 M_n$")
    ax.plot(x, fit_logM["slope"] * x + fit_logM["intercept"], "--", linewidth=1)
    ax.plot(n, S_n, "s", label=r"$S_n$")
    ax.plot(x, fit_S["slope"] * x + fit_S["intercept"], "--", linewidth=1)
    ax.plot(n, page, ":", label="Page average")
    ax.set_xlabel(r"number of qubits $n$")
    ax.set_ylabel("bits")
    ax.set_xlim(n.min() - 1, n.max() + 1)
    ax.set_ylim(0, max(log2_M.max(), page.max()) + 1)
    ax.legend(frameon=False, loc="upper left")
    save_figure(fig, args.outdir / "fig2_peak_scaling.pdf")

    print("log2(M_n) fit:", fit_logM)
    print("S_n fit:", fit_S)


if __name__ == "__main__":
    main()
