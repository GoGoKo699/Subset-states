#!/usr/bin/env python3
"""Figure 2: finite-size scaling of the retained Table-I estimates."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from scipy import stats

from subset_states.core import exact_page_entropy_bits
from subset_states.csv_plotting import plot_fig2
from subset_states.experiments import write_rows
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
    parser.add_argument("--table", type=Path, default=ROOT / "data" / "table_i_peaks.csv")
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig2")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    table = read_table_i(args.table)
    n = np.asarray([row["n"] for row in table], dtype=float)
    M_n = np.asarray([row["M_n"] for row in table], dtype=float)
    S_n = np.asarray([row["S_n"] for row in table], dtype=float)
    log2_M = np.log2(M_n)
    page = np.asarray([exact_page_entropy_bits(int(x)) for x in n])
    rows = [
        {
            "n": int(ni),
            "M_n": int(mi),
            "log2_M_n": float(logmi),
            "S_n": float(si),
            "Page_exact_bits": float(pagei),
            "Page_minus_S_n": float(pagei - si),
        }
        for ni, mi, si, logmi, pagei in zip(n, M_n, S_n, log2_M, page)
    ]
    write_rows(args.outdir / "fig2_table_with_page.csv", rows, rows[0].keys())
    fits = [
        {"quantity": "log2_M_n", **fit_line(n, log2_M)},
        {"quantity": "S_n", **fit_line(n, S_n)},
    ]
    write_rows(args.outdir / "fig2_linear_fit_summary.csv", fits, fits[0].keys())
    plot_fig2(args.outdir, args.outdir / "fig2_peak_scaling.pdf", write_fit_summary=False)


if __name__ == "__main__":
    main()
