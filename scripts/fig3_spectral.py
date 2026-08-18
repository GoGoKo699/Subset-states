#!/usr/bin/env python3
"""Figure 3: reduced-density-matrix bulk and isolated positive mode."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from subset_states.core import (
    matrix_from_support,
    mean_matrix_uniform_eigenvalue,
    random_subset,
    summary_stats,
)
from subset_states.csv_plotting import plot_fig3
from subset_states.experiments import write_rows
from subset_states.tables import PEAK_M


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=24)
    parser.add_argument("--m", type=int, default=None)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--bins", type=int, default=500)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig3")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    m = args.m if args.m is not None else PEAK_M[args.n]
    rng = np.random.default_rng(args.seed)
    support = random_subset(args.n, m, rng)
    matrix = matrix_from_support(args.n, support)
    rho = matrix.conj().T @ matrix
    eigvals = np.linalg.eigvalsh(rho).real
    eigvals[eigvals < 1e-15] = 0.0
    eigvals /= eigvals.sum()
    eigvals.sort()
    write_rows(
        args.outdir / "fig3_spectrum.csv",
        [{"index": i, "lambda": float(lam)} for i, lam in enumerate(eigvals)],
        ["index", "lambda"],
    )

    lambda0 = float(eigvals[-1])
    bulk = eigvals[:-1]
    stats = summary_stats(bulk[bulk > 0])
    write_rows(
        args.outdir / "fig3_spectrum_summary.csv",
        [
            {
                "n": args.n,
                "M": m,
                "lambda0": lambda0,
                "lambda_mean_matrix": mean_matrix_uniform_eigenvalue(args.n, m),
                "M_over_N": m / (1 << args.n),
                "bulk_nonzero_count": stats.count,
                "bulk_mean": stats.mean,
                "bulk_std": stats.std,
                "bulk_min": stats.minimum,
                "bulk_max": stats.maximum,
            }
        ],
        [
            "n",
            "M",
            "lambda0",
            "lambda_mean_matrix",
            "M_over_N",
            "bulk_nonzero_count",
            "bulk_mean",
            "bulk_std",
            "bulk_min",
            "bulk_max",
        ],
    )
    plot_fig3(args.outdir, args.outdir / "fig3_spectral_bulk.pdf", bins=args.bins)


if __name__ == "__main__":
    main()
