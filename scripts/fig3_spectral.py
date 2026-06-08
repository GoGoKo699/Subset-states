#!/usr/bin/env python3
"""Figure 3: spectrum of the reduced density matrix at peak support."""
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

from subset_states.core import matrix_from_support, random_subset, summary_stats
from subset_states.experiments import write_rows
from subset_states.plotting import apply_journal_style, save_figure
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
    omega = matrix_from_support(args.n, support)
    rho = omega.conj().T @ omega
    eigvals = np.linalg.eigvalsh(rho).real
    eigvals[eigvals < 1e-15] = 0.0
    eigvals = eigvals / eigvals.sum()
    eigvals.sort()

    rows = [{"index": i, "lambda": float(lam)} for i, lam in enumerate(eigvals)]
    write_rows(args.outdir / "fig3_spectrum.csv", rows, ["index", "lambda"])

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
            "M_over_N",
            "bulk_nonzero_count",
            "bulk_mean",
            "bulk_std",
            "bulk_min",
            "bulk_max",
        ],
    )

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    positive_bulk = bulk[bulk > 0]
    ax.hist(positive_bulk, bins=args.bins)
    ax.set_xlabel(r"eigenvalue $\lambda$ excluding $\lambda_0$")
    ax.set_ylabel("count")
    ax.set_xlim(0, positive_bulk.max() * 1.05)
    save_figure(fig, args.outdir / "fig3_spectral_bulk.pdf")
    print(f"lambda0={lambda0:.8g}, M/N={m / (1 << args.n):.8g}")


if __name__ == "__main__":
    main()
