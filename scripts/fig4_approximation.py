#!/usr/bin/env python3
"""Figure 4: numerical entropy and corrected fixed-cardinality approximations."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from subset_states.core import (
    dense_bulk_approximation,
    entropy_from_support,
    hypergeometric_occupancy_approximation,
    random_subset,
    summary_stats,
)
from subset_states.csv_plotting import plot_fig4
from subset_states.experiments import m_grid, write_rows


def compute(n: int, m_values: np.ndarray, samples: int, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for m in m_values:
        values = np.empty(samples, dtype=float)
        for sample in range(samples):
            values[sample] = entropy_from_support(n, random_subset(n, int(m), rng))
        stats = summary_stats(values)
        rows.append(
            {
                "m": int(m),
                "entropy_mean": stats.mean,
                "entropy_std": stats.std,
                "entropy_sem": stats.sem,
                "entropy_count": stats.count,
                "D_NM_hypergeometric": hypergeometric_occupancy_approximation(n, int(m)),
                "T_NM_mean_matrix": dense_bulk_approximation(n, int(m)),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=14)
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--zoom-stop", type=int, default=1600)
    parser.add_argument("--zoom-points", type=int, default=100)
    parser.add_argument("--zoom-samples", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig4")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    N = 1 << args.n
    zoom_stop = min(args.zoom_stop, N)
    rows = compute(args.n, m_grid(1, N, args.points), args.samples, args.seed)
    zoom_rows = compute(
        args.n,
        m_grid(1, zoom_stop, args.zoom_points),
        args.zoom_samples,
        args.seed + 1,
    )
    write_rows(args.outdir / "fig4_summary.csv", rows, rows[0].keys())
    write_rows(args.outdir / "fig4_zoom_summary.csv", zoom_rows, zoom_rows[0].keys())
    plot_fig4(args.outdir, args.outdir / "fig4_approximation.pdf", n=args.n)


if __name__ == "__main__":
    main()
