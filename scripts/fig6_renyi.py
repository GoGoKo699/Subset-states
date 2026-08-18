#!/usr/bin/env python3
"""Figure 6: Rényi entropies of random subset states."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from subset_states.core import summary_stats
from subset_states.csv_plotting import plot_fig6
from subset_states.experiments import m_grid, random_subset_renyi_samples, write_rows

ORDERS = (1.0, 2.0, np.inf)


def compute(n: int, m_values: np.ndarray, samples: int, seed: int) -> list[dict]:
    rows: list[dict] = []
    for offset, m in enumerate(m_values):
        values = random_subset_renyi_samples(n, int(m), samples, seed=seed + offset, orders=ORDERS)
        for order in ORDERS:
            stats = summary_stats(values[float(order)])
            rows.append(
                {
                    "m": int(m),
                    "order": "inf" if np.isinf(order) else order,
                    "mean": stats.mean,
                    "std": stats.std,
                    "sem": stats.sem,
                    "count": stats.count,
                    "minimum": stats.minimum,
                    "maximum": stats.maximum,
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=14)
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--zoom-stop", type=int, default=1600)
    parser.add_argument("--zoom-points", type=int, default=100)
    parser.add_argument("--zoom-samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig6")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    N = 1 << args.n
    zoom_stop = min(args.zoom_stop, N)
    rows = compute(args.n, m_grid(1, N, args.points), args.samples, args.seed)
    zoom = compute(
        args.n,
        m_grid(1, zoom_stop, args.zoom_points),
        args.zoom_samples,
        args.seed + 10_000,
    )
    write_rows(args.outdir / "fig6_renyi_summary.csv", rows, rows[0].keys())
    write_rows(args.outdir / "fig6_renyi_zoom_summary.csv", zoom, zoom[0].keys())
    plot_fig6(args.outdir, args.outdir / "fig6_renyi.pdf", n=args.n)


if __name__ == "__main__":
    main()
