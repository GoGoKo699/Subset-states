#!/usr/bin/env python3
"""Figure 6: Rényi entropies of random subset states."""
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

from subset_states.core import summary_stats
from subset_states.experiments import m_grid, random_subset_renyi_samples, write_rows
from subset_states.plotting import apply_journal_style, save_figure

ORDERS = (1.0, 2.0, np.inf)
ORDER_LABELS = {1.0: r"$S^{(1)}$", 2.0: r"$S^{(2)}$", np.inf: r"$S^{(\infty)}$"}


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


def rows_for_order(rows: list[dict], order: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    key = "inf" if np.isinf(order) else order
    filtered = [row for row in rows if row["order"] == key]
    return (
        np.array([row["m"] for row in filtered]),
        np.array([row["mean"] for row in filtered]),
        np.array([row["sem"] for row in filtered]),
    )


def draw(axis, rows: list[dict], *, markers: bool = True) -> None:
    styles = {1.0: "o-", 2.0: "s--", np.inf: "^:"}
    for order in ORDERS:
        m, mean, sem = rows_for_order(rows, order)
        fmt = styles[order] if markers else styles[order].replace("o", "").replace("s", "").replace("^", "")
        axis.errorbar(m, mean, yerr=sem, fmt=fmt, markersize=2.2, capsize=1.2, label=ORDER_LABELS[order])


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
    if args.zoom_stop > N:
        print(f"warning: --zoom-stop={args.zoom_stop} exceeds N=2**n={N}; using --zoom-stop={zoom_stop}.")

    rows = compute(args.n, m_grid(1, N, args.points), args.samples, args.seed)
    zoom_rows = compute(args.n, m_grid(1, zoom_stop, args.zoom_points), args.zoom_samples, args.seed + 10_000)
    write_rows(args.outdir / "fig6_renyi_summary.csv", rows, rows[0].keys())
    write_rows(args.outdir / "fig6_renyi_zoom_summary.csv", zoom_rows, zoom_rows[0].keys())

    apply_journal_style()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.7, 6.0), sharey=False)
    draw(ax1, rows)
    ax1.set_xlabel(r"support size $M$")
    ax1.set_ylabel("Rényi entropy")
    ax1.set_xlim(0, 1 << args.n)
    ax1.set_ylim(0, args.n / 2)
    ax1.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.42))

    draw(ax2, zoom_rows)
    ax2.set_xlabel(r"support size $M$")
    ax2.set_ylabel("Rényi entropy")
    ax2.set_xlim(0, zoom_stop)
    ax2.set_ylim(max(0, args.n / 2 - 3), args.n / 2)
    ax2.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.42))
    fig.subplots_adjust(hspace=0.55)
    save_figure(fig, args.outdir / "fig6_renyi.pdf")


if __name__ == "__main__":
    main()
