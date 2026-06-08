#!/usr/bin/env python3
"""Figure 4: numerical entropy vs. analytical approximations D_{N,M} and T_{N,M}."""
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

from subset_states.core import dnm_approximation, entropy_from_support, random_subset, summary_stats, tnm_approximation
from subset_states.experiments import m_grid, write_rows
from subset_states.plotting import apply_journal_style, make_inset_axis, save_figure


def compute(n: int, m_values: np.ndarray, samples: int, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for m in m_values:
        values = np.empty(samples, dtype=float)
        for sample in range(samples):
            support = random_subset(n, int(m), rng)
            values[sample] = entropy_from_support(n, support)
        stats = summary_stats(values)
        rows.append(
            {
                "m": int(m),
                "entropy_mean": stats.mean,
                "entropy_std": stats.std,
                "entropy_sem": stats.sem,
                "entropy_count": stats.count,
                "D_NM": dnm_approximation(n, int(m)),
                "T_NM": tnm_approximation(n, int(m)),
            }
        )
    return rows


def draw(axis, rows, label_data: bool) -> None:
    m = np.array([row["m"] for row in rows])
    mean = np.array([row["entropy_mean"] for row in rows])
    sem = np.array([row["entropy_sem"] for row in rows])
    D = np.array([row["D_NM"] for row in rows])
    T = np.array([row["T_NM"] for row in rows])
    axis.errorbar(m, mean, yerr=sem, fmt="o", markersize=2.2, capsize=1.2, label=r"numerical $S_{N,M}$" if label_data else None)
    axis.plot(m, D, ":", label=r"$D_{N,M}$" if label_data else None)
    axis.plot(m, T, "--", label=r"$T_{N,M}$" if label_data else None)


def plot(rows: list[dict], zoom_rows: list[dict], out_pdf: Path, n: int) -> None:
    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.8))
    draw(ax, rows, True)
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel("entropy")
    ax.set_xlim(0, 1 << n)
    ax.set_ylim(0, n / 2)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.35), ncol=3)

    inset = make_inset_axis(ax)
    draw(inset, zoom_rows, False)
    inset.set_xlim(0, max(row["m"] for row in zoom_rows) * 1.02)
    inset.set_ylim(max(0, n / 2 - 3), n / 2)
    inset.set_xlabel(r"$M$", labelpad=1)
    inset.set_ylabel(r"$S$", labelpad=1)
    save_figure(fig, out_pdf)


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
    if args.zoom_stop > N:
        print(f"warning: --zoom-stop={args.zoom_stop} exceeds N=2**n={N}; using --zoom-stop={zoom_stop}.")

    rows = compute(args.n, m_grid(1, N, args.points), args.samples, args.seed)
    zoom_rows = compute(args.n, m_grid(1, zoom_stop, args.zoom_points), args.zoom_samples, args.seed + 1)
    write_rows(args.outdir / "fig4_summary.csv", rows, rows[0].keys())
    write_rows(args.outdir / "fig4_zoom_summary.csv", zoom_rows, zoom_rows[0].keys())
    plot(rows, zoom_rows, args.outdir / "fig4_approximation.pdf", args.n)


if __name__ == "__main__":
    main()
