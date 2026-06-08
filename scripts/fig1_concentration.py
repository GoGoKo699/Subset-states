#!/usr/bin/env python3
"""Figure 1: concentration of random-subset entanglement over states and partitions."""
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

from subset_states.core import entropy_from_support, random_balanced_partition, random_subset, summary_stats
from subset_states.experiments import m_grid, write_rows
from subset_states.plotting import apply_journal_style, make_inset_axis, save_figure


def compute_curve(n: int, m_values: np.ndarray, samples: int, seed: int) -> tuple[list[dict], list[dict]]:
    rng = np.random.default_rng(seed)
    summary_rows: list[dict] = []
    raw_rows: list[dict] = []
    for m in m_values:
        state_values = np.empty(samples, dtype=float)
        partition_values = np.empty(samples, dtype=float)

        fixed_support = random_subset(n, int(m), rng)
        for sample in range(samples):
            support = random_subset(n, int(m), rng)
            state_values[sample] = entropy_from_support(n, support)
            right_bits = random_balanced_partition(n, rng)
            partition_values[sample] = entropy_from_support(n, fixed_support, right_bits)
            raw_rows.append(
                {
                    "m": int(m),
                    "sample": sample,
                    "entropy_multiple_states": state_values[sample],
                    "entropy_multiple_partitions": partition_values[sample],
                }
            )

        st = summary_stats(state_values)
        pt = summary_stats(partition_values)
        summary_rows.append(
            {
                "m": int(m),
                "states_count": st.count,
                "states_mean": st.mean,
                "states_std": st.std,
                "states_sem": st.sem,
                "states_min": st.minimum,
                "states_max": st.maximum,
                "partitions_count": pt.count,
                "partitions_mean": pt.mean,
                "partitions_std": pt.std,
                "partitions_sem": pt.sem,
                "partitions_min": pt.minimum,
                "partitions_max": pt.maximum,
            }
        )
    return summary_rows, raw_rows


def plot(summary: list[dict], zoom_summary: list[dict], out_pdf: Path, n: int) -> None:
    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.8))

    def draw(axis, rows, markersize: float = 2.5):
        m = np.array([row["m"] for row in rows])
        y_states = np.array([row["states_mean"] for row in rows])
        e_states = np.array([row["states_sem"] for row in rows])
        y_part = np.array([row["partitions_mean"] for row in rows])
        e_part = np.array([row["partitions_sem"] for row in rows])
        axis.errorbar(m, y_part, yerr=e_part, fmt="o-", markersize=markersize, capsize=1.5, label="fixed state, random partitions")
        axis.errorbar(m, y_states, yerr=e_states, fmt="s--", markersize=markersize, capsize=1.5, label="fixed partition, random states")

    draw(ax, summary)
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel(r"entropy $S_{N,M}$")
    ax.set_xlim(0, 1 << n)
    ax.set_ylim(0, n / 2)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.38), ncol=1, frameon=False)

    inset = make_inset_axis(ax)
    draw(inset, zoom_summary, markersize=2.0)
    inset.set_xlim(0, max(row["m"] for row in zoom_summary) * 1.02)
    inset.set_ylim(max(0, n / 2 - 3), n / 2)
    inset.set_xlabel(r"$M$", labelpad=1)
    inset.set_ylabel(r"$S$", labelpad=1)
    inset.legend().remove()

    save_figure(fig, out_pdf)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=14)
    parser.add_argument("--points", type=int, default=50)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--zoom-points", type=int, default=30)
    parser.add_argument("--zoom-samples", type=int, default=500)
    parser.add_argument("--zoom-stop", type=int, default=1600)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig1")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    N = 1 << args.n
    zoom_stop = min(args.zoom_stop, N)
    if args.zoom_stop > N:
        print(f"warning: --zoom-stop={args.zoom_stop} exceeds N=2**n={N}; using --zoom-stop={zoom_stop}.")

    m_values = m_grid(1, N, args.points, include_stop=True)
    zoom_values = m_grid(1, zoom_stop, args.zoom_points, include_stop=True)

    summary, raw = compute_curve(args.n, m_values, args.samples, args.seed)
    zoom_summary, zoom_raw = compute_curve(args.n, zoom_values, args.zoom_samples, args.seed + 1)

    fields = list(summary[0].keys())
    write_rows(args.outdir / "fig1_summary.csv", summary, fields)
    write_rows(args.outdir / "fig1_zoom_summary.csv", zoom_summary, fields)
    write_rows(args.outdir / "fig1_raw.csv", raw, list(raw[0].keys()))
    write_rows(args.outdir / "fig1_zoom_raw.csv", zoom_raw, list(zoom_raw[0].keys()))
    plot(summary, zoom_summary, args.outdir / "fig1_concentration.pdf", args.n)


if __name__ == "__main__":
    main()
