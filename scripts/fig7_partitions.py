#!/usr/bin/env python3
"""Figure 7: entropy over random bipartitions for subset vs. complex Haar states."""
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

from subset_states.core import exact_page_entropy_bits, summary_stats
from subset_states.experiments import (
    complex_haar_state,
    partition_entropy_samples_for_fixed_subset,
    partition_entropy_samples_for_state,
    write_rows,
)
from subset_states.plotting import apply_journal_style, save_figure
from subset_states.tables import PEAK_M


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=20)
    parser.add_argument("--m", type=int, default=None)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig7")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    m = args.m if args.m is not None else PEAK_M[args.n]

    subset_values = partition_entropy_samples_for_fixed_subset(args.n, m, args.samples, seed=args.seed)
    haar = complex_haar_state(args.n, seed=args.seed + 1)
    haar_values = partition_entropy_samples_for_state(args.n, haar, args.samples, seed=args.seed + 2)

    rows = []
    for idx, value in enumerate(subset_values):
        rows.append({"state_family": "subset", "sample": idx, "entropy": float(value)})
    for idx, value in enumerate(haar_values):
        rows.append({"state_family": "complex_haar", "sample": idx, "entropy": float(value)})
    write_rows(args.outdir / "fig7_partition_samples.csv", rows, rows[0].keys())

    subset_stats = summary_stats(subset_values)
    haar_stats = summary_stats(haar_values)
    write_rows(
        args.outdir / "fig7_partition_summary.csv",
        [
            {
                "state_family": "subset",
                "n": args.n,
                "M": m,
                "count": subset_stats.count,
                "mean": subset_stats.mean,
                "std": subset_stats.std,
                "sem": subset_stats.sem,
                "minimum": subset_stats.minimum,
                "maximum": subset_stats.maximum,
                "Page_exact_bits": exact_page_entropy_bits(args.n),
            },
            {
                "state_family": "complex_haar",
                "n": args.n,
                "M": "",
                "count": haar_stats.count,
                "mean": haar_stats.mean,
                "std": haar_stats.std,
                "sem": haar_stats.sem,
                "minimum": haar_stats.minimum,
                "maximum": haar_stats.maximum,
                "Page_exact_bits": exact_page_entropy_bits(args.n),
            },
        ],
        ["state_family", "n", "M", "count", "mean", "std", "sem", "minimum", "maximum", "Page_exact_bits"],
    )

    apply_journal_style()
    fig, axes = plt.subplots(2, 1, figsize=(4.7, 5.8))
    for ax, values, label in [
        (axes[0], subset_values, rf"subset, $M={m}$"),
        (axes[1], haar_values, "complex Haar"),
    ]:
        ax.hist(values, bins=60)
        ax.axvline(values.mean(), linestyle="--", linewidth=1, label="sample mean")
        ax.axvline(exact_page_entropy_bits(args.n), linestyle=":", linewidth=1, label="Page average")
        ax.set_xlabel("entropy")
        ax.set_ylabel("count")
        ax.set_title(label)
        ax.legend(frameon=False)
    fig.subplots_adjust(hspace=0.55)
    save_figure(fig, args.outdir / "fig7_partitions.pdf")


if __name__ == "__main__":
    main()
