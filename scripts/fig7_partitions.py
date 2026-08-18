#!/usr/bin/env python3
"""Figure 7: balanced-cut entropy distributions on one common entropy scale."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.core import exact_page_entropy_bits, summary_stats
from subset_states.csv_plotting import plot_fig7
from subset_states.experiments import (
    complex_haar_state,
    partition_entropy_samples_for_fixed_subset,
    partition_entropy_samples_for_state,
    write_rows,
)
from subset_states.tables import PEAK_M


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=20)
    parser.add_argument("--m", type=int, default=None)
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--bins", type=int, default=100)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig7")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    m = args.m if args.m is not None else PEAK_M[args.n]

    subset_values = partition_entropy_samples_for_fixed_subset(args.n, m, args.samples, seed=args.seed)
    haar = complex_haar_state(args.n, seed=args.seed + 1)
    haar_values = partition_entropy_samples_for_state(args.n, haar, args.samples, seed=args.seed + 2)
    rows = [
        *[
            {"state_family": "subset", "sample": idx, "entropy": float(value)}
            for idx, value in enumerate(subset_values)
        ],
        *[
            {"state_family": "complex_haar", "sample": idx, "entropy": float(value)}
            for idx, value in enumerate(haar_values)
        ],
    ]
    write_rows(args.outdir / "fig7_partition_samples.csv", rows, rows[0].keys())
    subset_stats, haar_stats = summary_stats(subset_values), summary_stats(haar_values)
    page = exact_page_entropy_bits(args.n)
    summary_rows = [
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
            "Page_exact_bits": page,
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
            "Page_exact_bits": page,
        },
    ]
    write_rows(args.outdir / "fig7_partition_summary.csv", summary_rows, summary_rows[0].keys())
    plot_fig7(args.outdir, args.outdir / "fig7_partitions.pdf", n=args.n, bins=args.bins)


if __name__ == "__main__":
    main()
