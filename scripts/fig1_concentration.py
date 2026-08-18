#!/usr/bin/env python3
"""Figure 1: concentration over supports and balanced bipartitions."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from subset_states.core import entropy_from_support, random_balanced_partition, random_subset, summary_stats
from subset_states.csv_plotting import plot_fig1
from subset_states.experiments import m_grid, write_rows


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
        st, pt = summary_stats(state_values), summary_stats(partition_values)
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
    summary, raw = compute_curve(args.n, m_grid(1, N, args.points), args.samples, args.seed)
    zoom_summary, zoom_raw = compute_curve(
        args.n,
        m_grid(1, zoom_stop, args.zoom_points),
        args.zoom_samples,
        args.seed + 1,
    )
    fields = list(summary[0].keys())
    write_rows(args.outdir / "fig1_summary.csv", summary, fields)
    write_rows(args.outdir / "fig1_zoom_summary.csv", zoom_summary, fields)
    write_rows(args.outdir / "fig1_raw.csv", raw, list(raw[0].keys()))
    write_rows(args.outdir / "fig1_zoom_raw.csv", zoom_raw, list(zoom_raw[0].keys()))
    plot_fig1(args.outdir, args.outdir / "fig1_concentration.pdf", n=args.n)


if __name__ == "__main__":
    main()
