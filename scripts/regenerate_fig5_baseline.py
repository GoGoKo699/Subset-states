#!/usr/bin/env python3
"""Regenerate the random position/Fourier baseline and arithmetic points.

The n=14 default is expensive. New CSVs go to generated/full/fig5.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.core import (
    entropy_from_state_vector,
    entropy_from_support,
    qft_state_from_support,
    summary_stats,
)
from subset_states.experiments import almost_prime_union_supports, m_grid, random_subset_qft_samples, write_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=14)
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--outdir", type=Path, default=ROOT / "generated" / "full" / "fig5")
    args = parser.parse_args()
    rows = []
    for offset, m in enumerate(m_grid(1, 1 << args.n, args.points)):
        position, fourier = random_subset_qft_samples(args.n, int(m), args.samples, seed=args.seed + offset)
        ps, fs = summary_stats(position), summary_stats(fourier)
        rows.append({
            "m": int(m), "sample_count": args.samples,
            "position_mean": ps.mean, "position_std": ps.std, "position_sem": ps.sem,
            "qft_mean": fs.mean, "qft_std": fs.std, "qft_sem": fs.sem,
            "delta_qft_minus_position": fs.mean - ps.mean,
        })
        if offset % 10 == 0:
            print(f"baseline {offset + 1}/{args.points}")
    unions = []
    for k, support in almost_prime_union_supports(args.n):
        unions.append({
            "k": k, "m": int(support.size),
            "position_entropy": entropy_from_support(args.n, support),
            "qft_entropy": entropy_from_state_vector(args.n, qft_state_from_support(args.n, support)),
        })
    write_rows(args.outdir / "fig5_random_qft_summary.csv", rows, rows[0].keys())
    write_rows(args.outdir / "fig5_almost_prime_unions.csv", unions, unions[0].keys())


if __name__ == "__main__":
    main()
