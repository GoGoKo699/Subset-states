#!/usr/bin/env python3
"""Screen random supports and retain the best mean balanced-cut entropy."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from subset_states.core import balanced_bipartitions, entropy_from_support, random_subset, summary_stats
from subset_states.tables import PEAK_M


def distribution(n: int, m: int, seed: int) -> np.ndarray:
    support = random_subset(n, m, np.random.default_rng(seed))
    return np.asarray([entropy_from_support(n, support, bits) for bits in balanced_bipartitions(n)], dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--m", type=int, default=None)
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--seed-stop", type=int, default=100)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    m = args.m if args.m is not None else PEAK_M[args.n]
    rows = []
    for seed in range(args.seed_start, args.seed_stop):
        stats = summary_stats(distribution(args.n, m, seed))
        rows.append({"seed": seed, "n": args.n, "M": m, "partition_count": stats.count, "mean": stats.mean, "std": stats.std, "minimum": stats.minimum, "maximum": stats.maximum})
    best = max(rows, key=lambda row: row["mean"])
    print(f"best seed in [{args.seed_start}, {args.seed_stop}): {best['seed']}, mean={best['mean']:.8f}")
    out = args.out or Path(f"best_of_random_candidates_n{args.n}_{args.seed_start}_{args.seed_stop}.csv")
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)

if __name__ == "__main__":
    main()
