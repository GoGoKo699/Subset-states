#!/usr/bin/env python3
"""Redraw Figure 7 from fig7_partition_samples.csv."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.csv_plotting import plot_fig7


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / "outputs" / "fig7")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--n", type=int, default=None)
    parser.add_argument("--bins", type=int, default=60)
    args = parser.parse_args()
    out = args.out if args.out is not None else args.data_dir / "fig7_partitions_from_csv.pdf"
    plot_fig7(args.data_dir, out, n=args.n, bins=args.bins)


if __name__ == "__main__":
    main()
