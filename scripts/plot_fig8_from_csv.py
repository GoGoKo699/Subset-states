#!/usr/bin/env python3
"""Redraw Figure 8 from fig8_samples_n*.csv."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.csv_plotting import plot_fig8


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / "outputs" / "fig8")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--n-values", type=int, nargs="*", default=None)
    parser.add_argument("--bins", type=int, default=120)
    args = parser.parse_args()
    out = args.out if args.out is not None else args.data_dir / "fig8_greedy_partitions_from_csv.pdf"
    plot_fig8(args.data_dir, out, n_values=args.n_values, bins=args.bins)


if __name__ == "__main__":
    main()
