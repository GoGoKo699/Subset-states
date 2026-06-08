#!/usr/bin/env python3
"""Redraw a local peak-verification plot from CSV files."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.csv_plotting import plot_peak_search


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--label", type=str, default="")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "outputs" / "peak_verification")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    suffix = f"_{args.label}" if args.label else ""
    samples = args.data_dir / f"peak_verification_n{args.n}{suffix}_samples.csv"
    fit = args.data_dir / f"peak_verification_n{args.n}{suffix}_fit.csv"
    out = args.out or (args.data_dir / f"peak_verification_n{args.n}{suffix}_from_csv.pdf")
    plot_peak_search(samples, fit, out)


if __name__ == "__main__":
    main()
