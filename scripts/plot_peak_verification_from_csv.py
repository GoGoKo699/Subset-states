#!/usr/bin/env python3
"""Redraw one local peak-verification plot from released CSV files."""
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
    parser.add_argument("--data-dir", type=Path, default=ROOT / "outputs" / "peak_verification")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    samples = args.data_dir / f"peak_verification_n{args.n}_samples.csv"
    fit = args.data_dir / f"peak_verification_n{args.n}_fit.csv"
    out = args.out or args.data_dir / f"peak_verification_n{args.n}_from_csv.pdf"
    plot_peak_search(samples, fit, out)

if __name__ == "__main__":
    main()
