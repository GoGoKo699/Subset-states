#!/usr/bin/env python3
"""Redraw Figure 5 from previously generated CSV files."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.csv_plotting import plot_fig5


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / "outputs" / "fig5")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--n", type=int, default=None, help="Optional; inferred from fig5_random_qft_summary.csv if omitted.")
    args = parser.parse_args()
    out = args.out if args.out is not None else args.data_dir / "fig5_qft_almost_primes_from_csv.pdf"
    plot_fig5(args.data_dir, out, n=args.n)


if __name__ == "__main__":
    main()
