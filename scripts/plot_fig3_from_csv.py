#!/usr/bin/env python3
"""Redraw Figure 3 from fig3_spectrum.csv."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.csv_plotting import plot_fig3


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / "outputs" / "fig3")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--bins", type=int, default=500)
    args = parser.parse_args()
    out = args.out if args.out is not None else args.data_dir / "fig3_spectral_bulk_from_csv.pdf"
    plot_fig3(args.data_dir, out, bins=args.bins)


if __name__ == "__main__":
    main()
