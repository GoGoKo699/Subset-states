#!/usr/bin/env python3
"""Redraw revised Figure 2 from existing CSV files."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.csv_plotting import plot_fig2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / "outputs" / "fig2")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    out = args.out if args.out is not None else args.data_dir / "fig2_peak_scaling.pdf"
    plot_fig2(args.data_dir, out)


if __name__ == "__main__":
    main()
