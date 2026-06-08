#!/usr/bin/env python3
"""Redraw Figure 2 from fig2_table_with_page.csv."""
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
    out = args.out if args.out is not None else args.data_dir / "fig2_peak_scaling_from_csv.pdf"
    plot_fig2(args.data_dir, out)


if __name__ == "__main__":
    main()
