#!/usr/bin/env python3
"""Figure 5: cardinality baseline and residue-matched residual deficits."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from subset_states.plotting import apply_journal_style, save_figure, COOL_PALETTE


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def values(rows: list[dict[str, str]], key: str) -> np.ndarray:
    return np.asarray([float(row[key]) for row in rows], dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig5")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    baseline = read_rows(args.data_dir / "fig5_random_qft_summary.csv")
    unions = read_rows(args.data_dir / "fig5_almost_prime_unions.csv")
    controls = read_rows(args.data_dir / "residue_matched_summary.csv")

    apply_journal_style()
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(7.25, 3.35), gridspec_kw={"width_ratios": [1.35, 1]}
    )

    m = values(baseline, "m")
    position_mean = values(baseline, "position_mean")
    position_std = values(baseline, "position_std")
    fourier_mean = values(baseline, "qft_mean")
    fourier_std = values(baseline, "qft_std")
    c0, c1 = COOL_PALETTE[0], COOL_PALETTE[2]
    ax.fill_between(m, position_mean - position_std, position_mean + position_std, color=c0, alpha=0.16, linewidth=0)
    ax.fill_between(m, fourier_mean - fourier_std, fourier_mean + fourier_std, color=c1, alpha=0.11, linewidth=0)
    ax.plot(m, position_mean, color=c0, linewidth=1.25, label="random supports")
    ax.plot(m, fourier_mean, color=c1, linewidth=1.1, linestyle="--", label="random supports after QFT")
    union_m = values(unions, "m")
    ax.plot(union_m, values(unions, "position_entropy"), "+", color=COOL_PALETTE[3], markersize=5.8, markeredgewidth=1.1, label=r"$U_{N,k}$")
    ax.plot(union_m, values(unions, "qft_entropy"), "x", color=COOL_PALETTE[5], markersize=4.8, markeredgewidth=1.0, label=r"QFT $U_{N,k}$")
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel("entropy (bits)")
    ax.set_xlim(0, 2**14)
    ax.set_ylim(0, 7)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2, columnspacing=0.9, handlelength=1.8)
    ax.text(0.02, 0.96, "(a)", transform=ax.transAxes, ha="left", va="top", fontweight="bold")

    labels = ["cardinality", "parity", r"mod $4$", r"mod $8$"]
    x = np.arange(4, dtype=float)
    colors = [COOL_PALETTE[0], COOL_PALETTE[1], COOL_PALETTE[2]]
    for idx, k in enumerate((1, 2, 3)):
        rows = sorted(
            [row for row in controls if int(row["k"]) == k],
            key=lambda row: int(row["matched_low_bits"]),
        )
        position = values(rows, "position_residual_deficit")
        fourier = values(rows, "fourier_residual_deficit")
        ax2.plot(x, position, "o-", color=colors[idx], linewidth=1.25, markersize=4.0)
        ax2.plot(x + 0.045, fourier, "o--", color=colors[idx], markerfacecolor="white", markeredgewidth=1.0, linewidth=1.0, markersize=3.8)
        ax2.text(-0.03, position[0] + (0.035 if k != 2 else -0.055), rf"$k={k}$", color=colors[idx], ha="right", va="center", fontsize=8)
    ax2.axhline(0, color="0.45", linewidth=0.7)
    ax2.set_xticks(x, labels, rotation=22, ha="right")
    ax2.set_ylabel(r"residual deficit $\Delta_t$ (bits)")
    ax2.set_xlabel("constraint retained in null ensemble", labelpad=2)
    ax2.set_xlim(-0.5, 3.35)
    ax2.set_ylim(0, 1.62)
    style_legend = [
        Line2D([0], [0], color="0.2", marker="o", linestyle="-", markersize=4, label="computational basis"),
        Line2D([0], [0], color="0.2", marker="o", markerfacecolor="white", linestyle="--", markersize=4, label="after QFT"),
    ]
    ax2.legend(handles=style_legend, frameon=False, loc="upper right", handlelength=2.0)
    ax2.text(0.02, 0.96, "(b)", transform=ax2.transAxes, ha="left", va="top", fontweight="bold")
    fig.subplots_adjust(left=0.085, right=0.99, top=0.97, bottom=0.26, wspace=0.32)
    save_figure(fig, args.outdir / "fig5_qft_residue_controls.pdf")


if __name__ == "__main__":
    main()
