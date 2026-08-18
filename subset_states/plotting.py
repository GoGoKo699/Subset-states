from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from cycler import cycler

# Cool, non-warm palette retained from the publication repository.
COOL_PALETTE = [
    "#25639C",  # blue
    "#008C95",  # teal
    "#6D5ACF",  # violet
    "#4B5563",  # slate
    "#4BA3C7",  # cyan-blue
    "#2E8B8B",  # blue-green
    "#415A77",  # steel/navy
    "#7C6BC4",  # muted purple
]


def make_inset_axis(ax, bounds=(0.57, 0.55, 0.39, 0.39)):
    """Create an inset axis using Matplotlib's built-in API."""

    return ax.inset_axes(bounds)


def apply_journal_style() -> None:
    """Conservative defaults suitable for the manuscript figures."""

    plt.rcParams.update(
        {
            "figure.figsize": (4.5, 3.6),
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.8,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "lines.linewidth": 1.15,
            "axes.prop_cycle": cycler(color=COOL_PALETTE),
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.bbox": "tight",
        }
    )


def save_figure(fig, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    if path.suffix.lower() == ".pdf":
        fig.savefig(path.with_suffix(".png"), dpi=300)
    plt.close(fig)
