from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from cycler import cycler

# A deliberately cool, non-warm palette: blue, teal, violet, slate, cyan, green-blue.
# Avoids orange/red/yellow so that all default plots have a cool visual tone.
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
    """Create an inset axis without importing mpl_toolkits.

    Some Ubuntu/Python installations mix user-site Matplotlib with the system
    mpl_toolkits package, which can break `mpl_toolkits.axes_grid1.inset_locator`
    with an AttributeError involving matplotlib._docstring.  Using the built-in
    Axes.inset_axes method avoids that dependency entirely.  The coordinates are
    axis-relative: (left, bottom, width, height).
    """

    return ax.inset_axes(bounds)


def apply_journal_style() -> None:
    """Small, conservative figure defaults suitable for journal submission.

    The colour cycle is intentionally cool-toned; it avoids the default orange/red
    entries that can dominate Matplotlib figures.
    """

    plt.rcParams.update(
        {
            "figure.figsize": (4.5, 3.6),
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
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
    # Store a PNG next to the PDF for quick inspection in repositories.
    if path.suffix.lower() == ".pdf":
        fig.savefig(path.with_suffix(".png"), dpi=300)
