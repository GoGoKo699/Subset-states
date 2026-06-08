#!/usr/bin/env python3
"""Figure 8 and greedy search for high-entanglement subset samples."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from subset_states.core import balanced_bipartitions, entropy_from_support, random_subset, summary_stats
from subset_states.experiments import random_subset_entropy_samples, write_rows
from subset_states.plotting import apply_journal_style, save_figure
from subset_states.tables import PEAK_M


def parse_seed_map(text: str) -> dict[int, int]:
    out: dict[int, int] = {}
    for item in text.split(","):
        if not item.strip():
            continue
        n_text, seed_text = item.split(":")
        out[int(n_text)] = int(seed_text)
    return out


def all_partition_distribution(n: int, m: int, subset_seed: int) -> np.ndarray:
    support = random_subset(n, m, np.random.default_rng(subset_seed))
    partitions = list(balanced_bipartitions(n))
    values = np.empty(len(partitions), dtype=float)
    for idx, right_bits in enumerate(partitions):
        values[idx] = entropy_from_support(n, support, right_bits)
    return values


def greedy_search(n: int, m: int, seed_start: int, seed_stop: int, outdir: Path) -> None:
    rows = []
    best_seed = None
    best_mean = -np.inf
    for seed in range(seed_start, seed_stop):
        values = all_partition_distribution(n, m, seed)
        stats = summary_stats(values)
        rows.append(
            {
                "seed": seed,
                "n": n,
                "M": m,
                "partition_count": stats.count,
                "mean": stats.mean,
                "std": stats.std,
                "minimum": stats.minimum,
                "maximum": stats.maximum,
            }
        )
        if stats.mean > best_mean:
            best_mean = stats.mean
            best_seed = seed
    write_rows(outdir / f"greedy_search_n{n}_{seed_start}_{seed_stop}.csv", rows, rows[0].keys())
    print(f"best seed in [{seed_start}, {seed_stop}): {best_seed}, mean={best_mean:.8f}")


def make_plot(n_values: list[int], seed_map: dict[int, int], background_samples: int, seed: int, outdir: Path) -> None:
    apply_journal_style()
    fig, axes = plt.subplots(len(n_values), 1, figsize=(4.8, 3.0 * len(n_values)))
    if len(n_values) == 1:
        axes = [axes]

    for panel, (ax, n) in enumerate(zip(axes, n_values)):
        m = PEAK_M[n]
        background = random_subset_entropy_samples(n, m, background_samples, seed=seed + panel)
        selected_seed = seed_map.get(n, seed + 10_000 + panel)
        selected = all_partition_distribution(n, m, selected_seed)

        raw_rows = []
        for idx, value in enumerate(background):
            raw_rows.append({"n": n, "M": m, "family": "random_subsets_natural_partition", "sample": idx, "entropy": float(value)})
        for idx, value in enumerate(selected):
            raw_rows.append({"n": n, "M": m, "family": f"selected_seed_{selected_seed}_all_partitions", "sample": idx, "entropy": float(value)})
        write_rows(outdir / f"fig8_samples_n{n}.csv", raw_rows, raw_rows[0].keys())

        lo = min(background.min(), selected.min())
        hi = max(background.max(), selected.max())
        bins = np.linspace(lo, hi, 120)
        ax.hist(background, bins=bins, alpha=0.65, label="random subsets, natural partition")
        ax.hist(selected, bins=bins, alpha=0.65, label=f"selected subset, all partitions")
        ax.axvline(background.mean(), linestyle="--", linewidth=1)
        ax.axvline(selected.mean(), linestyle=":", linewidth=1)
        ax.set_xlabel("entropy")
        ax.set_ylabel("count")
        ax.set_title(rf"$n={n}$, $M={m}$, selected seed={selected_seed}")
        ax.legend(frameon=False)

        b = summary_stats(background)
        s = summary_stats(selected)
        write_rows(
            outdir / f"fig8_summary_n{n}.csv",
            [
                {"family": "background", "n": n, "M": m, "seed": seed + panel, "count": b.count, "mean": b.mean, "std": b.std, "minimum": b.minimum, "maximum": b.maximum},
                {"family": "selected", "n": n, "M": m, "seed": selected_seed, "count": s.count, "mean": s.mean, "std": s.std, "minimum": s.minimum, "maximum": s.maximum},
            ],
            ["family", "n", "M", "seed", "count", "mean", "std", "minimum", "maximum"],
        )
    fig.subplots_adjust(hspace=0.55)
    save_figure(fig, outdir / "fig8_greedy_partitions.pdf")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["plot", "search"], default="plot")
    parser.add_argument("--n", type=int, default=10, help="n for --mode search")
    parser.add_argument("--m", type=int, default=None, help="M for --mode search; defaults to Table-I peak M")
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--seed-stop", type=int, default=100)
    parser.add_argument("--n-values", type=int, nargs="+", default=[10, 14], help="n values for --mode plot")
    parser.add_argument("--best-seeds", type=str, default="10:163765,14:541")
    parser.add_argument("--background-samples", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig8")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.mode == "search":
        m = args.m if args.m is not None else PEAK_M[args.n]
        greedy_search(args.n, m, args.seed_start, args.seed_stop, args.outdir)
    else:
        make_plot(args.n_values, parse_seed_map(args.best_seeds), args.background_samples, args.seed, args.outdir)


if __name__ == "__main__":
    main()
