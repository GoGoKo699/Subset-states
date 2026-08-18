#!/usr/bin/env python3
"""Local verification of the retained Table-I peak values.

This is a neighbourhood check, not a reconstruction of the historical global
search.  It samples support sizes near each recorded peak and fits a quadratic
in log2(M).
"""
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

from subset_states.core import summary_stats
from subset_states.experiments import random_subset_entropy_samples, write_rows
from subset_states.peak_fitting import local_m_grid, quadratic_peak_fit
from subset_states.plotting import apply_journal_style, save_figure
from subset_states.tables import read_table_i, select_table_rows


def read_schedule(path: Path) -> dict[int, dict[str, float | int]]:
    with path.open(newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    return {
        int(row["n"]): {
            "samples": int(row["samples"]),
            "points": int(row["points"]),
            "span": float(row["span"]),
        }
        for row in reader
    }


def plot_one(n: int, rows: list[dict], fit: dict, out_pdf: Path) -> None:
    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    m = np.asarray([row["m"] for row in rows], dtype=float)
    x = np.log2(m)
    mean = np.asarray([row["mean"] for row in rows], dtype=float)
    sem = np.asarray([row["sem"] for row in rows], dtype=float)
    ax.errorbar(x, mean, yerr=sem, fmt="o", markersize=3, capsize=1.5, label="sampled mean")
    xx = np.linspace(x.min(), x.max(), 300)
    yy = fit["quadratic_a"] * xx**2 + fit["quadratic_b"] * xx + fit["quadratic_c"]
    ax.plot(xx, yy, "--", linewidth=1, label="quadratic fit")
    ax.axvline(np.log2(float(fit["table_M_n"])), linestyle=":", linewidth=1, label=r"Table $M_n$")
    ax.axvline(float(fit["log2_M_peak"]), linestyle="-.", linewidth=1, label="fit maximum")
    ax.set_xlabel(r"$\log_2 M$")
    ax.set_ylabel(r"$S_{N,M}$")
    ax.set_title(fr"local peak check, $n={n}$")
    ax.legend(frameon=False)
    save_figure(fig, out_pdf)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=Path, default=ROOT / "data" / "table_i_peaks.csv")
    parser.add_argument("--schedule", type=Path, default=ROOT / "data" / "peak_verification_schedule.csv")
    parser.add_argument("--n-values", type=int, nargs="+", default=None)
    parser.add_argument("--all-table", action="store_true")
    parser.add_argument("--default-max-n", type=int, default=20)
    parser.add_argument("--samples", type=int, default=None)
    parser.add_argument("--points", type=int, default=None)
    parser.add_argument("--span", type=float, default=None)
    parser.add_argument("--seed", type=int, default=20250605)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "peak_verification")
    args = parser.parse_args()

    table = read_table_i(args.table)
    if args.all_table:
        n_values = [int(row["n"]) for row in table]
    elif args.n_values:
        n_values = args.n_values
    else:
        n_values = [int(row["n"]) for row in table if int(row["n"]) <= args.default_max_n]
    selected = select_table_rows(table, n_values)
    schedule = read_schedule(args.schedule)
    args.outdir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for index, row in enumerate(selected):
        n = int(row["n"])
        config = schedule.get(n, {"samples": 50, "points": 15, "span": 0.25})
        jobs.append(
            {
                "n": n,
                "table_M_n": int(row["M_n"]),
                "table_S_n": float(row["S_n"]),
                "samples": args.samples or int(config["samples"]),
                "points": args.points or int(config["points"]),
                "span": args.span if args.span is not None else float(config["span"]),
                "seed": args.seed + 10000 * index,
            }
        )
    for job in jobs:
        print(job)
    if args.dry_run:
        return

    summary = []
    for job in jobs:
        grid = local_m_grid(job["n"], job["table_M_n"], points=job["points"], span=job["span"])
        rows = []
        for offset, m in enumerate(grid):
            vals = random_subset_entropy_samples(job["n"], int(m), job["samples"], seed=job["seed"] + offset)
            stats = summary_stats(vals)
            rows.append({
                "n": job["n"], "m": int(m), "log2_m": float(np.log2(m)),
                "table_M_n": job["table_M_n"], "table_S_n": job["table_S_n"],
                "count": stats.count, "mean": stats.mean, "std": stats.std,
                "sem": stats.sem, "minimum": stats.minimum, "maximum": stats.maximum,
            })
        fit = quadratic_peak_fit(
            np.asarray([row["m"] for row in rows]),
            np.asarray([row["mean"] for row in rows]),
            np.asarray([row["sem"] for row in rows]),
        )
        fit.update(job)
        fit["M_peak_minus_table_M_n"] = float(fit["M_peak"] - job["table_M_n"])
        fit["S_peak_minus_table_S_n"] = float(fit["S_peak"] - job["table_S_n"])
        write_rows(args.outdir / f"peak_verification_n{job['n']}_samples.csv", rows, rows[0].keys())
        write_rows(args.outdir / f"peak_verification_n{job['n']}_fit.csv", [fit], fit.keys())
        plot_one(job["n"], rows, fit, args.outdir / f"peak_verification_n{job['n']}.pdf")
        summary.append(fit)
    write_rows(args.outdir / "peak_verification_summary.csv", summary, summary[0].keys())


if __name__ == "__main__":
    main()
