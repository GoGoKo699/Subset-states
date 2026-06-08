#!/usr/bin/env python3
"""Local verification of the Table-I peak scaling values.

This script is deliberately local. It does not perform a global search over all
support sizes. Instead, for each selected n it samples random subset states in a
window centered at the tabulated peak M_n, fits a quadratic in log2(M), and
reports whether the fitted local maximum is consistent with the table.

The workload is controlled by data/peak_verification_schedule.csv or by command
line overrides such as --samples, --points and --span.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from subset_states.core import summary_stats
from subset_states.experiments import random_subset_entropy_samples, write_rows
from subset_states.peak_fitting import local_m_grid, quadratic_peak_fit
from subset_states.plotting import apply_journal_style, save_figure
from subset_states.tables import read_table_i, select_table_rows


def read_schedule(path: Path | None) -> dict[int, dict[str, float | int]]:
    if path is None:
        return {}
    with path.open(newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        return {}
    reader = csv.DictReader(lines)
    required = {"n", "samples", "points"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")

    schedule: dict[int, dict[str, float | int]] = {}
    for row in reader:
        n = int(row["n"])
        schedule[n] = {
            "samples": int(row["samples"]),
            "points": int(row["points"]),
            "span": float(row.get("span", 0.25) or 0.25),
        }
        if row.get("half_width") not in (None, ""):
            schedule[n]["half_width"] = int(row["half_width"])
    return schedule


def selected_n_values(args: argparse.Namespace, table_rows: list[dict[str, int | float]]) -> list[int]:
    table_n = [int(row["n"]) for row in table_rows]
    if args.all_table:
        return table_n
    if args.n_values:
        return [int(n) for n in args.n_values]
    return [n for n in table_n if n <= args.default_max_n]


def row_parameters(n: int, schedule: dict[int, dict[str, float | int]], args: argparse.Namespace) -> tuple[int, int, float, int | None]:
    row = schedule.get(n, {})
    samples = int(args.samples if args.samples is not None else row.get("samples", 50))
    points = int(args.points if args.points is not None else row.get("points", 15))
    span = float(args.span if args.span is not None else row.get("span", 0.25))
    half_width = args.half_width if args.half_width is not None else row.get("half_width")
    half_width_int = None if half_width is None else int(half_width)
    return samples, points, span, half_width_int


def run_one(
    *,
    n: int,
    table_m: int,
    table_s: float,
    samples: int,
    points: int,
    span: float,
    half_width: int | None,
    seed: int,
    outdir: Path,
    label: str,
) -> dict:
    m_values = local_m_grid(n, table_m, points=points, span=span, half_width=half_width)
    rows: list[dict] = []

    for offset, m in enumerate(m_values):
        values = random_subset_entropy_samples(n, int(m), samples, seed=seed + offset)
        stats = summary_stats(values)
        rows.append(
            {
                "n": n,
                "m": int(m),
                "log2_m": float(np.log2(m)),
                "table_M_n": table_m,
                "table_S_n": table_s,
                "delta_m_from_table": int(m - table_m),
                "count": stats.count,
                "mean": stats.mean,
                "std": stats.std,
                "sem": stats.sem,
                "minimum": stats.minimum,
                "maximum": stats.maximum,
            }
        )

    suffix = f"_{label}" if label else ""
    sample_path = outdir / f"peak_verification_n{n}{suffix}_samples.csv"
    write_rows(sample_path, rows, rows[0].keys())

    m = np.asarray([row["m"] for row in rows], dtype=float)
    mean = np.asarray([row["mean"] for row in rows], dtype=float)
    sem = np.asarray([row["sem"] for row in rows], dtype=float)
    fit = quadratic_peak_fit(m, mean, sem)
    fit.update(
        {
            "n": n,
            "table_M_n": table_m,
            "table_S_n": table_s,
            "M_peak_minus_table_M_n": float(fit["M_peak"] - table_m),
            "S_peak_minus_table_S_n": float(fit["S_peak"] - table_s),
            "window_low": int(m.min()),
            "window_high": int(m.max()),
            "points_requested": points,
            "points_used": int(m.size),
            "samples_per_m": samples,
            "span": span,
            "half_width": "" if half_width is None else half_width,
            "seed": seed,
            "sample_csv": sample_path.name,
        }
    )
    fit_path = outdir / f"peak_verification_n{n}{suffix}_fit.csv"
    write_rows(fit_path, [fit], fit.keys())

    make_one_plot(n=n, rows=rows, fit=fit, out_pdf=outdir / f"peak_verification_n{n}{suffix}.pdf")
    return fit


def make_one_plot(*, n: int, rows: list[dict], fit: dict, out_pdf: Path) -> None:
    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))

    m = np.asarray([row["m"] for row in rows], dtype=float)
    x = np.log2(m)
    mean = np.asarray([row["mean"] for row in rows], dtype=float)
    sem = np.asarray([row["sem"] for row in rows], dtype=float)

    ax.errorbar(x, mean, yerr=sem, fmt="o", markersize=3, capsize=1.5, label="sampled mean")
    x_fit = np.linspace(x.min(), x.max(), 300)
    y_fit = fit["quadratic_a"] * x_fit**2 + fit["quadratic_b"] * x_fit + fit["quadratic_c"]
    ax.plot(x_fit, y_fit, "--", linewidth=1, label="quadratic fit")
    ax.axvline(np.log2(float(fit["table_M_n"])), linestyle=":", linewidth=1, label=r"Table $M_n$")
    ax.axvline(float(fit["log2_M_peak"]), linestyle="-.", linewidth=1, label="fit maximum")
    ax.axhline(float(fit["table_S_n"]), linestyle=":", linewidth=0.8)
    ax.set_xlabel(r"$\log_2 M$")
    ax.set_ylabel(r"$S_{N,M}$")
    ax.set_title(fr"local peak check, $n={n}$")
    ax.legend(frameon=False, loc="best")
    save_figure(fig, out_pdf)


def make_summary_plot(fit_rows: Iterable[dict], out_pdf: Path) -> None:
    rows = list(fit_rows)
    if not rows:
        return
    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    n = np.asarray([float(row["n"]) for row in rows])
    delta_m = np.asarray([float(row["M_peak_minus_table_M_n"]) for row in rows])
    delta_s = np.asarray([float(row["S_peak_minus_table_S_n"]) for row in rows])
    ax.axhline(0.0, linestyle=":", linewidth=1)
    ax.plot(n, delta_s, "o-", label=r"$S_{\rm fit}-S_n$ (bits)")
    ax.set_xlabel(r"number of qubits $n$")
    ax.set_ylabel(r"entropy difference")
    ax2 = ax.twinx()
    ax2.plot(n, delta_m, "s--", label=r"$M_{\rm fit}-M_n$")
    ax2.set_ylabel(r"support-size difference")
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles + handles2, labels + labels2, frameon=False, loc="best")
    save_figure(fig, out_pdf)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=Path, default=ROOT / "data" / "table_i_peaks.csv", help="CSV with columns n,M_n,S_n")
    parser.add_argument("--schedule", type=Path, default=ROOT / "data" / "peak_verification_schedule.csv", help="CSV with columns n,samples,points,span")
    parser.add_argument("--n-values", type=int, nargs="+", default=None, help="n values to verify, e.g. --n-values 10 12 14")
    parser.add_argument("--all-table", action="store_true", help="verify all n values in the table; can be expensive")
    parser.add_argument("--default-max-n", type=int, default=20, help="default n cutoff when neither --n-values nor --all-table is used")
    parser.add_argument("--samples", type=int, default=None, help="override samples per M for all selected n")
    parser.add_argument("--points", type=int, default=None, help="override number of M grid points for all selected n")
    parser.add_argument("--span", type=float, default=None, help="override fractional half-width around table M_n")
    parser.add_argument("--half-width", type=int, default=None, help="override span with an absolute half-width around table M_n")
    parser.add_argument("--seed", type=int, default=20250605)
    parser.add_argument("--label", type=str, default="", help="optional suffix for output files")
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "peak_verification")
    parser.add_argument("--dry-run", action="store_true", help="print selected jobs without sampling")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    table_rows_all = read_table_i(args.table)
    n_values = selected_n_values(args, table_rows_all)
    table_rows = select_table_rows(table_rows_all, n_values)
    schedule = read_schedule(args.schedule)

    print("Local peak-verification jobs:")
    jobs = []
    for index, row in enumerate(table_rows):
        n = int(row["n"])
        table_m = int(row["M_n"])
        table_s = float(row["S_n"])
        samples, points, span, half_width = row_parameters(n, schedule, args)
        seed = args.seed + 10000 * index
        jobs.append((n, table_m, table_s, samples, points, span, half_width, seed))
        width = f"half_width={half_width}" if half_width is not None else f"span={span:g}"
        print(f"  n={n:2d}, table M_n={table_m}, samples={samples}, points={points}, {width}, seed={seed}")

    if args.dry_run:
        return

    fit_rows = []
    for n, table_m, table_s, samples, points, span, half_width, seed in jobs:
        fit_rows.append(
            run_one(
                n=n,
                table_m=table_m,
                table_s=table_s,
                samples=samples,
                points=points,
                span=span,
                half_width=half_width,
                seed=seed,
                outdir=args.outdir,
                label=args.label,
            )
        )
        print(
            f"n={n:2d}: M_fit={fit_rows[-1]['M_peak']:.3f}, "
            f"S_fit={fit_rows[-1]['S_peak']:.6f}, "
            f"status={fit_rows[-1]['fit_status']}"
        )

    suffix = f"_{args.label}" if args.label else ""
    summary_path = args.outdir / f"peak_verification_summary{suffix}.csv"
    write_rows(summary_path, fit_rows, fit_rows[0].keys())
    make_summary_plot(fit_rows, args.outdir / f"peak_verification_summary{suffix}.pdf")
    print("written:", summary_path)


if __name__ == "__main__":
    main()
