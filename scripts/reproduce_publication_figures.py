#!/usr/bin/env python3
"""Redraw all seven figures from released CSV, or opt into new computations.

Default: generated/figures.  --smoke: small new runs in generated/smoke.
--full: original sampling defaults in generated/full (potentially expensive).
No mode writes into the released data/ or outputs/ trees.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

FIGURE_FILES = {
    1: "fig1/fig1_concentration.png",
    2: "fig2/fig2_peak_scaling.png",
    3: "fig3/fig3_spectral_bulk.png",
    4: "fig4/fig4_approximation.png",
    5: "fig5/fig5_qft_residue_controls.png",
    6: "fig6/fig6_renyi.png",
    7: "fig7/fig7_partitions.png",
}
CSV_INPUTS = {
    1: ("fig1_summary.csv", "fig1_zoom_summary.csv"),
    2: ("fig2_table_with_page.csv",),
    3: ("fig3_spectrum.csv", "fig3_spectrum_summary.csv"),
    4: ("fig4_summary.csv", "fig4_zoom_summary.csv"),
    5: ("fig5_random_qft_summary.csv", "fig5_almost_prime_unions.csv", "residue_matched_summary.csv"),
    6: ("fig6_renyi_summary.csv", "fig6_renyi_zoom_summary.csv"),
    7: ("fig7_partition_samples.csv", "fig7_partition_summary.csv"),
}


def validate_csv_inputs(released_root: Path, data_root: Path) -> None:
    """Fail before any drawing if even one required released table is absent."""
    missing = []
    for number, names in CSV_INPUTS.items():
        folder = data_root if number == 5 else released_root / f"fig{number}"
        missing.extend(str(folder / name) for name in names if not (folder / name).is_file())
    if missing:
        raise FileNotFoundError("missing released CSV inputs:\n" + "\n".join(missing))


def validate_output_root(outdir: Path, *protected: Path) -> None:
    """Reject destinations that could overwrite released inputs or figures."""
    target = outdir.resolve()
    for folder in protected:
        source = folder.resolve()
        if target == source or source in target.parents or target in source.parents:
            raise ValueError(f"output directory must be separate from released evidence: {folder}")


def build_commands(
    mode: str,
    outdir: Path,
    *,
    released_root: Path = ROOT / "outputs",
    data_root: Path = ROOT / "data",
) -> list[list[str]]:
    """Build explicit input/output commands so no child relies on output defaults."""
    commands: list[list[str]] = []

    def add(script: str, *args: str | Path) -> None:
        commands.append([sys.executable, str(ROOT / "scripts" / script), *map(str, args)])

    if mode == "redraw":
        for number, filename in FIGURE_FILES.items():
            if number == 5:
                add("fig5_qft_residue_controls.py", "--data-dir", data_root, "--outdir", outdir / "fig5")
            else:
                add(f"plot_fig{number}_from_csv.py", "--data-dir", released_root / f"fig{number}", "--out", outdir / filename)
        return commands
    if mode not in {"smoke", "full"}:
        raise ValueError(f"unknown figure mode: {mode}")

    small_grid = ["--n", "6", "--points", "5", "--samples", "2", "--zoom-points", "4", "--zoom-samples", "2", "--zoom-stop", "32"] if mode == "smoke" else []
    add("fig1_concentration.py", *small_grid, "--outdir", outdir / "fig1")
    add("fig2_peak_scaling.py", "--table", data_root / "table_i_peaks.csv", "--outdir", outdir / "fig2")
    spectrum = ["--n", "8", "--m", "40", "--bins", "20"] if mode == "smoke" else []
    add("fig3_spectral.py", *spectrum, "--outdir", outdir / "fig3")
    add("fig4_approximation.py", *small_grid, "--outdir", outdir / "fig4")
    baseline = ["--n", "6", "--points", "5", "--samples", "2"] if mode == "smoke" else []
    nulls = ["--n", "6", "--samples", "2"] if mode == "smoke" else []
    add("regenerate_fig5_baseline.py", *baseline, "--outdir", outdir / "fig5")
    add("run_residue_controls.py", *nulls, "--outdir", outdir / "fig5")
    add("fig5_qft_residue_controls.py", "--data-dir", outdir / "fig5", "--outdir", outdir / "fig5")
    add("fig6_renyi.py", *small_grid, "--outdir", outdir / "fig6")
    partitions = ["--n", "6", "--m", "20", "--samples", "3"] if mode == "smoke" else []
    add("fig7_partitions.py", *partitions, "--outdir", outdir / "fig7")
    return commands


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--full", action="store_true", help="run every figure computation with original defaults; expensive")
    modes.add_argument("--smoke", action="store_true", help="run reduced-size computations, including Figure 5 nulls")
    parser.add_argument("--outdir", type=Path, help="generated destination; must be separate from released evidence")
    parser.add_argument("--released-dir", type=Path, default=ROOT / "outputs", help="released figure CSV tree")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data", help="released Figure 5 CSV and peak-table inputs")
    args = parser.parse_args()
    args.released_dir = args.released_dir.resolve()
    args.data_dir = args.data_dir.resolve()
    mode = "full" if args.full else "smoke" if args.smoke else "redraw"
    outdir = (args.outdir or ROOT / "generated" / ("figures" if mode == "redraw" else mode)).resolve()
    try:
        validate_output_root(outdir, ROOT / "outputs", ROOT / "data", args.released_dir, args.data_dir)
        if mode == "redraw":
            validate_csv_inputs(args.released_dir, args.data_dir)
        elif not (args.data_dir / "table_i_peaks.csv").is_file():
            raise FileNotFoundError(f"missing peak-table input: {args.data_dir / 'table_i_peaks.csv'}")
    except (ValueError, FileNotFoundError) as exc:
        parser.error(str(exc))

    for command in build_commands(mode, outdir, released_root=args.released_dir, data_root=args.data_dir):
        print("+", " ".join(command), flush=True)
        subprocess.run(command, check=True, cwd=ROOT)
    missing = [str(outdir / rel) for rel in FIGURE_FILES.values() if not (outdir / rel).is_file()]
    if missing:
        raise SystemExit("missing generated figure files:\n" + "\n".join(missing))
    print(f"Generated all seven PNG figures in {outdir} ({mode}).")


if __name__ == "__main__":
    main()
