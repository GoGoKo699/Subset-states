#!/usr/bin/env python3
"""Convenience driver for publication figures.

Default mode regenerates the fast deterministic figures (2 and 5) and verifies
that all released vector figures are present.  ``--full`` launches every
sampling script with its publication defaults and can be computationally costly.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

FIGURE_FILES = {
    1: "fig1/fig1_concentration.pdf",
    2: "fig2/fig2_peak_scaling.pdf",
    3: "fig3/fig3_spectral_bulk.pdf",
    4: "fig4/fig4_approximation.pdf",
    5: "fig5/fig5_qft_residue_controls.pdf",
    6: "fig6/fig6_renyi.pdf",
    7: "fig7/fig7_partitions.pdf",
}


def run(script: str, *args: str) -> None:
    command = [sys.executable, str(ROOT / "scripts" / script), *args]
    print("+", " ".join(command))
    subprocess.run(command, check=True, cwd=ROOT)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true", help="run every publication computation; expensive")
    parser.add_argument("--smoke", action="store_true", help="run reduced-size computation checks")
    args = parser.parse_args()

    if args.full:
        for script in [
            "fig1_concentration.py", "fig2_peak_scaling.py", "fig3_spectral.py",
            "fig4_approximation.py", "regenerate_fig5_baseline.py",
            "run_residue_controls.py", "fig5_qft_residue_controls.py",
            "fig6_renyi.py", "fig7_partitions.py",
        ]:
            run(script)
    elif args.smoke:
        run("fig1_concentration.py", "--n", "6", "--points", "5", "--samples", "2", "--zoom-points", "4", "--zoom-samples", "2", "--zoom-stop", "32", "--outdir", str(ROOT / "outputs" / "smoke" / "fig1"))
        run("fig2_peak_scaling.py")
        run("fig3_spectral.py", "--n", "8", "--m", "40", "--bins", "20", "--outdir", str(ROOT / "outputs" / "smoke" / "fig3"))
        run("fig4_approximation.py", "--n", "6", "--points", "5", "--samples", "2", "--zoom-stop", "32", "--zoom-points", "4", "--zoom-samples", "2", "--outdir", str(ROOT / "outputs" / "smoke" / "fig4"))
        run("fig5_qft_residue_controls.py")
        run("fig6_renyi.py", "--n", "6", "--points", "5", "--samples", "2", "--zoom-stop", "32", "--zoom-points", "4", "--zoom-samples", "2", "--outdir", str(ROOT / "outputs" / "smoke" / "fig6"))
        run("fig7_partitions.py", "--n", "6", "--m", "20", "--samples", "3", "--outdir", str(ROOT / "outputs" / "smoke" / "fig7"))
    else:
        run("fig2_peak_scaling.py")
        run("fig5_qft_residue_controls.py")

    missing = [str(ROOT / "outputs" / rel) for rel in FIGURE_FILES.values() if not (ROOT / "outputs" / rel).exists()]
    if missing:
        raise SystemExit("missing released figure files:\n" + "\n".join(missing))
    print("All seven released figure PDFs are present.")


if __name__ == "__main__":
    main()
