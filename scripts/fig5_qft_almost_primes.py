#!/usr/bin/env python3
"""Figure 5: position/QFT entanglement and almost-prime union states."""
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

from subset_states.core import entropy_from_state_vector, entropy_from_support, qft_state_from_support
from subset_states.experiments import almost_prime_union_supports, m_grid, random_subset_qft_samples, write_rows
from subset_states.core import summary_stats
from subset_states.plotting import apply_journal_style, save_figure


def random_curve(n: int, m_values: np.ndarray, samples: int, seed: int) -> list[dict]:
    rows: list[dict] = []
    for offset, m in enumerate(m_values):
        position, fourier = random_subset_qft_samples(n, int(m), samples, seed=seed + offset)
        pos = summary_stats(position)
        qft = summary_stats(fourier)
        rows.append(
            {
                "m": int(m),
                "position_mean": pos.mean,
                "position_std": pos.std,
                "position_sem": pos.sem,
                "qft_mean": qft.mean,
                "qft_std": qft.std,
                "qft_sem": qft.sem,
                "delta_qft_minus_position": qft.mean - pos.mean,
            }
        )
    return rows


def almost_prime_rows(n: int) -> list[dict]:
    rows: list[dict] = []
    for k, support in almost_prime_union_supports(n):
        position = entropy_from_support(n, support)
        qft = entropy_from_state_vector(n, qft_state_from_support(n, support))
        rows.append(
            {
                "k": k,
                "m": int(support.size),
                "position_entropy": position,
                "qft_entropy": qft,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=14)
    parser.add_argument("--points", type=int, default=100)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20250604)
    parser.add_argument("--outdir", type=Path, default=ROOT / "outputs" / "fig5")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    random_rows = random_curve(args.n, m_grid(1, 1 << args.n, args.points), args.samples, args.seed)
    union_rows = almost_prime_rows(args.n)
    write_rows(args.outdir / "fig5_random_qft_summary.csv", random_rows, random_rows[0].keys())
    write_rows(args.outdir / "fig5_almost_prime_unions.csv", union_rows, union_rows[0].keys())

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    m = np.array([row["m"] for row in random_rows])
    pos = np.array([row["position_mean"] for row in random_rows])
    pos_sem = np.array([row["position_sem"] for row in random_rows])
    qft = np.array([row["qft_mean"] for row in random_rows])
    qft_sem = np.array([row["qft_sem"] for row in random_rows])
    ax.errorbar(m, pos, yerr=pos_sem, fmt="o", markersize=2, capsize=1.2, label=r"random subset")
    ax.errorbar(m, qft, yerr=qft_sem, fmt=".", markersize=2, capsize=1.2, label=r"random subset after QFT")
    union_m = np.array([row["m"] for row in union_rows])
    union_pos = np.array([row["position_entropy"] for row in union_rows])
    union_qft = np.array([row["qft_entropy"] for row in union_rows])
    ax.plot(union_m, union_pos, "+", markersize=6, label=r"$U_{N,k}$")
    ax.plot(union_m, union_qft, "x", markersize=5, label=r"QFT $U_{N,k}$")
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel("entropy")
    ax.set_xlim(0, 1 << args.n)
    ax.set_ylim(0, args.n / 2)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.45), ncol=2)
    save_figure(fig, args.outdir / "fig5_qft_almost_primes.pdf")

    last = union_rows[-1]
    print(
        "Largest almost-prime union: "
        f"k={last['k']}, M={last['m']} of N={1 << args.n}; "
        "0 and 1 are excluded under the standard Ω(k) definition."
    )


if __name__ == "__main__":
    main()
