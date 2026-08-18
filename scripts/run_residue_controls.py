#!/usr/bin/env python3
"""Generate cardinality- and residue-matched null ensembles for almost-prime supports.

The primary run uses n=14, k=1,2,3, t=0,1,2,3 matched low-order bits,
and paired position/Fourier entropy calculations.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states.core import (
    entropy_from_state_vector,
    entropy_from_support,
    qft_state_from_support,
)
from subset_states.experiments import almost_prime_union_supports


def write_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def shannon_bits(counts: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=np.float64)
    total = counts.sum()
    if total <= 0:
        return 0.0
    p = counts[counts > 0] / total
    return float(-np.sum(p * np.log2(p)))


def residue_counts(support: np.ndarray, modulus: int) -> np.ndarray:
    return np.bincount(support % modulus, minlength=modulus).astype(np.int64)


def sample_matched_support(
    n: int,
    counts: np.ndarray,
    modulus: int,
    rng: np.random.Generator,
) -> np.ndarray:
    N = 1 << n
    parts: list[np.ndarray] = []
    for residue, count in enumerate(counts):
        count = int(count)
        if count == 0:
            continue
        pool = np.arange(residue, N, modulus, dtype=np.int64)
        if count > pool.size:
            raise ValueError(
                f"cannot draw {count} labels from residue {residue} mod {modulus}; "
                f"pool has size {pool.size}"
            )
        parts.append(rng.choice(pool, size=count, replace=False).astype(np.int64))
    if not parts:
        raise ValueError("matched support would be empty")
    return np.sort(np.concatenate(parts))


def summarize(values: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(values, dtype=np.float64)
    std = float(values.std(ddof=1)) if values.size > 1 else 0.0
    q025, q25, q50, q75, q975 = np.quantile(values, [0.025, 0.25, 0.5, 0.75, 0.975])
    return {
        "count": int(values.size),
        "mean": float(values.mean()),
        "std": std,
        "sem": float(std / np.sqrt(values.size)) if values.size > 1 else 0.0,
        "q025": float(q025),
        "q25": float(q25),
        "median": float(q50),
        "q75": float(q75),
        "q975": float(q975),
        "minimum": float(values.min()),
        "maximum": float(values.max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=14)
    parser.add_argument("--k-values", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument("--matched-bits", type=int, nargs="+", default=[0, 1, 2, 3])
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026081701)
    parser.add_argument("--outdir", type=Path, default=ROOT / "data")
    args = parser.parse_args()

    if args.n <= 0 or args.n % 2:
        raise ValueError("n must be a positive even integer")
    if any(t < 0 or t > args.n // 2 for t in args.matched_bits):
        raise ValueError("matched-bits must lie between 0 and n/2")
    if args.samples < 2:
        raise ValueError("samples must be at least 2")

    supports = dict(almost_prime_union_supports(args.n))
    summary_rows: list[dict] = []
    raw_rows: list[dict] = []
    bound_rows: list[dict] = []
    half = args.n // 2

    structured: dict[int, tuple[float, float]] = {}
    for k in args.k_values:
        support = supports[k]
        structured[k] = (
            entropy_from_support(args.n, support),
            entropy_from_state_vector(args.n, qft_state_from_support(args.n, support)),
        )

    for k in args.k_values:
        support = supports[k]
        M = int(support.size)
        structured_position, structured_fourier = structured[k]
        for t in args.matched_bits:
            modulus = 1 << t
            counts = residue_counts(support, modulus)
            H_res = shannon_bits(counts)
            ceiling = float(half - t + H_res)
            forced_gap = float(half - ceiling)
            seed = int(args.seed + 1000 * k + 37 * t)
            rng = np.random.default_rng(seed)
            position = np.empty(args.samples, dtype=np.float64)
            fourier = np.empty(args.samples, dtype=np.float64)

            for sample in range(args.samples):
                matched = sample_matched_support(args.n, counts, modulus, rng)
                position[sample] = entropy_from_support(args.n, matched)
                fourier[sample] = entropy_from_state_vector(
                    args.n, qft_state_from_support(args.n, matched)
                )
                raw_rows.append(
                    {
                        "n": args.n,
                        "k": k,
                        "M": M,
                        "matched_low_bits": t,
                        "modulus": modulus,
                        "sample": sample,
                        "seed": seed,
                        "position_entropy": float(position[sample]),
                        "fourier_entropy": float(fourier[sample]),
                        "fourier_minus_position": float(fourier[sample] - position[sample]),
                    }
                )

            ps = summarize(position)
            fs = summarize(fourier)
            position_residual = float(ps["mean"] - structured_position)
            fourier_residual = float(fs["mean"] - structured_fourier)
            summary_rows.append(
                {
                    "n": args.n,
                    "k": k,
                    "M": M,
                    "matched_low_bits": t,
                    "modulus": modulus,
                    "samples": args.samples,
                    "seed": seed,
                    "residue_counts": json.dumps(counts.tolist(), separators=(",", ":")),
                    "residue_entropy_bits": H_res,
                    "entropy_ceiling_bits": ceiling,
                    "forced_gap_from_balanced_max_bits": forced_gap,
                    "structured_position_entropy": structured_position,
                    "structured_fourier_entropy": structured_fourier,
                    "null_position_mean": ps["mean"],
                    "null_position_std": ps["std"],
                    "null_position_sem": ps["sem"],
                    "null_position_q025": ps["q025"],
                    "null_position_median": ps["median"],
                    "null_position_q975": ps["q975"],
                    "null_position_min": ps["minimum"],
                    "null_position_max": ps["maximum"],
                    "position_residual_deficit": position_residual,
                    "position_standardized_deficit": float(position_residual / ps["std"]) if ps["std"] else float("inf"),
                    "position_null_count_at_or_below_structured": int(np.count_nonzero(position <= structured_position)),
                    "null_fourier_mean": fs["mean"],
                    "null_fourier_std": fs["std"],
                    "null_fourier_sem": fs["sem"],
                    "null_fourier_q025": fs["q025"],
                    "null_fourier_median": fs["median"],
                    "null_fourier_q975": fs["q975"],
                    "null_fourier_min": fs["minimum"],
                    "null_fourier_max": fs["maximum"],
                    "fourier_residual_deficit": fourier_residual,
                    "fourier_standardized_deficit": float(fourier_residual / fs["std"]) if fs["std"] else float("inf"),
                    "fourier_null_count_at_or_below_structured": int(np.count_nonzero(fourier <= structured_fourier)),
                    "paired_null_fourier_minus_position_mean": float(np.mean(fourier - position)),
                    "paired_null_fourier_minus_position_std": float(np.std(fourier - position, ddof=1)),
                }
            )
            bound_rows.append(
                {
                    "n": args.n,
                    "k": k,
                    "M": M,
                    "matched_low_bits": t,
                    "modulus": modulus,
                    "residue_counts": json.dumps(counts.tolist(), separators=(",", ":")),
                    "residue_probabilities": json.dumps((counts / M).tolist(), separators=(",", ":")),
                    "residue_entropy_bits": H_res,
                    "entropy_ceiling_bits": ceiling,
                    "forced_gap_from_balanced_max_bits": forced_gap,
                    "structured_position_entropy": structured_position,
                    "ceiling_minus_structured_entropy": float(ceiling - structured_position),
                }
            )
            print(
                f"k={k}, t={t}, M={M}: "
                f"position residual={position_residual:.6f}, "
                f"Fourier residual={fourier_residual:.6f}"
            )

    # Deficit reductions relative to the cardinality-matched reference.
    for k in args.k_values:
        base = next(
            row for row in summary_rows if row["k"] == k and row["matched_low_bits"] == 0
        )
        pos0 = float(base["position_residual_deficit"])
        four0 = float(base["fourier_residual_deficit"])
        for row in summary_rows:
            if row["k"] != k:
                continue
            row["cardinality_deficit_reduction_position"] = (
                float(1.0 - float(row["position_residual_deficit"]) / pos0)
                if pos0 != 0
                else float("nan")
            )
            row["cardinality_deficit_reduction_fourier"] = (
                float(1.0 - float(row["fourier_residual_deficit"]) / four0)
                if four0 != 0
                else float("nan")
            )

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_rows(args.outdir / "residue_matched_summary.csv", summary_rows)
    write_rows(args.outdir / "residue_matched_samples.csv", raw_rows)
    write_rows(args.outdir / "residue_entropy_bounds.csv", bound_rows)


if __name__ == "__main__":
    main()
