#!/usr/bin/env python3
"""Check retained residue-control evidence without rewriting it.

This validates the recorded n=14, k=1..3, t=0..3, 1000-sample experiment.
All raw-row statistics are checked; only a stated prefix is replayed from seeds.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subset_states import core
from subset_states.experiments import almost_prime_union_supports
from scripts.run_residue_controls import sample_matched_support


class ValidationError(ValueError):
    """An explicit acceptance failure, including when Python uses -O."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def close(actual: float, expected: float, label: str, *, atol: float = 2e-11) -> float:
    require(math.isfinite(actual) and math.isfinite(expected), f"{label}: non-finite value")
    error = abs(actual - expected)
    require(error <= atol, f"{label}: {actual:.16g} != {expected:.16g} (error {error:.3g})")
    return error


def read_rows(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        names = reader.fieldnames or []
        require(len(names) == len(set(names)), f"{path}: duplicate columns")
        require(required <= set(names), f"{path}: missing columns {sorted(required - set(names))}")
        rows = list(reader)
    require(bool(rows), f"{path}: no data rows")
    for index, row in enumerate(rows, 2):
        require(None not in row and all(value is not None for value in row.values()),
                f"{path}:{index}: malformed CSV row")
        for field in required:
            require(row[field] != "", f"{path}:{index}: empty {field}")
    return rows


def finite(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    require(math.isfinite(value), f"{field}: non-finite value")
    return value


def shannon(counts: np.ndarray) -> float:
    p = counts[counts > 0] / counts.sum()
    return float(-np.sum(p * np.log2(p)))


def independent_entropy(n: int, support: np.ndarray, *, fourier: bool = False) -> float:
    vector = np.zeros(1 << n, dtype=complex)
    vector[support] = 1 / math.sqrt(len(support))
    if fourier:
        # Opposite-sign QFT conjugates these real input amplitudes and preserves entropy.
        vector = np.fft.fft(vector, norm="ortho")
    matrix = vector.reshape(1 << (n // 2), -1)
    values = np.linalg.eigvalsh(matrix.conj().T @ matrix)
    require(float(values.min()) >= -2e-12, "independent density matrix is not positive")
    close(float(values.sum()), 1.0, "independent density matrix trace")
    values = values[values > 1e-14]
    return float(-np.sum(values * np.log2(values)))


def independent_almost_prime_supports(n: int) -> dict[int, np.ndarray]:
    # Trial division is independent of the production sieve.
    counts = np.zeros(1 << n, dtype=int)
    for number in range(2, 1 << n):
        remainder, factor = number, 2
        while factor * factor <= remainder:
            while remainder % factor == 0:
                counts[number] += 1
                remainder //= factor
            factor += 1
        if remainder > 1:
            counts[number] += 1
    return {k: np.flatnonzero((counts >= 1) & (counts <= k)) for k in (1, 2, 3)}


def residue_checks(root: Path = ROOT, *, replay_samples: int = 3) -> dict[str, float | int]:
    require(0 <= replay_samples <= 1000, "replay_samples must lie in [0, 1000]")
    metadata = {"n", "k", "M", "matched_low_bits", "modulus", "seed"}
    stat_names = {"mean", "std", "sem", "q025", "median", "q975", "min", "max"}
    summary_fields = metadata | {
        "samples", "residue_counts", "residue_entropy_bits", "entropy_ceiling_bits",
        "forced_gap_from_balanced_max_bits", "paired_null_fourier_minus_position_mean",
        "paired_null_fourier_minus_position_std",
    }
    for basis in ("position", "fourier"):
        summary_fields |= {f"null_{basis}_{stat}" for stat in stat_names}
        summary_fields |= {f"structured_{basis}_entropy", f"{basis}_residual_deficit",
                           f"{basis}_standardized_deficit", f"{basis}_null_count_at_or_below_structured",
                           f"cardinality_deficit_reduction_{basis}"}
    summary = read_rows(root / "data/residue_matched_summary.csv", summary_fields)
    raw = read_rows(root / "data/residue_matched_samples.csv", metadata | {
        "sample", "position_entropy", "fourier_entropy", "fourier_minus_position"})
    bounds = read_rows(root / "data/residue_entropy_bounds.csv", metadata - {"seed"} | {
        "residue_counts", "residue_probabilities", "residue_entropy_bits", "entropy_ceiling_bits",
        "forced_gap_from_balanced_max_bits", "structured_position_entropy", "ceiling_minus_structured_entropy"})
    key = lambda row: (int(row["n"]), int(row["k"]), int(row["matched_low_bits"]))
    expected = {(14, k, t) for k in (1, 2, 3) for t in range(4)}
    by_key = {key(row): row for row in summary}
    bound_by_key = {key(row): row for row in bounds}
    require(len(summary) == len(by_key) == 12 and set(by_key) == expected,
            "summary must contain each of the 12 recorded groups exactly once")
    require(len(bounds) == len(bound_by_key) == 12 and set(bound_by_key) == expected,
            "bounds must contain each of the 12 recorded groups exactly once")
    grouped: dict[tuple[int, int, int], list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        grouped[key(row)].append(row)
    require(set(grouped) == expected, "raw group keys disagree with the recorded experiment")

    supports = independent_almost_prime_supports(14)
    production_supports = dict(almost_prime_union_supports(14))
    structured = {}
    for k, support in supports.items():
        require(np.array_equal(support, production_supports[k]), f"k={k}: production support sieve")
        structured[k] = (independent_entropy(14, support), independent_entropy(14, support, fourier=True))
        close(core.entropy_from_support(14, support), structured[k][0], f"k={k}: production entropy")
        close(core.entropy_from_state_vector(14, core.qft_state_from_support(14, support)),
              structured[k][1], f"k={k}: production Fourier entropy")

    max_error, max_replay_error, replayed = 0.0, 0.0, 0
    for group in sorted(expected):
        n, k, t = group
        row, bound = by_key[group], bound_by_key[group]
        label = f"n={n}, k={k}, t={t}"
        support, q = supports[k], 1 << t
        target = np.bincount(support % q, minlength=q)
        count, seed = int(row["samples"]), int(row["seed"])
        require(count == 1000, f"{label}: expected 1000 recorded samples")
        require(seed == 2026081701 + 1000 * k + 37 * t, f"{label}: recorded seed schedule")
        require(int(row["M"]) == len(support) and int(row["modulus"]) == q, f"{label}: support metadata")
        require(json.loads(row["residue_counts"]) == target.tolist(), f"{label}: residue populations")
        group_rows = grouped[group]
        require(len(group_rows) == count, f"{label}: sample count mismatch")
        require(sorted(int(r["sample"]) for r in group_rows) == list(range(count)),
                f"{label}: duplicate or missing sample IDs")
        group_rows.sort(key=lambda r: int(r["sample"]))
        for sample_row in group_rows:
            require(all(int(sample_row[field]) == int(row[field]) for field in metadata),
                    f"{label}: raw metadata mismatch")
        position = np.array([finite(r, "position_entropy") for r in group_rows])
        fourier = np.array([finite(r, "fourier_entropy") for r in group_rows])
        differences = np.array([finite(r, "fourier_minus_position") for r in group_rows])
        require(np.max(np.abs(differences - (fourier - position))) < 2e-11,
                f"{label}: paired raw differences")
        H = shannon(target)
        ceiling = n / 2 - t + H
        close(core.residue_class_entropy_ceiling(n, support, t), ceiling, f"{label}: production ceiling")
        for field, value in {
            "residue_entropy_bits": H, "entropy_ceiling_bits": ceiling,
            "forced_gap_from_balanced_max_bits": n / 2 - ceiling,
            "structured_position_entropy": structured[k][0],
        }.items():
            max_error = max(max_error, close(finite(row, field), value, f"{label}: {field}"))
            max_error = max(max_error, close(finite(bound, field), value, f"{label}: bounds {field}"))
        require(all(int(bound[field]) == int(row[field]) for field in metadata - {"seed"}),
                f"{label}: bounds metadata")
        require(json.loads(bound["residue_counts"]) == target.tolist(), f"{label}: bounds populations")
        probabilities = np.asarray(json.loads(bound["residue_probabilities"]), dtype=float)
        require(probabilities.shape == target.shape and np.all(np.isfinite(probabilities)),
                f"{label}: bounds probabilities shape or finiteness")
        require(np.max(np.abs(probabilities - target / len(support))) < 2e-12,
                f"{label}: bounds probabilities")
        close(finite(bound, "ceiling_minus_structured_entropy"), ceiling - structured[k][0],
              f"{label}: ceiling slack")
        require(structured[k][0] <= ceiling + 2e-11 and float(position.max()) <= ceiling + 2e-11,
                f"{label}: position residue ceiling violated")
        for basis, values, entropy in zip(("position", "fourier"), (position, fourier), structured[k]):
            require(np.all((values >= -2e-11) & (values <= n / 2 + 2e-11)), f"{label}: entropy range")
            close(finite(row, f"structured_{basis}_entropy"), entropy, f"{label}: structured {basis}")
            sd = float(values.std(ddof=1))
            require(sd > 0, f"{label}: zero null variance")
            statistics = dict(mean=float(values.mean()), std=sd, sem=sd / math.sqrt(count),
                              q025=float(np.quantile(values, .025)), median=float(np.median(values)),
                              q975=float(np.quantile(values, .975)), min=float(values.min()), max=float(values.max()))
            for stat, value in statistics.items():
                field = f"null_{basis}_{stat}"
                max_error = max(max_error, close(finite(row, field), value, f"{label}: {field}"))
            residual = statistics["mean"] - entropy
            close(finite(row, f"{basis}_residual_deficit"), residual, f"{label}: {basis} residual")
            close(finite(row, f"{basis}_standardized_deficit"), residual / sd,
                  f"{label}: {basis} standardized deficit", atol=2e-8)
            below = int(np.count_nonzero(values <= entropy))
            require(int(row[f"{basis}_null_count_at_or_below_structured"]) == below,
                    f"{label}: {basis} tail count mismatch")
            require(below == 0, f"{label}: recorded all-samples-above-structured claim fails")
            baseline = finite(by_key[(n, k, 0)], f"{basis}_residual_deficit")
            require(baseline != 0, f"{label}: zero cardinality baseline")
            close(finite(row, f"cardinality_deficit_reduction_{basis}"), 1 - residual / baseline,
                  f"{label}: {basis} deficit reduction")
        close(finite(row, "paired_null_fourier_minus_position_mean"), float(differences.mean()), f"{label}: paired mean")
        close(finite(row, "paired_null_fourier_minus_position_std"), float(differences.std(ddof=1)), f"{label}: paired std")

        rng = np.random.default_rng(seed)
        for sample in range(replay_samples):
            sampled = sample_matched_support(n, target, q, rng)
            require(len(sampled) == len(support) and len(np.unique(sampled)) == len(sampled),
                    f"{label}: sampler cardinality or uniqueness")
            require(np.all((sampled >= 0) & (sampled < (1 << n))), f"{label}: sampler label range")
            require(np.array_equal(np.bincount(sampled % q, minlength=q), target), f"{label}: sampler populations")
            for basis, actual, recorded in (
                ("position", core.entropy_from_support(n, sampled), position[sample]),
                ("fourier", core.entropy_from_state_vector(n, core.qft_state_from_support(n, sampled)), fourier[sample]),
            ):
                max_replay_error = max(max_replay_error, close(actual, recorded, f"{label}, sample={sample}: replay {basis}"))
            replayed += 1
    return {"groups_checked": len(summary), "raw_rows_checked": len(raw),
            "max_summary_error": max_error, "seed_replayed_supports": replayed,
            "max_seed_replay_entropy_error": max_replay_error}


def markdown_report(title: str, sections: dict[str, dict[str, float | int]], scope: str) -> str:
    lines = [f"# {title}", "", "**Status: PASS for the checks listed below.**", "", scope, ""]
    for name, values in sections.items():
        lines.extend([f"## {name}", "", "| Check | Result |", "| --- | ---: |"])
        for key, value in values.items():
            label = key.replace("log2_M_n", r"$`\log_2 M_n`$").replace("S_n", "$`S_n`$")
            # Preserve mathematical subscripts while spelling out ordinary keys.
            label = label.replace("_slope", " slope").replace("_intercept", " intercept")
            if "$" not in label:
                label = label.replace("_", " ")
            number = f"{value:.12g}"
            if "e" in number:
                mantissa, exponent = number.split("e")
                number = rf"$`{mantissa}\times10^{{{int(exponent)}}}`$"
            lines.append(f"| {label} | {number} |")
        lines.append("")
    return "\n".join(lines)


def write_report(text: str, output: Path | None) -> None:
    if output is not None:
        destination = output.resolve()
        require(destination.suffix.lower() == ".md", "--output must name a Markdown (.md) file")
        require(not destination.exists(), "refusing to overwrite an existing file")
        for name in ("data", "outputs", "subset_states", "scripts", "tests"):
            require(not destination.is_relative_to(ROOT / name), "refusing to write inside source or evidence directories")
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("x", encoding="utf-8") as handle:
            handle.write(text)
    print(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="save the report to a new .md file")
    parser.add_argument("--replay-samples", type=int, default=3, help="seed-replay prefix per group (1..1000; default 3)")
    args = parser.parse_args()
    try:
        require(1 <= args.replay_samples <= 1000, "--replay-samples must lie in [1, 1000]")
        results = residue_checks(replay_samples=args.replay_samples)
        report = markdown_report("Residue-control validation", {"Recorded experiment": results},
            "Checks all 12,000 raw rows against summaries and bounds, recomputes the three structured states, "
            "and replays the stated number of matched supports. This is not a full rerun unless all 1,000 "
            "samples per group are requested; it does not establish an arithmetic mechanism or novelty.")
        write_report(report, args.output)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"# Residue-control validation\n\n**Status: FAIL**\n\n{exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
