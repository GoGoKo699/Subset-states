"""Manuscript tables used by the reproducibility scripts."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import numpy as np
from numpy.typing import NDArray

TABLE_I_PEAKS: tuple[tuple[int, int, float], ...] = (
    (10, 107, 4.072),
    (12, 276, 5.108),
    (14, 716, 6.143),
    (16, 1873, 7.176),
    (18, 4934, 8.196),
    (20, 13091, 9.215),
    (22, 34771, 10.231),
    (24, 93018, 11.242),
    (26, 250660, 12.251),
    (28, 672556, 13.258),
    (30, 1836685, 14.263),
)

PEAK_M: dict[int, int] = {n: m for n, m, _ in TABLE_I_PEAKS}
PEAK_S: dict[int, float] = {n: s for n, _, s in TABLE_I_PEAKS}


def table_i_array() -> NDArray[np.float64]:
    return np.asarray(TABLE_I_PEAKS, dtype=float)


def read_table_i(path: str | Path | None = None) -> list[dict[str, int | float]]:
    if path is None:
        return [{"n": n, "M_n": m, "S_n": s} for n, m, s in TABLE_I_PEAKS]
    path = Path(path)
    with path.open(newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    required = {"n", "M_n", "S_n"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return [
        {"n": int(row["n"]), "M_n": int(row["M_n"]), "S_n": float(row["S_n"])}
        for row in reader
    ]


def select_table_rows(
    rows: Iterable[dict[str, int | float]],
    n_values: Iterable[int] | None = None,
) -> list[dict[str, int | float]]:
    selected = None if n_values is None else {int(n) for n in n_values}
    out = [row for row in rows if selected is None or int(row["n"]) in selected]
    if not out:
        raise ValueError("No Table I rows matched the requested n values.")
    return out
