from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Iterable, Sequence
import csv

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .core import (
    entropy_from_state_vector,
    entropy_from_support,
    matrix_from_support,
    omega_prime_factors_sieve,
    qft_state_from_support,
    random_balanced_partition,
    random_subset,
    renyi_entropy_from_spectrum,
    spectrum_from_matrix,
    summary_stats,
)


def m_grid(start: int, stop: int, points: int, *, include_stop: bool = True) -> NDArray[np.int64]:
    """Return a sorted grid of unique integer support sizes."""

    if points <= 1:
        return np.array([start], dtype=np.int64)
    endpoint = include_stop
    grid = np.linspace(start, stop, points, endpoint=endpoint)
    grid = np.unique(np.rint(grid).astype(np.int64))
    return grid[(grid >= start) & (grid <= stop)]


def write_rows(path: str | Path, rows: Iterable[dict], fieldnames: Sequence[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_csv_columns(path: str | Path) -> dict[str, NDArray]:
    import csv

    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        data: dict[str, list] = {name: [] for name in reader.fieldnames or []}
        for row in reader:
            for key, value in row.items():
                try:
                    data[key].append(float(value))
                except (TypeError, ValueError):
                    data[key].append(value)
    return {key: np.asarray(value) for key, value in data.items()}


def random_subset_entropy_samples(
    n: int,
    m: int,
    samples: int,
    *,
    seed: int,
    order: float = 1.0,
    right_bits: Sequence[int] | None = None,
) -> NDArray[np.float64]:
    rng = np.random.default_rng(seed)
    values = np.empty(samples, dtype=np.float64)
    for idx in range(samples):
        support = random_subset(n, m, rng)
        values[idx] = entropy_from_support(n, support, right_bits, order=order)
    return values


def random_subset_renyi_samples(
    n: int,
    m: int,
    samples: int,
    *,
    seed: int,
    orders: Sequence[float] = (1.0, 2.0, np.inf),
) -> dict[float, NDArray[np.float64]]:
    rng = np.random.default_rng(seed)
    values = {float(order): np.empty(samples, dtype=np.float64) for order in orders}
    for idx in range(samples):
        support = random_subset(n, m, rng)
        eigvals = spectrum_from_matrix(matrix_from_support(n, support))
        for order in orders:
            values[float(order)][idx] = renyi_entropy_from_spectrum(eigvals, order=float(order))
    return values


def random_subset_qft_samples(
    n: int,
    m: int,
    samples: int,
    *,
    seed: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    rng = np.random.default_rng(seed)
    position = np.empty(samples, dtype=np.float64)
    fourier = np.empty(samples, dtype=np.float64)
    for idx in range(samples):
        support = random_subset(n, m, rng)
        position[idx] = entropy_from_support(n, support)
        fourier_state = qft_state_from_support(n, support)
        fourier[idx] = entropy_from_state_vector(n, fourier_state)
    return position, fourier


def partition_entropy_samples_for_fixed_subset(
    n: int,
    m: int,
    samples: int,
    *,
    seed: int,
) -> NDArray[np.float64]:
    rng = np.random.default_rng(seed)
    support = random_subset(n, m, rng)
    values = np.empty(samples, dtype=np.float64)
    for idx in range(samples):
        right_bits = random_balanced_partition(n, rng)
        values[idx] = entropy_from_support(n, support, right_bits)
    return values


def partition_entropy_samples_for_state(
    n: int,
    state: ArrayLike,
    samples: int,
    *,
    seed: int,
) -> NDArray[np.float64]:
    rng = np.random.default_rng(seed)
    values = np.empty(samples, dtype=np.float64)
    for idx in range(samples):
        right_bits = random_balanced_partition(n, rng)
        values[idx] = entropy_from_state_vector(n, state, right_bits)
    return values


def complex_haar_state(n: int, *, seed: int) -> NDArray[np.complex128]:
    rng = np.random.default_rng(seed)
    N = 1 << n
    state = rng.normal(size=N) + 1j * rng.normal(size=N)
    return (state / np.linalg.norm(state)).astype(np.complex128)


def almost_prime_union_supports(n: int) -> list[tuple[int, NDArray[np.int64]]]:
    """Return supports U_{N,k} for k=1,...,n-1.

    U_{N,k} contains integers below 2**n with 1 <= Ω(x) <= k. Thus U_{N,n-1}
    excludes both 0 and 1, as required by the usual definition of k-almost primes.
    """

    N = 1 << n
    omega = omega_prime_factors_sieve(N)
    supports: list[tuple[int, NDArray[np.int64]]] = []
    for k in range(1, n):
        support = np.flatnonzero((omega >= 1) & (omega <= k)).astype(np.int64)
        supports.append((k, support))
    return supports


def stats_row(prefix: str, values: ArrayLike) -> dict[str, float | int]:
    stats = summary_stats(values)
    return {f"{prefix}_{key}": value for key, value in asdict(stats).items()}
