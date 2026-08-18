"""Local quadratic fits used to verify Table-I peak values."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def local_m_grid(
    n: int,
    center: int,
    *,
    points: int,
    span: float | None = 0.25,
    half_width: int | None = None,
) -> NDArray[np.int64]:
    if n <= 0:
        raise ValueError("n must be positive.")
    N = 1 << n
    if not (1 <= center <= N):
        raise ValueError(f"center must satisfy 1 <= center <= 2**n={N}; got {center}.")
    if points < 3:
        raise ValueError("points must be at least 3 for a quadratic fit.")
    if half_width is not None:
        width = int(half_width)
        if width <= 0:
            raise ValueError("half_width must be positive.")
        low, high = max(1, center - width), min(N, center + width)
    else:
        span = 0.25 if span is None else span
        if not (0 < span <= 1):
            raise ValueError("span must satisfy 0 < span <= 1.")
        low = max(1, int(round(center * (1.0 - span))))
        high = min(N, int(round(center * (1.0 + span))))
    grid = np.unique(np.rint(np.linspace(low, high, points)).astype(np.int64))
    grid = np.unique(np.concatenate([grid, np.asarray([center], dtype=np.int64)]))
    grid = grid[(grid >= 1) & (grid <= N)]
    if grid.size < 3:
        raise ValueError("local grid contains fewer than 3 distinct support sizes.")
    return grid


def quadratic_peak_fit(
    m: ArrayLike,
    mean: ArrayLike,
    sem: ArrayLike | None = None,
) -> dict[str, float | str | bool]:
    m_arr = np.asarray(m, dtype=np.float64)
    y = np.asarray(mean, dtype=np.float64)
    if m_arr.ndim != 1 or y.ndim != 1 or m_arr.size != y.size:
        raise ValueError("m and mean must be one-dimensional arrays of the same length.")
    if m_arr.size < 3 or np.any(m_arr <= 0):
        raise ValueError("at least three positive support sizes are required.")
    x = np.log2(m_arr)
    weights = None
    if sem is not None:
        err = np.asarray(sem, dtype=np.float64)
        positive = err[np.isfinite(err) & (err > 0)]
        if positive.size:
            weights = 1.0 / np.maximum(err, positive.min())
    coeff = np.polyfit(x, y, deg=2, w=weights)
    a, b, c = [float(v) for v in coeff]
    best = int(np.argmax(y))
    interior = False
    if a < 0:
        x_peak = -b / (2.0 * a)
        interior = bool(x.min() <= x_peak <= x.max())
    if a >= 0 or not interior:
        x_peak = float(x[best])
        S_peak = float(y[best])
        status = "best_grid_point_reported"
    else:
        S_peak = c - b * b / (4.0 * a)
        status = "concave_quadratic"
    return {
        "fit_status": status,
        "fit_variable": "log2_M",
        "quadratic_a": a,
        "quadratic_b": b,
        "quadratic_c": c,
        "log2_M_peak": float(x_peak),
        "log2_M_peak_sem": float("nan"),
        "M_peak": float(2.0**x_peak),
        "M_peak_sem": float("nan"),
        "S_peak": float(S_peak),
        "S_peak_sem": float("nan"),
        "best_grid_M": float(m_arr[best]),
        "best_grid_S": float(y[best]),
        "peak_inside_window": interior,
    }
