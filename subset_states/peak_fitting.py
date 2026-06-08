"""Local quadratic fits used to verify Table-I peak values."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def local_m_grid(n: int, center: int, *, points: int, span: float | None = 0.25, half_width: int | None = None) -> NDArray[np.int64]:
    """Return a sorted integer grid of support sizes around a tabulated peak.

    Parameters
    ----------
    n:
        Number of qubits. The allowed support sizes are 1 <= M <= 2**n.
    center:
        Tabulated peak support size M_n.
    points:
        Number of grid points requested. The center is always included.
    span:
        Fractional half-width of the grid, used when ``half_width`` is not given.
        For example, span=0.20 gives approximately [0.8 center, 1.2 center].
    half_width:
        Optional absolute half-width around center. If supplied, it overrides
        ``span``.
    """

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
        low = max(1, center - width)
        high = min(N, center + width)
    else:
        if span is None:
            span = 0.25
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


def quadratic_peak_fit(m: ArrayLike, mean: ArrayLike, sem: ArrayLike | None = None) -> dict[str, float | str | bool]:
    """Fit mean entropy as a quadratic function of log2(M).

    The returned dictionary reports the quadratic coefficients, the estimated
    peak, and a conservative status string. If the fitted parabola is not
    concave or its maximum falls outside the sampled window, the best observed
    grid point is reported instead of a fitted interior maximum.
    """

    m_arr = np.asarray(m, dtype=np.float64)
    y = np.asarray(mean, dtype=np.float64)
    if m_arr.ndim != 1 or y.ndim != 1 or m_arr.size != y.size:
        raise ValueError("m and mean must be one-dimensional arrays of the same length.")
    if m_arr.size < 3:
        raise ValueError("at least 3 grid points are required for a quadratic fit.")
    if np.any(m_arr <= 0):
        raise ValueError("all support sizes M must be positive.")

    x = np.log2(m_arr)

    weights = None
    if sem is not None:
        err = np.asarray(sem, dtype=np.float64)
        if err.shape != y.shape:
            raise ValueError("sem must have the same shape as mean.")
        positive = err[np.isfinite(err) & (err > 0)]
        if positive.size:
            floor = positive.min()
            weights = 1.0 / np.maximum(err, floor)

    coeff = np.polyfit(x, y, deg=2, w=weights)
    a, b, c = [float(v) for v in coeff]

    cov = None
    if x.size > 3:
        try:
            _, cov = np.polyfit(x, y, deg=2, w=weights, cov=True)
        except Exception:
            cov = None

    best = int(np.argmax(y))
    status = "concave_quadratic"
    interior_peak = False

    if a < 0:
        x_peak = -b / (2.0 * a)
        interior_peak = bool(x.min() <= x_peak <= x.max())
        if interior_peak:
            S_peak = c - b * b / (4.0 * a)
        else:
            status = "concave_quadratic_peak_outside_window_best_grid_point_reported"
            x_peak = float(x[best])
            S_peak = float(y[best])
    else:
        status = "non_concave_quadratic_best_grid_point_reported"
        x_peak = float(x[best])
        S_peak = float(y[best])

    x_sem = float("nan")
    S_sem = float("nan")
    if cov is not None and a < 0 and interior_peak:
        grad_x = np.asarray([b / (2.0 * a * a), -1.0 / (2.0 * a), 0.0])
        grad_S = np.asarray([b * b / (4.0 * a * a), -b / (2.0 * a), 1.0])
        x_var = float(grad_x @ cov @ grad_x)
        S_var = float(grad_S @ cov @ grad_S)
        x_sem = float(np.sqrt(max(x_var, 0.0)))
        S_sem = float(np.sqrt(max(S_var, 0.0)))
    elif sem is not None:
        err = np.asarray(sem, dtype=np.float64)
        if np.isfinite(err[best]):
            S_sem = float(err[best])

    M_peak = float(2.0**x_peak)
    M_sem = float(np.log(2.0) * M_peak * x_sem) if np.isfinite(x_sem) else float("nan")

    return {
        "fit_status": status,
        "fit_variable": "log2_M",
        "quadratic_a": a,
        "quadratic_b": b,
        "quadratic_c": c,
        "log2_M_peak": float(x_peak),
        "log2_M_peak_sem": x_sem,
        "M_peak": M_peak,
        "M_peak_sem": M_sem,
        "S_peak": float(S_peak),
        "S_peak_sem": S_sem,
        "best_grid_M": float(m_arr[best]),
        "best_grid_S": float(y[best]),
        "peak_inside_window": bool(interior_peak),
    }
