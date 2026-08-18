from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from .core import (
    dense_bulk_approximation,
    exact_page_entropy_bits,
    hypergeometric_occupancy_approximation,
    mean_matrix_uniform_eigenvalue,
)
from .experiments import write_rows
from .plotting import COOL_PALETTE, apply_journal_style, save_figure


def _coerce(value: str):
    if value is None:
        return value
    text = str(value).strip()
    if text == "":
        return ""
    try:
        return float(text)
    except ValueError:
        return text


def read_csv_rows(path: str | Path) -> list[dict]:
    path = Path(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return [{key: _coerce(value) for key, value in row.items()} for row in reader]


def _arr(rows: Sequence[dict], key: str) -> np.ndarray:
    return np.asarray([float(row[key]) for row in rows], dtype=float)


def _str_arr(rows: Sequence[dict], key: str) -> np.ndarray:
    return np.asarray([str(row[key]) for row in rows])


def _infer_n_from_m(rows: Sequence[dict]) -> int:
    mmax = int(max(float(row["m"]) for row in rows))
    n = int(round(math.log2(mmax)))
    return n if (1 << n) == mmax else int(math.ceil(math.log2(mmax)))


def _fit_line(x: np.ndarray, y: np.ndarray) -> dict:
    result = stats.linregress(x, y)
    residuals = y - (result.slope * x + result.intercept)
    return {
        "slope": float(result.slope),
        "intercept": float(result.intercept),
        "r_value": float(result.rvalue),
        "r_squared": float(result.rvalue**2),
        "p_value": float(result.pvalue),
        "slope_stderr": float(result.stderr),
        "intercept_stderr": float(result.intercept_stderr),
        "residual_std": float(residuals.std(ddof=2)) if len(x) > 2 else float("nan"),
    }


def plot_fig1(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None) -> None:
    """Plot the two averaging procedures with explicit entropy units."""

    data_dir = Path(data_dir)
    summary = read_csv_rows(data_dir / "fig1_summary.csv")
    zoom = read_csv_rows(data_dir / "fig1_zoom_summary.csv")
    if n is None:
        n = _infer_n_from_m(summary)

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.8))

    def draw(axis, rows, markersize: float = 2.5):
        m = _arr(rows, "m")
        y_states = _arr(rows, "states_mean")
        e_states = _arr(rows, "states_sem")
        y_part = _arr(rows, "partitions_mean")
        e_part = _arr(rows, "partitions_sem")
        axis.errorbar(
            m,
            y_part,
            yerr=e_part,
            fmt="o-",
            markersize=markersize,
            capsize=1.5,
            label="fixed support, random balanced cuts",
        )
        axis.errorbar(
            m,
            y_states,
            yerr=e_states,
            fmt="s--",
            markersize=markersize,
            capsize=1.5,
            label="fixed cut, random supports",
        )

    draw(ax, summary)
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel("mean entropy (bits)")
    ax.set_xlim(0, 1 << n)
    ax.set_ylim(0, n / 2)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.38), ncol=1, frameon=False)

    inset = ax.inset_axes((0.59, 0.59, 0.38, 0.38))
    draw(inset, zoom, markersize=2.0)
    inset.set_xlim(0, max(float(row["m"]) for row in zoom) * 1.02)
    inset.set_ylim(max(0, n / 2 - 3), n / 2)
    inset.set_xlabel(r"$M$", labelpad=1)
    inset.set_ylabel("bits", labelpad=1)
    inset.legend().remove()
    save_figure(fig, out_pdf)


def plot_fig2(data_dir: str | Path, out_pdf: str | Path, *, write_fit_summary: bool = True) -> None:
    """Plot Table-I scaling, retaining all n=10,...,30 points."""

    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig2_table_with_page.csv")
    n = _arr(rows, "n")
    M_n = _arr(rows, "M_n")
    S_n = _arr(rows, "S_n")
    log2_M = _arr(rows, "log2_M_n") if "log2_M_n" in rows[0] else np.log2(M_n)
    page = (
        _arr(rows, "Page_exact_bits")
        if "Page_exact_bits" in rows[0]
        else np.array([exact_page_entropy_bits(int(x)) for x in n])
    )

    fit_logM = _fit_line(n, log2_M)
    fit_S = _fit_line(n, S_n)
    if write_fit_summary:
        fit_rows = [{"quantity": "log2_M_n", **fit_logM}, {"quantity": "S_n", **fit_S}]
        write_rows(data_dir / "fig2_linear_fit_summary_from_csv.csv", fit_rows, fit_rows[0].keys())

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    x = np.linspace(n.min() - 1, n.max() + 1, 200)

    points_m = ax.plot(n, log2_M, "o", label=r"$\log_2\widehat M_n$")[0]
    ax.plot(
        x,
        fit_logM["slope"] * x + fit_logM["intercept"],
        "--",
        linewidth=1,
        color=points_m.get_color(),
    )
    points_s = ax.plot(n, S_n, "s", label=r"$\widehat S_n$")[0]
    ax.plot(
        x,
        fit_S["slope"] * x + fit_S["intercept"],
        "--",
        linewidth=1,
        color=points_s.get_color(),
    )
    ax.plot(n, page, ":", label="Page mean")
    ax.set_xlabel(r"number of qubits $n$")
    ax.set_ylabel("value (bits)")
    ax.set_xlim(n.min() - 1, n.max() + 1)
    ax.set_xticks([10, 15, 20, 25, 30])
    ax.set_ylim(0, max(log2_M.max(), page.max()) + 1)
    ax.legend(frameon=False, loc="upper left")
    save_figure(fig, out_pdf)


def plot_fig3(data_dir: str | Path, out_pdf: str | Path, *, bins: int = 500) -> None:
    """Plot the spectral bulk and compare three isolated-mode estimates."""

    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig3_spectrum.csv")
    eigvals = _arr(rows, "lambda")
    eigvals.sort()
    lambda0 = float(eigvals[-1])
    positive_bulk = eigvals[:-1]
    positive_bulk = positive_bulk[positive_bulk > 0]
    if positive_bulk.size == 0:
        raise ValueError("fig3_spectrum.csv contains no positive bulk eigenvalues.")

    summary_path = data_dir / "fig3_spectrum_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(
            "fig3_spectrum_summary.csv is required for the fixed-cardinality mean-mode comparison."
        )
    summary = read_csv_rows(summary_path)[0]
    n = int(float(summary["n"]))
    m = int(float(summary["M"]))
    lambda_mf = mean_matrix_uniform_eigenvalue(n, m)
    lambda_leading = m / (1 << n)

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.9, 3.6))
    ax.hist(positive_bulk, bins=bins, color=COOL_PALETTE[0])
    ax.set_xlabel(r"bulk eigenvalue $\lambda$")
    ax.set_ylabel("count")
    ax.set_xlim(0, positive_bulk.max() * 1.05)
    ax.set_title(rf"$n={n}$, $M={m:,}$", loc="left", pad=3)

    inset = ax.inset_axes((0.52, 0.52, 0.45, 0.41))
    values = np.asarray([lambda0, lambda_mf, lambda_leading])
    y = np.asarray([2, 1, 0])
    labels = [r"observed $\lambda_0$", r"mode of $\mathbb{E}[\rho]$", r"leading $M/N$"]
    for value, yi, color in zip(values, y, COOL_PALETTE[:3]):
        inset.plot(value, yi, "o", color=color, markersize=4)
        inset.axvline(value, color=color, linewidth=0.8, alpha=0.45)
    pad = max(values.max() - values.min(), 1e-6) * 0.18
    inset.set_xlim(values.min() - pad, values.max() + pad)
    inset.set_ylim(-0.55, 2.55)
    inset.set_yticks(y, labels)
    inset.set_xlabel("isolated eigenvalue", labelpad=1)
    inset.tick_params(axis="y", labelsize=6.8)
    inset.tick_params(axis="x", labelsize=6.8)
    save_figure(fig, out_pdf)


def _draw_fig4(axis, rows: Sequence[dict], n: int, label_data: bool) -> None:
    m = _arr(rows, "m").astype(int)
    mean = _arr(rows, "entropy_mean")
    sem = _arr(rows, "entropy_sem")
    sparse = np.asarray([hypergeometric_occupancy_approximation(n, int(x)) for x in m])
    dense = np.asarray([dense_bulk_approximation(n, int(x)) for x in m])
    axis.errorbar(
        m,
        mean,
        yerr=sem,
        fmt="o",
        markersize=2.2,
        capsize=1.2,
        label="numerical mean" if label_data else None,
    )
    axis.plot(m, sparse, ":", label=r"diagonal occupancy $D_{N,M}$" if label_data else None)
    axis.plot(m, dense, "--", label=r"dense bulk $T_{N,M}$" if label_data else None)


def plot_fig4(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None) -> None:
    """Plot numerical data with corrected fixed-cardinality approximations."""

    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig4_summary.csv")
    zoom = read_csv_rows(data_dir / "fig4_zoom_summary.csv")
    if n is None:
        n = _infer_n_from_m(rows)

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.8))
    _draw_fig4(ax, rows, n, True)
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel("entropy (bits)")
    ax.set_xlim(0, 1 << n)
    ax.set_ylim(0, n / 2)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.35), ncol=3)

    inset = ax.inset_axes((0.59, 0.59, 0.38, 0.38))
    _draw_fig4(inset, zoom, n, False)
    inset.set_xlim(0, max(float(row["m"]) for row in zoom) * 1.02)
    inset.set_ylim(max(0, n / 2 - 3), n / 2)
    inset.set_xlabel(r"$M$", labelpad=1)
    inset.set_ylabel("bits", labelpad=1)
    save_figure(fig, out_pdf)


def plot_fig5(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None) -> None:
    """Compare deterministic arithmetic supports with the random null distribution.

    The shaded regions are one ensemble standard deviation.  They represent the
    spread of the random-support null distribution, not uncertainty in its mean.
    """

    data_dir = Path(data_dir)
    random_rows = read_csv_rows(data_dir / "fig5_random_qft_summary.csv")
    union_rows = read_csv_rows(data_dir / "fig5_almost_prime_unions.csv")
    if n is None:
        n = _infer_n_from_m(random_rows)

    m = _arr(random_rows, "m")
    pos = _arr(random_rows, "position_mean")
    pos_std = _arr(random_rows, "position_std")
    qft = _arr(random_rows, "qft_mean")
    qft_std = _arr(random_rows, "qft_std")
    union_m = _arr(union_rows, "m")
    union_pos = _arr(union_rows, "position_entropy")
    union_qft = _arr(union_rows, "qft_entropy")

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.7))
    line_pos = ax.plot(m, pos, "-", linewidth=1.2, label=r"random supports: mean $\pm$ 1 SD")[0]
    ax.fill_between(
        m,
        np.maximum(0.0, pos - pos_std),
        np.minimum(n / 2, pos + pos_std),
        color=line_pos.get_color(),
        alpha=0.16,
        linewidth=0,
    )
    line_qft = ax.plot(
        m,
        qft,
        "--",
        linewidth=1.2,
        label=r"after QFT: mean $\pm$ 1 SD",
    )[0]
    ax.fill_between(
        m,
        np.maximum(0.0, qft - qft_std),
        np.minimum(n / 2, qft + qft_std),
        color=line_qft.get_color(),
        alpha=0.14,
        linewidth=0,
    )
    ax.plot(union_m, union_pos, "+", markersize=6, label=r"$U_{N,k}$")
    ax.plot(union_m, union_qft, "x", markersize=5, label=r"QFT $U_{N,k}$")
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel("entropy (bits)")
    ax.set_xlim(0, 1 << n)
    ax.set_ylim(0, n / 2)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.46), ncol=2)
    save_figure(fig, out_pdf)


def _order_key(value) -> float:
    text = str(value).strip().lower()
    if text in {"inf", "infinity", "np.inf"}:
        return float("inf")
    return float(value)


def _rows_for_order(rows: Sequence[dict], order: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    filtered = [
        row
        for row in rows
        if _order_key(row["order"]) == order
        or (np.isinf(order) and np.isinf(_order_key(row["order"])))
    ]
    return _arr(filtered, "m"), _arr(filtered, "mean"), _arr(filtered, "sem")


def _draw_fig6(axis, rows: Sequence[dict]) -> list:
    orders = [1.0, 2.0, np.inf]
    labels = {1.0: r"$S^{(1)}$", 2.0: r"$S^{(2)}$", np.inf: r"$S^{(\infty)}$"}
    styles = {1.0: "o-", 2.0: "s--", np.inf: "^:"}
    handles = []
    for order in orders:
        m, mean, sem = _rows_for_order(rows, order)
        handle = axis.errorbar(
            m,
            mean,
            yerr=sem,
            fmt=styles[order],
            markersize=2.2,
            capsize=1.2,
            label=labels[order],
        )
        handles.append(handle)
    return handles


def plot_fig6(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None) -> None:
    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig6_renyi_summary.csv")
    zoom = read_csv_rows(data_dir / "fig6_renyi_zoom_summary.csv")
    if n is None:
        n = _infer_n_from_m(rows)
    zoom_stop = int(max(float(row["m"]) for row in zoom))

    apply_journal_style()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.7, 5.7), sharey=False)
    handles = _draw_fig6(ax1, rows)
    ax1.set_xlabel(r"support size $M$")
    ax1.set_ylabel("Rényi entropy (bits)")
    ax1.set_xlim(0, 1 << n)
    ax1.set_ylim(0, n / 2)
    ax1.text(0.02, 0.94, "full range", transform=ax1.transAxes, va="top", fontsize=8)

    _draw_fig6(ax2, zoom)
    ax2.set_xlabel(r"support size $M$")
    ax2.set_ylabel("Rényi entropy (bits)")
    ax2.set_xlim(0, zoom_stop)
    ax2.set_ylim(max(0, n / 2 - 3), n / 2)
    ax2.text(0.02, 0.94, "peak region", transform=ax2.transAxes, va="top", fontsize=8)

    labels = [r"$S^{(1)}$", r"$S^{(2)}$", r"$S^{(\infty)}$"]
    fig.legend(handles=handles, labels=labels, frameon=False, ncol=3, loc="lower center")
    fig.subplots_adjust(hspace=0.45, bottom=0.12)
    save_figure(fig, out_pdf)


def _families(rows: Sequence[dict]) -> list[str]:
    key = "state_family" if "state_family" in rows[0] else "family"
    seen: list[str] = []
    for fam in _str_arr(rows, key):
        if str(fam) not in seen:
            seen.append(str(fam))
    return seen


def plot_fig7(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None, bins: int = 100) -> None:
    """Plot both bipartition distributions on a common entropy scale."""

    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig7_partition_samples.csv")
    summary_path = data_dir / "fig7_partition_summary.csv"
    summary = read_csv_rows(summary_path) if summary_path.exists() else []
    if n is None and summary and summary[0].get("n", "") != "":
        n = int(float(summary[0]["n"]))
    n = 20 if n is None else n
    page = (
        float(summary[0]["Page_exact_bits"])
        if summary and "Page_exact_bits" in summary[0]
        else exact_page_entropy_bits(n)
    )
    m_subset = None
    for row in summary:
        if str(row.get("state_family", "")) == "subset" and row.get("M", "") != "":
            m_subset = int(float(row["M"]))

    key = "state_family" if "state_family" in rows[0] else "family"
    preferred = [fam for fam in ("subset", "complex_haar") if fam in _families(rows)]
    families = preferred or _families(rows)
    arrays = {
        family: np.asarray(
            [float(row["entropy"]) for row in rows if str(row[key]) == family],
            dtype=float,
        )
        for family in families
    }
    global_min = min(values.min() for values in arrays.values())
    global_max = max(values.max() for values in arrays.values())
    padding = 0.04 * (global_max - global_min)
    lo, hi = global_min - padding, global_max + padding
    bin_edges = np.linspace(lo, hi, bins + 1)

    apply_journal_style()
    fig, axes = plt.subplots(len(families), 1, figsize=(4.7, 2.65 * len(families)), sharex=True)
    if len(families) == 1:
        axes = [axes]
    for index, (ax, family) in enumerate(zip(axes, families)):
        values = arrays[family]
        ax.hist(values, bins=bin_edges)
        ax.axvline(values.mean(), linestyle="--", linewidth=1, label="sample mean")
        ax.axvline(page, linestyle=":", linewidth=1, label="Page mean")
        if family == "complex_haar":
            title = "complex Haar state"
        elif family == "subset" and m_subset is not None:
            title = rf"peak-support subset state, $M={m_subset}$"
        else:
            title = family.replace("_", " ")
        ax.set_ylabel("count")
        ax.set_title(title)
        ax.set_xlim(lo, hi)
        if index == 0:
            ax.legend(frameon=False, loc="upper right")
    axes[-1].set_xlabel("entropy (bits)")
    fig.subplots_adjust(hspace=0.30)
    save_figure(fig, out_pdf)


def _n_from_fig8_path(path: Path) -> int:
    match = re.search(r"fig8_samples_n(\d+)\.csv$", path.name)
    if not match:
        raise ValueError(f"cannot infer n from {path}")
    return int(match.group(1))


def plot_fig8(
    data_dir: str | Path,
    out_pdf: str | Path,
    *,
    n_values: Sequence[int] | None = None,
    bins: int = 120,
) -> None:
    """Repository-only best-of-random-candidate screening figure.

    Histograms are normalized because the two sample pools have different sizes.
    This exploratory plot is intentionally not part of the revised manuscript.
    """

    data_dir = Path(data_dir)
    files = (
        [data_dir / f"fig8_samples_n{int(n)}.csv" for n in n_values]
        if n_values
        else sorted(data_dir.glob("fig8_samples_n*.csv"), key=_n_from_fig8_path)
    )
    if not files:
        raise FileNotFoundError(f"no fig8_samples_n*.csv files found in {data_dir}")

    apply_journal_style()
    fig, axes = plt.subplots(len(files), 1, figsize=(4.8, 3.0 * len(files)))
    if len(files) == 1:
        axes = [axes]
    for ax, path in zip(axes, files):
        rows = read_csv_rows(path)
        n = _n_from_fig8_path(path)
        key = "family" if "family" in rows[0] else "state_family"
        arrays = []
        for family in _families(rows):
            values = np.asarray(
                [float(row["entropy"]) for row in rows if str(row[key]) == family],
                dtype=float,
            )
            arrays.append((family, values))
        lo = min(values.min() for _, values in arrays)
        hi = max(values.max() for _, values in arrays)
        bin_edges = np.linspace(lo, hi, bins)
        for family, values in arrays:
            label = (
                "random supports, natural cut"
                if family.startswith("random")
                else "selected random candidate, all balanced cuts"
            )
            ax.hist(values, bins=bin_edges, density=True, alpha=0.55, label=label)
            ax.axvline(values.mean(), linestyle="--" if family.startswith("random") else ":", linewidth=1)
        M = int(float(rows[0].get("M", rows[0].get("m", 0))))
        ax.set_xlabel("entropy (bits)")
        ax.set_ylabel("density")
        ax.set_title(rf"$n={n}$, $M={M}$")
        ax.legend(frameon=False)
    fig.subplots_adjust(hspace=0.55)
    save_figure(fig, out_pdf)


def plot_peak_search(samples_csv: str | Path, fit_csv: str | Path, out_pdf: str | Path) -> None:
    samples = read_csv_rows(samples_csv)
    fit_rows = read_csv_rows(fit_csv)
    if not fit_rows:
        raise ValueError(f"no fit row found in {fit_csv}")
    fit = fit_rows[0]
    m = _arr(samples, "m")
    mean = _arr(samples, "mean")
    sem = _arr(samples, "sem")
    x = np.log2(m.astype(float))
    x_fit = np.linspace(x.min(), x.max(), 300)
    y_fit = (
        float(fit["quadratic_a"]) * x_fit**2
        + float(fit["quadratic_b"]) * x_fit
        + float(fit["quadratic_c"])
    )

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    ax.errorbar(x, mean, yerr=sem, fmt="o", capsize=1.5, label="samples")
    ax.plot(x_fit, y_fit, "--", label="quadratic fit")
    ax.axvline(float(fit["log2_M_peak"]), linestyle=":", linewidth=1, label="estimated peak")
    ax.set_xlabel(r"$\log_2 M$")
    ax.set_ylabel("entropy (bits)")
    ax.legend(frameon=False)
    save_figure(fig, out_pdf)
