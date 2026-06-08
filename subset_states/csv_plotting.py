from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from .core import exact_page_entropy_bits
from .experiments import write_rows
from .plotting import apply_journal_style, make_inset_axis, save_figure


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
    if (1 << n) == mmax:
        return n
    return int(math.ceil(math.log2(mmax)))


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
        axis.errorbar(m, y_part, yerr=e_part, fmt="o-", markersize=markersize, capsize=1.5, label="fixed state, random partitions")
        axis.errorbar(m, y_states, yerr=e_states, fmt="s--", markersize=markersize, capsize=1.5, label="fixed partition, random states")

    draw(ax, summary)
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel(r"entropy $S_{N,M}$")
    ax.set_xlim(0, 1 << n)
    ax.set_ylim(0, n / 2)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.38), ncol=1, frameon=False)

    inset = ax.inset_axes((0.59, 0.59, 0.38, 0.38))
    draw(inset, zoom, markersize=2.0)
    inset.set_xlim(0, max(float(row["m"]) for row in zoom) * 1.02)
    inset.set_ylim(max(0, n / 2 - 3), n / 2)
    inset.set_xlabel(r"$M$", labelpad=1)
    inset.set_ylabel(r"$S$", labelpad=1)
    inset.legend().remove()
    save_figure(fig, out_pdf)


def plot_fig2(data_dir: str | Path, out_pdf: str | Path, *, write_fit_summary: bool = True) -> None:
    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig2_table_with_page.csv")
    n = _arr(rows, "n")
    M_n = _arr(rows, "M_n")
    S_n = _arr(rows, "S_n")
    log2_M = _arr(rows, "log2_M_n") if "log2_M_n" in rows[0] else np.log2(M_n)
    page = _arr(rows, "Page_exact_bits") if "Page_exact_bits" in rows[0] else np.array([exact_page_entropy_bits(int(x)) for x in n])

    fit_logM = _fit_line(n, log2_M)
    fit_S = _fit_line(n, S_n)
    if write_fit_summary:
        fit_rows = [{"quantity": "log2_M_n", **fit_logM}, {"quantity": "S_n", **fit_S}]
        write_rows(data_dir / "fig2_linear_fit_summary_from_csv.csv", fit_rows, fit_rows[0].keys())

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    x = np.linspace(n.min() - 1, n.max() + 1, 200)
    ax.plot(n, log2_M, "o", label=r"$\log_2 M_n$")
    ax.plot(x, fit_logM["slope"] * x + fit_logM["intercept"], "--", linewidth=1)
    ax.plot(n, S_n, "s", label=r"$S_n$")
    ax.plot(x, fit_S["slope"] * x + fit_S["intercept"], "--", linewidth=1)
    ax.plot(n, page, ":", label="Page average")
    ax.set_xlabel(r"number of qubits $n$")
    ax.set_ylabel("bits")
    ax.set_xlim(n.min() - 1, n.max() + 1)
    ax.set_ylim(0, max(log2_M.max(), page.max()) + 1)
    ax.legend(frameon=False, loc="upper left")
    save_figure(fig, out_pdf)


def plot_fig3(data_dir: str | Path, out_pdf: str | Path, *, bins: int = 500) -> None:
    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig3_spectrum.csv")
    eigvals = _arr(rows, "lambda")
    eigvals.sort()
    bulk = eigvals[:-1]
    positive_bulk = bulk[bulk > 0]
    if positive_bulk.size == 0:
        raise ValueError("fig3_spectrum.csv contains no positive bulk eigenvalues.")

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    ax.hist(positive_bulk, bins=bins)
    ax.set_xlabel(r"eigenvalue $\lambda$ excluding $\lambda_0$")
    ax.set_ylabel("count")
    ax.set_xlim(0, positive_bulk.max() * 1.05)
    save_figure(fig, out_pdf)


def _draw_fig4(axis, rows, label_data: bool) -> None:
    m = _arr(rows, "m")
    mean = _arr(rows, "entropy_mean")
    sem = _arr(rows, "entropy_sem")
    D = _arr(rows, "D_NM")
    T = _arr(rows, "T_NM")
    axis.errorbar(m, mean, yerr=sem, fmt="o", markersize=2.2, capsize=1.2, label=r"numerical $S_{N,M}$" if label_data else None)
    axis.plot(m, D, ":", label=r"$D_{N,M}$" if label_data else None)
    axis.plot(m, T, "--", label=r"$T_{N,M}$" if label_data else None)


def plot_fig4(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None) -> None:
    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig4_summary.csv")
    zoom = read_csv_rows(data_dir / "fig4_zoom_summary.csv")
    if n is None:
        n = _infer_n_from_m(rows)

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.8))
    _draw_fig4(ax, rows, True)
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel("entropy")
    ax.set_xlim(0, 1 << n)
    ax.set_ylim(0, n / 2)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.35), ncol=3)

    inset = ax.inset_axes((0.59, 0.59, 0.38, 0.38))
    _draw_fig4(inset, zoom, False)
    inset.set_xlim(0, max(float(row["m"]) for row in zoom) * 1.02)
    inset.set_ylim(max(0, n / 2 - 3), n / 2)
    inset.set_xlabel(r"$M$", labelpad=1)
    inset.set_ylabel(r"$S$", labelpad=1)
    save_figure(fig, out_pdf)


def plot_fig5(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None) -> None:
    data_dir = Path(data_dir)
    random_rows = read_csv_rows(data_dir / "fig5_random_qft_summary.csv")
    union_rows = read_csv_rows(data_dir / "fig5_almost_prime_unions.csv")
    if n is None:
        n = _infer_n_from_m(random_rows)

    m = _arr(random_rows, "m")
    pos = _arr(random_rows, "position_mean")
    pos_sem = _arr(random_rows, "position_sem")
    qft = _arr(random_rows, "qft_mean")
    qft_sem = _arr(random_rows, "qft_sem")
    union_m = _arr(union_rows, "m")
    union_pos = _arr(union_rows, "position_entropy")
    union_qft = _arr(union_rows, "qft_entropy")

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    ax.errorbar(m, pos, yerr=pos_sem, fmt="o", markersize=2, capsize=1.2, label=r"random subset")
    ax.errorbar(m, qft, yerr=qft_sem, fmt=".", markersize=2, capsize=1.2, label=r"random subset after QFT")
    ax.plot(union_m, union_pos, "+", markersize=6, label=r"$U_{N,k}$")
    ax.plot(union_m, union_qft, "x", markersize=5, label=r"QFT $U_{N,k}$")
    ax.set_xlabel(r"support size $M$")
    ax.set_ylabel("entropy")
    ax.set_xlim(0, 1 << n)
    ax.set_ylim(0, n / 2)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.45), ncol=2)
    save_figure(fig, out_pdf)


def _order_key(value) -> float:
    text = str(value).strip().lower()
    if text in {"inf", "infinity", "np.inf"}:
        return float("inf")
    return float(value)


def _rows_for_order(rows: Sequence[dict], order: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    filtered = [row for row in rows if (_order_key(row["order"]) == order or (np.isinf(order) and np.isinf(_order_key(row["order"]))))]
    return _arr(filtered, "m"), _arr(filtered, "mean"), _arr(filtered, "sem")


def _draw_fig6(axis, rows: Sequence[dict], markers: bool = True) -> None:
    orders = [1.0, 2.0, np.inf]
    labels = {1.0: r"$S^{(1)}$", 2.0: r"$S^{(2)}$", np.inf: r"$S^{(\infty)}$"}
    styles = {1.0: "o-", 2.0: "s--", np.inf: "^:"}
    for order in orders:
        m, mean, sem = _rows_for_order(rows, order)
        fmt = styles[order] if markers else styles[order].replace("o", "").replace("s", "").replace("^", "")
        axis.errorbar(m, mean, yerr=sem, fmt=fmt, markersize=2.2, capsize=1.2, label=labels[order])


def plot_fig6(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None) -> None:
    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig6_renyi_summary.csv")
    zoom = read_csv_rows(data_dir / "fig6_renyi_zoom_summary.csv")
    if n is None:
        n = _infer_n_from_m(rows)
    zoom_stop = int(max(float(row["m"]) for row in zoom))

    apply_journal_style()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.7, 6.0), sharey=False)
    _draw_fig6(ax1, rows)
    ax1.set_xlabel(r"support size $M$")
    ax1.set_ylabel("Rényi entropy")
    ax1.set_xlim(0, 1 << n)
    ax1.set_ylim(0, n / 2)
    ax1.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.42))

    _draw_fig6(ax2, zoom)
    ax2.set_xlabel(r"support size $M$")
    ax2.set_ylabel("Rényi entropy")
    ax2.set_xlim(0, zoom_stop)
    ax2.set_ylim(max(0, n / 2 - 3), n / 2)
    ax2.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.42))
    fig.subplots_adjust(hspace=0.55)
    save_figure(fig, out_pdf)


def _families(rows: Sequence[dict]) -> list[str]:
    seen = []
    for fam in _str_arr(rows, "state_family" if "state_family" in rows[0] else "family"):
        if fam not in seen:
            seen.append(str(fam))
    return seen


def plot_fig7(data_dir: str | Path, out_pdf: str | Path, *, n: int | None = None, bins: int = 60) -> None:
    data_dir = Path(data_dir)
    rows = read_csv_rows(data_dir / "fig7_partition_samples.csv")
    summary_path = data_dir / "fig7_partition_summary.csv"
    page = None
    if summary_path.exists():
        summary = read_csv_rows(summary_path)
        if n is None and summary and summary[0].get("n", "") != "":
            n = int(float(summary[0]["n"]))
        if summary and "Page_exact_bits" in summary[0]:
            page = float(summary[0]["Page_exact_bits"])
    if n is None:
        n = 20
    if page is None:
        page = exact_page_entropy_bits(n)

    families = _families(rows)
    apply_journal_style()
    fig, axes = plt.subplots(len(families), 1, figsize=(4.7, 2.9 * len(families)))
    if len(families) == 1:
        axes = [axes]
    for ax, family in zip(axes, families):
        values = np.asarray([float(row["entropy"]) for row in rows if str(row["state_family"]) == family], dtype=float)
        ax.hist(values, bins=bins)
        ax.axvline(values.mean(), linestyle="--", linewidth=1, label="sample mean")
        ax.axvline(page, linestyle=":", linewidth=1, label="Page average")
        label = "complex Haar" if family == "complex_haar" else family.replace("_", " ")
        ax.set_xlabel("entropy")
        ax.set_ylabel("count")
        ax.set_title(label)
        ax.legend(frameon=False)
    fig.subplots_adjust(hspace=0.55)
    save_figure(fig, out_pdf)


def _n_from_fig8_path(path: Path) -> int:
    match = re.search(r"fig8_samples_n(\d+)\.csv$", path.name)
    if not match:
        raise ValueError(f"cannot infer n from {path}")
    return int(match.group(1))


def plot_fig8(data_dir: str | Path, out_pdf: str | Path, *, n_values: Sequence[int] | None = None, bins: int = 120) -> None:
    data_dir = Path(data_dir)
    if n_values:
        files = [data_dir / f"fig8_samples_n{int(n)}.csv" for n in n_values]
    else:
        files = sorted(data_dir.glob("fig8_samples_n*.csv"), key=_n_from_fig8_path)
    if not files:
        raise FileNotFoundError(f"no fig8_samples_n*.csv files found in {data_dir}")

    apply_journal_style()
    fig, axes = plt.subplots(len(files), 1, figsize=(4.8, 3.0 * len(files)))
    if len(files) == 1:
        axes = [axes]
    for ax, path in zip(axes, files):
        rows = read_csv_rows(path)
        n = _n_from_fig8_path(path)
        families = _families(rows)
        arrays = []
        for family in families:
            key = "family" if "family" in rows[0] else "state_family"
            values = np.asarray([float(row["entropy"]) for row in rows if str(row[key]) == family], dtype=float)
            arrays.append((family, values))
        lo = min(values.min() for _, values in arrays)
        hi = max(values.max() for _, values in arrays)
        bin_edges = np.linspace(lo, hi, bins)
        for family, values in arrays:
            label = "random subsets, natural partition" if family.startswith("random") else "selected subset, all partitions"
            ax.hist(values, bins=bin_edges, alpha=0.65, label=label)
            ax.axvline(values.mean(), linestyle="--" if family.startswith("random") else ":", linewidth=1)
        M = int(float(rows[0].get("M", rows[0].get("m", 0))))
        ax.set_xlabel("entropy")
        ax.set_ylabel("count")
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
    y_fit = float(fit["quadratic_a"]) * x_fit**2 + float(fit["quadratic_b"]) * x_fit + float(fit["quadratic_c"])

    apply_journal_style()
    fig, ax = plt.subplots(figsize=(4.7, 3.6))
    ax.errorbar(x, mean, yerr=sem, fmt="o", capsize=1.5, label="samples")
    ax.plot(x_fit, y_fit, "--", label="quadratic fit")
    ax.axvline(float(fit["log2_M_peak"]), linestyle=":", linewidth=1, label="estimated peak")
    ax.set_xlabel(r"$\log_2 M$")
    ax.set_ylabel(r"$S_{N,M}$")
    ax.legend(frameon=False)
    save_figure(fig, out_pdf)
