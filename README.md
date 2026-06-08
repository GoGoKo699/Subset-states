# Subset-state entanglement code

This repository contains the numerical scripts used for the manuscript **Support-size entanglement trajectories of random subset states**.

The code has two roles:

1. regenerate the manuscript figures from numerical sampling; and
2. provide a transparent, local verification script for the Table-I peak values.

The peak values in Table I are treated as production data. The verification script samples in a window around the tabulated peak support size \(M_n\), then fits a quadratic curve in \(\log_2 M\). It is not a global rediscovery of the peaks.

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Run the smoke test:

```bash
python3 scripts/run_smoke_tests.py
```

Expected output:

```text
smoke tests passed
```

## Important data files

```text
data/table_i_peaks.csv
```

contains the production Table-I values:

```text
n, M_n, S_n
10, 107, 4.072
...
30, 1836685, 14.263
```

```text
data/peak_verification_schedule.csv
```

contains a modifiable schedule for the local peak checks. Edit `samples`, `points`, and `span` there, or override them from the command line.

## Regenerate figures from computation

Main figures:

```bash
python3 scripts/fig1_concentration.py
python3 scripts/fig2_peak_scaling.py
python3 scripts/fig3_spectral.py
python3 scripts/fig4_approximation.py
python3 scripts/fig5_qft_almost_primes.py
```

Appendix figures:

```bash
python3 scripts/fig6_renyi.py
python3 scripts/fig7_partitions.py
python3 scripts/fig8_greedy_partitions.py
```

The scripts write CSV data and PDF/PNG figures under `outputs/`.

## Redraw figures from existing CSV files

These scripts do not recompute entropy values:

```bash
python3 scripts/plot_fig1_from_csv.py
python3 scripts/plot_fig2_from_csv.py
python3 scripts/plot_fig3_from_csv.py
python3 scripts/plot_fig4_from_csv.py
python3 scripts/plot_fig5_from_csv.py
python3 scripts/plot_fig6_from_csv.py
python3 scripts/plot_fig7_from_csv.py
python3 scripts/plot_fig8_from_csv.py
```

## Figure 2 and Table-I scaling

Figure 2 is generated directly from `data/table_i_peaks.csv`:

```bash
python3 scripts/fig2_peak_scaling.py
```

Outputs:

```text
outputs/fig2/fig2_peak_scaling.pdf
outputs/fig2/fig2_peak_scaling.png
outputs/fig2/fig2_table_with_page.csv
outputs/fig2/fig2_linear_fit_summary.csv
```

The linear-fit summary contains the coefficients for \(S_n\) and \(\log_2 M_n\).

## Local verification of Table-I peaks

Use:

```bash
python3 scripts/peak_scaling_verification.py --n-values 10 12 14
```

This samples support sizes in a local window around the tabulated peak \(M_n\), computes the average entropy at each sampled support size, and fits a quadratic curve in \(\log_2 M\).

For a larger but still manageable run:

```bash
python3 scripts/peak_scaling_verification.py --n-values 10 12 14 16 18 20
```

To use one uniform setting for all selected rows:

```bash
python3 scripts/peak_scaling_verification.py \
  --n-values 10 12 14 16 \
  --samples 500 \
  --points 31 \
  --span 0.25 \
  --label check
```

To inspect the planned workload without running it:

```bash
python3 scripts/peak_scaling_verification.py --all-table --dry-run
```

Outputs:

```text
outputs/peak_verification/peak_verification_n*_samples.csv
outputs/peak_verification/peak_verification_n*_fit.csv
outputs/peak_verification/peak_verification_n*.pdf
outputs/peak_verification/peak_verification_summary.csv
outputs/peak_verification/peak_verification_summary.pdf
```

To redraw one verification plot from CSV:

```bash
python3 scripts/plot_peak_verification_from_csv.py --n 14
```

## Colour palette

The default palette is deliberately cool-toned and avoids orange/red/yellow. To modify it, edit:

```text
subset_states/plotting.py
```

The palette is defined in `COOL_PALETTE`.

## Reproducibility notes

- All stochastic scripts accept a `--seed` argument.
- The local peak-verification script is designed to verify the tabulated values in a neighbourhood of the known peaks, not to repeat the original high-cost global search.
- Large values such as `n=26,28,30` can be expensive with the dense exact entropy calculation. Use the schedule file to lower `samples` and `points`, or run those checks on a machine with sufficient memory and time.


### Matplotlib note

The figure scripts do not import `mpl_toolkits.axes_grid1`; inset panels use `Axes.inset_axes` through `subset_states.plotting.make_inset_axis`.  This avoids a common Ubuntu issue where user-site Matplotlib and system `mpl_toolkits` are mixed.  For the cleanest environment, run the code inside a Python virtual environment.
