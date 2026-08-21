# Support-size entanglement trajectories of random subset states

Code, data, released figures, and validation for the manuscript
**“Support-size entanglement trajectories of random subset states.”**
The manuscript updates the existing preprint record [arXiv:2501.06292](https://arxiv.org/abs/2501.06292).

Equal-positive-amplitude subset states exhibit a rise–peak–fall entanglement
trajectory as their computational-basis support grows. This repository contains
the fixed-cardinality ensemble calculations, exact-moment checks, numerical
figure scripts, and constrained reference ensembles used in the paper.

## Key results represented in this repository

- Exact ensemble-mean reduced state for uniformly random supports of fixed size.
- Exact average purity across a balanced bipartition, yielding a rigorous
  power-law support window, an optimal interior balance, and a dense-side
  boundary.
- Numerical rise–peak–fall trajectories and retained peak estimates through
  $n=30$.
- Exact hypergeometric diagonal-entropy and residue-class entropy bounds.
- Cardinality-, parity-, mod-4-, and mod-8-matched reference ensembles for
  almost-prime supports, evaluated before and after a quantum Fourier transform.

## Analytic scale hierarchy

Let $N=2^n$ and write a power-law support size as $M=cN^\gamma$, with
$c>0$. The exact fixed-cardinality average-purity formula has the leading
structure

```math
\overline{P}_{N,M}
=
2N^{-1/2}
+c^{-1}N^{-\gamma}
+c^2N^{2\gamma-2}
+\text{lower-order terms}.
```

**Near-maximal window.** Every fixed exponent
$\frac{1}{2}<\gamma<\frac{3}{4}$ gives

```math
\overline{S}_{N,M}
\geq
\frac{n}{2}-1-o(1).
```

**Purity-optimal interior scale.** The exponent
$\gamma=\frac{2}{3}$ uniquely balances the sparse
$N^{-\gamma}$ correction and the dense rectangle
$N^{2\gamma-2}$ correction. Optimizing the remaining prefactor gives

```math
M
=
2^{-1/3}N^{\frac{2}{3}}+O(1).
```

**Dense-side boundary.** The exponent $\gamma=\frac{3}{4}$ is the point at
which the rectangle contribution enters at the same $N^{-1/2}$ order as the
balanced-cut background.

The numerical maximum of the mean von Neumann entropy is a separate quantity.
Over $n=10,\ldots,30$, the retained estimates give

```math
\log_2 \widehat{M}_n
=
0.703541\,n-0.357734.
```

The pointwise effective exponent is

```math
\gamma_n^{\mathrm{eff}}
=
\frac{\log_2 \widehat{M}_n}{n},
\qquad
\gamma_{10}^{\mathrm{eff}}\approx 0.674,
\qquad
\gamma_{30}^{\mathrm{eff}}\approx 0.694.
```

It increases over the retained range. This motivates comparison with the
$\frac{3}{4}$ dense-side boundary but does not establish an asymptotic peak law.

## Repository map

```text
subset_states/     reusable scientific routines
scripts/           figure, validation, and local-verification entry points
data/              Table I and released matched-null datasets
outputs/            manuscript-ready figures and generated outputs
tests/              lightweight and exhaustive small-system regression tests
validation/         final independent validation reports
exploratory/        historical experiments not used in the manuscript
docs/               audit and release documentation
```

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 scripts/run_smoke_tests.py
```

Run the independent publication validation:

```bash
python3 scripts/final_scientific_validation.py
```

Regenerate the fast deterministic figures and verify that all seven released
PDFs are present:

```bash
python3 scripts/reproduce_publication_figures.py
```

A reduced computational smoke run is available through:

```bash
python3 scripts/reproduce_publication_figures.py --smoke
```

The full publication computations can be launched with `--full`, but several
runs are intentionally expensive.

## Publication figures

| Figure | Scientific content | Main script | Released output |
|---|---|---|---|
| 1 | Support-size trajectory and cut comparison | `scripts/fig1_concentration.py` | `outputs/fig1/fig1_concentration.pdf` |
| 2 | Peak entropy and support-size scaling | `scripts/fig2_peak_scaling.py` | `outputs/fig2/fig2_peak_scaling.pdf` |
| 3 | Reduced-state spectrum near the peak | `scripts/fig3_spectral.py` | `outputs/fig3/fig3_spectral_bulk.pdf` |
| 4 | Sparse and dense approximations | `scripts/fig4_approximation.py` | `outputs/fig4/fig4_approximation.pdf` |
| 5 | Cardinality and residue controls | `scripts/fig5_qft_residue_controls.py` | `outputs/fig5/fig5_qft_residue_controls.pdf` |
| 6 | Rényi trajectories | `scripts/fig6_renyi.py` | `outputs/fig6/fig6_renyi.pdf` |
| 7 | Entropy over balanced cuts | `scripts/fig7_partitions.py` | `outputs/fig7/fig7_partitions.pdf` |

Figure 5 can be redrawn exactly from the released CSV files. The matched-null
raw dataset contains 12,000 paired computational/Fourier samples.

## Table I and local verification

`data/table_i_peaks.csv` stores the retained peak estimates used in the paper.
The command

```bash
python3 scripts/peak_scaling_verification.py --n-values 10 12 14
```

performs a local neighbourhood check around selected tabulated values. It is not
a reconstruction of the original global search. Full technical provenance and
the robustness comparison are documented in [PROVENANCE.md](PROVENANCE.md).

## Figure 5 constrained ensembles

The released files are:

```text
data/fig5_random_qft_summary.csv
data/fig5_almost_prime_unions.csv
data/residue_matched_summary.csv
data/residue_matched_samples.csv
data/residue_entropy_bounds.csv
```

Regenerate the null samples with:

```bash
python3 scripts/run_residue_controls.py
python3 scripts/fig5_qft_residue_controls.py
```

## Reproducibility levels

1. **Fast:** unit tests, exact small-system checks, Figure 2, and Figure 5 redraw.
2. **Standard:** regenerate selected numerical figures with reduced or custom
   workloads.
3. **Computational:** use publication defaults for the full sampling runs and
   local high-dimensional peak checks.

All stochastic scripts expose seeds. The tested software environment is recorded
in `requirements-tested.txt`.

## Citation

Use the citation metadata in `CITATION.cff`. The preferred paper citation is the
updated arXiv record and, once available, the journal version.

## License

The repository is released under the MIT License. See `LICENSE`.
