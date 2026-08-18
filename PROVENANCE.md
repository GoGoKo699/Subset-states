# Numerical provenance and robustness notes

This document records the evidence hierarchy behind the publication release.
It is intentionally more technical than the repository landing page.

## Table I peak estimates

`data/table_i_peaks.csv` contains the retained production estimates
$(\widehat{M}_n,\widehat{S}_n)$ for even $n=10,\ldots,30$.

The complete intermediate global-search records for
$n=22,24,26,28,30$ were not retained. Repeating those searches with the
released dense exact-entropy implementation is computationally expensive. The
recorded point estimates are therefore preserved as historical production data,
not reconstructed raw data.

The repository supplies `scripts/peak_scaling_verification.py`, which samples a
local window around a recorded $\widehat{M}_n$ and fits a quadratic in
$\log_2 M$. This verifies neighbourhood consistency; it is not a global
rediscovery procedure.

## Robustness of the finite-size regressions

Using all eleven retained rows, $n=10,\ldots,30$,

$$
\widehat{S}_n=0.509300\,n-0.990091,
$$

$$
\log_2\widehat{M}_n=0.703541\,n-0.357734.
$$

Restricting the fit to the computationally accessible range
$n=10,\ldots,20$ gives

$$
\widehat{S}_n=0.514457\,n-1.065190,
$$

$$
\log_2\widehat{M}_n=0.693450\,n-0.210116.
$$

The slope changes are approximately 1.01% for peak entropy and 1.43% for peak
support size. The unreconstructed larger-size rows extend the observed trend but
do not create or materially change it. Neither regression is claimed as an
asymptotic theorem.

## Exact analytical validation

The release exposes code for:

- the exact fixed-cardinality ensemble-mean reduced state;
- the exact ensemble-average purity;
- the hypergeometric diagonal-entropy curve;
- the residue-class entropy ceiling
  $S(A)\leq n/2-t+H_2(\mathbf{p}^{(t)})$.

`tests/test_exact_purity_and_residue.py` and
`scripts/final_scientific_validation.py` independently check these formulas by
exhaustively enumerating all 65,535 nonempty supports at $n=4$.

## Residue-matched reference ensembles

For each almost-prime support $U_{N,k}$ with $k=1,2,3$, the released
Figure 5 analysis uses four nested null ensembles:

| matched low bits $t$ | retained information |
|---:|---|
| 0 | support cardinality only |
| 1 | even/odd populations |
| 2 | populations modulo 4 |
| 3 | populations modulo 8 |

For every $(k,t)$, 1,000 supports are independently sampled while preserving
the complete residue-count vector. The same support is evaluated in the
computational and Fourier bases, giving paired observations. Seeds are recorded
in `data/residue_matched_summary.csv`; all 12,000 paired samples are in
`data/residue_matched_samples.csv`.

The residual deficit after matching a constraint means “not explained by the
matched information under uniform placement.” It is not claimed to be a unique
fingerprint of primality or to isolate all higher-order arithmetic correlations.

## Figure files and generated data

The manuscript-ready vector PDFs in `outputs/fig1/` through `outputs/fig7/` are
the files used for the final manuscript candidate. Existing historical output
CSVs in the GitHub repository may be retained when applying this release as an
overlay. Figure 5’s final raw and summary data are included directly in this
release because the constrained-null calculation is new and central to the
revised manuscript.

## Historical best-of-random-candidates experiment

The former “greedy” Figure 8 experiment is not part of the manuscript. It
screened independent random candidates and selected the one with the largest
mean balanced-cut entropy; it was not a greedy construction. The code is kept
under `exploratory/best_of_random_candidates/` with corrected terminology.
