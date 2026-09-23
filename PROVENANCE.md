# Evidence and provenance

This document describes the evidence present in the repository, rather than treating earlier review verdicts as evidence. The cleanup started from commit [4cb99f1](https://github.com/GoGoKo699/Subset-states/tree/4cb99f109194124f500d75dd065ff43581a2682d), dated 21 August 2026.

## Evidence hierarchy

| Material | What can be checked | Remaining limitation |
|---|---|---|
| Exact mean-state, purity, diagonal-entropy and residue formulas | Complete proofs in the [manuscript](PAPER.md), companion [research notes](docs/RESEARCH.md), independent enumeration, production-code comparisons | Enumeration alone is not a proof for arbitrary size |
| `data/table_i_peaks.csv` | Transcription of Table I in the [2025 preprint](https://arxiv.org/html/2501.06292v1), arithmetic regressions, agreement with released derived tables | No complete original global-search records, raw peak samples, search uncertainty or full execution provenance |
| `outputs/peak_verification/` | Local grids and quadratic fits for $n = 10, 12, 14, 16, 18, 20$ | These are neighborhood checks around given centers, not global rediscovery; files named `samples` hold aggregate statistics at each $M$ |
| `outputs/fig1/` | Released individual entropy samples and summaries | Support labels and a full original execution manifest are absent |
| `outputs/fig3/` | One $n = 24$ reduced-state spectrum | One sample does not establish typicality or a random-matrix limit |
| `outputs/fig4/`, `outputs/fig6/` | Summary curves and plotting inputs | Complete per-state records are not released here |
| `data/residue_matched_*` | 12,000 paired entropy observations, 12 group summaries, seeds, and residue populations | Raw support labels are absent; a limited seed replay checks implementation consistency, not every stored observation |
| `outputs/fig7/` | Partition entropy samples for two selected $n = 20$ states | Does not establish a uniform statement about all supports or all cuts |

CSV files and retained PNGs are unchanged by the cleanup. Their SHA-256 checksums are recorded in [the evidence manifest](validation/evidence_sha256.json), alongside the original license. This freezes the inherited evidence; it does not certify the original data-generation process.

## Retained peaks and finite fits

All eleven Table I rows, for even $n = 10$ through $30$, appear in the public 2025 preprint. Their publication establishes the source of the numbers, not independent reproducibility of the global searches. In particular, the complete intermediate searches at $n = 22, 24, 26, 28, 30$ were not retained.

The unweighted regressions recomputed from the CSV are:

| Range | Entropy fit | Log-support fit |
|---|---|---|
| $n = 10$–$30$ | $\widehat{S} = 0.509300n - 0.990091$ | $\log_2 \widehat{M} = 0.703541n - 0.357734$ |
| $n = 10$–$20$ | $\widehat{S} = 0.514457n - 1.065190$ | $\log_2 \widehat{M} = 0.693450n - 0.210116$ |

The slopes differ by roughly 1.01% and 1.43%, respectively, relative to the full-range slopes. This is a descriptive sensitivity check. It supplies neither missing uncertainty estimates nor a proof that the historical estimates are unbiased. Regression standard errors describe scatter of the rounded table entries, not uncertainty in the original search.

The effective ratio $\log_2 \widehat{M}/n$ grows from about 0.674 to 0.694 across the table. These finite observations establish no limiting exponent. In particular, the purity crossover at exponent $3/4$ does not predict convergence of the entropy maximizer to that exponent. An entropy slope above $1/2$ also cannot persist asymptotically, since $S \le n/2$.

## Residue-matched comparisons

The structured supports contain labels $x$ with $1 \le \Omega(x) \le k$, where $\Omega$ counts prime factors with multiplicity and $k = 1, 2, 3$. For $n = 14$, each support is compared with four ensembles matching populations modulo $2^t$, for $t = 0, 1, 2, 3$. Each group contains 1,000 samples, and the same sampled support is evaluated before and after the QFT.

A remaining deficit means that uniform placement conditional on those counts does not reproduce the structured entropy. It does not identify the cause or a unique fingerprint of primality. Zero sampled null entropies below a structured value is a finite rank observation, not a zero tail probability. The exact computational-basis residue ceiling is not automatically a Fourier-basis ceiling.

## Historical material

The old manuscript-readiness review, bibliography audit, duplicated text validation reports, and PDF copies of PNG figures were removed from the active tree. They remain available at the starting commit above. No manuscript source or complete revised bibliography was present in that tree during the initial cleanup.

On 23 September 2026, the complete revised manuscript package dated 21 August 2026 was recovered from the project files. Its text, mathematics, table, seven figures and captions, four appendices, declarations, and 33 references are now integrated in [PAPER.md](PAPER.md). The [integration record](paper/INTEGRATION.md) identifies the source snapshots and explains the scientific annotations. This recovery supplies the complete written manuscript; it does not supply the missing original peak-search records.

The PNG images under `paper/figures/` are direct renderings of that manuscript's seven supplied figures, with separate checksums in [the manuscript manifest](validation/manuscript_manifest.json). They preserve the manuscript presentation, including its revised Figure 2. The released CSV/PNG evidence under `data/` and `outputs/`, and its original checksum manifest, remain unchanged.

The historical best-of-random-candidates experiment remains under [exploratory](exploratory/best_of_random_candidates/README.md). It selects among independent random candidates; it is not a greedy construction.
