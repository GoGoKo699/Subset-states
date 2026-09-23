# Sanity check and cleanup — 22 September 2026

**Assessment: the core exact results survive this audit; the repository does not establish a complete publication-ready contribution.** The remaining high-$n$ provenance and novelty questions are substantive. They are now visible in the landing page and research notes.

The baseline was commit [4cb99f109194124f500d75dd065ff43581a2682d](https://github.com/GoGoKo699/Subset-states/tree/4cb99f109194124f500d75dd065ff43581a2682d). Independent mathematical, numerical-code, validator, and primary-literature reviews were used. Prior automated PASS reports were treated as historical claims to inspect.

## Findings and actions

| Finding | Action |
|---|---|
| Exact mean reduced state, average purity, entropy lower bound, hypergeometric diagonal entropy, and residue ceiling check out | Added derivations and explicit assumptions in [research notes](RESEARCH.md), including the exact cubic equation for the purity minimizer |
| Old validation hardcoded Table I and mainly checked copies of formulas | Validator now reads actual CSVs, compares derived tables, and checks production routines against independent calculations |
| Old acceptance assertions disappear under optimized Python | Replaced with explicit failures; mutation tests exercise both `-O` and `-OO` |
| QFT routine accepted duplicate and negative support labels; for example, duplicate labels could yield a state with squared norm $1/2$ | Validate labels before indexing or integer conversion; added rejection and normalization regressions |
| State-vector and support constructors used inconsistent subsystem basis ordering | Aligned ordering and documented the complex partial-trace convention; valid historical entropy spectra are invariant under the corrected row/column permutations |
| Invalid Rényi orders/spectra and extreme orders could silently misbehave | Added input checks and stable entropy evaluation, with targeted regressions |
| Default figure “reproduction” redrew only two figures and checked that old PDFs existed; smoke mode overwrote some released outputs | Redraw all seven from CSV; put all default fresh outputs under `generated/`; test destination isolation and missing-input failures |
| Repository cited a 2026 renamed arXiv paper that is not the current public record | Corrected the citation to *Arithmetic sequences as quantum states* (2025); documented the later title as a working title |
| Old reports claimed submission readiness and audited absent manuscript/build files | Removed those reports from the active tree; linked the baseline commit for history; replaced them with scoped Markdown documentation |
| No TeX sources were present, but documentation contained TeX notation and the repository retained 14 PDF figures | Kept prose mathematics in Markdown, removed PDF duplicates, and made plot defaults PNG without an external TeX dependency |

## What was verified

The initial ten tests passed before editing. All 39 tests in the updated integrated suite passed. Exact/data checks and figure workflows are recorded in [the scientific report](../validation/SCIENTIFIC_SANITY_CHECKS.md) and [reproducibility instructions](REPRODUCIBILITY.md).

The mathematical audit rederived the identities and separately checked integer/rational enumeration at $n = 4$. The maintained validator exhaustively checks 65,550 supports in total across $n = 2$ and $n = 4$, including comparisons with the production coefficient matrices, entropies, mean matrices, purity formula, diagonal expectation, residue ceilings, and entropy lower bound. QFT amplitudes receive an explicit small-system DFT check.

The data check reconstructs all 12,000 released residue observations into their summaries and checks residue bounds. It replays three seeded supports per group, 36 in total. That is a limited seed replay, not a regeneration of every raw observation.

All seven CSV-based figure redraws and all seven reduced-run outputs completed successfully. The complete original-size simulations and missing global peak searches were not run. SHA-256 comparison against the baseline confirms all 50 retained CSV/PNG evidence files and the original license are unchanged.

## What remains unresolved

- Historical peak values are publicly traceable to the 2025 preprint, but the complete global-search records and uncertainty estimates are absent. Later local checks cannot supply them.
- The exact purity minimizer does not locate the mean von Neumann entropy maximum. No asymptotic entropy-peak exponent is proved.
- The dense-bulk approximation lacks a uniform error theorem. Checking its implementation or algebra does not validate its accuracy throughout the parameter range.
- A single spectrum and two selected states do not prove concentration or universality.
- Existing sparse-state and fixed-size subset-state literature overlaps the broad framing. The targeted [literature check](REFERENCES.md) does not establish novelty of the revised contribution.

The cleanup makes the project inspectable and reproducible within the stated limits. It does not replace peer review or recover missing experiments.
