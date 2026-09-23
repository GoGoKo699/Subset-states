# Changelog

## Markdown cleanup — 2026-09-22

- Replaced manuscript-readiness reports with a scoped sanity check, research proofs, reproducibility instructions, and primary references in Markdown.
- Corrected the public preprint title and year; aligned the existing rc2 software version labels without creating a release tag.
- Fixed support validation, coefficient-matrix basis ordering, and invalid Rényi inputs.
- Reworked scientific validation to read released data and compare production routines, with explicit failures that survive Python optimization.
- Made all seven CSV figure redraws and reduced simulations write under `generated/`; changed default graphics to PNG.
- Removed PDF duplicates and obsolete text reports. Preserved all released CSV/PNG evidence and the license byte-for-byte, with checksums.


## 1.0.0-rc2 — 2026-08-21

- Added the general power-law support window
  1/2 < γ < 3/4 implied by the exact average-purity formula.
- Distinguished the purity-optimal interior balance γ = 2/3 from the
  dense-side boundary γ = 3/4.
- Clarified that the finite-size peak-support regression is a separate numerical
  observation and not an asymptotic consequence of the purity bound.
- Standardized ordinary Shannon-entropy notation from H₂ to H.
- Kept all seven released figures, Table I values, and numerical datasets
  unchanged.

## 1.0.0-rc1 — 2026-08-18

- Added exact average-purity and residue-class ceiling implementations.
- Added the final cardinality- and residue-matched Figure 5 workflow and data.
- Replaced the old almost-prime plot with the two-panel controlled comparison.
- Added publication-level README, provenance record, citation metadata, tests,
  validation, and release figure files.
- Reclassified the former “greedy” experiment as a best-of-random-candidates
  exploratory calculation outside the manuscript.
