# Changelog

## Post-release metadata — 2026-09-23

- Aligned package and citation version labels to the published v1.0.0 release and recorded its publication date.
- Added release navigation and a citation example tied to the exact release commit. The published tag and its archived files are unchanged; that snapshot retains the earlier rc2 metadata labels.

## 1.0.0 — 2026-09-23

[Published release](https://github.com/GoGoKo699/Subset-states/releases/tag/v1.0.0), commit [`d3729b4`](https://github.com/GoGoKo699/Subset-states/commit/d3729b4d4c5e6f4441da469d626c2d998500cfb8).

- Integrated the complete revised manuscript in Markdown, including all 88 numbered equations, two abstract displays, seven figures, four appendices, and 33 references.
- Typeset manuscript and supporting-document formulas with GitHub-native mathematics, with the rendered desktop layout checked on GitHub.
- Added reader and LLM navigation, arXiv citation guidance, and manuscript completeness checks.
- Passed 39 unit/regression tests, the scoped scientific validator, and repository checks; preserved all 51 frozen evidence/license checksums.

## Markdown cleanup — 2026-09-22

- Replaced manuscript-readiness reports with a scoped sanity check, research proofs, reproducibility instructions, and primary references in Markdown.
- Corrected the public preprint title and year; aligned the existing rc2 software version labels without creating a release tag.
- Fixed support validation, coefficient-matrix basis ordering, and invalid Rényi inputs.
- Reworked scientific validation to read released data and compare production routines, with explicit failures that survive Python optimization.
- Made all seven CSV figure redraws and reduced simulations write under `generated/`; changed default graphics to PNG.
- Removed PDF duplicates and obsolete text reports. Preserved all released CSV/PNG evidence and the license byte-for-byte, with checksums.


## 1.0.0-rc2 — 2026-08-21

- Added the general power-law support window
  $`\frac{1}{2} < \gamma < \frac{3}{4}`$ implied by the exact average-purity formula.
- Distinguished the purity-optimal interior balance $`\gamma = 2/3`$ from the
  dense-side boundary $`\gamma = 3/4`$.
- Clarified that the finite-size peak-support regression is a separate numerical
  observation and not an asymptotic consequence of the purity bound.
- Standardized ordinary Shannon-entropy notation from $`H_2`$ to $`H`$.
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
