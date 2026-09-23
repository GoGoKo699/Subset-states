# Final repository review — 23 September 2026

**Result: the checks below passed.** This review covers the repository's scientific consistency, manuscript completeness, rendered typography, and citation guidance. It does not establish novelty or reconstruct missing historical numerical searches.

## Science and mathematical fidelity

An independent review found no blocking mathematical error in the exact mean reduced state, average purity, asymptotic expansions and entropy bounds, or residue-ceiling proof. A separate comparison of the typographic revision against the preceding Markdown edition checked all 88 numbered equations, both abstract displays, and the inline mathematics. No changed exponent, missing term, altered inequality, misplaced index, or fraction-scope error was found.

The final edits distinguish the released sampled trajectories from retained historical peak estimates for n = 10–30 in the abstract, introduction, peak discussion, and conclusion. The text no longer suggests that an upward finite-size trend identifies the 3/4 purity boundary as an entropy-peak limit. Equation (12) explicitly states its real-coefficient assumption. These changes leave the reported numerical values and mathematical results intact.

## Rendered layout

The manuscript was inspected in GitHub's actual Markdown preview, including the mean-state formulas highlighted in the reported screenshot. All 90 mathematical displays now use numbered or unnumbered equation blocks with Unicode and supported subscript/superscript tags. Inline mathematics uses matching typography. Long formulas break at logical operations, with explicit grouping around fractions.

At the inspected desktop viewport of 1363 × 936 pixels, the manuscript reading column measured 929 pixels. A full-document DOM check found no horizontal overflow in the paragraphs, equation blocks, table, or images. All seven figures loaded; GitHub retained 307 subscript and 158 superscript elements, and the manuscript contained no code blocks. The equation and reference anchors remain intact.

This records the actual desktop inspection; it is not a claim of testing every browser, viewport, or mobile device. The manuscript remains readable without a mathematics renderer or a document compiler.

## Automated verification

| Check | Result |
|---|---|
| Unit and regression suite | 39 tests passed |
| Exhaustive small systems | 65,550 nonempty supports at n = 2 and n = 4 checked |
| Exact purity vs enumeration | Maximum absolute discrepancy 1.12 × 10⁻¹⁵ or less |
| Retained peak table and derived fits | All 11 rows consistent; maximum discrepancy 3.56 × 10⁻¹⁵ or less |
| Residue-control summaries | All 12,000 observations across 12 groups checked |
| Seed replay | 36 sampled supports reproduced |
| Dense-bulk implementation | 16,660 parameter pairs checked against the stated formula |
| Manuscript inventory | 88 equation numbers, 33 references, seven figures, all original section/statement anchors, and 11 table rows retained |
| Repository integrity | All 51 frozen evidence/license checksums unchanged; local paths and fragment links valid; no prohibited document sources or PDF artifacts |
| Citation/package metadata | YAML and TOML parsed; public-paper title, authors, year, identifier, and DOI matched the official arXiv record |

The dense-bulk check verifies implementation of the stated approximation, not its accuracy as a model of the ensemble. Small-system enumeration and seeded replay have the scope stated in the [scientific validation report](SCIENTIFIC_SANITY_CHECKS.md).

Run the checks with:

```bash
python scripts/run_smoke_tests.py
python scripts/final_scientific_validation.py
python scripts/check_repository.py
```

## Discovery and citation

[llms.md](../llms.md) explains relevant research questions and links readers or retrieval systems to the appropriate manuscript, code, data, and provenance. The README and package metadata use descriptive research keywords. These additions aid interpretation and relevance assessment; they do not guarantee search indexing, ranking, or recommendations.

[CITING.md](../docs/CITING.md) and [CITATION.cff](../CITATION.cff) direct users to Ruge Lin, Germán Sierra, and José I. Latorre, *Arithmetic sequences as quantum states*, arXiv:2501.06292 [quant-ph] (2025), [doi:10.48550/arXiv.2501.06292](https://doi.org/10.48550/arXiv.2501.06292). The [official record](https://arxiv.org/abs/2501.06292), checked on this date, lists v1 submitted on 10 January 2025. The later repository manuscript has a different working title; users should also identify the repository commit for revised results, code, or data.

The remaining scientific and provenance questions are recorded in [SANITY_CHECK.md](../docs/SANITY_CHECK.md) and [PROVENANCE.md](../PROVENANCE.md).
