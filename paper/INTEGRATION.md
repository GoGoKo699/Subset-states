# Manuscript integration record

[PAPER.md](../PAPER.md) is the complete Markdown edition of *Support-size entanglement trajectories of random subset states*, by Ruge Lin, Germán Sierra, and José I. Latorre. It incorporates the recovered revised manuscript dated 21 August 2026, rather than reconstructing a paper from the repository's shorter research notes.

## Source and scope

The source package is `Subset_States_revised_manuscript.zip`. Its manuscript source and 24-page rendered manuscript match the separately supplied `Subset_States_main_revised.tex` and `Subset_States_main_revised.pdf` byte for byte. The package also supplies the complete bibliography and seven figure files. It is the latest complete manuscript package recovered from the project files during this integration, on 23 September 2026.

The [machine-readable manifest](../validation/manuscript_manifest.json) records the source filenames, sizes, SHA-256 hashes, section destinations, equation labels, figure hashes, and reference keys. These identify the conversion input without introducing a document compiler or retaining source-format artifacts in the repository. The Markdown edition is the maintained paper; no external source package is needed to read it.

| Source material | Markdown destination | Coverage |
|---|---|---|
| Title, authors, affiliations, abstract, keywords | [Front matter](../PAPER.md) | Full text, including two unnumbered abstract formulas |
| Main text | [Sections 1–7](../PAPER.md#sec-introduction) | All sections and subsections |
| Mathematical statements and arguments | [Section 3](../PAPER.md#sec-exact-moments), [Section 6](../PAPER.md#sec-structured), appendices C–D | Both propositions, theorem, corollary, remark, and complete proofs/derivations |
| Numbered mathematical displays | Equations (1)–(88) | Original numbering, including unlabelled displays and individually numbered aligned rows |
| Peak estimates | [Table 1](../PAPER.md#tab-peak-estimates) | All 11 rows, n = 10–30, checked against the released CSV |
| Figures | Figures 1–7 | All panels, legends, original captions, and sampling qualifications |
| Supporting analysis | [Appendices A–D](../PAPER.md#app-renyi) | Rényi trajectories, balanced-cut distributions, combinatorial derivations, residue proof and sampler |
| Declarations | Data availability, acknowledgements, conflict of interest | Preserved author statements, with a separate data-availability qualification |
| Bibliography | [References 1–33](../PAPER.md#references) | All cited entries, with supplied identifiers and links |

## Representation

All prose and mathematics are ordinary Markdown with Unicode or explicit plain-text notation. Fenced `text` blocks keep multiline formulas readable without a mathematics renderer. Parentheses and braces group indices and exponents; `binom(a,b)` denotes a binomial coefficient, `⊗` a tensor product, and `⊕` a direct sum. Equation numbers and reference numbers have stable local anchors. The original section and appendix order is retained.

The supplied figure files were rendered once with Poppler 26.05.0 at 180 dots per inch into [figures](figures). They are manuscript assets, not new simulations. Their separate checksums preserve the exact imported presentation; Figure 2 comes from the revised manuscript and differs from the older released rendering. Existing CSV and PNG evidence under `data/` and `outputs/` has not been replaced. The existing Python figure scripts remain available to redraw plots from released data, without a document build dependency.

## Scientific and editorial qualifications

The conversion preserves the source's scientific content and adds explicitly marked repository audit notes. A faithful transcription alone cannot establish reproducibility, novelty, or publication readiness.

- The working manuscript title is not substituted into the metadata for the public 2025 preprint, *Arithmetic sequences as quantum states*, arXiv:2501.06292v1. [CITATION.cff](../CITATION.cff) continues to identify the verified public version.
- Historical peak values, fitted lines, and captions are retained. Statements about a numerical peak describe those finite estimates. Complete original global-search records, especially for n = 22–30, remain unavailable; current neighborhood checks do not replace them. The source's protocol description is qualified accordingly.
- The purity-minimizing scale and the dense-side purity crossover do not locate the von Neumann entropy maximizer. The dense separated-mode ansatz remains heuristic; numerical agreement is not a limiting-spectrum theorem. Ensemble means, individual supports, and distributions over cuts remain distinct.
- The supplied acknowledgements and conflict-of-interest statement are retained as author declarations. The data-availability statement is accompanied by the repository's actual evidence limitations.
- All 33 source references are included. Their inclusion is separate from the narrower [verified literature check](../docs/REFERENCES.md), and does not certify novelty or independently verify every bibliography entry.

The [research notes](../docs/RESEARCH.md) remain a concise companion and retain additional audit derivations, including the exact finite-size purity-minimization equation. The [sanity check](../docs/SANITY_CHECK.md) and [provenance record](../PROVENANCE.md) remain part of the paper's supporting repository documentation.

## Completeness checks

The conversion was reviewed against the source by section, including inline quantities, all numbered and unnumbered displays, figure captions, table cells, declarations, and bibliography entries. Particular attention was given to the appendix remainder terms and to the distinction between standard deviations and standard errors.

An independent comparison found no substantive omissions or transcription errors. All 58 citation occurrences retain their source order, and all 76 original cross-reference labels have mapped destinations. Rendering the Markdown produced the expected 90 formula blocks, seven images, one table, and 15 subsections. The existing 39 tests and scientific validator also passed after integration; all 51 frozen evidence/license checksums remained unchanged.

Run `python scripts/check_repository.py` to check the manuscript inventory and all local fragment links, the 88 equation markers, 33 reference entries, seven figure assets and their hashes, and the table against `data/table_i_peaks.csv`, alongside the existing Markdown policy and frozen-evidence checks. This structural check prevents missing objects and broken references; it does not replace mathematical review or prove semantic equivalence of edited prose.
