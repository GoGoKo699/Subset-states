# Repository working rules

- Keep research prose and audit reports in Markdown. Use plain-text or Unicode mathematics; do not introduce TeX sources, BibTeX, PDF deliverables, or a TeX build dependency.
- Preserve the existing license and author attribution.
- Maintain `PAPER.md` as the complete manuscript, including its proofs, appendices, figure captions, declarations, and bibliography. Keep the manuscript inventory and source record in sync with intentional structural changes; use the research notes as a companion, not a replacement for paper content.
- Keep proved finite-ensemble identities, asymptotic consequences, heuristic approximations, finite numerical observations, and historical estimates explicitly distinguished.
- Never infer a von Neumann entropy maximum from the purity minimum or its dense-side crossover. Never label local peak checks as a global search.
- Treat `data/` CSV files and `outputs/` CSV/PNG files as released evidence. Default computations write under `generated/`. Evidence changes require an explicit reason and an updated provenance record; never silently refresh the checksum manifest.
- Keep public preprint metadata tied to a verified public version, rather than an unpublished working title.
- Before delivery run the unit tests, `scripts/final_scientific_validation.py`, `scripts/check_repository.py`, and a CSV figure redraw if plotting changed. For computational changes use a small smoke run; do not launch high-n searches without a concrete research need.
- A validator PASS is scoped to its listed checks. It is not a novelty or submission-readiness verdict.
