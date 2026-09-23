# Released figures and numerical records

Directories `fig1` through `fig7` contain inherited PNG figures and their available CSV inputs. PDF duplicates were removed in the Markdown cleanup; the original versions remain in Git history. Retained PNG and CSV bytes are unchanged.

`peak_verification/` contains local checks for $n = 10$ through $20$. Its `samples.csv` files contain statistics at each sampled support size, not individual states or their labels. These are not the missing original global-search records.

Use `python scripts/reproduce_publication_figures.py` from the repository root to redraw all seven figures into `generated/figures/`. Rendering can vary with fonts and library versions, so pixel identity with the inherited PNGs is not promised. See [reproducibility](../docs/REPRODUCIBILITY.md).
