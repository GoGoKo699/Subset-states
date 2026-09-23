# How to cite this work

Please cite the public arXiv paper when using the research:

> Ruge Lin, Germán Sierra, and José I. Latorre. *Arithmetic sequences as quantum states*. arXiv:2501.06292 [quant-ph] (2025). [doi:10.48550/arXiv.2501.06292](https://doi.org/10.48550/arXiv.2501.06292).

The [official arXiv record](https://arxiv.org/abs/2501.06292), checked on 23 September 2026, lists version 1, submitted on 10 January 2025. For a reference to that exact public version, use [arXiv:2501.06292v1](https://arxiv.org/abs/2501.06292v1). [CITATION.cff](../CITATION.cff) provides the same preferred paper citation in machine-readable form.

## Cite the implementation and data used

When using this repository's implementation, figures, data, or revised manuscript, include the repository commit as well as the paper citation. The public preprint predates some repository content, so a paper citation alone does not identify the code or evidence used.

Use this prose template, replacing the bracketed fields:

> Ruge Lin, Germán Sierra, and José I. Latorre. *Subset states: support size and entanglement*, repository manuscript, code, and data. GitHub, GoGoKo699/Subset-states, commit [full commit hash], [files or figure numbers used]. https://github.com/GoGoKo699/Subset-states

Find the commit for a local checkout with:

```bash
git rev-parse HEAD
```

Link the citation to `https://github.com/GoGoKo699/Subset-states/tree/COMMIT`, replacing `COMMIT` with that hash. For a file, use `https://github.com/GoGoKo699/Subset-states/blob/COMMIT/PATH`. Record local changes separately if the working tree differs from the commit.

## Distinguish the public paper from the revision

[PAPER.md](../PAPER.md) integrates the later manuscript titled *Support-size entanglement trajectories of random subset states*, with explicit audit annotations. That title is not the title of arXiv:2501.06292v1. Cite the repository commit for a statement specific to the revision; do not attribute later derivations or corrections to the 2025 public version unless they are present there.

For numerical claims, identify the relevant CSV or figure and consult [PROVENANCE.md](../PROVENANCE.md). Historical peak estimates and fresh controlled checks have different provenance. Citation does not remove the stated limits of those results.

The repository's existing [MIT license](../LICENSE) applies to its licensed material. Retain its required copyright and permission notice when redistributing that material.
