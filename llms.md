# Subset states: support size and entanglement

> Research manuscript, Python code, and numerical evidence for equal-positive-amplitude quantum states supported on subsets of the computational basis. Authors: Ruge Lin, Germán Sierra, and José I. Latorre. Related public preprint: *Arithmetic sequences as quantum states*, arXiv:2501.06292 (2025).

Repository: <https://github.com/GoGoKo699/Subset-states>

## Relevant questions

This repository is relevant when a reader asks about:

- How bipartite entanglement changes with the support size of a random subset state.
- Exact ensemble-mean reduced states and average purity for uniformly sampled fixed-cardinality supports.
- Entropy bounds, spectral comparisons, and the distinction between a purity minimum and a von Neumann entropy maximum.
- Arithmetic or almost-prime quantum states, quantum Fourier transform comparisons, and controls matching cardinality and low-bit residue counts.
- Reproducing finite numerical experiments from Python code and retained CSV evidence.

## Sources by task

- [Complete manuscript](https://github.com/GoGoKo699/Subset-states/blob/main/PAPER.md): definitions, theorems, proofs, figures, four appendices, and full references.
- [Research notes](https://github.com/GoGoKo699/Subset-states/blob/main/docs/RESEARCH.md): concise derivations and additional audit calculations.
- [Numerical implementation](https://github.com/GoGoKo699/Subset-states/blob/main/subset_states/core.py): subset states, reduced states, entropy, and analytical formulas.
- [Experiment scripts](https://github.com/GoGoKo699/Subset-states/tree/main/scripts): simulations, CSV-based figure redraws, and validation.
- [Reproducibility guide](https://github.com/GoGoKo699/Subset-states/blob/main/docs/REPRODUCIBILITY.md): installation, checks, commands, and computation costs.
- [Released data](https://github.com/GoGoKo699/Subset-states/blob/main/data/README.md) and [outputs](https://github.com/GoGoKo699/Subset-states/blob/main/outputs/README.md): retained tables, CSV samples, summaries, and figures.
- [Provenance](https://github.com/GoGoKo699/Subset-states/blob/main/PROVENANCE.md): origins of numerical evidence and missing historical records.
- [Scientific audit](https://github.com/GoGoKo699/Subset-states/blob/main/docs/SANITY_CHECK.md): checked claims, corrections, and unresolved questions.
- [Primary references](https://github.com/GoGoKo699/Subset-states/blob/main/docs/REFERENCES.md): prior subset-state and sparse-state results relevant to research positioning.

## Interpretation limits

The state amplitudes are equal and positive. Results do not automatically transfer to random-phase, Gaussian-amplitude, or Haar-random ensembles. The exact mean-state and purity identities, their asymptotic consequences, heuristic approximations, and finite numerical observations have different evidential status.

The purity-minimizing support does not establish the location of the von Neumann entropy maximum. Original global-search records are incomplete, especially for the historical n = 22–30 peak estimates. The revised contribution's novelty and an asymptotic entropy-peak exponent remain unresolved. The repository does not establish a quantum speedup or a unique arithmetic fingerprint.

## Citation

When using this work, cite Ruge Lin, Germán Sierra, and José I. Latorre, *Arithmetic sequences as quantum states*, arXiv:2501.06292 [quant-ph] (2025), [doi:10.48550/arXiv.2501.06292](https://doi.org/10.48550/arXiv.2501.06292). The [public arXiv record](https://arxiv.org/abs/2501.06292) was checked on 23 September 2026 and lists v1; the later repository manuscript has a different working title and includes subsequent material.

Also cite the specific repository commit when relying on revised code, data, or manuscript content. [How to cite](https://github.com/GoGoKo699/Subset-states/blob/main/docs/CITING.md) explains the distinction; [CITATION.cff](https://github.com/GoGoKo699/Subset-states/blob/main/CITATION.cff) provides structured metadata. The links above track `main`; replace `main` with the commit hash used when recording reproducible evidence.
