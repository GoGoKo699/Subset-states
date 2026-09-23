# References and research positioning

Primary records checked on 22 September 2026. This is a targeted literature check, not a completed novelty assessment. The recovered manuscript's [full 33-entry bibliography](../PAPER.md#references) is maintained separately; its inclusion does not extend this verification to every cited work.

| Source | Role in this project |
|---|---|
| Ruge Lin, Germán Sierra, José I. Latorre, [Arithmetic sequences as quantum states](https://arxiv.org/abs/2501.06292v1), arXiv:2501.06292v1 (2025) | Public predecessor. Its [Table I](https://arxiv.org/html/2501.06292v1) contains the retained $`n = 10`$–$`30`$ peak values. The currently listed title is not the later working title used in old repository metadata. |
| Fernando Granha Jeronimo, Nir Magrafta, Pei Wu, [Pseudorandom and Pseudoentangled States from Subset States](https://arxiv.org/abs/2312.15285), arXiv:2312.15285 (2023; revised 2024) | Prior analysis of the same equal-positive-amplitude, fixed-size random subset ensemble. The multiple-copy calculation includes falling-factorial inclusion probabilities; see Equation 3.2 in the [full text](https://arxiv.org/html/2312.15285v2). |
| Tudor Giurgica-Tiron, Adam Bouland, [Pseudorandomness from Subset States](https://arxiv.org/abs/2312.09206), arXiv:2312.09206 (2023) | Independent prior subset-state moment and indistinguishability analysis. Equal-amplitude subset states and their phase-free pseudorandomness are not introduced by this repository. |
| Giuseppe De Tomasi, Ivan M. Khaymovich, [Multifractality meets entanglement: relation for non-ergodic extended states](https://arxiv.org/abs/2001.03173), [Physical Review Letters 124, 200602 (2020)](https://doi.org/10.1103/PhysRevLett.124.200602) | Sparse random states can reach Page-scale entanglement despite occupying a vanishing fraction of the basis. Their Gaussian nonzero amplitudes distinguish that ensemble from the equal-positive-amplitude model here. |
| Don N. Page, [Average Entropy of a Subsystem](https://arxiv.org/abs/gr-qc/9305007), [Physical Review Letters 71, 1291 (1993)](https://doi.org/10.1103/PhysRevLett.71.1291) | Haar-random reference mean and its asymptotic behavior; a benchmark for a different ensemble. |

## What the comparison establishes

The broad fact that sparse support can coexist with high entanglement is established prior work. Fixed-cardinality subset moments are also established objects of study. The counting arguments in the research notes should therefore be read as direct derivations and checks for this repository, without a claim that moment counting itself is new.

The equal-positive-amplitude assumption matters: full support is a product state, and the coherent dense-side contribution cannot be transferred unchanged from random-phase or Gaussian-amplitude ensembles.

It remains necessary to compare the exact balanced-cut purity formula, its finite-size minimization, and the residue ceilings directly against prior equations before asserting a new publishable theorem. This cleanup makes no claim of quantum speedup, computational pseudorandomness, a unique primality fingerprint, or a settled asymptotic entropy-peak exponent.

## Citation metadata

`CITATION.cff` names the existing public preprint with its actual title and year. Its software title describes this repository separately. Cite the commit used for revised code or data. A working title in a README is not evidence of a new arXiv version or a journal publication.
