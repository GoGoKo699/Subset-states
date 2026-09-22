# Research notes: equal-positive-amplitude subset states

These notes state the mathematical model, derive the exact results represented in the code, and separate them from approximations and numerical evidence. They are the maintained research text for this repository. Entropies and logarithms use bits unless `ln` is written.

## 1. Model and endpoints

Let n be a positive even integer, N = 2ⁿ, and d = 2^(n/2). Choose a support A uniformly among all M-element subsets of the N computational-basis labels, with 1 ≤ M ≤ N. Define

```text
|A⟩ = (1/√M) ∑ over x in A of |x⟩.
```

All nonzero amplitudes are real and positive. Random phases or Gaussian nonzero amplitudes give different ensembles. A fixed balanced cut identifies each label with a pair (a,b) in a d-by-d grid. Write X for its binary incidence matrix and C = X/√M. The right reduced state is ρ = XᵀX/M. For a general complex coefficient matrix with entries C[a,b], the corresponding formula is ρ = CᵀC*, where * denotes entrywise conjugation. Its eigenvalues are the squared singular values of C.

The von Neumann entropy is S = −∑ λ log₂ λ, with zero terms omitted. Rényi-2 entropy is −log₂ Tr(ρ²); min-entropy is −log₂ λmax. Always 0 ≤ S ≤ log₂ d = n/2.

A singleton support is a computational-basis product state. Full support is |+⟩ raised to the n-fold tensor product and is also unentangled. A permutation-matrix support of size d is maximally entangled. Thus support size alone does not specify an individual state's entropy. The curves here describe ensembles indexed by M, not a physical time evolution or a nested growth process.

Uniform support sampling is invariant under bit permutations, so every fixed balanced cut has the same ensemble distribution. This does not imply that all cuts of an individual support have the same entropy.

## 2. Exact mean reduced state

For k distinct grid cells, the probability that all are occupied is

```text
pₖ = (M)ₖ / (N)ₖ,
(z)ₖ = z(z−1)…(z−k+1),
pₖ = 0 when k > M.
```

Each diagonal entry of ρ sums d occupied-cell indicators divided by M. Each off-diagonal entry sums d products of two distinct indicators divided by M. Therefore

```text
E[ρ] = (1/d − β) I + β J,
β = (M−1) / [d(N−1)],
```

where J is the all-ones matrix. Its uniform-vector eigenvalue is 1/d + (d−1)β. Each of its other d−1 eigenvalues equals 1/d − β.

This is the spectrum of the mean matrix. It is not the mean of the ordered sample eigenvalues. For example, M = 1 gives E[ρ] = I/d even though every sample has zero entropy. Concavity gives E[S(ρ)] ≤ S(E[ρ]), not equality.

Implementation: `fixed_cardinality_mean_spectrum` in [core.py](../subset_states/core.py).

## 3. Exact average purity

Expand Tr(ρ²) as the sum of X[a,b] X[a,c] X[a′,c] X[a′,b] over all four indices, divided by M². There are three cases:

| Index pattern | Number of terms | Distinct occupied cells |
|---|---:|---:|
| a = a′ and b = c | N | 1 |
| Exactly one equality | 2N(d−1) | 2 |
| Neither equality | N(d−1)² | 4 |

Substituting the inclusion probabilities gives

```text
P̄ = E[Tr(ρ²)]
   = 1/M
     + 2(d−1)(M−1) / [M(N−1)]
     + (d−1)²(M−1)(M−2)(M−3) / [M(N−1)(N−2)(N−3)].
```

The formula covers n ≥ 2 even and every valid M, including both product-state endpoints. The last contribution counts rectangles in the grid; no independence approximation is used.

For each sample S ≥ −log₂ Tr(ρ²). Since −log₂ is convex, Jensen's inequality yields the ensemble bound

```text
E[S] ≥ E[−log₂ Tr(ρ²)] ≥ −log₂ P̄.
```

This is a lower bound on mean entropy. A concentration claim or an individual-state guarantee requires an additional argument.

Implementation: `fixed_cardinality_average_purity` and `average_entropy_lower_bound_from_purity`.

## 4. Three distinct support scales

For fixed c > 0 and a fixed exponent 0 < γ < 1, take integer M asymptotic to cN^γ. Expanding the exact purity gives

```text
P̄ = 2N^(−1/2) + c^(−1)N^(−γ) + c²N^(2γ−2)
     + terms smaller than their respective leading contributions.
```

For every 1/2 < γ < 3/4, the last two displayed contributions are smaller than N^(−1/2), so

```text
P̄ = [2 + o(1)] N^(−1/2),
E[S] ≥ n/2 − 1 − o(1).
```

The deficit stays bounded by approximately one bit. This is not a statement that the absolute deficit vanishes, nor that the ensemble attains the Page mean.

The smallest correction to this purity background occurs when the sparse and rectangle terms balance, giving γ = 2/3 and c = 2^(−1/3). The O(1) correction to the resulting minimizer can be checked directly. Set

```text
a = 2/(d+1),
b = (d−1)/[(d+1)(N−2)(N−3)].

P̄(M) = bM² − 6bM + a + 11b + (1−a−6b)/M.
```

For N ≥ 16, the last numerator is positive. The second derivative is positive, and the derivative vanishes exactly when

```text
M²(M−3) = N(N−5)/2.
```

The continuous minimizer is 2^(−1/3)N^(2/3) + 1 + O(N^(−1/3)); comparing the two neighboring integers gives an integer minimizer of the form 2^(−1/3)N^(2/3) + O(1). At N = 4, direct evaluation gives the minimum at M = 3.

At γ = 3/4 instead,

```text
P̄ = (2+c²)N^(−1/2) + o(N^(−1/2)).
```

This is the crossover where the rectangle term becomes comparable with the balanced-cut background. It is not a theorem about the von Neumann entropy maximum. The average-purity minimizer, the maximizer of average Rényi-2 entropy, and the maximizer of average von Neumann entropy are different optimization problems.

## 5. Exact diagonal and residue bounds

The occupancy W of any column is hypergeometric with population size N, d marked cells, and M draws. The diagonal probability of that column is W/M. Hence the mean diagonal entropy is exactly

```text
E[H(diag ρ)] = −d E[(W/M) log₂(W/M)].
```

Deleting off-diagonal coherences cannot reduce entropy, so E[S] ≤ E[H(diag ρ)]. The function named `hypergeometric_occupancy_approximation` computes this exact diagonal expectation; treating it as an approximation to von Neumann entropy is appropriate only when coherences are sufficiently small.

For the natural cut, the right subsystem consists of the low n/2 bits. Let pᵣ be the fraction of occupied labels congruent to r modulo 2ᵗ, with 0 ≤ t ≤ n/2. Dephase the t low bits. The resulting blocks have total weights pᵣ and dimension at most d/2ᵗ. Their entropy is at most H(p) + log₂(d/2ᵗ), while dephasing cannot decrease entropy. Therefore every such support satisfies

```text
S ≤ n/2 − t + H(p).
```

This is a cut-specific computational-basis statement. Applying a global QFT changes the state; the input residue populations do not automatically bound its output entanglement in the same way.

## 6. Approximations and reference ensembles

The dense-bulk ansatz uses the uniform-mode eigenvalue λ of E[ρ] as an estimate for a separated sample eigenvalue, then assigns the remaining weight a Page-like bulk:

```text
T = h₂(λ) + (1−λ)[log₂ d − 1/(2 ln 2)],
h₂(λ) = −λ log₂ λ − (1−λ) log₂(1−λ).
```

The mean-state formula is exact. The identification of its eigenvalue with a typical separated eigenvalue, the bulk shape, and this entropy expression are approximations. Algebraically equivalent formulas and one spectral plot do not establish their accuracy uniformly in M or n.

The Page mean for a complex Haar-random pure state across a d-by-d cut is

```text
S_Page = [H_(d²) − H_d − (d−1)/(2d)] / ln 2,
```

where H_j is the j-th harmonic number. It is a reference average for a different ensemble, not an entropy ceiling. The ceiling is log₂ d. Related work and ensemble distinctions are in [references](REFERENCES.md).

## 7. Numerical evidence

The retained peak table fits log₂ M̂ = 0.703541n − 0.357734 over eleven even sizes n = 10–30. This is a descriptive regression of historical point estimates. The original global-search records are incomplete, especially for n = 22–30. Local checks near the supplied centers cannot recover global optimality or supply the original search uncertainty. See [provenance](../PROVENANCE.md).

| Figure | What is shown | Scope |
|---|---|---|
| [1](../outputs/fig1/fig1_concentration.png) | Support-size curves and cut sampling | Finite ensemble samples |
| [2](../outputs/fig2/fig2_peak_scaling.png) | Historical peak table and regressions | No asymptotic exponent established |
| [3](../outputs/fig3/fig3_spectral_bulk.png) | Spectrum and a bulk comparison | One sampled n = 24 state |
| [4](../outputs/fig4/fig4_approximation.png) | Sparse and dense approximations | Numerical comparison, not a uniform error theorem |
| [5](../outputs/fig5/fig5_qft_residue_controls.png) | Almost-prime supports and matched ensembles | Paired computational/Fourier comparisons at n = 14 |
| [6](../outputs/fig6/fig6_renyi.png) | Several Rényi entropy curves | Different entropy orders need not peak together |
| [7](../outputs/fig7/fig7_partitions.png) | Entropy over sampled balanced cuts | Two selected n = 20 states |

The QFT uses the positive-exponent convention, implemented by NumPy's orthonormal inverse FFT. Changing its sign conjugates the output for these real inputs and leaves bipartite entropies unchanged. This convention check says nothing about implementing a physical QFT at large n.

## 8. Open scientific questions

1. Determine the asymptotic maximizer of the mean von Neumann entropy, with error control; the purity calculation alone does not do this.
2. Establish a quantitative dense-bulk approximation for this fixed-cardinality, positive-amplitude ensemble, with a stated regime and finite-size errors.
3. Compare the exact formulas and their consequences directly with existing fixed-size subset-state moment calculations before making a novelty claim.
4. Recover or replace the missing global-search evidence if the high-n numerical peak law is retained as a central contribution.

The present cleanup strengthens the repository's reliability and boundaries. It does not settle these research questions.
