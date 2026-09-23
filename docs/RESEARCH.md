# Research notes: equal-positive-amplitude subset states

These notes are a concise companion to the [complete manuscript](../PAPER.md). They state the mathematical model, derive the exact results represented in the code, and separate them from approximations and numerical evidence. They also retain additional audit derivations, including finite-size purity minimization. Entropies and logarithms use bits unless $`\ln`$ is written.

## 1. Model and endpoints

Let $`n`$ be a positive even integer, $`N = 2^n`$, and $`d = 2^{n/2}`$. Choose a support $`A`$ uniformly among all $`M`$-element subsets of the $`N`$ computational-basis labels, with $`1 \le M \le N`$. Define

```math
\lvert A\rangle = \frac{1}{\sqrt M}\sum_{x\in A}\lvert x\rangle.
```

All nonzero amplitudes are real and positive. Random phases or Gaussian nonzero amplitudes give different ensembles. A fixed balanced cut identifies each label with a pair $`(a,b)`$ in a $`d`$-by-$`d`$ grid. Write $`X`$ for its binary incidence matrix and $`C = X/\sqrt M`$. The right reduced state is $`\rho = X^{\mathsf T}X/M`$. For a general complex coefficient matrix with entries $`C_{a,b}`$, the corresponding formula is $`\rho = C^{\mathsf T}C^*`$, where $`*`$ denotes entrywise conjugation. Its eigenvalues are the squared singular values of $`C`$.

The von Neumann entropy is $`S = -\sum_i\lambda_i\log_2\lambda_i`$, with zero terms omitted. Rényi-2 entropy is $`-\log_2\mathrm{Tr}(\rho^2)`$; min-entropy is $`-\log_2\lambda_{\max}`$. Always $`0\le S\le\log_2 d = n/2`$.

A singleton support is a computational-basis product state. Full support is $`\lvert+\rangle`$ raised to the $`n`$-fold tensor product and is also unentangled. A permutation-matrix support of size $`d`$ is maximally entangled. Thus support size alone does not specify an individual state's entropy. The curves here describe ensembles indexed by $`M`$, not a physical time evolution or a nested growth process.

Uniform support sampling is invariant under bit permutations, so every fixed balanced cut has the same ensemble distribution. This does not imply that all cuts of an individual support have the same entropy.

## 2. Exact mean reduced state

For $`k`$ distinct grid cells, the probability that all are occupied is

```math
\begin{aligned}
p_k &= \frac{(M)_k}{(N)_k},\\
(z)_k &= z(z-1)\cdots(z-k+1),\\
p_k &= 0\qquad\text{when }k>M.
\end{aligned}
```

Each diagonal entry of $`\rho`$ sums $`d`$ occupied-cell indicators divided by $`M`$. Each off-diagonal entry sums $`d`$ products of two distinct indicators divided by $`M`$. Therefore

```math
\begin{aligned}
\mathbb E[\rho] &= \left(\frac1d-\beta\right)I+\beta J,\\
\beta &= \frac{M-1}{d(N-1)},
\end{aligned}
```

where $`J`$ is the all-ones matrix. Its uniform-vector eigenvalue is $`1/d+(d-1)\beta`$. Each of its other $`d-1`$ eigenvalues equals $`1/d-\beta`$.

This is the spectrum of the mean matrix. It is not the mean of the ordered sample eigenvalues. For example, $`M=1`$ gives $`\mathbb E[\rho]=I/d`$ even though every sample has zero entropy. Concavity gives $`\mathbb E[S(\rho)]\le S(\mathbb E[\rho])`$, not equality.

Implementation: `fixed_cardinality_mean_spectrum` in [core.py](../subset_states/core.py).

## 3. Exact average purity

Expand $`\mathrm{Tr}(\rho^2)`$ as the sum of $`X_{a,b}X_{a,c}X_{a',c}X_{a',b}`$ over all four indices, divided by $`M^2`$. There are three cases:

| Index pattern | Number of terms | Distinct occupied cells |
|---|---:|---:|
| $`a=a'`$ and $`b=c`$ | $`N`$ | $`1`$ |
| Exactly one equality | $`2N(d-1)`$ | $`2`$ |
| Neither equality | $`N(d-1)^2`$ | $`4`$ |

Substituting the inclusion probabilities gives

```math
\begin{aligned}
\overline P &= \mathbb E[\mathrm{Tr}(\rho^2)]\\
&= \frac1M+\frac{2(d-1)(M-1)}{M(N-1)}\\
&\quad+\frac{(d-1)^2(M-1)(M-2)(M-3)}{M(N-1)(N-2)(N-3)}.
\end{aligned}
```

The formula covers $`n\ge2`$ even and every valid $`M`$, including both product-state endpoints. The last contribution counts rectangles in the grid; no independence approximation is used.

For each sample $`S\ge-\log_2\mathrm{Tr}(\rho^2)`$. Since $`-\log_2`$ is convex, Jensen's inequality yields the ensemble bound

```math
\mathbb E[S]
\ge\mathbb E[-\log_2\mathrm{Tr}(\rho^2)]
\ge-\log_2\overline P.
```

This is a lower bound on mean entropy. A concentration claim or an individual-state guarantee requires an additional argument.

Implementation: `fixed_cardinality_average_purity` and `average_entropy_lower_bound_from_purity`.

## 4. Three distinct support scales

For fixed $`c>0`$ and a fixed exponent $`0<\gamma<1`$, take integer $`M`$ asymptotic to $`cN^\gamma`$. Expanding the exact purity gives

```math
\overline P = 2N^{-1/2}+c^{-1}N^{-\gamma}+c^2N^{2\gamma-2}
\; +\; \cdots,
```

where the omitted terms are smaller than their respective leading contributions.

For every $`1/2<\gamma<3/4`$, the last two displayed contributions are smaller than $`N^{-1/2}`$, so

```math
\begin{aligned}
\overline P &= [2+o(1)]N^{-1/2},\\
\mathbb E[S] &\ge\frac n2-1-o(1).
\end{aligned}
```

The deficit stays bounded by approximately one bit. This is not a statement that the absolute deficit vanishes, nor that the ensemble attains the Page mean.

The smallest correction to this purity background occurs when the sparse and rectangle terms balance, giving $`\gamma=2/3`$ and $`c=2^{-1/3}`$. The $`O(1)`$ correction to the resulting minimizer can be checked directly. Set

```math
\begin{aligned}
a &= \frac2{d+1},\\
b &= \frac{d-1}{(d+1)(N-2)(N-3)},\\[4pt]
\overline P(M) &= bM^2-6bM+a+11b+\frac{1-a-6b}{M}.
\end{aligned}
```

For $`N\ge16`$, the last numerator is positive. The second derivative is positive, and the derivative vanishes exactly when

```math
M^2(M-3)=\frac{N(N-5)}2.
```

The continuous minimizer is $`2^{-1/3}N^{2/3}+1+O(N^{-1/3})`$; comparing the two neighboring integers gives an integer minimizer of the form $`2^{-1/3}N^{2/3}+O(1)`$. At $`N=4`$, direct evaluation gives the minimum at $`M=3`$.

At $`\gamma=3/4`$ instead,

```math
\overline P=(2+c^2)N^{-1/2}+o(N^{-1/2}).
```

This is the crossover where the rectangle term becomes comparable with the balanced-cut background. It is not a theorem about the von Neumann entropy maximum. The average-purity minimizer, the maximizer of average Rényi-2 entropy, and the maximizer of average von Neumann entropy are different optimization problems.

## 5. Exact diagonal and residue bounds

The occupancy $`W`$ of any column is hypergeometric with population size $`N`$, $`d`$ marked cells, and $`M`$ draws. The diagonal probability of that column is $`W/M`$. Hence the mean diagonal entropy is exactly

```math
\mathbb E[H(\mathrm{diag}\rho)]
=-d\,\mathbb E\!\left[\frac WM\log_2\!\left(\frac WM\right)\right].
```

Deleting off-diagonal coherences cannot reduce entropy, so $`\mathbb E[S]\le\mathbb E[H(\mathrm{diag}\rho)]`$. The function named `hypergeometric_occupancy_approximation` computes this exact diagonal expectation; treating it as an approximation to von Neumann entropy is appropriate only when coherences are sufficiently small.

For the natural cut, the right subsystem consists of the low $`n/2`$ bits. Let $`p_r`$ be the fraction of occupied labels congruent to $`r`$ modulo $`2^t`$, with $`0\le t\le n/2`$. Dephase the $`t`$ low bits. The resulting blocks have total weights $`p_r`$ and dimension at most $`d/2^t`$. Their entropy is at most $`H(p)+\log_2(d/2^t)`$, while dephasing cannot decrease entropy. Therefore every such support satisfies

```math
S\le\frac n2-t+H(p).
```

This is a cut-specific computational-basis statement. Applying a global QFT changes the state; the input residue populations do not automatically bound its output entanglement in the same way.

## 6. Approximations and reference ensembles

The dense-bulk ansatz uses the uniform-mode eigenvalue $`\lambda`$ of $`\mathbb E[\rho]`$ as an estimate for a separated sample eigenvalue, then assigns the remaining weight a Page-like bulk:

```math
\begin{aligned}
T &= h_2(\lambda)+(1-\lambda)\left[\log_2 d-\frac1{2\ln2}\right],\\
h_2(\lambda) &= -\lambda\log_2\lambda-(1-\lambda)\log_2(1-\lambda).
\end{aligned}
```

The mean-state formula is exact. The identification of its eigenvalue with a typical separated eigenvalue, the bulk shape, and this entropy expression are approximations. Algebraically equivalent formulas and one spectral plot do not establish their accuracy uniformly in $`M`$ or $`n`$.

The Page mean for a complex Haar-random pure state across a $`d`$-by-$`d`$ cut is

```math
S_{\mathrm{Page}}=\frac{H_{d^2}-H_d-\dfrac{d-1}{2d}}{\ln2},
```

where $`H_j`$ is the $`j`$-th harmonic number. It is a reference average for a different ensemble, not an entropy ceiling. The ceiling is $`\log_2 d`$. Related work and ensemble distinctions are in [references](REFERENCES.md).

## 7. Numerical evidence

The retained peak table fits $`\log_2\widehat M=0.703541n-0.357734`$ over eleven even sizes $`n=10`$–$`30`$. This is a descriptive regression of historical point estimates. The original global-search records are incomplete, especially for $`n=22`$–$`30`$. Local checks near the supplied centers cannot recover global optimality or supply the original search uncertainty. See [provenance](../PROVENANCE.md).

| Figure | What is shown | Scope |
|---|---|---|
| [1](../outputs/fig1/fig1_concentration.png) | Support-size curves and cut sampling | Finite ensemble samples |
| [2](../outputs/fig2/fig2_peak_scaling.png) | Historical peak table and regressions | No asymptotic exponent established |
| [3](../outputs/fig3/fig3_spectral_bulk.png) | Spectrum and a bulk comparison | One sampled $`n=24`$ state |
| [4](../outputs/fig4/fig4_approximation.png) | Sparse and dense approximations | Numerical comparison, not a uniform error theorem |
| [5](../outputs/fig5/fig5_qft_residue_controls.png) | Almost-prime supports and matched ensembles | Paired computational/Fourier comparisons at $`n=14`$ |
| [6](../outputs/fig6/fig6_renyi.png) | Several Rényi entropy curves | Different entropy orders need not peak together |
| [7](../outputs/fig7/fig7_partitions.png) | Entropy over sampled balanced cuts | Two selected $`n=20`$ states |

The QFT uses the positive-exponent convention, implemented by NumPy's orthonormal inverse FFT. Changing its sign conjugates the output for these real inputs and leaves bipartite entropies unchanged. This convention check says nothing about implementing a physical QFT at large $`n`$.

## 8. Open scientific questions

1. Determine the asymptotic maximizer of the mean von Neumann entropy, with error control; the purity calculation alone does not do this.
2. Establish a quantitative dense-bulk approximation for this fixed-cardinality, positive-amplitude ensemble, with a stated regime and finite-size errors.
3. Compare the exact formulas and their consequences directly with existing fixed-size subset-state moment calculations before making a novelty claim.
4. Recover or replace the missing global-search evidence if the high-$`n`$ numerical peak law is retained as a central contribution.

The present cleanup strengthens the repository's reliability and boundaries. It does not settle these research questions.
