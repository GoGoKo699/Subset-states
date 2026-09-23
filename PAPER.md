# Support-size entanglement trajectories of random subset states

**Ruge Lin**  
Thrust of Artificial Intelligence, Information Hub,  
The Hong Kong University of Science and Technology (Guangzhou), China

**Germán Sierra**  
Instituto de Física Teórica UAM-CSIC, Universidad Autónoma de Madrid, Spain

**José I. Latorre**  
Centre for Quantum Technologies, National University of Singapore, Singapore

> **Repository edition — 23 September 2026.** This is the complete Markdown conversion of the revised manuscript dated 21 August 2026, with its original section, equation, figure, table, and reference numbering. The [integration record](paper/INTEGRATION.md) documents source identity and coverage. Explicit repository audit notes qualify historical numerical claims and data availability. In particular, the n = 10–30 peak estimates are retained historical values; the original global searches have not been independently reconstructed. Novelty and the asymptotic entropy-maximizing support remain open, as detailed in the [sanity check](docs/SANITY_CHECK.md).

<a id="abstract"></a>

## Abstract

Equal-amplitude subset states are constrained by small Schmidt rank when their support is sparse, but approach a coherent product direction when their support becomes dense. We study this competition for uniformly random supports of fixed cardinality M in an n-qubit computational basis of dimension N = 2<sup>n</sup>. We derive the exact ensemble-mean reduced state and an exact formula for the average purity across a balanced bipartition. At the support scale that optimally balances the two purity corrections,

<a id="abstract-purity-scale"></a>

> M = 2<sup>−1/3</sup> N<sup>2/3</sup> + O(1),

the fixed-cardinality ensemble satisfies

<a id="abstract-entropy-bound"></a>

> S̄<sub>N,M</sub> ≥ n/2 − 1 − o(1),

so a vanishing fraction of the basis can produce average entanglement within asymptotically one bit of the maximum. The released sampled trajectories exhibit a rise–peak–fall pattern: the von Neumann entropy vanishes at single and full support and reaches a broad intermediate maximum. Retained historical estimates of the peak mean entropy approach the balanced Page value as n increases from 10 to 30. Over this range, the estimated peak support grows exponentially while occupying a decreasing fraction of the basis. An exact hypergeometric diagonal-entropy bound and a dense separated-mode ansatz explain the sparse and dense sides of the trajectory. As an application, we use cardinality- and residue-matched random ensembles to separate support-size, low-order congruence, and residual arrangement effects in almost-prime subset states, before and after a quantum Fourier transform.

**Keywords:** subset states; entanglement entropy; random quantum states; fixed-cardinality ensembles; random matrices; almost-prime states

## Contents

- [Abstract](#abstract)
- [1. Introduction](#sec-introduction)
- [2. Subset states and the fixed-cardinality ensemble](#sec-model)
- [3. Exact moments and compressed-support entanglement](#sec-exact-moments)
- [4. Support-size entanglement trajectory](#sec-trajectory)
- [5. Sparse-to-dense spectral mechanism](#sec-spectral-mechanism)
- [6. Structured supports and Fourier comparison](#sec-structured)
- [7. Conclusions](#sec-conclusions)
- [Data availability](#data-availability), [acknowledgements](#acknowledgements), and [conflict of interest](#conflict-of-interest)
- [Appendix A. Rényi-entropy trajectories](#app-renyi)
- [Appendix B. Entropy distributions over balanced cuts](#app-partitions)
- [Appendix C. Combinatorial derivations](#app-combinatorial-derivations)
- [Appendix D. Residue-class bound and constrained ensembles](#app-residue-proof)
- [References](#references)

**Notation.** Formulas use Unicode mathematics with HTML subscripts and superscripts in Markdown. Parentheses and brackets make the scope of inline fractions explicit; binom(a,b) is a binomial coefficient. A bar denotes an ensemble mean, a hat a numerical estimate. Logarithms are base 2 unless written ln. Equation references link directly to their numbered displays.

<a id="sec-introduction"></a>

## 1. Introduction

At support size M = 1, an equal-amplitude subset state is a computational-basis product state. At full support M = 2<sup>n</sup>, it is again a product state, |+⟩<sup>⊗n</sup>. Between these endpoints, however, a support of size 2<sup>n/2</sup> can be arranged as a permutation matrix across a balanced cut and produce exactly maximal entanglement. This elementary contrast poses the question studied here: for a support chosen uniformly at random, how does bipartite entanglement depend on the number of occupied basis labels?

A subset state has the form

<a id="eq-1"></a>

> **(1)** |A⟩ = (1/√|A|) ∑<sub>x ∈ A</sub> |x⟩,

so all randomness is carried by the support A; there are no random phases or amplitude fluctuations. Subset states were introduced as a restricted witness family in quantum complexity theory [[1]](#ref-1) and have recently become central examples in pseudorandomness and pseudoentanglement [[2]](#ref-2), [[3]](#ref-3), [[4]](#ref-4), [[5]](#ref-5). Those works ask whether suitable ensembles can be computationally or information-theoretically difficult to distinguish from Haar-random states. Our question is complementary and structural: before any indistinguishability criterion is imposed, what entanglement is generated by support placement alone?

For unrestricted random pure states, typical balanced subsystems are nearly maximally entangled. Exact mean formulas were obtained through purity and entropy calculations [[6]](#ref-6), [[7]](#ref-7), [[8]](#ref-8), and concentration results explain why high entanglement is generic in large Hilbert spaces [[9]](#ref-9). The associated reduced states are standard induced random matrices, with well-developed spectral descriptions [[10]](#ref-10), [[11]](#ref-11), [[12]](#ref-12). Support-constrained states need not follow the same reference. In particular, sparse random states can reach volume-law or Page-like entanglement while occupying only a vanishing part of Hilbert space, with the answer depending on the support and coefficient statistics [[13]](#ref-13). Equal-positive-amplitude subset states add a further coherent constraint: as the support becomes dense, constructive overlaps drive the state toward the uniform product vector. This produces a descending branch absent from models with independent random occupied amplitudes.

We first compute the ensemble-mean reduced state and the exact average purity by counting occupied rectangles in the bipartite support matrix. For a power-law support M = cN<sup>γ</sup>, the resulting terms define a near-maximal window 1/2 < γ < 3/4: γ = 2/3 is the optimal interior balance, while γ = 3/4 is the dense-side boundary. We then examine the finite-size von Neumann trajectory numerically. The released sampled trajectories show a broad intermediate maximum, and retained historical peak estimates for even n = 10, …, 30 approach the balanced Page benchmark while occupying a shrinking fraction of the basis. The exact first moment, a hypergeometric occupancy bound, and a separated-mode random-matrix ansatz explain the sparse rise and dense decline.

A second purpose of the fixed-cardinality ensemble is to provide a controlled reference for structured supports. Typical entanglement depends on the sector or manifold from which random states are drawn; particle-number and symmetry constraints can substantially modify the appropriate Page-type benchmark [[14]](#ref-14), [[15]](#ref-15), [[16]](#ref-16), [[17]](#ref-17), [[18]](#ref-18). We apply the same principle to almost-prime supports. A cardinality-matched ensemble first removes the dominant effect of support size, and residue-matched ensembles then remove elementary congruence populations visible in the binary encoding. This hierarchy turns the random trajectory into a diagnostic that separates support-size, low-order congruence, and residual arrangement effects. The application connects with quantum states built from primes and almost primes [[19]](#ref-19), [[20]](#ref-20), [[21]](#ref-21), while keeping the random subset-state problem as the main subject of the paper.

<a id="sec-model"></a>

## 2. Subset states and the fixed-cardinality ensemble

For any positive integer K, write

<a id="eq-2"></a>

> **(2)** [K]<sub>0</sub> := {0, …, K − 1}.

Let N = 2<sup>n</sup>, with n even, and let

<a id="eq-3"></a>

> **(3)** A ⊆ [N]<sub>0</sub>, |A| = M.

The associated equal-positive-amplitude subset state is

<a id="eq-4"></a>

> **(4)** |A⟩ = (1/√M) ∑<sub>x ∈ A</sub> |x⟩.

The two extreme cardinalities are product states:

<a id="eq-5"></a>

> **(5)** M = 1 ⇒ |A⟩ = |x⟩,

and

<a id="eq-6"></a>

> **(6)** M = N ⇒ |A⟩ = (1/√N) ∑<sub>x = 0</sub><sup>N−1</sup> |x⟩ = |+⟩<sup>⊗n</sup>.

Set

<a id="eq-7"></a>

> **(7)** d = 2<sup>n/2</sup> = √N

and use the balanced decomposition

<a id="eq-8"></a>

> **(8)** ℋ<sup>(n)</sup> = ℋ<sub>L</sub> ⊗ ℋ<sub>R</sub>, dim ℋ<sub>L</sub> = dim ℋ<sub>R</sub> = d.

Unless stated otherwise, the right subsystem contains the n/2 least significant qubits. Each computational label is then identified with a pair x = (a,b) ∈ [d]<sub>0</sub> × [d]<sub>0</sub>. Define the binary support matrix X<sub>A</sub> ∈ {0,1}<sup>d×d</sup> by

<a id="eq-9"></a>

> **(9)** (X<sub>A</sub>)<sub>ab</sub> = 1[(a,b) ∈ A]

and the normalized coefficient matrix

<a id="eq-10"></a>

> **(10)** C<sub>A</sub> = X<sub>A</sub> / √M.

The state is

<a id="eq-11"></a>

> **(11)** |A⟩ = ∑<sub>a,b = 0</sub><sup>d−1</sup> (C<sub>A</sub>)<sub>ab</sub> |a⟩<sub>L</sub> |b⟩<sub>R</sub>,

Since C<sub>A</sub> has real entries, its right reduced density matrix is

<a id="eq-12"></a>

> **(12)** ρ<sub>A</sub><sup>R</sup> = C<sub>A</sub><sup>†</sup> C<sub>A</sub> = (X<sub>A</sub><sup>†</sup> X<sub>A</sub>) / M.

We measure entanglement in bits,

<a id="eq-13"></a>

> **(13)** S(A) := −Tr(ρ<sub>A</sub><sup>R</sup> log₂ ρ<sub>A</sub><sup>R</sup>).

The nonzero spectrum is the same if the other half is traced out [[22]](#ref-22), [[23]](#ref-23).

The support constraint itself does not prevent maximal balanced entanglement. If M = d and X<sub>A</sub> is a permutation matrix, then

<a id="eq-14"></a>

> **(14)** ρ<sub>A</sub><sup>R</sup> = I<sub>d</sub> / d, S(A) = log₂ d = n/2.

This permutation-support family includes the rainbow pairing familiar from inhomogeneous spin-chain constructions [[24]](#ref-24), [[25]](#ref-25). Together with [Eq. (6)](#eq-6), it shows that support size cannot impose a monotone entanglement law.

For the random ensemble, A is sampled uniformly from the binom(N,M) supports of cardinality M. We write

<a id="eq-15"></a>

> **(15)** S̄<sub>N,M</sub> := 𝔼<sub>A</sub> S(A).

The uniform support law is invariant under qubit permutations, so the ensemble mean is the same for every balanced choice of left and right qubits, although a fixed support can have cut-dependent entropy.

Numerical estimates of the peak support and peak mean entropy are denoted by M̂<sub>n</sub> and Ŝ<sub>n</sub>. The hats distinguish them from exact maximizers of S̄<sub>N,M</sub>.

<a id="sec-exact-moments"></a>

## 3. Exact moments and compressed-support entanglement

Although fixing the support cardinality correlates the entries of X<sub>A</sub>, the first moment of the reduced state and the average purity can be evaluated exactly. Multi-copy moment operators of fixed-size random subset states have also been analyzed in recent pseudorandomness work [[5]](#ref-5), [[4]](#ref-4). Here the balanced bipartition and support-matrix geometry lead to especially direct closed formulas for the ensemble-mean reduced state and the ensemble-average purity. For a cell i = (a,b), write X<sub>i</sub> = 1[i ∈ A]. If i₁, …, i<sub>r</sub> are distinct cells, uniform sampling without replacement gives

<a id="eq-16"></a>

> **(16)** 𝔼<sub>A</sub>[∏<sub>j = 1</sub><sup>r</sup> X<sub>iⱼ</sub>] = (M)<sub>r</sub> / (N)<sub>r</sub>,<br>
> (x)<sub>r</sub> = x(x − 1) ⋯ (x − r + 1).

This factorial-moment identity is the only combinatorial input needed below.

### 3.1. Mean reduced state

<a id="prop-mean-reduced-state"></a>

**Proposition 3.1 (Mean reduced state).** For every 1 ≤ M ≤ N, the ensemble-mean reduced state is

<a id="eq-17"></a>

> **(17)** 𝔼<sub>A</sub> ρ<sub>A</sub><sup>R</sup> = (1/d − β<sub>N,M</sub>) I<sub>d</sub> + β<sub>N,M</sub> J<sub>d</sub>,<br>
> β<sub>N,M</sub> = (M − 1) / [d(N − 1)],

where J<sub>d</sub> is the all-ones matrix. Its uniform-mode eigenvalue is

<a id="eq-18"></a>

> **(18)** λ<sub>mf</sub>(N,M) = 1/d + [(d − 1)(M − 1)] / [d(N − 1)],

and its orthogonal eigenvalue is

<a id="eq-19"></a>

> **(19)** λ<sub>⊥</sub>(N,M) = 1/d − (M − 1) / [d(N − 1)],

with multiplicity d−1.

**Proof.** From [Eq. (12)](#eq-12),

<a id="eq-20"></a>

> **(20)** (ρ<sub>A</sub><sup>R</sup>)<sub>bb′</sub> = (1/M) ∑<sub>a = 0</sub><sup>d−1</sup> X<sub>ab</sub> X<sub>ab′</sub>.

For b = b′, each summand has expectation M/N, hence 𝔼<sub>A</sub>(ρ<sub>A</sub><sup>R</sup>)<sub>bb</sub> = 1/d. For b ≠ b′, the two cells are distinct and [Eq. (16)](#eq-16) gives

<a id="eq-21"></a>

> **(21)** 𝔼<sub>A</sub>(ρ<sub>A</sub><sup>R</sup>)<sub>bb′</sub> = (d/M) (M)<sub>2</sub> / (N)<sub>2</sub><br>
> = (M − 1) / [d(N − 1)].

[Equation (17)](#eq-17) follows. A matrix of the form αI<sub>d</sub> + βJ<sub>d</sub> has eigenvalue α + dβ on the uniform vector and α on its orthogonal complement. ∎

The subscript “mf” emphasizes that λ<sub>mf</sub> is the uniform eigenvalue of the mean matrix, not an identity for 𝔼<sub>A</sub> λ<sub>max</sub>(ρ<sub>A</sub><sup>R</sup>). At M = 1, for example, every individual reduced state is pure, whereas 𝔼<sub>A</sub> ρ<sub>A</sub><sup>R</sup> = I<sub>d</sub>/d. The quantity becomes a useful finite-size mean-field scale in the regime where sampled spectra develop an isolated leading coherent mode.

### 3.2. Average purity

Purity is the second Rényi moment and is a standard diagnostic of random bipartite states [[6]](#ref-6), [[10]](#ref-10). Define

<a id="eq-22"></a>

> **(22)** P(A) := Tr[(ρ<sub>A</sub><sup>R</sup>)²], P̄<sub>N,M</sub> := 𝔼<sub>A</sub> P(A).

<a id="thm-average-purity"></a>

**Theorem 3.2 (Average purity).** For N = d² ≥ 4 and 1 ≤ M ≤ N,

<a id="eq-23"></a>

> **(23)** P̄<sub>N,M</sub> = 1/M<br>
> + 2(d − 1)(M − 1) / [M(N − 1)]<br>
> + [(d − 1)²(M − 1)(M − 2)(M − 3)] / [M(N − 1)(N − 2)(N − 3)].

The counting mechanism is visible directly from

<a id="eq-24"></a>

> **(24)** P(A) = (1/M²) ∑<sub>a,a′,b,b′</sub> X<sub>ab</sub> X<sub>ab′</sub> X<sub>a′b′</sub> X<sub>a′b</sub>.

The four factors are the corners of an ordered rectangle. Depending on whether a = a′ and b = b′, they contain one, two, or four distinct cells. The respective numbers of ordered index quadruples are

<a id="eq-25"></a>

> **(25)** N, 2N(d − 1), N(d − 1)².

Applying [Eq. (16)](#eq-16) to these three cases gives [Eq. (23)](#eq-23). The complete count, including the small-M cases, is given in [Appendix C](#app-combinatorial-derivations). The endpoint checks are immediate: P̄<sub>N,1</sub> = P̄<sub>N,N</sub> = 1.

To expose the relevant powers before choosing a particular scale, let M<sub>N</sub> = cN<sup>γ</sup> + O(1) with fixed c > 0 and 0 < γ < 1. Expanding [Eq. (23)](#eq-23) gives

<a id="eq-26"></a>

> **(26)** P̄<sub>N,M<sub>N</sub></sub> = 2N<sup>−1/2</sup> + c<sup>−1</sup>N<sup>−γ</sup> + c²N<sup>2γ−2</sup><br>
> + o(N<sup>−1/2</sup> + N<sup>−γ</sup> + N<sup>2γ−2</sup>).

The three displayed terms are, respectively, the balanced-cut background, the sparse 1/M contribution, and the four-distinct-cell rectangle contribution. The sparse term is comparable with the background at γ = 1/2, while the rectangle term is comparable with it at γ = 3/4.

<a id="cor-compressed-support-entropy"></a>

**Corollary 3.3 (Power-law support window and optimal balance).** Let M<sub>N</sub> = cN<sup>γ</sup> + O(1) with fixed c > 0 and 0 < γ < 1.

If 1/2 < γ < 3/4, then

<a id="eq-27"></a>

> **(27)** P̄<sub>N,M<sub>N</sub></sub> = 2N<sup>−1/2</sup> (1 + o(1)),

and therefore

<a id="eq-28"></a>

> **(28)** S̄<sub>N,M<sub>N</sub></sub> ≥ n/2 − 1 − o(1).

Within this interval, the slower of the two support-dependent corrections in [Eq. (26)](#eq-26) decays fastest at γ = 2/3. More precisely, for M<sub>N</sub> = cN<sup>2/3</sup> + O(1),

<a id="eq-29"></a>

> **(29)** P̄<sub>N,M<sub>N</sub></sub> = 2N<sup>−1/2</sup> + (c<sup>−1</sup> + c²) N<sup>−2/3</sup> + O(N<sup>−1</sup>),

and the coefficient c<sup>−1</sup> + c² is minimized at c = 2<sup>−1/3</sup>. Hence, for M<sub>N</sub> = 2<sup>−1/3</sup>N<sup>2/3</sup> + O(1),

<a id="eq-30"></a>

> **(30)** S̄<sub>N,M<sub>N</sub></sub> ≥ n/2 − 1 − [3 / (2<sup>5/3</sup> ln 2)] N<sup>−1/6</sup><br>
> + O(N<sup>−1/3</sup>).

At the dense-side boundary M<sub>N</sub> = cN<sup>3/4</sup> + O(1),

<a id="eq-31"></a>

> **(31)** P̄<sub>N,M<sub>N</sub></sub> = (2 + c²) N<sup>−1/2</sup> + c<sup>−1</sup>N<sup>−3/4</sup> + O(N<sup>−1</sup>),

so

<a id="eq-32"></a>

> **(32)** S̄<sub>N,M<sub>N</sub></sub> ≥ n/2 − log₂(2 + c²) − o(1).

**Proof.** [Equation (26)](#eq-26) follows by inserting M<sub>N</sub> = cN<sup>γ</sup> + O(1) and d = N<sup>1/2</sup> into [Eq. (23)](#eq-23). For 1/2 < γ < 3/4, both support-dependent terms are o(N<sup>−1/2</sup>). For every density matrix,

<a id="eq-33"></a>

> **(33)** S(A) ≥ S<sup>(2)</sup>(A) = −log₂ P(A).

Since −log₂ x is convex, Jensen's inequality yields

<a id="eq-34"></a>

> **(34)** S̄<sub>N,M</sub> ≥ 𝔼<sub>A</sub>[−log₂ P(A)] ≥ −log₂ P̄<sub>N,M</sub>.

This proves [Eq. (28)](#eq-28). The decay rate of the larger support-dependent correction is governed by min{γ, 2−2γ}, which is uniquely maximized at γ = 2/3. Substitution of γ = 2/3 and minimization over c give [Eqs. (29)](#eq-29) and [(30)](#eq-30). Substitution of γ = 3/4 gives [Eqs. (31)](#eq-31) and [(32)](#eq-32). Detailed remainders are collected in [Appendix C](#app-combinatorial-derivations). ∎

<a id="remark-optimal-balance"></a>

**Remark 3.4.** At the optimal interior balance, M<sub>N</sub> = 2<sup>(2/3)n + O(1)</sup>, so a vanishing fraction of the computational basis already gives average balanced entanglement within asymptotically one bit of the maximum. The exponent 2/3 uniquely balances the sparse and rectangle corrections; 3/4 instead marks the dense-side boundary, beyond which this purity lower bound loses the maximal balanced slope. These statements do not locate the maximizer of S̄<sub>N,M</sub>: minimizing average purity and maximizing average von Neumann entropy are different problems. They also concern ensemble averages; concentration is separate.

<a id="sec-trajectory"></a>

## 4. Support-size entanglement trajectory

The endpoint and permutation-support examples in [Sec. 2](#sec-model) show that the same family contains product states at M = 1 and M = N and maximally entangled states at M = d. They do not determine the behavior of a typical uniformly random support. Two further exact relations delimit the possible trajectory. First, the Schmidt rank cannot exceed either M or d, so

<a id="eq-35"></a>

> **(35)** S(A) ≤ log₂(rank C<sub>A</sub>) ≤ min{log₂ M, n/2}.

Second, every subset state has overlap

<a id="eq-36"></a>

> **(36)** |⟨+|<sup>⊗n</sup>|A⟩|² = M/N

with the full-support product state. Sparse supports are therefore constrained by rank, whereas dense supports acquire increasing weight in the coherent uniform direction.

<a id="subsec-numerical-protocol"></a>

### 4.1. Numerical protocol

Unless stated otherwise, the numerical calculations use the natural balanced cut, with the low-order n/2 bits assigned to the right subsystem. A support is drawn uniformly without replacement from [N]<sub>0</sub>. If σ<sub>j</sub>(C<sub>A</sub>) are the singular values of the coefficient matrix, then the nonzero reduced-state eigenvalues are

<a id="eq-37"></a>

> **(37)** λ<sub>j</sub> = σ<sub>j</sub>(C<sub>A</sub>)²,

and the entropy is evaluated as

<a id="eq-38"></a>

> **(38)** S(A) = −∑<sub>j: λⱼ > 0</sub> λ<sub>j</sub> log₂ λ<sub>j</sub>.

Using singular values avoids constructing both reduced matrices and is stable for the positive semidefinite spectrum.

At fixed (N,M), independent supports A₁,…,A<sub>R</sub> give the Monte Carlo estimator

<a id="eq-39"></a>

> **(39)** S̄̂<sub>N,M</sub> = (1/R) ∑<sub>r=1</sub><sup>R</sup> S(A<sub>r</sub>).

Error bars on sampled mean trajectories show the standard error s<sub>R</sub>/√R, where s<sub>R</sub> is the sample standard deviation. Figure captions state when the average is instead over balanced cuts of one fixed support.

For comparison with unrestricted random pure states, the exact Page mean for a balanced d × d bipartition is

<a id="eq-40"></a>

> **(40)** S<sub>Page</sub>(n) = [ℋ<sub>d²</sub> − ℋ<sub>d</sub> − (d−1)/(2d)] / ln 2<br>
> = n/2 − 1/(2 ln 2) + O(2<sup>−n</sup>),

where ℋ<sub>k</sub> := ∑<sub>j = 1</sub><sup>k</sup> (1/j) [[7]](#ref-7), [[8]](#ref-8).

Peak estimates are obtained in two stages. A coarse scan first identifies a candidate maximum region. Mean entropies in a local support-size window are then fitted by a quadratic function of x = log₂ M; the interior maximum defines M̂<sub>n</sub> and Ŝ<sub>n</sub>. This procedure estimates a characteristic peak scale rather than an exact integer maximizer. The support space is combinatorial, and the trajectory is broad enough near its maximum that unit-level changes in M̂<sub>n</sub> have no separate physical meaning.

> **Repository audit note.** This paragraph records the manuscript’s stated peak-estimation protocol. The released peak table contains historical point estimates. The original complete scan records, fit windows, and uncertainty estimates are not available in the repository, and the n = 22–30 peak searches have not been independently replayed. The retained data permit verification of the displayed table and finite-range fits; local peak checks do not establish a global maximum. See [reproducibility limits](docs/REPRODUCIBILITY.md).

### 4.2. Rise, broad maximum, and return to the product state

[Figure 1](#fig-trajectory) compares two diagnostics at n = 14. The first fixes the natural balanced cut and averages over independently sampled supports. The second draws one support at each M and averages its entropy over random balanced cuts. The fixed-cardinality ensemble mean is cut-independent, but the second procedure is conditional on one support and therefore tests how a sampled support behaves across cuts.

<a id="fig-trajectory"></a>

![Figure 1: Mean balanced-cut entropy at n = 14](paper/figures/fig1_concentration.png)

**Figure 1.** Mean balanced-cut entropy of random subset states at n = 14. The solid curve fixes one random support at each M and averages over random balanced cuts; the dashed curve fixes the natural cut and averages over random supports. Full-range points use 100 samples and inset points use 500 samples. Error bars are standard errors of the displayed means. The two procedures recover the same global rise–peak–fall trajectory, with small support-specific finite-size offsets. The inset resolves the broad maximum.

Both procedures exhibit the same global pattern. The entropy rises rapidly from zero, reaches a broad maximum at an intermediate support size, and then falls smoothly to zero at full support. The close curves show that support cardinality provides a strong baseline in these samples. They do not imply that all supports, or all cuts of one support, have identical entropy.

The inset also shows why the peak position should be read as a scale rather than as a sharply distinguished integer. On the fixed-cut curve, the sampled means for 552 ≤ M ≤ 993 remain within approximately 1.5 × 10⁻² bits of the largest displayed mean. The independently generated trajectory therefore identifies the same peak region as the n = 14 estimate in [Table 1](#tab-peak-estimates), while displaying the flatness of the local profile.

The descending branch has a specifically coherent origin. Increasing M does not merely increase the number of populated basis vectors: by [Eq. (36)](#eq-36), it also increases the fidelity with |+⟩<sup>⊗n</sup>. Near full support, this positive direction dominates and the state returns to a product state. [Section 5](#sec-spectral-mechanism) turns this observation into a spectral picture involving an isolated uniform mode and a residual bulk.

<a id="subsec-peak-scaling"></a>

### 4.3. Estimated peak and finite-size scaling

For each even n in the studied range, [Table 1](#tab-peak-estimates) lists the retained historical estimates of the peak location M̂<sub>n</sub> and peak mean entropy Ŝ<sub>n</sub>. The protocol above describes the manuscript’s stated estimation procedure; the original complete scan and fit records are unavailable. Hats are retained because the entries are numerical estimates of the ensemble peak, not exact maximizers of S̄<sub>N,M</sub>.

<a id="tab-peak-estimates"></a>

| n | M̂<sub>n</sub> | Ŝ<sub>n</sub> (bits) |
|---:|---:|---:|
| 10 | 107 | 4.072 |
| 12 | 276 | 5.108 |
| 14 | 716 | 6.143 |
| 16 | 1873 | 7.176 |
| 18 | 4934 | 8.196 |
| 20 | 13091 | 9.215 |
| 22 | 34771 | 10.231 |
| 24 | 93018 | 11.242 |
| 26 | 250660 | 12.251 |
| 28 | 672556 | 13.258 |
| 30 | 1836685 | 14.263 |

**Table 1.** Estimated support size M̂<sub>n</sub> at the maximum of the random fixed-cardinality entropy trajectory and the corresponding estimated mean entropy Ŝ<sub>n</sub>. All rows are used in the finite-size fits in [Figure 2](#fig-peak-scaling).

[Figure 2](#fig-peak-scaling) compares these estimates with the exact balanced Page mean. Least-squares fits over all eleven rows give

<a id="eq-41"></a>

> **(41)** Ŝ<sub>n</sub> ≃ 0.5093 n − 0.9901,<br>
> log₂ M̂<sub>n</sub> ≃ 0.70354 n − 0.35773.

<a id="fig-peak-scaling"></a>

![Figure 2: Estimated peak support size and peak entropy](paper/figures/fig2_peak_scaling.png)

**Figure 2.** Estimated peak support size and peak entropy for even n = 10,…,30. Dashed lines are the least-squares finite-size fits to log₂ M̂<sub>n</sub> and Ŝ<sub>n</sub> in [Eq. (41)](#eq-41); the dotted curve is the exact Page mean for a balanced bipartition. The common vertical coordinate reports base-2 values: Ŝ<sub>n</sub> and the Page mean are entropies in bits, whereas log₂ M̂<sub>n</sub> is the base-2 logarithm of the support size. The fits summarize the studied range and are not asserted as asymptotic laws.

Two finite-size observations are robust at the scale of the data. First, the peak entropy tracks the Page mean increasingly closely: the difference S<sub>Page</sub>(n) − Ŝ<sub>n</sub> decreases from about 0.21 bits at n = 10 to about 0.016 bits at n = 30. Second, the estimated peak occupies a shrinking fraction of the computational basis. The ratio M̂<sub>n</sub>/N decreases from approximately 10.4% at n = 10 to 0.171% at n = 30, while the second fit in [Eq. (41)](#eq-41) summarizes its exponential growth in absolute size.

The slopes in [Eq. (41)](#eq-41) are effective finite-range slopes. The Page curve has asymptotic slope 1/2, so the fitted value 0.5093 is not evidence for entropy growth above the maximal balanced scaling. Likewise, 0.70354 is not claimed as a limiting support exponent. A second finite-size diagnostic is the pointwise ratio

<a id="eq-42"></a>

> **(42)** γ<sub>n</sub><sup>eff</sup> := (log₂ M̂<sub>n</sub>)/n.

Across the retained rows it increases monotonically from approximately 0.674 at n = 10 to 0.694 at n = 30. Both these ratios and the fitted slope lie between the correction-balanced interior value 2/3 and the dense-side boundary 3/4 identified in [Corollary 3.3](#cor-compressed-support-entropy). These finite-size trends do not determine an asymptotic exponent: a finite intercept and other subleading corrections can produce the same behavior. The purity scales and the finite-size location of the von Neumann maximum remain different quantities and need not coincide.

The empirical conclusion is nevertheless strong: equal-positive-amplitude states can attain near-Page balanced entanglement while occupying an exponentially compressed support. This entanglement statement is logically distinct from pseudorandomness. The scale M ≍ N<sup>2/3</sup> from [Corollary 3.3](#cor-compressed-support-entropy) lies within fixed-size regimes where polynomial-copy moment averages of random subset states are known to approach their Haar counterparts [[5]](#ref-5), [[4]](#ref-4). The present work isolates the balanced-entanglement consequence and its support-size trajectory; it does not rederive those indistinguishability results.

<a id="sec-spectral-mechanism"></a>

## 5. Sparse-to-dense spectral mechanism

The exact moments above do not determine the von Neumann entropy trajectory, but they identify the two structures that govern its opposite sides. Sparse supports are limited by occupancy and rank. Dense positive supports reinforce the uniform mode of the reduced state. Between them, a broad residual spectrum carries most of the entropy.

<a id="subsec-dense-description"></a>

### 5.1. Separated mode and dense-side ansatz

[Figure 3](#fig-spectrum) shows one sampled reduced-state spectrum near the estimated peak. One eigenvalue is separated from a broad bulk. Induced random density matrices have closely related bulk descriptions [[10]](#ref-10), [[11]](#ref-11), [[12]](#ref-12), but the positive fixed-cardinality ensemble also has the coherent uniform direction identified exactly in [Proposition 3.1](#prop-mean-reduced-state).

<a id="fig-spectrum"></a>

![Figure 3: Reduced-state spectrum at n = 24](paper/figures/fig3_spectral_bulk.png)

**Figure 3.** Reduced-state spectrum for one random subset state at n = 24 and M = 93018, near the estimated entropy peak. The histogram contains the bulk eigenvalues after removing the largest eigenvalue λ₀ := λ<sub>max</sub>(ρ<sub>A</sub><sup>R</sup>). The inset compares the observed λ₀ with the exact uniform eigenvalue λ<sub>mf</sub> of the ensemble-mean reduced state and with the leading approximation M/N. The mean-matrix mode is a finite-size mean-field scale, not an exact formula for the largest eigenvalue of each sample.

To model the dense side, we use λ<sub>mf</sub> as the separated mode and approximate the remaining trace by a balanced random-matrix bulk. If the normalized residual bulk has the Page entropy log₂ d − 1/(2 ln 2), then the entropy decomposes into the binary mixing entropy of the separated mode and the rescaled bulk entropy. Writing h₂(x) := −x log₂ x − (1−x) log₂(1−x), this gives the fixed-mode ansatz

<a id="eq-43"></a>

> **(43)** T<sub>N,M</sub> = h₂(λ<sub>mf</sub>)<br>
> + (1−λ<sub>mf</sub>)[log₂ d − 1/(2 ln 2)],

with the usual endpoint convention and λ<sub>mf</sub> = λ<sub>mf</sub>(N,M) from [Eq. (18)](#eq-18). This is a dense-bulk ansatz, not an exact moment formula. It uses the exact first moment of the fixed-cardinality ensemble but still assumes that the residual spectrum behaves like a normalized random bulk [[7]](#ref-7), [[26]](#ref-26).

The mean uniform mode also gives a spectral interpretation of the 3/4 boundary. Since N = d², [Eq. (18)](#eq-18) can be written exactly as

<a id="eq-44"></a>

> **(44)** λ<sub>mf</sub>(N,M) = (M+d)/[d(d+1)].

In the dense range M ≫ d, this is asymptotic to M/N. For M = cN<sup>γ</sup>, its squared weight therefore scales as c²N<sup>2γ−2</sup>, the same power as the rectangle term in [Eq. (26)](#eq-26). At γ = 3/4 this contribution is of order N<sup>−1/2</sup> and becomes comparable with the balanced random-matrix background at the level of the second moment. This identifies the upper boundary without asserting that the observed largest eigenvalue equals the mean-matrix mode in every sample.

<a id="subsec-sparse-description"></a>

### 5.2. Hypergeometric occupancy and the sparse side

For a fixed right-basis label b, let

<a id="eq-45"></a>

> **(45)** W<sub>b</sub> := ∑<sub>a=0</sub><sup>d−1</sup> X<sub>ab</sub>

be the number of occupied cells in column b. Because the support contains exactly M cells sampled from the N = d² available cells, W<sub>b</sub> is hypergeometric rather than binomial:

<a id="eq-46"></a>

> **(46)** Pr(W<sub>b</sub> = w) = [binom(d,w) binom(N−d,M−w)] / binom(N,M),

where

<a id="eq-47"></a>

> **(47)** max(0, M−N+d) ≤ w ≤ min(d,M).

Since (ρ<sub>A</sub><sup>R</sup>)<sub>bb</sub> = W<sub>b</sub>/M, exchangeability of the columns gives the exact mean diagonal entropy

<a id="eq-48"></a>

> **(48)** D<sub>N,M</sub> := −d ∑<sub>w</sub> Pr(W<sub>b</sub> = w) (w/M) log₂(w/M),

where the w = 0 summand is zero.

Let Δ denote dephasing in the right computational basis. Dephasing does not decrease von Neumann entropy [[22]](#ref-22), so

<a id="eq-49"></a>

> **(49)** S(A) ≤ S(Δ(ρ<sub>A</sub><sup>R</sup>)) = H(diag ρ<sub>A</sub><sup>R</sup>).

Averaging yields the exact upper bound

<a id="eq-50"></a>

> **(50)** S̄<sub>N,M</sub> ≤ D<sub>N,M</sub>.

The quantity D<sub>N,M</sub> is exact as a diagonal-entropy average. It approximates the actual entanglement only in the sparse regime, where column overlaps and the associated off-diagonal coherences remain small.

[Figure 4](#fig-approximation) displays the numerical mean together with the exact diagonal upper bound and the dense ansatz. The diagonal curve describes the initial sparse rise but becomes loose after off-diagonal overlaps proliferate. The dense ansatz captures the descending side semiquantitatively. Their complementary behavior supports a crossover interpretation rather than a single formula valid over the full support range.

<a id="fig-approximation"></a>

![Figure 4: Sparse and dense descriptions at n = 14](paper/figures/fig4_approximation.png)

**Figure 4.** Mean balanced entropy of random subset states at n = 14 compared with the exact expected diagonal entropy D<sub>N,M</sub> and the dense-bulk ansatz T<sub>N,M</sub>. The quantity D<sub>N,M</sub> is an upper bound on the mean von Neumann entropy and is informative on the sparse side; T<sub>N,M</sub> uses the exact mean-matrix uniform mode and captures the descending dense-side trend semiquantitatively. Error bars on the numerical points show the standard error of the sampled mean. The inset resolves the peak region.

The mechanism is therefore simple. At small M, increasing the support populates more row and column sectors and enlarges the effective Schmidt support. At large M, positive overlaps reinforce the uniform coherent mode, which tends continuously to the product-state eigenvalue 1 at full support. The entropy maximum lies between these regimes: a broad bulk has formed, but the coherent direction has not yet become dominant. This interpretation does not by itself locate the exact maximizer.

<a id="sec-structured"></a>

## 6. Structured supports: almost-prime states and Fourier comparison

A random-state benchmark is informative only after the constraints known in advance have been specified. Typical entanglement changes when states are drawn from invariant subspaces or restricted manifolds [[14]](#ref-14), [[16]](#ref-16); particle-number, conserved-charge, and non-Abelian symmetry sectors provide concrete examples [[15]](#ref-15), [[17]](#ref-17), [[18]](#ref-18). The same principle applies to subset states. Support cardinality already has a leading entanglement effect, as also occurs in sparse random-state ensembles [[13]](#ref-13). A systematic departure from a matched random ensemble can therefore reveal organization absent from that ensemble, although it does not by itself identify a unique cause [[27]](#ref-27).

### 6.1. Residue-class entropy ceiling

Low-order residue classes have a direct entanglement consequence for the natural balanced cut. Let ℓ := n/2, choose 0 ≤ t ≤ ℓ, and partition a support A according to its residues modulo 2<sup>t</sup>:

<a id="eq-51"></a>

> **(51)** A<sub>r</sub><sup>(t)</sup> := {x ∈ A : x ≡ r (mod 2<sup>t</sup>)},<br>
> p<sub>r</sub><sup>(t)</sup> := |A<sub>r</sub><sup>(t)</sup>|/|A|.

The t residue bits are the t least significant qubits and lie in the right subsystem.

<a id="prop-residue-ceiling"></a>

**Proposition 6.1 (Residue-class entropy ceiling).** For every subset state |A⟩ and every 0 ≤ t ≤ ℓ,

<a id="eq-52"></a>

> **(52)** S(A) ≤ ℓ−t + H(p<sup>(t)</sup>),<br>
> H(p) := −∑<sub>r</sub> p<sub>r</sub> log₂ p<sub>r</sub>.

Consequently, if A occupies at most q residue classes modulo 2<sup>t</sup>, then

<a id="eq-53"></a>

> **(53)** S(A) ≤ ℓ−t + log₂ q.

If a single residue class is occupied, the corresponding t qubits factor from the state and S(A) ≤ ℓ−t.

The proof is given in [Appendix D](#app-residue-proof). Its main step is to dephase the t residue qubits. The dephased reduced state is a direct sum whose entropy is the Shannon entropy of the residue populations plus the average within-sector entropy, and each sector has only ℓ−t remaining right-subsystem qubits.

The parity case t = 1 is especially transparent. If a fraction ε of the support is even, then

<a id="eq-54"></a>

> **(54)** S(A) ≤ ℓ−1 + h₂(ε),

where h₂ is the binary entropy. A support confined to odd labels has one right-subsystem qubit fixed and therefore loses at least one bit of possible balanced entropy.

### 6.2. Nested matched ensembles for almost-prime supports

Let Ω(x) denote the number of prime factors of x, counted with multiplicity, and define

<a id="eq-55"></a>

> **(55)** U<sub>N,k</sub> := {x ∈ {2,…,N−1} : 1 ≤ Ω(x) ≤ k},<br>
> M<sub>N,k</sub> := |U<sub>N,k</sub>|.

The corresponding subset state is the equal superposition over U<sub>N,k</sub>. The case k = 1 is the Prime state, and increasing k successively adds almost primes. Prime and almost-prime quantum states were introduced and studied in Refs. [[19]](#ref-19), [[20]](#ref-20), [[21]](#ref-21); the number-theoretic terminology follows the classical almost-prime literature [[28]](#ref-28), [[29]](#ref-29).

For a deterministic support B, define the residue-count vector

<a id="eq-56"></a>

> **(56)** M<sub>t</sub>(B) := (|B<sub>0</sub><sup>(t)</sup>|, …, |B<sub>2ᵗ−1</sub><sup>(t)</sup>|)

and let ℰ<sub>t</sub>(B) be the uniform ensemble of supports A satisfying

<a id="eq-57"></a>

> **(57)** M<sub>t</sub>(A) = M<sub>t</sub>(B).

The case t = 0 matches only cardinality. The cases t = 1,2,3 additionally match parity, residues modulo 4, or residues modulo 8. Define

<a id="eq-58"></a>

> **(58)** μ<sub>t</sub>(B) := E<sub>A ∼ ℰₜ(B)</sub> S(A),<br>
> Δ<sub>t</sub>(B) := μ<sub>t</sub>(B) − S(B).

Since μ₀(B) = S̄<sub>N,|B|</sub>, a positive Δ₀(B) means that B is less entangled than a uniformly random support of the same size. This cardinality match fixes the Hilbert-space dimension, support size, amplitude magnitudes, and bipartition, but not residue populations or higher-order arrangements. After the Fourier transform we use the paired definitions

<a id="eq-59"></a>

> **(59)** μ<sub>t</sub><sup>F</sup>(B) := E<sub>A ∼ ℰₜ(B)</sub> S<sub>F</sub>(A),<br>
> Δ<sub>t</sub><sup>F</sup>(B) := μ<sub>t</sub><sup>F</sup>(B) − S<sub>F</sub>(B).

These quantities define a hierarchy of conditional reference ensembles, not a universal scalar measure of structure.

The Prime support illustrates why the first refinement matters. At n = 14 it contains M = 1900 labels, of which only the prime 2 is even. [Equation (54)](#eq-54) therefore gives

<a id="eq-60"></a>

> **(60)** S(U<sub>N,1</sub>) ≤ 6.00649 bits.

Parity alone forces a deficit of at least 0.99351 bits from the balanced maximum of 7 bits. The observed entropy, 4.89935 bits, is lower still. This ceiling concerns the guaranteed gap from the algebraic maximum. The matched-null calculation below asks a distinct conditional question: how much of the deficit from a cardinality-matched random mean remains when the parity counts are held fixed?

For comparison in a global basis, let

<a id="eq-61"></a>

> **(61)** F<sub>N</sub>|x⟩ = (1/√N) ∑<sub>y=0</sub><sup>N−1</sup> exp(2πixy/N)|y⟩

be the positive-exponent quantum Fourier transform, and write S<sub>F</sub>(A) := S(F<sub>N</sub>|A⟩). Fourier observables of the Prime state expose residue-class biases [[21]](#ref-21); more generally, the QFT is not a local basis change and its entangling structure depends on qubit ordering [[30]](#ref-30). We therefore apply the same constrained ensembles before and after F<sub>N</sub>, without assuming entanglement invariance.

<a id="fig-qft-residue-controls"></a>

![Figure 5: Almost-prime cardinality and residue controls](paper/figures/fig5_qft_residue_controls.png)

**Figure 5.** Cardinality- and residue-matched diagnostics for almost-prime subset states at n = 14. (a) The random fixed-cardinality curves are means over 100 supports at each displayed M; shaded regions show one ensemble standard deviation. The + and × markers are the deterministic states |U<sub>N,k</sub>⟩ and F<sub>N</sub>|U<sub>N,k</sub>⟩. (b) Residual deficits Δ<sub>t</sub> and Δ<sub>t</sub><sup>F</sup> for k = 1,2,3 after matching successively the cardinality, parity, mod-4, and mod-8 populations. Each conditional mean uses 1000 independently sampled supports with the exact residue counts of the corresponding arithmetic support. Solid filled markers are computational-basis deficits; dashed open markers are Fourier-basis deficits.

[Figure 5](#fig-qft-residue-controls)(a) shows large cardinality-matched deficits. For k = 1,2,3, the computational-basis values lie below the matched random means by

<a id="eq-62"></a>

> **(62)** 1.09168, 1.45165, 1.24688 bits,

respectively. Panel (b) shows what remains after elementary congruence information is retained by the null ensemble. Matching parity for k = 1 leaves 0.12488 bits, reducing the original deficit by 88.6%. Matching residues modulo 4 for k = 2 leaves 0.23020 bits and reduces it by 84.1%. Matching residues modulo 8 for k = 3 leaves 0.11150 bits and reduces it by 91.1%.

The Fourier-space hierarchy is nearly parallel. The corresponding residuals are 0.11444, 0.21223, and 0.10399 bits, so matching the congruence populations reduces the Fourier-space cardinality deficits by 89.4%, 85.2%, and 91.6%. In each of these six strongest-control comparisons, the structured entropy lies below all 1000 sampled null values. This is a finite-sample rank statement, not a Gaussian-tail estimate.

The application therefore separates three levels of information. Cardinality matching removes the dominant support-size dependence; residue matching reduces most of the low-k deficits associated with elementary binary congruence populations; and a smaller, clearly resolved residual remains beyond the chosen mod-2<sup>t</sup> counts. That residual records arrangement information absent from the matched statistics. It is not claimed to be a unique fingerprint of primality or to isolate every higher-order arithmetic correlation.

<a id="sec-conclusions"></a>

## 7. Conclusions

Random equal-positive-amplitude subset states interpolate between product states at single and full support and a broad intermediate regime of high balanced entanglement. For the uniform fixed-cardinality ensemble, we obtained the exact mean reduced state and average purity. If M = cN<sup>γ</sup>, every fixed 1/2 < γ < 3/4 gives S̄<sub>N,M</sub> ≥ n/2 − 1 − o(1); γ = 2/3 uniquely balances the sparse and rectangle corrections, while γ = 3/4 marks the dense-side boundary. This is a support-size guarantee, not a theorem locating the von Neumann entropy maximum.

The released numerical trajectories rise from zero to a broad maximum and return to zero at full support. The retained historical peak estimates for even n = 10,…,30 approach the balanced Page mean. Across those estimates, the peak support occupies a decreasing basis fraction, while its effective exponent increases over the studied range and remains between 2/3 and 3/4. The hypergeometric diagonal-entropy bound and separated-mode ansatz explain the opposite sides of the curve: sparse supports gain Schmidt sectors, whereas dense positive supports reinforce the uniform coherent mode.

> **Repository audit note.** The n = 10–30 statement above refers to the retained historical peak estimates. Full original entropy trajectories over that entire size range are unavailable; the released full curves cover selected sizes.

The almost-prime application uses this random trajectory as a controlled reference. Cardinality matching removes the leading support-size effect, and parity, mod-4, or mod-8 matching removes most of the remaining low-order congruence deficit before and after the quantum Fourier transform. The residuals record arrangement information beyond the chosen controls without defining a unique arithmetic fingerprint.

Two questions remain immediate. The first is the asymptotic location and width of the von Neumann entropy maximum, including whether it is related to the 3/4 dense-side boundary. The second is concentration: the exact formulas control ensemble averages, while fluctuations determine when near-maximal entanglement is typical for individual supports and across many cuts.

## Data availability

The code, processed figure data, and numerical provenance records supporting this study are openly available in the accompanying GitHub repository [[31]](#ref-31).

> **Repository audit note.** The repository includes the released processed tables and figures, available raw residue-control samples, source code, and provenance and validation records. Historical entropy-peak scan records and full raw support samples for all manuscript figures are not included. The [provenance record](PROVENANCE.md) and [reproducibility guide](docs/REPRODUCIBILITY.md) identify these limits and distinguish replayed evidence from retained historical estimates.

## Acknowledgements

The authors thank S. Carrazza for valuable support with numerical simulations. A substantial part of this work was carried out while R.L. and J.I.L. were affiliated with the Quantum Research Centre, Technology Innovation Institute, United Arab Emirates, and while R.L. was affiliated with the Departament de Física Quàntica i Astrofísica and Institut de Ciències del Cosmos, Universitat de Barcelona, Spain. G.S. acknowledges financial support through the Spanish MINECO grant PID2021-127726NB-I00, the CSIC Research Platform on Quantum Technologies PTI-001, and the QUANTUM ENIA project Quantum Spain through the RTRP-Next Generation within the framework of the Digital Spain 2026 Agenda.

## Conflict of interest

The authors declare no competing interests.

<a id="app-renyi"></a>

## Appendix A. Rényi-entropy trajectories

For α > 0, α ≠ 1, the Rényi entropy of a support A is

<a id="eq-63"></a>

> **(63)** S<sup>(α)</sup>(A) := [log₂ Tr[(ρ<sub>A</sub><sup>R</sup>)<sup>α</sup>]] / (1 − α),

with the von Neumann entropy recovered as α → 1 [[32]](#ref-32). The min-entropy is

<a id="eq-64"></a>

> **(64)** S<sup>(∞)</sup>(A) = −log₂ λ<sub>max</sub>(ρ<sub>A</sub><sup>R</sup>).

We write S̄<sub>N,M</sub><sup>(α)</sup> for the fixed-cardinality ensemble mean. For every individual reduced state,

<a id="eq-65"></a>

> **(65)** S<sup>(1)</sup>(A) ≥ S<sup>(2)</sup>(A) ≥ S<sup>(∞)</sup>(A),

and the same ordering is preserved after averaging.

[Figure 6](#fig-renyi) shows that the rise–peak–fall trajectory is not specific to the von Neumann entropy. All three orders vanish at M = 1 and M = N and attain an intermediate maximum. The peak moves toward smaller supports as α increases. In the refined n = 14 scan, the largest sampled means occur at

<a id="eq-66"></a>

> **(66)** M ≈ 712, 518, 259

for α = 1, 2, ∞, respectively. The corresponding mean entropies are approximately 6.144, 5.766, and 4.484 bits. Larger Rényi orders give more weight to the largest eigenvalue, so they respond earlier to the emergence of the leading coherent mode discussed in [Section 5.1](#subsec-dense-description).

<a id="fig-renyi"></a>

![Mean Rényi entropies for random subset states at n = 14, with a full-range panel and an enlarged peak-region panel.](paper/figures/fig6_renyi.png)

**Figure 6.** Mean Rényi entropies of orders 1, 2, and ∞ for random subset states at n = 14. The upper panel displays the full support range; the lower panel resolves the peak region. Full-range points use 200 random supports per M, and peak-region points use 1000. Error bars are standard errors of the displayed means. Increasing the order lowers the entropy and shifts the sampled maximum toward smaller support sizes.

<a id="app-partitions"></a>

## Appendix B. Entropy distributions over balanced cuts

The ensemble mean S̄<sub>N,M</sub> is independent of the chosen balanced cut, but a fixed support can have cut-dependent entanglement. To examine this dependence directly, we fix one random subset state at the tabulated n = 20 peak scale, M = 13091, and sample 1000 balanced cuts. Complementary cuts are identified by fixing one qubit on the right subsystem. For comparison, we fix one complex Haar-random state and evaluate it on an independently sampled set of 1000 balanced cuts. Haar vectors are generated by normalizing complex Gaussian vectors, a standard realization of Haar measure [[23]](#ref-23), [[33]](#ref-33). This experiment concerns two selected states; it is not an ensemble average over supports or Haar vectors.

[Figure 7](#fig-partition) places both distributions on one entropy axis. For the subset state, the sample mean and standard deviation are

<a id="eq-67"></a>

> **(67)** 9.21520 bits, 2.62 × 10<sup>−3</sup> bits,

with sampled range [9.20546, 9.22237] bits. For the Haar state they are

<a id="eq-68"></a>

> **(68)** 9.27871 bits, 7.09 × 10<sup>−4</sup> bits,

with range [9.27643, 9.28123] bits. The exact balanced Page mean at n = 20 is 9.27865 bits. Thus the selected Haar state is tightly concentrated around the Page benchmark, consistent with generic-entanglement concentration [[9]](#ref-9), whereas the selected subset state has a mean lower by about 6.35 × 10<sup>−2</sup> bits and a cut-to-cut standard deviation about 3.7 times larger.

<a id="fig-partition"></a>

![Entropy distributions over 1000 sampled balanced cuts for one subset state and one complex Haar state, shown on the same horizontal scale.](paper/figures/fig7_partitions.png)

**Figure 7.** Entropy distributions over 1000 sampled balanced cuts for one peak-support subset state with n = 20 and M = 13091 (upper panel) and one complex Haar state (lower panel). Both panels use the same horizontal scale. Dashed lines mark the sample means and the dotted line marks the exact Page mean. The figure describes the cut distribution of two selected states; it does not assert maximal entanglement for every cut or a concentration theorem for the subset-state ensemble.

The two sampled distributions do not overlap. Nevertheless, the subset state remains highly entangled across every sampled cut: its lowest sampled entropy is more than 9.20 bits out of the maximum 10 bits. The comparison clarifies the meaning of the near-Page peak in [Table 1](#tab-peak-estimates). Together with the trajectory data, this selected-state example illustrates how random positive-amplitude subset states near the peak can have volume-law entanglement that is uniformly high over many cuts while retaining spectral and cut-dependent structure distinct from that of a Haar-random state.

<a id="app-combinatorial-derivations"></a>

## Appendix C. Combinatorial derivations for the fixed-cardinality ensemble

This appendix supplies the detailed count behind [Theorem 3.2](#thm-average-purity), the asymptotic expansion used in [Corollary 3.3](#cor-compressed-support-entropy), and the sparse diagonal-entropy formula.

### C.1. Indicator moments under sampling without replacement

Let A be uniform among all M-element subsets of [N]₀, and let X<sub>i</sub> = 1[i ∈ A]. For r distinct cells i₁, …, i<sub>r</sub>, the event i₁, …, i<sub>r</sub> ∈ A leaves M − r support elements to be chosen from the remaining N − r cells. Hence

<a id="eq-69"></a>

> **(69)** E<sub>A</sub>[∏<sub>j = 1</sub><sup>r</sup> X<sub>iⱼ</sub>] = Pr(i₁, …, i<sub>r</sub> ∈ A)<br>
> = binom(N − r, M − r) / binom(N, M)<br>
> = (M)<sub>r</sub> / (N)<sub>r</sub>.

Here binom(u, v) denotes the binomial coefficient. Repeated indicator factors may be removed because X<sub>i</sub><sup>k</sup> = X<sub>i</sub>.

### C.2. Rectangle count for the average purity

Starting from [Eq. (12)](#eq-12),

<a id="eq-70"></a>

> **(70)** Tr[(ρ<sub>A</sub><sup>R</sup>)²]<br>
> = ∑<sub>b,b′</sub> (ρ<sub>A</sub><sup>R</sup>)<sub>bb′</sub> (ρ<sub>A</sub><sup>R</sup>)<sub>b′b</sub><br>
> = (1/M²) ∑<sub>a,a′,b,b′</sub> X<sub>ab</sub> X<sub>ab′</sub> X<sub>a′b′</sub> X<sub>a′b</sub>.

For every ordered quadruple (a, a′, b, b′), the factors correspond to the four corners

<a id="eq-71"></a>

> **(71)** (a, b), (a, b′), (a′, b′), (a′, b).

There are three possible coincidence patterns.

**One distinct cell.** When a = a′ and b = b′, all factors coincide. There are d² = N ordered quadruples, each contributing expectation (M)₁/(N)₁.

**Two distinct cells.** If a = a′ and b ≠ b′, the two cells are (a, b) and (a, b′); there are d²(d − 1) = N(d − 1) ordered choices. The case a ≠ a′ and b = b′ gives the same count. Together these contribute 2N(d − 1) terms with expectation (M)₂/(N)₂.

**Four distinct cells.** When a ≠ a′ and b ≠ b′, all four corners are distinct. There are d²(d − 1)² = N(d − 1)² ordered rectangles, each with expectation (M)₄/(N)₄.

Combining the three cases gives

<a id="eq-72"></a>

> **(72)** P̄<sub>N,M</sub> = [N (M)₁/(N)₁<br>
> + 2N(d − 1) (M)₂/(N)₂<br>
> + N(d − 1)² (M)₄/(N)₄] / M².

Cancelling common factors gives [Eq. (23)](#eq-23). The falling-factorial form also covers M < 4 automatically because (M)₄ = 0.

### C.3. Power-law window and asymptotic entropy bounds

Let M = cN<sup>γ</sup> + O(1) with fixed c > 0, 0 < γ < 1, and d = N<sup>1/2</sup>. The three terms in [Eq. (23)](#eq-23) satisfy

<a id="eq-73"></a>

> **(73)** 1/M = c<sup>−1</sup>N<sup>−γ</sup> + O(N<sup>−2γ</sup>),

<a id="eq-74"></a>

> **(74)** [2(d − 1)(M − 1)] / [M(N − 1)]<br>
> = [2/(d + 1)] (1 − 1/M)<br>
> = 2N<sup>−1/2</sup> + O(N<sup>−1</sup> + N<sup>−γ−1/2</sup>),

<a id="eq-75"></a>

> **(75)** [(d − 1)²(M − 1)(M − 2)(M − 3)]<br>
> ÷ [M(N − 1)(N − 2)(N − 3)]<br>
> = c²N<sup>2γ−2</sup><br>
> + O(N<sup>2γ−5/2</sup> + N<sup>γ−2</sup>).

Consequently,

<a id="eq-76"></a>

> **(76)** P̄<sub>N,M</sub> = 2N<sup>−1/2</sup> + c<sup>−1</sup>N<sup>−γ</sup> + c²N<sup>2γ−2</sup><br>
> + O(N<sup>−1</sup> + N<sup>−2γ</sup> + N<sup>−γ−1/2</sup><br>
> + N<sup>2γ−5/2</sup> + N<sup>γ−2</sup>),

which implies [Eq. (26)](#eq-26).

[Equation (76)](#eq-76) gives the three regimes in [Corollary 3.3](#cor-compressed-support-entropy): at γ = 1/2 the sparse term joins the background, throughout 1/2 < γ < 3/4 both support-dependent terms are o(N<sup>−1/2</sup>), and at γ = 3/4 the rectangle term joins the background. Within the interior window, the slower correction is governed by min{γ, 2 − 2γ} and is optimized by

<a id="eq-77"></a>

> **(77)** γ = 2 − 2γ = 2/3.

At this exponent the remaining coefficient is c<sup>−1</sup> + c²; the condition −c<sup>−2</sup> + 2c = 0 gives c = 2<sup>−1/3</sup>. Substitution into [Eq. (76)](#eq-76), followed by [Eq. (34)](#eq-34) and log₂(1 + x) = x/(ln 2) + O(x²), yields [Eqs. (29)](#eq-29)–[(30)](#eq-30).

### C.4. Hypergeometric diagonal entropy

A fixed column contains d of the N = d² support cells. Drawing an M-element support without replacement therefore gives

<a id="eq-78"></a>

> **(78)** W<sub>b</sub> ∼ Hypergeometric(N, d, M),

which is [Eq. (46)](#eq-46). Because the d columns are exchangeable,

<a id="eq-79"></a>

> **(79)** E<sub>A</sub> H(diag ρ<sub>A</sub><sup>R</sup>)<br>
> = −∑<sub>b = 0</sub><sup>d−1</sup> E<sub>A</sub>[(W<sub>b</sub>/M) log₂(W<sub>b</sub>/M)]<br>
> = −d ∑<sub>w</sub> Pr(W<sub>b</sub> = w) (w/M) log₂(w/M)<br>
> = D<sub>N,M</sub>.

Computational-basis dephasing is a unital channel and therefore does not decrease von Neumann entropy. Hence

<a id="eq-80"></a>

> **(80)** S(ρ<sub>A</sub><sup>R</sup>) ≤ S(Δ(ρ<sub>A</sub><sup>R</sup>)) = H(diag ρ<sub>A</sub><sup>R</sup>),

and averaging proves [Eq. (50)](#eq-50).

<a id="app-residue-proof"></a>

## Appendix D. Residue-class bound and constrained ensembles

### D.1. Proof of Proposition 6.1

With ℓ = n/2, let the right subsystem factor as

<a id="eq-81"></a>

> **(81)** ℋ<sub>R</sub> = ℋ<sub>R′</sub> ⊗ ℋ<sub>Q</sub>,<br>
> dim ℋ<sub>Q</sub> = 2<sup>t</sup>,<br>
> dim ℋ<sub>R′</sub> = 2<sup>ℓ−t</sup>,

where Q contains the t least significant qubits. For every occupied residue r, strip those bits from the labels in A<sub>r</sub><sup>(t)</sup> and define the normalized branch state |ψ<sub>r</sub>⟩ ∈ ℋ<sub>L</sub> ⊗ ℋ<sub>R′</sub>. Then

<a id="eq-82"></a>

> **(82)** |A⟩ = ∑<sub>r : pᵣ⁽ᵗ⁾ > 0</sub> √(p<sub>r</sub><sup>(t)</sup>) |ψ<sub>r</sub>⟩ |r⟩<sub>Q</sub>.

Let ρ<sub>R</sub> = Tr<sub>L</sub> |A⟩⟨A| and dephase Q in its computational basis. The dephasing channel is a uniform mixture of commuting phase unitaries, so concavity and unitary invariance of the von Neumann entropy imply

<a id="eq-83"></a>

> **(83)** S(ρ<sub>R</sub>) ≤ S(𝒟<sub>Q</sub>(ρ<sub>R</sub>)).

Writing

<a id="eq-84"></a>

> **(84)** ρ<sub>R′</sub><sup>(r)</sup> := Tr<sub>L</sub> |ψ<sub>r</sub>⟩⟨ψ<sub>r</sub>|,

the dephased state is block diagonal:

<a id="eq-85"></a>

> **(85)** 𝒟<sub>Q</sub>(ρ<sub>R</sub>) = ⨁<sub>r</sub> p<sub>r</sub><sup>(t)</sup> ρ<sub>R′</sub><sup>(r)</sup>.

The entropy of a direct sum gives

<a id="eq-86"></a>

> **(86)** S(𝒟<sub>Q</sub>(ρ<sub>R</sub>))<br>
> = H(p<sup>(t)</sup>) + ∑<sub>r</sub> p<sub>r</sub><sup>(t)</sup> S(ρ<sub>R′</sub><sup>(r)</sup>)

<a id="eq-87"></a>

> **(87)** S(𝒟<sub>Q</sub>(ρ<sub>R</sub>)) ≤ H(p<sup>(t)</sup>) + ℓ − t,

because each ρ<sub>R′</sub><sup>(r)</sup> acts on 2<sup>ℓ−t</sup> dimensions. This proves [Eq. (52)](#eq-52) in [Proposition 6.1](#prop-residue-ceiling). If only q residues are occupied, H(p<sup>(t)</sup>) ≤ log₂ q, proving [Eq. (53)](#eq-53). For a single occupied residue, [Eq. (82)](#eq-82) contains one term, so the t qubits in Q factor exactly from the remaining state.

### D.2. Conditional support sampler

For a deterministic support B, the ensemble ℰ<sub>t</sub>(B) is sampled by partitioning [N]₀ into the 2<sup>t</sup> residue classes and choosing M<sub>r</sub><sup>(t)</sup>(B) labels uniformly without replacement from class r, independently across classes. This produces the uniform distribution over all supports with the prescribed count vector because every admissible support has probability

<a id="eq-88"></a>

> **(88)** ∏<sub>r = 0</sub><sup>2ᵗ−1</sup> [binom(N/2<sup>t</sup>, M<sub>r</sub><sup>(t)</sup>(B))]<sup>−1</sup>.

The calculations in [Figure 5(b)](#fig-qft-residue-controls) use 1000 supports for each pair (k, t) with k = 1, 2, 3 and t = 0, 1, 2, 3. Position- and Fourier-basis entropies are evaluated on the same sampled supports, giving a paired comparison.

<a id="references"></a>

## References

<a id="ref-1"></a>

**[1]** Alex Bredariol Grilo, Iordanis Kerenidis, and Jamie Sikora. “QMA with subset state witnesses.” *Chicago Journal of Theoretical Computer Science* **2016**(4), 4:1–4:17 (2016). DOI: [10.4086/cjtcs.2016.004](https://doi.org/10.4086/cjtcs.2016.004).

<a id="ref-2"></a>

**[2]** Zhengfeng Ji, Yi-Kai Liu, and Fang Song. “Pseudorandom Quantum States.” In *Advances in Cryptology – CRYPTO 2018, Part III*, edited by Hovav Shacham and Alexandra Boldyreva, *Lecture Notes in Computer Science* **10993**, pp. 126–152. Springer, Cham (2018). DOI: [10.1007/978-3-319-96878-0_5](https://doi.org/10.1007/978-3-319-96878-0_5).

<a id="ref-3"></a>

**[3]** Scott Aaronson, Adam Bouland, Bill Fefferman, Soumik Ghosh, Umesh Vazirani, Chenyi Zhang, and Zixin Zhou. “Quantum Pseudoentanglement.” In *15th Innovations in Theoretical Computer Science Conference (ITCS 2024)*, edited by Venkatesan Guruswami, *Leibniz International Proceedings in Informatics (LIPIcs)* **287**, 2:1–2:21. Schloss Dagstuhl – Leibniz-Zentrum für Informatik, Dagstuhl, Germany (2024). DOI: [10.4230/LIPIcs.ITCS.2024.2](https://doi.org/10.4230/LIPIcs.ITCS.2024.2).

<a id="ref-4"></a>

**[4]** Fernando Granha Jeronimo, Nir Magrafta, and Pei Wu. *Pseudorandom and pseudoentangled states from subset states* (2023). [arXiv:2312.15285 [quant-ph]](https://arxiv.org/abs/2312.15285).

<a id="ref-5"></a>

**[5]** Tudor Giurgica-Tiron and Adam Bouland. *Pseudorandomness from subset states* (2023). [arXiv:2312.09206 [quant-ph]](https://arxiv.org/abs/2312.09206).

<a id="ref-6"></a>

**[6]** Elihu Lubkin. “Entropy of an n-system from its correlation with a k-reservoir.” *Journal of Mathematical Physics* **19**(5), pp. 1028–1031 (1978). DOI: [10.1063/1.523763](https://doi.org/10.1063/1.523763).

<a id="ref-7"></a>

**[7]** Don N. Page. “Average entropy of a subsystem.” *Physical Review Letters* **71**(9), pp. 1291–1294 (1993). DOI: [10.1103/PhysRevLett.71.1291](https://doi.org/10.1103/PhysRevLett.71.1291).

<a id="ref-8"></a>

**[8]** Siddhartha Sen. “Average entropy of a quantum subsystem.” *Physical Review Letters* **77**(1), pp. 1–3 (1996). DOI: [10.1103/PhysRevLett.77.1](https://doi.org/10.1103/PhysRevLett.77.1).

<a id="ref-9"></a>

**[9]** Patrick Hayden, Debbie W. Leung, and Andreas Winter. “Aspects of Generic Entanglement.” *Communications in Mathematical Physics* **265**(1), pp. 95–117 (2006). DOI: [10.1007/s00220-006-1535-6](https://doi.org/10.1007/s00220-006-1535-6).

<a id="ref-10"></a>

**[10]** Karol Życzkowski and Hans-Jürgen Sommers. “Induced Measures in the Space of Mixed Quantum States.” *Journal of Physics A: Mathematical and General* **34**(35), pp. 7111–7125 (2001). DOI: [10.1088/0305-4470/34/35/335](https://doi.org/10.1088/0305-4470/34/35/335).

<a id="ref-11"></a>

**[11]** Ion Nechita. “Asymptotics of random density matrices.” *Annales Henri Poincaré* **8**(8), pp. 1521–1538 (2007). DOI: [10.1007/s00023-007-0345-5](https://doi.org/10.1007/s00023-007-0345-5).

<a id="ref-12"></a>

**[12]** Santosh Kumar and Akhilesh Pandey. “Entanglement in Random Pure States: Spectral Density and Average von Neumann Entropy.” *Journal of Physics A: Mathematical and Theoretical* **44**(44), 445301 (2011). DOI: [10.1088/1751-8113/44/44/445301](https://doi.org/10.1088/1751-8113/44/44/445301).

<a id="ref-13"></a>

**[13]** Giuseppe De Tomasi and Ivan M. Khaymovich. “Multifractality Meets Entanglement: Relation for Nonergodic Extended States.” *Physical Review Letters* **124**(20), 200602 (2020). DOI: [10.1103/PhysRevLett.124.200602](https://doi.org/10.1103/PhysRevLett.124.200602).

<a id="ref-14"></a>

**[14]** Yoshifumi Nakata and Mio Murao. “Generic Entanglement Entropy for Quantum States with Symmetry.” *Entropy* **22**(6), 684 (2020). DOI: [10.3390/e22060684](https://doi.org/10.3390/e22060684).

<a id="ref-15"></a>

**[15]** Lev Vidmar and Marcos Rigol. “Entanglement Entropy of Eigenstates of Quantum Chaotic Hamiltonians.” *Physical Review Letters* **119**(22), 220603 (2017). DOI: [10.1103/PhysRevLett.119.220603](https://doi.org/10.1103/PhysRevLett.119.220603).

<a id="ref-16"></a>

**[16]** Eugenio Bianchi, Lucas Hackl, Mario Kieburg, Marcos Rigol, and Lev Vidmar. “Volume-Law Entanglement Entropy of Typical Pure Quantum States.” *PRX Quantum* **3**(3), 030201 (2022). DOI: [10.1103/PRXQuantum.3.030201](https://doi.org/10.1103/PRXQuantum.3.030201).

<a id="ref-17"></a>

**[17]** Sara Murciano, Pasquale Calabrese, and Lorenzo Piroli. “Symmetry-resolved Page curves.” *Physical Review D* **106**(4), 046015 (2022). DOI: [10.1103/PhysRevD.106.046015](https://doi.org/10.1103/PhysRevD.106.046015).

<a id="ref-18"></a>

**[18]** Rohit Patil, Lucas Hackl, George R. Fagan, and Marcos Rigol. “Average pure-state entanglement entropy in spin systems with SU(2) symmetry.” *Physical Review B* **108**(24), 245101 (2023). DOI: [10.1103/PhysRevB.108.245101](https://doi.org/10.1103/PhysRevB.108.245101).

<a id="ref-19"></a>

**[19]** José Ignacio Latorre and Germán Sierra. “Quantum computation of prime number functions.” *Quantum Information & Computation* **14**(7&8), pp. 577–588 (2014). DOI: [10.26421/QIC14.7-8-3](https://doi.org/10.26421/QIC14.7-8-3).

<a id="ref-20"></a>

**[20]** José Ignacio Latorre and Germán Sierra. “There is entanglement in the primes.” *Quantum Information & Computation* **15**(7&8), pp. 622–659 (2015). DOI: [10.26421/QIC15.7-8-6](https://doi.org/10.26421/QIC15.7-8-6).

<a id="ref-21"></a>

**[21]** Diego García-Martín, Eduard Ribas, Stefano Carrazza, José I. Latorre, and Germán Sierra. “The Prime state and its quantum relatives.” *Quantum* **4**, 371 (2020). DOI: [10.22331/q-2020-12-11-371](https://doi.org/10.22331/q-2020-12-11-371).

<a id="ref-22"></a>

**[22]** Michael A. Nielsen and Isaac L. Chuang. *Quantum Computation and Quantum Information*. 10th anniversary edition. Cambridge University Press, Cambridge (2010). DOI: [10.1017/CBO9780511976667](https://doi.org/10.1017/CBO9780511976667).

<a id="ref-23"></a>

**[23]** Ingemar Bengtsson and Karol Życzkowski. *Geometry of Quantum States: An Introduction to Quantum Entanglement*. Second edition. Cambridge University Press, Cambridge (2017). DOI: [10.1017/9781139207010](https://doi.org/10.1017/9781139207010).

<a id="ref-24"></a>

**[24]** G. Vitagliano, A. Riera, and J. I. Latorre. “Volume-law scaling for the entanglement entropy in spin-1/2 chains.” *New Journal of Physics* **12**(11), 113049 (2010). DOI: [10.1088/1367-2630/12/11/113049](https://doi.org/10.1088/1367-2630/12/11/113049).

<a id="ref-25"></a>

**[25]** Giovanni Ramírez, Javier Rodríguez-Laguna, and Germán Sierra. “From conformal to volume law for the entanglement entropy in exponentially deformed critical spin 1/2 chains.” *Journal of Statistical Mechanics: Theory and Experiment* **2014**(10), P10004 (2014). DOI: [10.1088/1742-5468/2014/10/P10004](https://doi.org/10.1088/1742-5468/2014/10/P10004).

<a id="ref-26"></a>

**[26]** Ruge Lin. “Entanglement trajectory and its boundary.” *Quantum* **8**, 1282 (2024). DOI: [10.22331/q-2024-03-14-1282](https://doi.org/10.22331/q-2024-03-14-1282).

<a id="ref-27"></a>

**[27]** Masudul Haque, Paul A. McClarty, and Ivan M. Khaymovich. “Entanglement of midspectrum eigenstates of chaotic many-body systems: Reasons for deviation from random ensembles.” *Physical Review E* **105**(1), 014109 (2022). DOI: [10.1103/PhysRevE.105.014109](https://doi.org/10.1103/PhysRevE.105.014109).

<a id="ref-28"></a>

**[28]** A. A. Rényi. “On the representation of an even number as the sum of a single prime and single almost-prime number.” *Izvestiya Akademii Nauk SSSR. Seriya Matematicheskaya* **12**(1), pp. 57–78 (1948). [Math-Net.Ru](https://www.mathnet.ru/eng/im3018).

<a id="ref-29"></a>

**[29]** D. Roger Heath-Brown. “Almost-primes in arithmetic progressions and short intervals.” *Mathematical Proceedings of the Cambridge Philosophical Society* **83**(3), pp. 357–375 (1978). DOI: [10.1017/S0305004100054657](https://doi.org/10.1017/S0305004100054657).

<a id="ref-30"></a>

**[30]** Jielun Chen, E. M. Stoudenmire, and Steven R. White. “Quantum Fourier Transform Has Small Entanglement.” *PRX Quantum* **4**(4), 040318 (2023). DOI: [10.1103/PRXQuantum.4.040318](https://doi.org/10.1103/PRXQuantum.4.040318).

<a id="ref-31"></a>

**[31]** Ruge Lin. *Data and code for support-size entanglement trajectories of random subset states* (2026). GitHub repository: [GoGoKo699/Subset-states](https://github.com/GoGoKo699/Subset-states).

<a id="ref-32"></a>

**[32]** Alfréd Rényi. “On Measures of Entropy and Information.” In *Proceedings of the Fourth Berkeley Symposium on Mathematical Statistics and Probability*, vol. **1**, pp. 547–561. University of California Press, Berkeley, CA (1961). [Project Euclid](https://projecteuclid.org/euclid.bsmsp/1200512181).

<a id="ref-33"></a>

**[33]** Antonio Anna Mele. “Introduction to Haar Measure Tools in Quantum Information: A Beginner's Tutorial.” *Quantum* **8**, 1340 (2024). DOI: [10.22331/q-2024-05-08-1340](https://doi.org/10.22331/q-2024-05-08-1340).
