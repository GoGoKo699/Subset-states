# Final target-reader and submission-level review

**Manuscript:** *Support-size entanglement trajectories of random subset states*  
**Target:** *Journal of Physics A: Mathematical and Theoretical*  
**Review date:** 17 August 2026

## Verdict

The manuscript is scientifically coherent, internally consistent, readable as a quantum-information paper, and ready for author-level confirmation and repository release. I found no remaining mathematical, numerical, citation, or LaTeX defect that blocks submission.

The main result hierarchy is now clear:

1. exact ensemble mean reduced state;
2. exact average-purity theorem;
3. rigorous near-maximal average-entanglement guarantee at support size \(2^{(2/3)n+O(1)}\);
4. finite-size rise--peak--fall von Neumann trajectory through \(n=30\);
5. sparse hypergeometric and dense separated-mode explanations;
6. cardinality- and residue-matched application to almost-prime supports.

The paper no longer depends on the historical arithmetic framing for its identity. The almost-prime analysis functions as a controlled application of the general subset-state framework.

## Substantive refinements made in the final pass

### Prior-work positioning

The exact-moment section now acknowledges recent multi-copy moment analyses of fixed-size random subset states. The manuscript states precisely that its contribution is the direct balanced-bipartition formula for the ensemble-mean reduced state and average purity, rather than implying that fixed-size subset-state moments were previously unstudied.

The finite-size scaling section now relates the rigorous \(N^{2/3}\) support scale to known pseudorandom-subset-state regimes while explicitly distinguishing entanglement from indistinguishability.

### Claim calibration

The abstract now says that the estimated peak support grows exponentially while occupying a decreasing fraction of the basis **over the studied range**. It does not convert the finite-size fit into an asymptotic law.

Single-instance results are consistently described as sampled or selected states. Figure 3 is one sampled spectrum; Figure 7 concerns two selected states. No concentration theorem is inferred from those examples.

### Dense-side interpretation

The dense approximation is written in the transparent form

\[
T_{N,M}=h_2(\lambda_{\mathrm{mf}})
 +(1-\lambda_{\mathrm{mf}})
 \left(\log_2 d-\frac{1}{2\ln2}\right),
\]

which separates the binary mixing entropy of the coherent mode from the rescaled Page-like bulk entropy. This is algebraically identical to the previous expression.

### Residue-matched application

Fourier-space matched means and deficits are now explicitly defined:

\[
\mu_t^F(B)=\mathbb E_{A\sim\mathcal E_t(B)}S_F(A),
\qquad
\Delta_t^F(B)=\mu_t^F(B)-S_F(B).
\]

The manuscript uses **deficit reduction** rather than causal language such as “fraction explained.” It also distinguishes two different statements:

- the exact parity ceiling gives a guaranteed gap from the algebraic balanced maximum;
- the parity-matched null ensemble measures the residual gap from a cardinality-matched random mean.

### Readability

The long Section 3 heading was manually broken between “support” and “entanglement” to avoid an unattractive hyphenated section title. No substantive compression or expansion was needed elsewhere.

## Independent scientific validation

The final validation script exhaustively enumerated all 65,535 nonempty supports at \(n=4\).

| Check | Maximum discrepancy or violation |
|---|---:|
| Ensemble-mean reduced matrix | \(2.07\times10^{-14}\) |
| Exact average-purity formula | \(8.88\times10^{-16}\) |
| Residue ceiling, \(t=1\) | 0 |
| Residue ceiling, \(t=2\) | 0 |
| Independent entropy implementations | \(5.55\times10^{-16}\) bits |

The Table I regressions were reconstructed exactly:

\[
\widehat S_n=0.5093n-0.9900909\ldots,
\]

\[
\log_2\widehat M_n
=0.703540951\ldots n-0.357734368\ldots.
\]

The Page gap and support fraction decrease monotonically across all eleven retained rows. The old and revised dense-ansatz formulas agree to \(2.66\times10^{-15}\) bits over the complete \(n=14\) support range.

The residue-control summary reproduces the displayed residual deficits, deficit reductions, Prime-state parity ceiling, and all six finite-sample rank statements. Full numerical output is in `FINAL_SCIENTIFIC_VALIDATION.txt`.

## Bibliography audit

- 33 bibliography entries;
- 33 entries cited;
- no missing keys;
- no unused keys;
- no duplicate keys;
- no Biber warning or error.

The citations are generous but role-specific. Prior work is cited for terminology, random-state benchmarks, constrained ensembles, subset-state moments, sparse supports, prime-state structure, and QFT entanglement without obscuring which propositions and theorems are new here.

## LaTeX and PDF validation

The clean source compiles with

```text
pdflatex main.tex
biber main
pdflatex main.tex
pdflatex main.tex
```

Final status:

- 24 A4 pages;
- approximately 5,581 body words by `texcount`;
- 189 abstract prose words, excluding its two displayed formulas;
- no unresolved citation or reference;
- no LaTeX/package warning;
- no overfull or underfull box;
- all fonts embedded;
- PDF openable and unencrypted;
- not a scanned PDF;
- all pages rendered and visually inspected;
- no clipping, overlap, black box, or broken glyph.

## Journal fit

The paper is best submitted as a regular **Paper** in the journal's quantum-mechanics and quantum-information section, with additional relevance to physical combinatorics, asymptotics, number theory, and random matrices. The manuscript is shorter in words than a typical long JPhysA Paper, but its 24-page review layout contains seven figures, formal results, and four proof/diagnostic appendices. Padding would reduce readability and is not recommended.

## Author and release confirmations still required

These are administrative or authorial decisions, not scientific defects:

1. Confirm the author order and each current affiliation.
2. Confirm the corresponding author, current email addresses, and submission-system ORCID information.
3. Confirm the acknowledgement wording and all funder/grant identifiers.
4. Confirm the conflict-of-interest declaration.
5. Make the GitHub repository public before submission, create a release tag, and preferably archive the release with a persistent DOI.
6. Once the repository is public, change the data-availability sentence from future tense to a stable “openly available” statement and update the software citation with the release identifier.
7. Upload this manuscript as the revised version of arXiv:2501.06292 and ensure the arXiv metadata match the new title and abstract.

## Final assessment

After the seven confirmations above, the manuscript is suitable for submission-level packaging. No further scientific expansion is recommended before peer review.
