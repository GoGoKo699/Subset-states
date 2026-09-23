# Scientific sanity checks

**Status: PASS for the checks listed below.**

Exhaustive checks cover every nonempty support at $n=2$ and $n=4$ against production routines. Table I checks read the actual retained CSV and its derived outputs; they do not rerun the peak searches or validate an asymptotic law. The residue checks reconstruct all 12,000-row statistics and replay three seeded supports per group (36 total). Formula agreement does not establish the dense-bulk ansatz's approximation accuracy. Novelty remains a separate question.

## Retained Table I consistency

| Check | Result |
| --- | ---: |
| table rows checked | 11 |
| $\log_2 M_n$ slope | 0.703540951388 |
| $\log_2 M_n$ intercept | -0.357734368148 |
| $S_n$ slope | 0.5093 |
| $S_n$ intercept | -0.990090909091 |
| max derived table or fit error | $3.5527136788\times10^{-15}$ |

## Exhaustive small systems

| Check | Result |
| --- | ---: |
| supports enumerated | 65550 |
| max coefficient matrix error | 0 |
| max entropy error | $1.7763568394\times10^{-15}$ |
| max residue ceiling error | 0 |
| max mean spectrum error | $2.06779038336\times10^{-14}$ |
| max average purity error | $1.11022302463\times10^{-15}$ |
| max diagonal entropy error | $2.68673971959\times10^{-14}$ |
| max entropy lower bound error | $2.44249065418\times10^{-15}$ |
| max qft amplitude error | $2.77555756156\times10^{-16}$ |
| max qft entropy error | $3.33066907388\times10^{-16}$ |
| max entropy bound violation | $3.20342650381\times10^{-16}$ |

## Dense-bulk formula implementation

| Check | Result |
| --- | ---: |
| parameter pairs checked | 16660 |
| max formula error | $2.6645352591\times10^{-15}$ |

## Residue-control evidence

| Check | Result |
| --- | ---: |
| groups checked | 12 |
| raw rows checked | 12000 |
| max summary error | $3.5527136788\times10^{-15}$ |
| seed replayed supports | 36 |
| max seed replay entropy error | $2.6645352591\times10^{-15}$ |
