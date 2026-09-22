# Scientific sanity checks

**Status: PASS for the checks listed below.**

Exhaustive checks cover every nonempty support at n=2 and n=4 against production routines. Table I checks read the actual retained CSV and its derived outputs; they do not rerun the peak searches or validate an asymptotic law. The residue checks reconstruct all 12,000-row statistics and replay three seeded supports per group (36 total). Formula agreement does not establish the dense-bulk ansatz's approximation accuracy. Novelty remains a separate question.

## Retained Table I consistency

| Check | Result |
| --- | ---: |
| table rows checked | 11 |
| log2 M n slope | 0.703540951388 |
| log2 M n intercept | -0.357734368148 |
| S n slope | 0.5093 |
| S n intercept | -0.990090909091 |
| max derived table or fit error | 3.5527136788e-15 |

## Exhaustive small systems

| Check | Result |
| --- | ---: |
| supports enumerated | 65550 |
| max coefficient matrix error | 0 |
| max entropy error | 1.7763568394e-15 |
| max residue ceiling error | 0 |
| max mean spectrum error | 2.06779038336e-14 |
| max average purity error | 1.11022302463e-15 |
| max diagonal entropy error | 2.68673971959e-14 |
| max entropy lower bound error | 2.44249065418e-15 |
| max qft amplitude error | 2.77555756156e-16 |
| max qft entropy error | 3.33066907388e-16 |
| max entropy bound violation | 3.20342650381e-16 |

## Dense-bulk formula implementation

| Check | Result |
| --- | ---: |
| parameter pairs checked | 16660 |
| max formula error | 2.6645352591e-15 |

## Residue-control evidence

| Check | Result |
| --- | ---: |
| groups checked | 12 |
| raw rows checked | 12000 |
| max summary error | 3.5527136788e-15 |
| seed replayed supports | 36 |
| max seed replay entropy error | 2.6645352591e-15 |
