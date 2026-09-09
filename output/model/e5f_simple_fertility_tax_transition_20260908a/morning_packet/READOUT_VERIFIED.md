# Verified rebated property-tax transition

All 22 full transition dates and four smoke dates passed the recorded numerical gates. Both paths reproduce the verified 2023 impacts and their 2023/2027 smoke solutions. The comparison raises annual property tax from 1% to 2%, with equal household rebates in both cases.

| Year | Births/HH % | Total births % | Households % | Price % | Young ownership pp | Young parent ownership pp | Young parent rooms % |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2023 | 0.5093 | 0.5093 | 0.0000 | -4.8016 | 1.2360 | 0.1065 | -3.4767 |
| 2027 | 0.6076 | 0.6076 | -0.0000 | -6.0979 | 2.1823 | 0.0645 | -4.1287 |
| 2031 | 0.7084 | 0.7084 | 0.0000 | -6.9462 | 2.9659 | 0.3898 | -4.4433 |
| 2035 | 0.7921 | 0.7921 | 0.0000 | -7.5112 | 3.3830 | 1.1429 | -4.4725 |
| 2039 | 0.8586 | 0.8586 | 0.0000 | -7.8698 | 3.2284 | 1.2939 | -4.2853 |
| 2043 | 0.9268 | 0.9497 | 0.0227 | -8.0086 | 2.4025 | 0.6682 | -4.3637 |
| 2047 | 0.9435 | 0.9940 | 0.0500 | -7.9666 | 2.1987 | 0.8964 | -4.3337 |
| 2051 | 0.8913 | 0.9736 | 0.0816 | -7.8329 | 2.3303 | 1.5597 | -4.2301 |
| 2055 | 0.9315 | 1.0490 | 0.1164 | -7.7011 | 2.1447 | 1.5090 | -4.3157 |
| 2059 | 0.9805 | 1.1370 | 0.1550 | -7.6276 | 2.0468 | 1.5331 | -4.2871 |
| 2063 | 1.0497 | 1.2513 | 0.1995 | -7.5911 | 2.1541 | 1.7457 | -4.1838 |

Independent checks verified 324 local artifact hashes and all 11 dated effect rows. The largest recomputation discrepancy is 0. Maximum market residual is 1.50843e-05; absolute fiscal residual 2.42821e-05; mass-accounting residual 1.11022e-15.

Births are top-code adjusted four-year flows. Housing demand includes occupied rental and owner rooms. Population means household decision units, not resident persons. No household-mass difference occurs before 2043, as required by the twenty-year birth-entry lag.

This is a finite temporary-equilibrium diagnostic: M=0, retention 1, births/2.1, unchanged supply rule with elasticity 0.63, and current prices perceived permanent each date. It does not certify a perfect-foresight transition, a terminal equilibrium, or production adoption.

Young-parent group means are conditional on dependent-child status and can change through composition. The household-formation bridge and annual-age alignment remain outstanding. Remote checkpoint hashing/full policy-array comparisons are verified by the completed collector; local checkpoints were deliberately not downloaded.

Evidence: verification.json and ../results/comparison_receipt.json; full dated tables and gate packets under ../results/full/.
