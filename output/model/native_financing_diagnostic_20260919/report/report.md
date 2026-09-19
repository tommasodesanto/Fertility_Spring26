# Native financing diagnostic

Fixed-price partial-equilibrium comparison using the saved `g_pre` for every arm.

First births equal the exact loss of `n=0` mass between `g_pre` and `g_post_fertility`, before housing and current-tenure choices. Total births are the native `births` scalar. The report does not infer an adjusted top-code total unless native accounting supplies it. Tenure and location probability checks sum in float64 and allow only the stored array dtype's single machine epsilon in addition to the reporting arithmetic tolerance of 2e-11; they never renormalize probabilities. Fertility normalization is not separately audited.

## Comparison

| Case | Births / initial HH | First births / initial HH | First births age 26–38 / initial HH | Ownership | Mean physical rooms | Saving boundary (realized) |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 0.11559 | 0.0506498 | 0.02013 | 0.566411 | 6.63987 | 0 |
| mortgage_only | 0.116388 | 0.0511798 | 0.0202911 | 0.672661 | 6.75717 | 0 |
| unsecured_only | 0.123774 | 0.0575803 | 0.0222071 | 0.579136 | 6.77715 | 0 |
| both | 0.124516 | 0.0580773 | 0.0223529 | 0.687434 | 6.91302 | 0 |

All grouped birth outcomes use initial states in `g_pre`; current-tenure redistribution is not reclassified as initial tenure. Renter rooms use realized renter mass and `hR_pol`; owner rooms use realized owner mass and `H_own`. Saving-boundary censoring uses realized current mass and both saving-grid boundaries.

Supplementary figures: `native_financing_supplement.png` and `.pdf`.
