# Saved-policy forward decomposition

The one bounded execution completed in 121.964 seconds with zero household or GE solves, no Slurm job and no arrays collected locally. Peak resident memory was 1,860,472 KiB (about 1.77 GiB), above the preparation estimate but below 16 GiB. The 180-second cap was respected.

Baseline reconstruction reproduced the retained PRE array exactly, including its SHA-256. Credit reconstruction re-nested into its saved POST distribution with L1 error 1.14e-15. Feasibility projection was zero, public parameters remained unchanged, every retained native birth-order flow replayed, and age 40–44 uniform-clock childlessness and topcode-adjusted fertility reconciled. Grouped flows add to the saved scalar decomposition within floating-point precision. `forward_result.json` and `forward_completion_receipt.json` retain the full checks.

At identical baseline states, expanded-credit policy increases raw birth flows at every fecund age. The composition term is essentially zero at age18 and negative at every later fecund age. This is a policy-first accounting identity, not a decomposition of direct utility and continuation utility.

| Age cell starts at | Policy change | Composition change | Total change |
|---|---:|---:|---:|
| 18 | +0.942 | -0.000 | +0.942 |
| 22 | +0.857 | -1.530 | -0.673 |
| 26 | +0.656 | -2.283 | -1.627 |
| 30 | +0.452 | -2.139 | -1.688 |
| 34 | +0.257 | -1.599 | -1.341 |
| 38 | +0.131 | -1.147 | -1.016 |
| 42 | +0.061 | -0.715 | -0.654 |

All flow changes in the tables are raw birth events per 1,000 normalized households per four-year period.

| Inherited net financial assets | Policy change | Composition change | Total change |
|---|---:|---:|---:|
| Nonnegative | +3.283 | -29.937 | -26.654 |
| Negative | +0.072 | +20.524 | +20.596 |

Negative net financial assets include net mortgage positions; this table does not identify unsecured borrowing separately. The share of all-age PRE households with negative net financial assets rises from 25.9801% to 45.4221%. Membership changes endogenously, so the table does not establish a causal effect of debt on fertility.

| Children ever born before the period | Children at home before the period | Policy change | Composition change | Total change |
|---|---|---:|---:|---:|
| 0 | 0 | +2.776 | -5.737 | -2.961 |
| 1 | 0 | +0.222 | -0.996 | -0.774 |
| 1 | 1 | +0.224 | -1.678 | -1.454 |
| 2 | 0 | +0.042 | -0.174 | -0.132 |
| 2 | 1 | +0.057 | -0.485 | -0.428 |
| 2 | 2 | +0.035 | -0.344 | -0.309 |

Birth-order contributions use the exact native operator with risk pools snapshotted before births, fecundability applied and no within-period chaining. Origin family states are masked before that operator, then age and wealth contributions are summed. All saved probabilities, including endpoints, are used; the earlier logit-value exclusions do not apply. No occupied mass is deliberately dropped, and no projection occurs.

Completed fertility remains the native topcode-adjusted birth-children/entry ratio; the third flow above counts entry into the literal3+ state. The childlessness moment is the uniform-birth-time age 40–44 projection, not the all-age childless share. The original support-limited native-credit interpretation remains: this does not certify all unoccupied alternatives or grid convergence. The reverse decomposition ordering was not computed.

Source freeze and exact numerical/mathematical proposal are retained in `FORWARD_REVIEW.md`. The original failed analytical pass and its22 zero-weight sum discrepancies remain preserved separately.
