# Quarter-rule overlap at the selected 80% fit

This is a bounded household diagnostic for the freshly postchecked quarter-saving chain 54, loss 51.5560360491. It uses the exact saved pre-fertility distribution, attempt probabilities, age fecundity, value functions, prices, and selected parameters. The selected runtime authenticates with the same two-file frozen read-only overlay used by the local postcheck. The diagnostic makes no equilibrium solve, no Bellman recursion, and no calibration proposal.

For each childless renter state, the [native financial-access map](../../buyer_diagnostics/financial_access.py) identifies whether an owner product first becomes financially feasible when the financed share rises from 80% to 100%. Separately, the native renter saving kernel recomputes current-period renter values when the six-room cap is relaxed to 100 rooms, holding the saved next-age values and all prices fixed. It does this for both the childless wait family state and the parent family state reached after a successful first birth. The difference between their cap shadows measures whether the cap particularly penalizes parenting at the *same* starting state. This is a renter-only, one-period value comparison; it is not a full birth-value or policy-effect decomposition.

The weights are pre-fertility childless renter mass times age fecundity times $p(1-p)/\kappa_f$, where $p$ is the saved attempt probability and $\kappa_f=0.1276143394$. These weights are proportional to the local response of first-birth flow to an **attempt-versus-wait value** difference. A change in the value conditional on a *successful birth*, such as a parent housing benefit, enters the attempt value with another factor of age fecundity. Actual first-birth flow uses mass times fecundity times $p$ instead.

| Childless renter states | Share of renter-origin first births | Share of local birth responsiveness | Mean parent-minus-childless cap shadow, responsiveness weighted |
|---|---:|---:|---:|
| Owner access opens at 100% | 7.411% | 14.360% | 0.00000662 |
| Already owner feasible at 80% | 92.589% | 85.640% | 0.01215627 |

Among the newly financially reached states, zero responsiveness weight has a differential cap shadow exceeding $0.1\kappa_f=0.0127614$. In the states already owner feasible at 80%, 38.298% of responsiveness weight exceeds that threshold. The native financial-access share exactly reproduces the previously collected selected v4 result to displayed precision.

With the additional fecundity factor appropriate to a successful-birth branch value, newly reached states account for 14.497% of local response weight; their mean differential cap shadow is 0.00000678, versus 0.01213036 among states already owner feasible at 80%. This confirms that the overlap comparison is insensitive to that age reweighting.

## Who is reached, and what rooms would they rent?

The following mass shares use **all fertile-age childless renters before the fertility choice**. Birth shares use actual first-birth flow; response shares use the local attempt-probability derivative defined above. Income and wealth are in the model's native units, and their means are response weighted. Parent rooms describe the renter branch after a successful first birth, with the same saved continuation values. Counts are positive-response state-grid cells, not households.

| Financial access at the same starting state | Renter mass | Actual renter-origin birth flow | Local response weight | Mean income | Mean wealth | Birth flow at ages 18–22 |
|---|---:|---:|---:|---:|---:|---:|
| Opens at 100% | 50.315% | 7.411% | 14.360% | 1.935 | 0.023 | 73.820% |
| Already feasible at 80% | 47.326% | 92.589% | 85.640% | 3.586 | 0.514 | 60.008% |
| Still infeasible at 100% | 2.359% | 0% | 0% | — | — | — |

| Financial access | Parent rooms, cap 6 | Parent rooms, cap 100 | Positive-response cells at cap 6 | Response weight at cap 6 | Largest parent cap value gain | Childless wait rooms, cap 6 → 100 |
|---|---:|---:|---:|---:|---:|---:|
| Opens at 100% | 5.692 | 5.704 | 1 of 83 | 23.029% | 0.00002874 | 3.845 → 3.845 |
| Already feasible at 80% | 5.994 | 7.986 | 2,621 of 2,948 | 96.543% | 0.04317 | 5.228 → 5.601 |

For newly reached states, response-weighted wealth has median zero and 90th percentile 0.163; for already feasible states, the corresponding values are 0.302 and 1.140. The four-room owner product becomes financially feasible in every newly reached birth-flow state. Six-, eight- and ten-room owner products also become feasible for 100.000%, 99.9999% and 99.897% of that group's birth flow, respectively; the two-room product is never feasible under the parent housing floor. These are overlapping product-access shares, not purchase choices. The three access groups exhaust fertile renter mass, and the first two exhaust renter-origin births and local response weight to numerical precision.

### Fixed-continuation rental-cap sweep

Each row recomputes the parent and childless-wait renter branches with the indicated room cap, measured against the **same 100-room slack benchmark**. These figures use the successful-birth response weights (including the additional fecundity factor); the saved fertility probabilities, origin distribution, owner financial-access masks, prices and next-age values remain fixed. They are one-period diagnostics, not policy or equilibrium effects.

| Rental cap | Newly reached: parent / wait rooms | Newly reached: parent-minus-wait cap shadow | Already feasible: parent / wait rooms | Already feasible: parent-minus-wait cap shadow |
|---:|---:|---:|---:|---:|
| 6.0 | 5.696 / 3.847 | 0.00000678 | 5.994 / 5.221 | 0.012130 |
| 5.5 | 5.435 / 3.847 | 0.001270 | 5.499 / 5.061 | 0.020524 |
| 5.0 | 4.992 / 3.847 | 0.008411 | 5.000 / 4.811 | 0.033533 |
| 4.5 | 4.499 / 3.847 | 0.026348 | 4.500 / 4.448 | 0.053221 |

At the fitted six-room limit, newly reached parent renter demand is largely below the cap. Lower limits increasingly constrain that group's parent branch while leaving its childless wait rooms near 3.85. The group's extra parent cap penalty consequently grows, although it remains below the already-feasible group's penalty throughout this sweep. This shows that the weak overlap finding is specific to the six-room cap and saved fit; it does not prove that the mechanism is absent under other rental limits.

**Interpretation.** The newly reached renters have lower income and wealth, and their parent renter-room choice barely moves when the six-room limit is relaxed. The already owner-feasible group tends to choose six parent renter rooms and raises that choice by about two rooms when allowed. The diagnostic therefore gives a concrete reason for the weak three-way overlap at this fitted point: the credit change reaches low-resource states while the rental cap's larger parenting penalty lies in higher-resource states. It does not establish whether those newly reached households want to buy, whether price movements offset a direct gain, or how a permanent change to the cap would work.

The native reconstructed childless renter-room policy differs from the saved policy by at most 0.005835 room across the audited fertile-age cells. This checks a policy output but does not bound error in the renter value difference. The script's cap relaxation uses 100 rooms; the largest relaxed choice with positive responsiveness weight is 92.593 parent rooms or 90.768 childless rooms, so 100 is nonbinding in that set. The result should not be treated as a new calibration specification. Full precision, source hashes, weights and group moments are in [result.json](result.json). Reproduce with one core:

```bash
OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_overlap_v1/run.py
```
