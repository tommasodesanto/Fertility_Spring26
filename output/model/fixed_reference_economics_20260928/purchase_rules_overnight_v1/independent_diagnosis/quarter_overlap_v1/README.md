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

**Interpretation.** In this quarter fit, the rental cap matters particularly for parenting mostly where ownership was already financially feasible. Financing newly reaches a different set of households whose parent-versus-childless renter cap penalty is very small in this diagnostic. This weakens the proposed three-way overlap mechanism for the quarter rule at this fitted point. It does not establish whether those newly reached households want to buy, whether price movements offset a direct gain, or how a permanent change to the cap would work.

The native reconstructed childless renter-room policy differs from the saved policy by at most 0.005835 room across the audited fertile-age cells. This checks a policy output but does not bound error in the renter value difference. The script's cap relaxation uses 100 rooms; the largest relaxed choice with positive responsiveness weight is 92.593 parent rooms or 90.768 childless rooms, so 100 is nonbinding in that set. The result should not be treated as a new calibration specification. Full precision, source hashes, weights and group moments are in [result.json](result.json). Reproduce with one core:

```bash
OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_overlap_v1/run.py
```
