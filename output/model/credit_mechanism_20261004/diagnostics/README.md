# Common-state credit diagnostic

**Verified October 4, 2026.** This packet reads the retained heterogeneous,
fixed-price \(\phi=0.8\) and \(1\) cases. It changes no economic input, target,
solver, manuscript or acceptance gate. The recovered \(0.8\) versus \(0.95\)
scratch tests saved summaries, not paired policy arrays. Their summaries cannot
supply this value decomposition.

**New matched diagnostics:** one \(\phi=0.95\) case and one \(+10\%\) price
case completed serially, reusing the exactly matched cached baseline. Each wrote
all 87 numerical arrays, executed inputs and 17 standard graphs. Neither is GE.
Three-plus adjusted births change −0.242% under credit and −6.703% under price.

| Change | Raw births | Common-policy term | Exposure term | Ownership, ages 18–29 | Rooms, ages 18–29 |
|---|---:|---:|---:|---:|---:|
| \(\phi:0.8\to0.95\) | −0.250% | +0.156% | −0.406% | +10.112 pp | −0.548% |
| Price \(+10\%\), \(\phi=0.8\) | −6.434% | −7.417% | +0.983% | −3.818 pp | −8.605% |

The successful-birth-minus-wait gap, weighted by baseline local susceptibility
\(G\pi_j^2a(1-a)/\kappa\), changes by +0.000602 utility units under credit and
−0.067731 under price. Credit has full susceptibility support; price's unsupported
weight is \(1.4\times10^{-41}\), disclosed. This confirms a strong housing-price
channel alongside weak mortgage-credit sensitivity without changing preferences
or the income gradient. Price jointly changes rents, purchase costs, collateral
and inherited owner housing wealth. It does not isolate the family housing floor.
The price case lowers first/second/third birth flows 5.015%/6.753%/8.861%.
Complete tables/plots are in `phi095_decomposition/` and `price110_decomposition/`;
small ownership/room/gap tables are `phi095_small_table.json` and
`price110_small_table.json`.

## Findings

The aggregate credit response combines a small policy response with a larger
change in the households exposed to each policy. “Prospective parents are rich”
does not by itself explain this result: low-income, low-wealth households also
receive less credit value on the successful-birth branch than on the wait branch.
Correcting the income gradient may matter, but these comparisons do not establish
that it would create a strong fertility-credit interaction.

Here a birth flow is the mass attempting an upward birth times its conception
probability, before current tenure choice. These are household flows in a
normalized stationary cross-section, not national births or a dated transition.

| Flow accounting | Baseline | Relaxed | Change | Policy at baseline exposure | Composition |
|---|---:|---:|---:|---:|---:|
| Explicit first, second and third births | 0.115288 | 0.114280 | −0.874% | +0.166% | −1.040% |
| Three-plus adjusted child units | 0.129640 | 0.128521 | −0.863% | +0.171% | −1.034% |

Percent contributions divide by the corresponding baseline flow. The saved
three-plus weight is \(w_3=3.602359422009\); the exact implemented correction is
\(B^{adj}=B+(w_3-3)B_3\), from `production/engine/adult_entry.py:14`.
First-, second- and third-birth raw flows decline 0.891%, 0.913%, and 0.772%.
The symmetric Shapley decomposition attributes +0.190% to policy and −1.065%
to composition; the conclusion does not depend on the one-sided allocation.

Beginning debt exposure among households below the top child-count state grows
from 0.112887 to 0.178808 across fertile ages. Debt-cell composition adds
0.008351 births, while nonnegative-asset composition removes 0.009550. Beginning
renter composition contributes −0.009474, offset principally by small-owner
states. This establishes movement across debt/tenure/wealth states; it does not
uniquely assign the effect to borrowing rather than endogenous tenure, saving,
birth history, and current-income transitions. Age 18 has identical entry
composition and a negative policy effect. From age 26 onward the policy effect
is positive, but composition remains negative.

The first-birth table uses common baseline prebirth weights. Income groups are
pooled fertile-age childless-exposure terciles of current annual gross earnings,
with boundaries 0.338391 and 0.732667 in original gross-earnings units. They
differ from the scratch comparison's fixed productivity-state groups and CPS
family money-income bins. Beginning assets retain the same normalization.

| Income group; assets | Baseline birth rate | Common-policy response, pp | Wait credit gain | Successful-birth credit gain | Change in success-minus-wait gap | Gap-support share |
|---|---:|---:|---:|---:|---:|---:|
| Low; \(b\leq0\) | 0.065% | −0.000196 | 0.065652 | 0.025974 | −0.039679 | 92.509% |
| Middle; \(b\leq0\) | 14.032% | −0.030215 | 0.014257 | 0.013440 | −0.000816 | 100% |
| High; \(b\leq0\) | 59.792% | +0.138349 | 0.004967 | 0.006727 | +0.001759 | 100% |

Value gains are cardinal lifetime utility, **not earnings equivalents**. The
lowest-income success-minus-wait gap becomes more negative despite positive
credit value in both branches. Literal zero attempt probabilities prevent exact
logit inversion on 7.491% of that group's exposure, or 2.199% of all first-birth
exposure. All such mass remains in the birth-rate and flow decompositions;
conditional value means explicitly report their support. No aggregate consumption
or saving claim is made for the \(\phi=1\) arm, whose native reporter fails.

## Exact measurement and checks

The saved beginning distribution is postbirth and pre-tenure. For current
children-ever-born \(n\), children-at-home \(m\), and state \(s=(b,t,i,j,z)\),
the native independent-count map is inverted as
\[
 G_{nm}(s)=\frac{g^{post}_{nm}(s)-G_{n-1,m-1}(s)h_{n-1,m-1}(s)}{1-h_{nm}(s)},
 \qquad h=\pi_j a.
\]
Absent inflows are zero; the top count has no upward attempt. Native forward
code snapshots all at-risk pools before any births, preventing two births in a
period (`production/engine/distribution.py:711`). Reconstructed native flows
match within \(4\times10^{-17}\), the postbirth map within \(2.2\times10^{-19}\).
Signed negative roundoff totals approximately \(1.1\times10^{-17}\), retained
and disclosed. Newly occupied \(\phi=1\) debt states with exactly zero baseline
exposure carry 0.041521 mass and remain in composition.

At each supported state, \(\Delta^{try}=\kappa\log(a/(1-a))\) is the expected
try-minus-wait value gap, **not** success-minus-wait. Inclusive value \(I\)
gives \(V^{wait}=I+\kappa\log(1-a)\). With unchanged conception probability
and birth cost, credit changes the successful-birth relative value by
\(d\Delta^{try}/\pi_j\); its absolute gain equals
\(dV^{wait}+d\Delta^{try}/\pi_j\).
These identities follow `production/engine/household.py:972` and are checked.

The soft purchase screen is \(Rb+y+S\geq(1-\phi)Ph\), with net sale receipts
\(S=0\) for renters. Current after-tax income enters the purchase test. The
owner saving floor remains \(-\phi Ph\), with separate stayer/death constraints.
See `household.py:894` and `kernels.py:311`; \(\phi=1\) removes the nominal
down payment but does not make ownership costless. Gross earnings follow
`shared.py:175`: \(y z/[4(1-\tau_{pay})]\), not raw budget income.

## Reproduction and bounded follow-up

Run `extract_common_states.py` with the bundled Python and all BLAS/OpenMP/Numba
thread caps set to one. It writes `flow_decomposition.csv`,
`age_debt_tenure_decomposition.csv`, `common_state_groups.csv`, the compact
income/wealth CSV, two **supplemental** PNGs, and hash/check receipt. All 34
retained standard graphs are linked in `receipt.json`; their locations and
validity are indexed by the original credit-distribution manifest.

`paired_phi095.py` and `supervise_phi095.py` enforce the reviewed contracts.
The first strict inventory failure and macOS address-space-limit failure are
retained; the latter occurred before any solve. The lead approved exact reached
runtime identity, excluding only unused calibration/workflow/parameter-file
frontends and three tests, individually named with reasons in the receipts.
Every other archived hash, bound input, grid and baseline price matched.
External `ps` polling at 0.25 seconds enforced an 8 GiB **RSS** cap and 600-second
deadline; each complete case took under 13 seconds and peaked near 1.54 GiB.
Before/after imports pin 21 production modules and record package versions.
Both new cases pass native distribution/probability and aggregate-policy checks.
The dated-budget/purchase audit is **unexecuted**. An independent absolute-path
attempt identified the missing frozen-authentication source:
`calibration_archive/model_legacy_20261003/intergen_housing_fertility_howard_test/__init__.py`,
already deleted before this investigation. These are initialization failures,
not measured budget violations. [Audit receipts](../audit/README.md). No native
gate was relaxed, and no parent-targeted loan or strict-cash case was run.

The zero-solve [joint-budget check](joint_budget/README.md) verifies simultaneous
consumption and saving within housing alternatives. Its age-22 renter extension
audits every positive-probability owner option over 249 susceptible states,
covering 21.5% of all first-birth susceptibility. Current borrowing relief
improves waiting more than successful birth on average in this slice; this
does not replace the full dated-budget audit or identify the causes of the
permanent policy response.
