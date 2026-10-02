# Quarter selected fit: one-date fixed-price financing diagnostic

**Status (October 2): completed on the authorized one-core local fallback.**
Torch staging was blocked by intermittent DNS resolution and missing Kerberos
authentication (`klist -s` reported no ticket). The exact zero-Bellman
preflight and two-branch smoke passed locally. The production run made exactly
two native Bellman calls in 9.007 seconds; the watchdog recorded peak resident
memory 1.773 GiB, below its 24-GiB limit and 15-minute wall budget.

The native 80% control reproduces every saved fertility attempt probability
exactly (maximum absolute gap zero) and the accepted first-birth flow to
floating-point precision: 0.04974551614189719 versus 0.04974551614189717.
At the same prices, initial distribution and next-age values, today's 100%
financing raises first-birth flow to **0.049837104956647164**, a **0.184115%**
increase. The fertile-childless hazard rises from 0.238028027 to 0.238466271,
or 0.043824 percentage points. Of the absolute flow increase 0.000091589,
0.000083322 is renter-origin and 0.000008267 is owner-origin. The largest
age contribution is at age 18: +0.000096647; age 22 offsets part of it at
−0.000008236. All seven ages and six tenure products are in
[`local_run/production/result.json`](local_run/production/result.json).

The accepted quarter temporary GE date-zero policy flow is
0.0495313221352156, **0.430580% below** its matched control. The fixed-price
response and accepted dated-transition response have opposite signs. To test
how much of the difference is present at the accepted current prices, a second
two-call factorial used the [source-checked date-zero price receipt](../value_decomposition_receipt.json)
(SHA-256 `cca4cb630891aeca7b46cabe8b3d81ea950ae4efa5137b4aa1cbabe40c8a8534`).
It reports the accepted quarter 48-date temporary path's asset-price change
+0.290240368% and **separate** renter unit-price change +2.009844881%.
Applied to authenticated baseline price 0.6744838540900874 and rent
0.12159554040311076, these give 0.6764414785119185 and 0.12403942214723254.
The receipt records SHA-256 digests for the original accepted control and
policy date-zero packets; those packets are not copied into this folder.

| First-birth flow on saved initial households and next-age values | 80% financing | 100% financing |
|---|---:|---:|
| Baseline asset price and rent | 0.04974551614189719 | 0.049837104956647164 |
| Accepted temporary-path date-zero asset price and rent | 0.049443479745558655 | 0.04953224457049744 |

At 100% financing, moving both current prices to the accepted temporary-path
values lowers flow by 0.000304860386. The financing change at baseline prices
raises flow by 0.000091588815. Their sum gives −0.000213271571, compared with
the accepted GE date-zero change −0.000214194007. The remaining
−0.000000922435 is 0.001854 percentage points of baseline flow. This small
remainder includes anticipated future price and fiscal paths, updated
continuation values, and any other dated-transition differences. The factorial
shows that the **joint change in current asset price and rent** nearly
reproduces the accepted date-zero sign reversal with fixed continuation values;
it is not a structural attribution to one price or a general-equilibrium
counterfactual. Full age-by-tenure outputs for all four cells and SHA pins are
in [`local_run/factorial/result.json`](local_run/factorial/result.json).

Two final one-period solves at 100% financing cross the current asset price
and rent while keeping the same saved continuation value and initial
distribution. The exact-loop zero-Bellman smoke passed first. The two native
calls took 9.571 seconds, with peak observed resident memory 1.720 GiB.

| Current asset price, rent | First-birth flow at 100% financing |
|---|---:|
| Baseline asset price, baseline rent | 0.049837104956647164 |
| Accepted asset price, baseline rent | 0.04982435483434676 |
| Baseline asset price, accepted rent | 0.0495489664992902 |
| Accepted asset price, accepted rent | 0.04953224457049744 |

Changing only the asset price at baseline rent lowers flow by 0.000012750122;
changing only rent at baseline asset price lowers it by 0.000288138457.
The two-order Shapley accounting averages each price's marginal effect over
both orders: asset price contributes **−0.000014736026** (4.834% of the joint
price effect), and rent contributes **−0.000290124361** (95.166%). They sum
exactly to the joint current-price effect −0.000304860386. These crossed
prices are household calculations, not market-clearing interventions or
separately identified causal equilibrium price effects. Full precision,
age-by-tenure flows, and source pins are in
[`local_run/shapley/result.json`](local_run/shapley/result.json).

The direct fixed-price birth gain is concentrated at the youngest model age:
the age-18 cell contributes +0.000096647, slightly more than the total
+0.000091589 because older cells partly offset it. At 100% financing and
baseline asset price, raising only rent to the accepted value lowers age-18
birth flow by 0.000099490 and lowers total first-birth flow by 0.000288138.
Of that total rent-only reduction, 96.4% originates with households that were
renters before the fertility choice. These are matched-origin accounting
results; they do not by themselves identify the parent-versus-childless
housing-expenditure difference inside each household's value comparison.

The reference is the fresh postchecked experimental quarter-saving 80% fit,
local chain 54/case0094_nm, loss 51.5560360491. The driver authenticates its
frozen two-file source overlay, selected receipt, price, wealth grid, saved
value function, pre-fertility distribution, attempt probabilities, and shared
arrays. It tests the saved first-birth flow against the accepted matched
control 0.04974551614189717. These are the factual inputs for this diagnostic;
the fit is not an adopted paper baseline.

The only economic change is **experimental**: current-date financed share
`phi` rises from 0.8 to 1.0 in all four family entries, including existing
owner stayers under the native credit rule. Earnings, initial wealth and income
distributions, timing, transfers and floors, preferences, and targets remain at
the selected quarter fit. Asset price, renter unit price, next-age value
function, and current pre-fertility distribution are fixed. The native
one-period Bellman solver recomputes today's full household policy twice, once
at 80% and once at 100%. The 80% call must reproduce the saved attempt
probabilities within 1e-6 at every grid cell and aggregate first-birth flow
within 1e-8 before the 100% result is interpreted. No market, distribution,
or fertility-normalization solve is made.

For each childless state, first-birth flow is current pre-fertility mass times
age fecundity times the attempt probability. The driver reports total flow,
fertile-childless hazard, and age-by-tenure flows and risks. Its comparison
flow for the accepted 48-date temporary GE path is 0.0495313221352156.
The crossed-price calculation shows that the accepted current asset price and
rent almost reproduce the date-zero GE flow with fixed continuation values.
The remaining gap combines future anticipated equilibrium paths and other
transition feedback. This numerical decomposition does not identify separate
causal equilibrium effects of either price.

Under the model's no-arbitrage rule, $r_t=u p_t+(p_t-p_{t+1})$, where $u$ is
the stationary housing user cost. The baseline has $u=0.180279394$. The
accepted policy date-zero asset price and renter unit price above imply a
date-one asset price of 0.674350516, nearly back at the baseline 0.674483854.
Of the date-zero rent increase 0.002443882, the increase in $u p_t$ accounts
for 0.000352919 (14.4%) and the current-minus-next asset-price term accounts
for 0.002090962 (85.6%). This is price-path accounting, not evidence that
renter housing demand increased: renter and owner rooms clear one aggregate
housing-services market. The local receipt does not save both tenure-specific
room quantities at the two dates.

This motivates a duration test without claiming its outcome. At the accepted
date-zero asset price and fixed saved continuation value, the two crossed-rent
cells give a secant first-birth response of -0.119527 per unit of rent. Linear
interpolation places the rent increase that exactly exhausts the direct
credit gain at approximately 0.000660, or 0.54% of baseline rent. The accepted
temporary rent increase is 0.002444, or 2.01%. If the asset price instead
remained flat at its accepted date-zero policy value into the next period,
the no-arbitrage formula would imply a rent increase of only 0.000353 (0.29%)
at that same current asset price. This is a fixed-state arithmetic comparison,
not a forecast: a longer credit policy would change the equilibrium asset
price path, future values, and household decisions, all of which require a
separate accepted transition run.

Reproduce locally with the existing source tree. The `--smoke` pass exercises
the same control and policy branches using the saved fertility array and makes
zero native calls. Use the watchdog launcher for production; it limits the run
to 900 seconds, one numerical thread and 24 GiB observed resident memory.

```bash
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_fixedprice_v1/run.py --preflight --out /private/tmp/quarter_fixedprice_preflight
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_fixedprice_v1/run.py --smoke --out /private/tmp/quarter_fixedprice_smoke
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_fixedprice_v1/launch_local.py
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_fixedprice_v1/run.py --smoke --price-factorial --out /private/tmp/quarter_fixedprice_factorial_smoke
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_fixedprice_v1/launch_local.py --price-factorial
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_fixedprice_v1/run.py --smoke --price-shapley --out /private/tmp/quarter_fixedprice_shapley_smoke
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/independent_diagnosis/quarter_fixedprice_v1/launch_local.py --price-shapley
```

`progress.json` updates after authentication, after the control, and after the
policy; `result.json` records the full comparison. The local production folder
also contains the watchdog status and run log. The authenticated source is the
same selected quarter chain54 frozen overlay; source and selected-array SHA-256
digests are in the result's `source_receipt`.
