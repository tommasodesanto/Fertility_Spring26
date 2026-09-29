# Fixed-calibration economic investigations

**2007 stationary reference — block0506, September 28 verified export**

The author authorized all four proposed investigations on September 28: saved
birth/constraint anatomy, occupied housing/tenure/retirement anatomy, a 10%
fixed-price housing-cost shock, and a 10% housing-supply-intercept expansion.
He also requested a credit experiment with impact, a heuristic/partial
transition and a new steady state, and a transition with fixed physical housing.

Reference identity is permanently pinned by
`../fertility_identification_20260928/fixed_reference_manifest.json`, SHA256
`147f9e2cb20f66350f1ceaa16cb41f822041ec869676ef5d5b9d04f16e4190d4`.
The authoritative primary export remains unchanged. Common-primary loss is
19.581310760138322. All preferences, including
`psi_child=0.1355551166583114`, remain fixed. No target revisions, recalibration,
fertility renormalization, or automatic adoption of another chat's candidate.

## Author decisions and ownership

- **Credit experiment:** remove artificial borrowing and down-payment limits;
  retain lifetime solvency and repayment. Setting LTV to one or choosing a very
  negative arbitrary debt floor alone does not establish this experiment.
- **Fixed-housing transition:** hold physical housing stock/supply fixed; prices
  and rents clear markets. Household housing and tenure choices remain free.
- The separate user-owned chat **Prepare transitions for frozen block0506**,
  id `01a0ea09-e255-7f13-8962-40ec1c8d4c1b`, owns the historical September10–14
  reconstruction and isolated transition implementation, including the
  fixed-physical-stock comparison. Its workspace is
  `../fixed_reference_transition_20260928/`. It receives verified credit rules
  and terminal endpoints from this economic-analysis chat.
- This economic-analysis chat owns credit-rule implementation, fixed-price
  impact, the new steady state, elasticities, economic interpretation and the
  supply-shock comparison. The author explicitly corrected this division;
  the transition chat's withdrawn credit stage launched no numerical jobs.
  Calibration improvement remains in its separate chat.

## Current priorities and useful output

The author prioritizes understanding elasticities and financial constraints.
Choose the economic question and useful output before expanding computation
or reporting. Default delivery is a compact comparison table and a short
mechanism explanation; broad overviews and long assembled reports require a
specific purpose. Retain the underlying standard diagnostics and complete
fit/parameter tables, with links, without assembling them into another PDF.

The next comparison should answer two questions: does removing artificial
credit limits change fertility at given prices, and how does it change the
fertility response to housing costs? Compare the calibrated borrowing rule
with solvency-only borrowing, keeping preferences (including psi), earnings,
entry endowments, fiscal inputs and other primitives fixed. Evaluate impact
responses on identical inherited occupied states. Report first and subsequent
birth responses separately; report cohort completed fertility and first-birth
timing separately from impact flows. Inspect low-financial-wealth households
and inherited tenure to locate the response, without calling a binding-share
correlation causal. Market-clearing endpoints remain a separate stage.

Existing evidence measures a permanent +10% house-price and implied-rent
change. Its log-change elasticities are finite-change responses, not local
derivatives or rent-only elasticities. The paired credit comparison and local
step-size checks remain uncomputed. No new numerical budget or job is launched
by this output specification.

## Saved-state packet

`analyze_saved_state.py` completed a read-only Torch extraction in
`saved_anatomy_v2/` (job **18738188**, 30 seconds):
zero solves; first/subsequent births by occupied age/wealth/income/tenure;
realized housing choices; conditional-policy/occupied-mass comparisons;
net financial wealth, gross housing value and net worth through retirement.
The checkpoint does not identify a separate mortgage balance or housing equity.
Retain all
17 standard diagnostic plots. New figures are supplemental. The existing
authenticated `../fertility_identification_20260928/measurement_audit_v1/`
owns the empirical motherhood/conditional-child-count decomposition and B15
income verification; reuse those results instead of rebuilding them.

The first extraction (18737458) failed an unchanged `2e-10` accounting gate
because it normalized stored float32 choice probabilities before casting to
float64. Version 2 matches the native realization convention exactly; its
independent transaction replay agrees to machine precision. The failed remote
source/output remains preserved.

`audit_saved_constraints.py` produced `constraints_v1/` (job **18738587**,
32 seconds, zero solves; five-minute cap, one CPU, 24 GiB). It measures occupied
native saving-floor boundaries using separate buyer/renter and owner-stayer
policies. The within-branch shares are 25.112% for renters, 4.722% for buyers,
and 4.538% for owner stayers. Purchase exclusion incidence is unmeasured;
these shares do not establish whether down payments are quantitatively weak.

## Fixed-price experiment contract

Executed budget: at most three lifecycle solves, sequentially on
Torch with one computational thread and 16 GiB. Two exact reference-price
controls must pass before the single shock. Maximum 20 minutes including
preparation/reporting, maximum six minutes per lifecycle case, no retries or
deadline extensions. Existing observed selected solve time is about 133 seconds;
three solves therefore imply roughly seven minutes before audit/report overhead.
If a control fails, stop before the shock and retain the failure evidence.
The driver and plan were pinned before dispatch. Job **18738157** completed
exactly three solves in **376.76 seconds**, in `fixed_price_v1/`.
Executed driver SHA256:
`96d6923a252f57bc4d8c44fd6479b13f48ba217d74edf8ef629d120428b03b44`.
Its immutable bytes and plan remain on Torch in
`sources/fixed_price_v1/`. Each of the two controls matched all **113** numeric
arrays, 14 fit rows, 31 parameter values and 17 PNG hashes exactly. The single
shock passed the same applicable household/accounting checks.

The sole economic change is a permanent 10% increase in the housing asset
price and its implied rent under the unchanged user-cost mapping. Hold
preferences, earnings, interest rate, taxes, pension benefit, entry distribution,
survival, timing, housing menus, credit rules and supply primitives fixed.
Reoptimize household choices. Evaluate immediate behavioral differences on
the exact baseline occupied pre-choice distribution. A separately reported
normalized-cohort distribution is a partial-equilibrium diagnostic, not a
renewed demographic steady state or an equilibrium transition. Report market,
fiscal, estate and replacement residuals rather than normalizing them away.
Pass the original applicable household budget, value, probability, feasibility,
purchase accounting and operator gates; market clearing is not imposed on the
explicitly prescribed-price experiment. Keep complete 14-row target-fit and
31-row parameter tables and the standard 17 plots for each reported solution.

| Object | Reference | Immediate response | Normalized cohort |
|---|---:|---:|---:|
| Births per household per four-year period | 0.115419 | 0.110588 | 0.109497 |
| First births per household per period | 0.050164 | 0.045956 | 0.048126 |
| Ownership, all households | 0.668165 | 0.658504 | 0.646425 |
| Rooms per household | 5.847942 | 5.664364 | 5.453161 |
| Nonhousing consumption per household | 2.341419 | 2.383004 | 2.317723 |
| Completed fertility | 2.099998 | Not a one-date object | 1.988120 |

Immediate births fall 4.186%, with 87.10% of the decline from first births.
The normalized cohort has a 5.328% adult-entry replacement gap. Housing excess
supply is 0.545477 on impact and 0.756680 in the cohort calculation. These are
prescribed-price mechanisms, not equilibrium or historical-transition results.
`comparison.json`, all case receipts and full tables are in `fixed_price_v1/`.

Post-run authentication confirms the checkpoint and all 30 authenticated export
files are unchanged. Nine new below-tolerance inherited-state logs followed
the saved reference logging path. Only those current-job files were identified
by their job window, price and case receipt, then moved to their respective
new case folders with unchanged hashes. The relocation and post-run reference
authentication receipts are in `fixed_price_v1/`. The current local driver
includes an output-only logging-path correction for future use; it is not the
immutable version-1 driver that produced these results. No completed case was
restarted, and no feasibility tolerance or population projection was changed.

## Readout and regeneration

The reviewed report is `../../pdf/fixed_reference_economics_block0506.pdf`:
30 pages, full 14-row fit/31-row parameter tables for reference and shock,
the unchanged 17 standard plots for each, and five supplemental views.
PDF SHA256: `71a857c8415ca9771e3db0f5d8a771cb148c43d4caea55c8048240c5d4be7a46`.
Torch rendering job **18739272** took eight seconds; all 30 pages were visually
inspected. Existing standard-plot layouts and their dense legends are retained.

`build_readout.py --findings findings.json --output NEW_REPORT.pdf` assembles
the complete report inside the authenticated Torch mount, with no solves.
The exact rendering sources and findings are pinned in `sources/readout_v1/`.
To regenerate the standard graph set from an authenticated loaded case packet,
use `run_e5f_independent_numerical_audit.standard_diagnostics(packet, new_output,
validate_production_young=False)` in the frozen runtime. Do not replace the
actual saved parameters with the ancestry setup's parameters.

## Remaining equilibrium work

The 10% supply-intercept experiment preserves preferences and other primitives
and endogenizes prices. Population/adult-entry, fiscal and estate closure must
be reconciled with the new transition machinery before claiming a cleared
demographic equilibrium. The calibrated reference uses normalized population;
continuing that numerical normalization after fertility changes does not
establish demographic renewal. Historical 2007-to-2023 transition estimation
and new policy transitions from the 2007 reference are separate objects.

The separate transition chat's preparation note now writes the intended closed
endpoint conditions explicitly: endogenous prices satisfy renewal at fixed
preferences, PAYGO determines the pension, and housing clearing determines
population. Replacement in a genuine closed stationary endpoint is an
endogenous accounting condition, not permission to reset `psi_child`.
For a pure 10% supply-intercept increase, scaling population and all aggregate
flows by 1.1 at unchanged prices/policies is an algebraic candidate under the
present level-linear closure. Native operator/scaling verification, uniqueness,
stability and the transition remain unestablished; no supply result is claimed.

All numerical work, imports, tests and rendering belong on Torch. Keep large
checkpoints there, use isolated versioned experiment sources, preserve receipts,
and leave shared active model code and other chats' work untouched.
