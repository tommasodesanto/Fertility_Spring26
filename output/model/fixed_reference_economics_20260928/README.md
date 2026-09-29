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

The author-requested [Claude visual storyboard and prototypes](slide_inputs_v1/claude_visuals_v1/README.md)
now include occupied-state contributions to additional first births and a
matched-grid comparison of fixed-price cohort outcomes with the credit GE
endpoint. The lead independently checked the plotted values and actual images;
use reviewed_v5, with plain-language captions and unchanged numbers. No new model solves ran.
Further economic discussion prioritizes occupied policy maps, timing versus
completed family size, and heterogeneous responses before selecting new shocks.

The [three-price elasticity comparison](elasticity_v1/recovery_v1/README.md)
is complete and authenticated. Recovery **18815133** used exactly two new solves
and finished in 17 minutes 37 seconds within the approved 35-minute limit.
Immediate-birth price elasticities are **−0.437** with baseline borrowing limits
and **−0.469** with lifetime repayment only; completed-cohort-fertility slopes
are **−0.536** and **−0.555**. Removing artificial limits raises fertility levels
but does not attenuate the local price response at the common reference price.
House price and mapped rent move together with preferences fixed. These are
prescribed-price responses, separate from the GE endpoints below.
The lead checked the [supplemental figure](slide_inputs_v1/recovery_render_v2/actual_output/price_response.pdf)
and [compact table](slide_inputs_v1/recovery_render_v2/actual_output/local_elasticities.csv).
Full 14-row fits, 31 parameters, 17 standard plots and source/checkpoint evidence
are linked in the readout. The monitor is paused; ±2% robustness remains uncomputed.
The frozen calibration and completed supply/borrowing GE evidence are preserved.

The [supplemental housing-supply figure](slide_inputs_v1/rendered_output_v2/supplemental_housing_supply.pdf)
and [borrowing table](slide_inputs_v1/rendered_output_v2/borrowing_comparison.tex)
are ready for slide integration. Their [full-precision CSV](slide_inputs_v1/rendered_output_v2/borrowing_comparison.csv)
and [input hashes](slide_inputs_v1/rendered_output_v2/manifest.json) accompany them.
These are 2007 stationary and prescribed-price comparisons; none is a 2023 result.

The [supply replay](supply_v1/README.md) passed on Torch (job 18801003, zero
lifecycle solves). Removing artificial borrowing limits raises stationary
household population by **2.434% with fixed physical housing**, versus **5.000%
with elastic supply**. Prices and mapped rents rise **4.007% in both endpoints**.
Renewal determines their common price, while supply determines the population
accommodated. A pure 10% supply-intercept expansion instead raises population
10% at unchanged reference prices and household behavior. Native population,
birth-entry queue, housing, PAYGO and estate checks pass within the existing
tolerances; the provisional estate counterparty and grid-convergence caveats
remain. Full 14-row fit and 31-row parameter tables are linked in that README.

The [closed stationary credit GE](credit_ge_v1/README.md) is now complete and
exactly repeated: household population +5.000%, house prices and implied rents
+4.007%, with child benefit and all other economic primitives fixed. Replacement
fertility returns through endogenous prices, not benefit normalization. The
full fit/parameter tables and 17 standard plots are retained. This is an endpoint;
the borrowing transition and matched local elasticity comparison remain separate.

Current borrowing implementation and bounded Torch checks: [credit_v1](credit_v1/README.md).
The first matched-grid credit comparison is complete: immediate births +6.04%,
first births +11.59%, completed fertility 2.1008 to 2.1482, and mean first-birth
age 25.927 to 25.334. Repayment is retained and all preferences, including psi,
are fixed. These are fixed-price results, not equilibrium or transition results.
The compact readout links the full tables and unchanged 17-plot packets.

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
step-size checks remain uncomputed.

September 29 follow-through: the author authorized completing the borrowing
and housing-supply comparisons for a small set of slide figures and tables.
The new `elasticity_v1/` preparation owns price factors 0.98, 0.99, 1, 1.01
and 1.02 under both credit regimes. Zero-solve preflight **18800762** passed
support and repayment recurrence checks but exposed a numerical problem:
the fixed 262-node grid tightens saving floors by up to 0.449 away from the
reference price, compared with 0.001601 on price-specific boundary grids.
The proposed ten-case solve was therefore **not launched**. Its sources and
receipt are retained.

The revised design uses one common union of the boundary grids at all five
prices. Zero-solve preflight 18801218 passed: the common grid has 602 nodes,
preserves inherited and entrant atoms, and adds no price-dependent floor
tightening relative to the candidate-specific grids. Main job **18801318**
and continuation **18803216** have terminated. Four q0 controls and both −1%
price cases passed, but the reference +1% case hit its ten-minute cap during
post-processing. The symmetric elasticity comparison is incomplete; see the
[failure status, plan and original deadline](elasticity_v1/README.md).
The revised pre-launch budget is twelve solves: two reference-credit controls,
two solvency-only controls (one exact repeat of each), then eight price shocks;
600 seconds per case and 5400 seconds total, one CPU and 16 GiB on Torch.
The first controls must also establish a feasible time estimate for the
remaining cases. This is a new numerical design before production launch,
not an extension or restart of a completed search. A failed gate stops the
run, without automatic retries. The common inherited occupied distribution
is fixed for impact comparisons. Cohort completed fertility is reported
separately. Price and implied rent move together; these are prescribed-price
responses, not an equilibrium transition.

The completed `supply_v1/` replay verified that supply enters clearing and
reporting, but not household optimization at given prices. Each supply case
changes its explicitly named supply object while retaining the corresponding
household preferences, credit rule and fiscal inputs. It replays actual birth
renewal, the 16/20 entry queue, scaled housing clearing and PAYGO. These
stationary endpoints do not establish a transition path. Estimated 2023
comparisons remain with the transition chat until its transition is available.

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

### Exact baseline credit rules (September 29 source check)

Let $b$ denote net financial assets (negative means debt), $b'$ next-period
assets, $pH$ a home's value, $y$ current after-tax period income, and
$R=1.08243216$ the four-year gross return. These are the frozen baseline rules:

- **Renters:** $b'\geq s_{j+1}\min(b,0)$. The new unsecured credit line is
  zero (`lambda_d=0`); inherited debt can roll over. The next-age multiplier
  equals one through age 42, then 0.8, 0.6, 0.4, 0.2 at ages 46–58, and zero
  from age 62. The multiplier applies to then-current debt each period.
- **Buyers:** $b'\geq-0.8pH$. Purchase entry additionally requires
  $b+S+y/R\geq0.2pH$, where $S$ is sale proceeds from an old home after the
  6% selling cost (zero for an initial renter). Thus current income counts
  toward the down payment. The budget, consumption feasibility and grid
  support still apply; the entry inequality alone is not sufficient.
- **Owners keeping their home:** $b'\geq\min(b,-0.8pH)$. They can borrow up to
  80% of current home value or retain existing greater debt without increasing
  its principal. Interest remains payable. Where death is possible, the
  separate bound $b'\geq-0.94pH$ protects net estates. There is no scheduled
  owner amortization in this reference.

The credit experiment removes all three artificial limits, including purchase
entry, while preserving lifetime repayment and nonnegative estates. It does
not isolate the down-payment channel. All cases retain finite-grid limitations.

Read-only authentication matched local and frozen Torch solver, parameters and
kernels to the original source manifest: respective SHA256
`b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1`,
`66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464`,
`639c9a21797dbc9f2a0e9a891f283c115353c2edfcb89c959a7fe9f32b86ca27`.
Decisive source: `code/model/intergen_eqscale_seq_optimized/solver.py`
lines 127–187, 2697–2713 and 3410–3417; `parameters.py` lines 718–785;
`kernels.py` purchase checks and owner floor at lines 1273–1275.
Serialized settings come from the immutable reference manifest above. No
model solve or economic change was made for this clarification.

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

The verified supply and credit stationary endpoints retain fixed preferences.
Endogenous prices satisfy renewal, actual PAYGO balances, and housing clearing
determines household population. The native operator and queue checks are now
complete; uniqueness, stability and the adjustment path remain unestablished.
Replacement is an endogenous stationary accounting condition, not permission
to reset `psi_child`. Historical 2007-to-2023 estimation and subsequent policy
transitions require the separate transition chat's new computed path.

All numerical work, imports, tests and rendering belong on Torch. Keep large
checkpoints there, use isolated versioned experiment sources, preserve receipts,
and leave shared active model code and other chats' work untouched.
