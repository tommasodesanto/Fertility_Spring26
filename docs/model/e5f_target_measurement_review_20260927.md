# Measurement review: first-birth housing and bequests

27 September 2026. Read-only source/receipt review; no model solves, empirical reruns, target changes, weighting changes, or Google document edits. The current planning context supplied by the coordinating lead is 14 rows (10 scored, three validation, one normalization), with 11 fitted parameters. Historical frozen-contract classifications below are evidence about measurement, not a replacement for that new contract.

## Decision

Both targets have an acknowledged model-to-data approximation. Neither review establishes a new numerical units bug that requires discarding the author-approved first pass. Neither approximation should be described as exact estimator matching in the next calibration report. Two concrete provenance corrections should be included in the next versioned contract: the first-birth builder and its calendar-year clock; and the existence of the bequest builder. Keep historical receipts unchanged.

The first-birth gap can be investigated with a small additional observer calculation using existing policies. Exact reproduction of its empirical estimator requires more than renaming the current observer. The bequest gap cannot be fixed uniquely by an observer-only scaling: the relevant marital/recipient mapping is absent from the model and is itself partly mechanical in the empirical construction.

## 1. First birth: 1.465 additional rooms

### What the data estimate

The authoritative A2h receipt reports **1.465293280235685 rooms**, with clustered SE **0.0502703781254777**; the author-selected calibration value is **1.465**. It contains 117,853 estimation observations, 9,310 person clusters, 3,302 treated individuals and 6,008 control individuals. The baseline mean is 4.687270020033742 rooms. The count is of individual observations; A2h is not the alternative one-woman-per-household-year design.

The executable builder is `code/data/psid_followup_mar2026/sa_rooms_first_birth_v2.do`, arm `A2h`. It uses current adults age 18+, positive longitudinal individual weights, observed age/education/rooms, biological first-birth history, and confirmed-childless controls. Treated adults must have been reference person or spouse in the baseline window; all eligible controls are retained. Births before the first adult observation are excluded. Official year-specific rooms are matched to the interview year with the stated non-room codes excluded.

The estimator is a Sun–Abraham interacted-cohort event study, with individual and calendar-year fixed effects, age and education covariates, individual probability weights, and standard errors clustered by person ID. Its event clock is explicitly `K = year - f_c_y`. With baseline −3/−2, the headline `Wp3` is **calendar years +3/+4 since first birth**, relative to **calendar years −3/−2**. It is not the third and fourth post-birth interviews. The frozen `final_readout/target_provenance.json` incorrectly says “interviews” and names the older `sa_rooms_first_birth_household_aligned_v1.do` builder. The reported estimate itself agrees with the A2h receipt.

Sources: builder lines 97–114, 152–196; `code/data/psid_followup_mar2026/output/sa_rooms_first_birth_v2/A2h/fit_receipt.csv` and `run_receipt.json`; frozen provenance row `first_birth_rooms`.

### What the model observes

`code/model/tools/e5f_initial_housing_observer.py`, `_stationary_birth_diagnostic`, calls the dated first-birth branch functions in `run_e5f_transition_calibration.py`. These select successful first births from the pre-fertility risk set, create equal-mass treated/control branches, and advance both one **four-year model period**. They use the saving, income, tenure, location, survival and child-aging laws, without the aggregate Census age bridge. The treated branch can have continuation births; the control remains childless. The observed response is the destination mean-housing difference between these branches, using uncapped rooms.

This is a matched-state model contrast. It does not reproduce the empirical pre-birth housing path, the Sun–Abraham cohort/event weights, or the baseline head/spouse selection. The state space has no household relationship-status counterpart. A four-year horizon is reasonably close to the empirical +3/+4 calendar-year window, but does not by itself reproduce the −3/−2 reference normalization or the estimator. The observer already records the non-flat prepath and weighting caveats. A non-flat empirical prepath is a substantive interpretive limitation, not something the observer can remove mechanically.

### Bounded next check, without changing the objective

Using one retained solution, report treated and control housing separately at origin and destination, their masses, and the implied change-in-gap as a **supplemental diagnostic** beside the existing destination gap. This can show how much contemporaneous housing adjustment or anticipation the current contrast embeds. It must not be relabeled the exact empirical −3/−2 to +3/+4 event-study coefficient. Full matching would require a defined synthetic panel/event clock and the empirical regression applied to it, with explicit treatment of cohort weights and unavailable relationship status. Do not silently substitute that observer, move the horizon, or change 1.465.

## 2. Bequests: annual child-directed flow / aggregate wealth = 0.00729102347

### What the data estimate

The adopted 2007 SCF receipt divides **$674.489508 billion** in annual child-directed flow by **$92.509578 trillion** in aggregate signed net worth, both in 2022 dollars. The ratio is **0.007291023472616158**, or **0.7291023473% annually**. It weights positive `NETWORTH - TRUSTS` by head mortality, uses a 25% child share for legally married heads and 75% otherwise, and maps top-coded age 95+ to mortality at 95. There are 4,417 primary economic units, with five implicates each and weights divided by five.

This is a mortality-weighted proxy, not observed annual inheritance receipts. The adoption receipt explicitly defers ownership allocation when one spouse dies, child eligibility/share refinement, and the oldest-age mapping. The available-data calculation is marked complete. In particular, SCF current-roster children do not identify lifetime offspring, so the child shares are applied mechanically, including units whose eventual child recipients are unknown. These acknowledged limitations are not newly pending empirical work.

The frozen provenance says no executable builder was retained. A builder now exists at `code/data/scf/build_bequest_flow_2007.py`; its source implements the scenario battery. This review did not execute it or verify a source fingerprint against the historical receipt. The next contract should cite and pin the actual builder, retaining the adoption receipt as the authority for the selected row. The builder's generic diagnostic header does not overturn the subsequent explicit adoption.

Sources: `output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/bequest_flow_2007/calculation/verified_receipt.json`, its `published_method_reconciliation.json`, the SCF builder, and frozen provenance row `bequest_wealth`.

### What the model observes

`_wealth_diagnostics` in `e5f_initial_housing_observer.py` delegates to `add_aggregate_wealth_bequest_flow_moments` in `code/model/intergen_eqscale_seq_optimized/solver.py`. The numerator sums positive post-saving estates over **all** child-history/dependent-child states and weights them by model death probabilities. It annualizes once by dividing by the four-year period length. The denominator is signed beginning-period net worth for all living households. The native DUE branch separately incorporates stayer saving; that accounting adjustment does not allocate estates to recipients.

With the estate receiver inactive, housing enters at gross value. The empirical ratio instead excludes trusts from the estate base and allocates only a fraction to children. The model has no marital or spouse-ownership state and no trust asset. Household mortality, forced terminal death and the model's 18–85 age coverage also differ from the SCF head-mortality construction with a 95+ top code. The shared annual-flow/stock units are coherent; the recipient and population definitions differ. The estate-financing ledger that funds entrants and burns the remainder is a separate accounting object and does not resolve this mismatch.

The target has no sampling SE in this contract. Its inherited synthetic 5% calibration scale is a weighting convention, not empirical precision. Do not interpret a large weighted contribution as precise evidence for the child-recipient mapping.

### Bounded next check, without changing the objective

From one saved solution, decompose gross positive death estates by children ever born and age, and report a separate net-of-selling-cost estate total. Show 25% and 75% recipient-share scalings only as explicitly arbitrary bounds/sensitivities, never as an adopted target mapping. Existing empirical scenario tables can separately show how much the adopted child-directed numerator differs from its all-estates counterpart. These are small accounting observations, not new equilibria.

Excluding model households with no children ever born is not an exact correction: the selected empirical row applies mechanical child shares without lifetime-child eligibility. A justified common recipient definition requires an explicit mapping decision. Until then, retaining the approved first-pass proxy is transparent; claiming exact child-directed measurement would be false.

## 3. Income context and launch consequence

The adopted B15 earnings specification remains the single persistent four-year process: persistence 0.7345934906 and innovation SD 0.4838308245. Its external estimation receipt is `output/model/native_financing_diagnostic_20260919/specification_followup/earnings_entry_battery_v1/single_process_external_estimate.json`. This review neither re-estimates it nor reopens the income architecture; its four-year block-mean estimates should not be presented as a separately identified annual persistent/transitory decomposition.

No active source, target value, weight, or scientific gate was changed. The new contract should preserve these two mapping warnings and correct the provenance text before claiming measurement review complete. The diagnostics above are recommended future checks, not results already computed and not automatic grounds to suspend an explicitly provisional, author-approved calibration. Any later observer substitution needs a new contract, a fixed-parameter old/new comparison, and explicit disclosure before the resulting losses are compared.
