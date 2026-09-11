# Bounded beta diagnosis — September 11, 2026

Read-only diagnosis from frozen-source hashes and saved calibration outputs; no model solve, cluster submission, target change or bound change. The selected initial stationary candidate has loss 158.07639074681484 and annual discount factor 0.9957960668523341. This is a working minimum-distance candidate, not an identified final preference estimate.

## Verified implementation

The actual initial driver inherits the normalized old checkpoint and applies all nine coordinates through `bind_parenthood_utility` (`tmp/e5f_matched_pf/code/model/tools/run_e5f_initial_revision_probe.py:75` and `:93`). The annual beta coordinate is explicitly raised to the fourth power (`e5f_parenthood_utility.py:144`) and the rho aliases are updated (`:154`). Thus beta4 = 0.983290008859785 and the annual pure discount rate 1/beta−1 = 0.4221681%. This is high patience; it is not an annual/four-year reporting confusion.

The checked local parameters, solver, local_panel, parenthood utility, initial driver, sequential configuration, wealth observer and utils files match SHA256 pins in the working contract's frozen 70abd4a8 observation snapshot. This establishes identity for inspected files, not the complete local tree.

`intergen_eqscale_seq_optimized/parameters.py:168` defines q = 1.04^4−1 and `:730` scales annual labor income and pensions by four. `solver.py:1075` sets R_gross=1+q; the Markov household branch at `:2523` budgets R_gross*b + period income. `solver.py:6127` divides period net labor income by four and by (1−payroll tax) to construct annual gross earnings; `:6136` annualizes death-estate flow by four. There is no obvious unit mismatch in these inspected formulas. A subsequent SHA256-verified read of the selected checkpoint establishes actual q=0.08243216 and R_gross=1.08243216 per four-year period, equivalent to an annual return of2%, rather than the4% source default. Actual beta4*R4=1.064345. The source-default calculation is not the runtime economy. See current_candidate_transition/runtime_beta_check.json.

The Markov Bellman continuation uses survival*expected next value+(1−survival)*bequest value (`solver.py:2499–2513`), then applies beta once in the saving objective/kernel. Terminal continuation is the bequest value. The KFE also survival-weights advancing mass (`:4891`), and the wealth observer counts terminal deaths with probability one (`:6116`). The subsequent verified checkpoint read confirms use_age_survival=True, with survival1 for the first12 transitions and then0.9391263,0.9184976,0.8849522,0.8300468. The wealth grid has120points with b_max=30. This verifies the serialized settings, not their empirical adequacy.

## What saved probes establish

The transformed beta coordinate is x=log(−log(beta_annual)); increasing x lowers beta (`initial_calibration_contract/extended_refinement/run_search.py:27`). Each coordinate probe re-normalizes psi_child to completed fertility 2.1 and re-solves housing and the balanced pension, so it is a conditional equilibrium response, not a partial derivative holding every other economic object fixed.

The previous selected beta was 0.9989441232 at loss278.708079. The latest selected beta is lower, 0.9957960669 at loss158.076391. Beta is not at its 0.9995 ceiling and the saved 1%-of-range near-bound flag is False. The last saved coordinate probes favor the higher beta: beta0.9952743905 gives loss160.185508, while beta0.9950820097 gives loss163.161957. Rounds2–3 instead favor the lower-beta probe; rounds0,1,4,5 favor the higher. The objective does not uniformly drive beta upward.

At the last derivative center r4_joint_10 (not at the final selected point), lowering beta lowers wealth/annual earnings, ownership30–55, and recent-parent ownership. These are below target at the final candidate, so their economic direction supports high patience. Lowering beta also lowers average rooms and the first-birth room response, which overshoot their targets; those directions oppose high patience. Therefore a generic claim that all housing targets require high beta is false.

A sharper concern is the old-age wealth/income p90/p50 row. Its last beta derivative is14.355 weighted residual units per unit x, versus only0.021–0.048 in rounds0–4. It explains97.1% of the last beta column's squared norm. The observer constructs wealth/income values on the asset/tenure grid and takes p90/p50 (`e5f_initial_housing_observer.py:177–214`); `utils.py:225–245` returns the first grid-supported value whose cumulative weight crosses each quantile. A percentile crossing is therefore a plausible explanation for the derivative jump. The saved Jacobian alone does not prove the jump's exact source; retrieve the two beta probes' p50/p90 numerators and local cumulative weights before interpreting this sensitivity economically.

The last 12x9 weighted Jacobian has rank9 and condition2467.12 in transformed coordinates. Its weakest right-singular vector is dominated by theta1 (coefficient−0.966), with theta0−0.219 and beta-coordinate−0.119. Projecting the beta column onto only the two bequest columns leaves norm2.298 versus original14.570; this supports a local discounting/bequest substitution concern. It is not proof of identification: the column is dominated by a possibly nonsmooth percentile response. The Jacobian center and final selected residuals differ, so this is not a final-selected-point gradient. The exact saved beta-pair contribution decomposition, collected afterward, is below.

## Exact saved beta-pair decomposition (follow-up)

The lead collected the three frozen scored outputs at `current_candidate_transition/beta_saved_scores/{r4_joint_10,r5_d0_-1,r5_d0_+1}.json`. All three share the approved objective fingerprint; each complete loss recomputes exactly to floating-point tolerance; all eight non-beta structural coordinates agree. The finite difference of the12 weighted residuals reproduces the saved last beta Jacobian column. This is stronger evidence than applying that Jacobian to a different selected point.

**Yes: the old-age dispersion movement favors higher beta, and accounts for all of this pair's net improvement.** Raising annual beta from0.9950820097 to0.9952743905 lowers loss from163.161957 to160.185508 (−2.976449). The old-age wealth/income p90/p50 falls from4.410197 to4.233969 toward target3.515935; its contribution falls by3.016451. The other11 rows together increase loss by0.040002. Thus the saving and ownership gains alone do not explain the direction of this particular beta comparison. The97.1% figure above concerns the derivative column's squared norm; it is a different statistic from this exact loss-contribution attribution (101.34% of the net improvement).

The local tradeoff is economically clear: higher beta improves wealth/earnings (−0.232837 loss), overall ownership (−0.259160) and recent-parent ownership (−0.129698), but worsens the first-birth rooms overshoot (+0.535487), mean rooms (+0.145520), mean first-birth age (+0.114787), and bequest/wealth (+0.039610). The remaining fertility/family rows partly offset those costs. This is the last pair around beta0.995179, not a conclusion about every earlier stage or the selected beta0.995796.

The files contain the p90/p50 ratio, not separate p50/p90 values or cumulative weights, so they establish the dispersion movement but do not prove that a discrete quantile crossing caused it. Raw observer numerators/denominators or checkpoint distributions remain needed for that distinction. In this exact pair psi_child stays fixed at0.1787090653 because both initial guesses pass the normalization gate; completed fertility is2.10006145 at higher beta and2.10026235 at lower beta. The controller attempted normalization as required, but these particular probes did not need an intercept update.

| Scored row | Model, lower beta | Model, higher beta | Loss change, higher minus lower |
|---|---:|---:|---:|
| cps_childlessness | 0.2051592294 | 0.205140391 | -0.009198611 |
| cps_exactly_one | 0.2045783556 | 0.2046460414 | -0.032995169 |
| nchs_mean_age | 26.10168631 | 26.1049173 | +0.114787290 |
| nchs_share30 | 0.2306227727 | 0.2307531023 | -0.067190433 |
| wealth_earnings | 5.397983201 | 5.418767497 | -0.232837295 |
| bequest_wealth | 0.008276240705 | 0.008268970485 | +0.039610265 |
| old_dispersion | 4.410197053 | 4.233968658 | -3.016450750 |
| mean_rooms | 6.330234787 | 6.330973373 | +0.145520470 |
| ownership_30_55 | 0.5443594921 | 0.5448936011 | -0.259159701 |
| first_birth_rooms | 1.143640028 | 1.148212244 | +0.535486682 |
| family_rooms | 0.1871203269 | 0.1878387282 | -0.064323851 |
| recent_parent_ownership | 0.1490077914 | 0.1491814656 | -0.129697864 |

## Saved percentile components now recovered

The three raw early_measurement.json outputs provide separate percentiles, now
verified against their scored ratios and saved with source hashes in
current_candidate_transition/beta_percentile_components.json. From lower to
higher beta, p50 rises4.098095536→4.269338274 (about4.18%), while p90 rises only
18.073408858→18.076244443 (about0.016%). The ratio change is therefore driven
by the median. The median at the derivative center is4.164790385. This narrows
the next diagnostic to the weighted CDF and support near its50%crossing; it
still does not prove that the movement is a discretization artifact. No moment,
weight, parameter or grid has been changed.

## Specific follow-up, without changing the target contract

The saved center/pair scores have now been checked above. Next retrieve their raw early wealth-observer details and separate p50/p90 values; the exact contribution changes are already established. Read runtime beta, q, R_gross, period_years, flow-scaling flag, survival schedule, pension and bequest specification from the selected checkpoint, and verify their frozen-source provenance.

Then, if a numerical diagnostic is authorized, use one bounded selected-point panel: exact replay plus beta perturbations at x±0.005 and x±0.02 (five candidate evaluations, at most40stationary solves under the existing eight-solve normalization cap, each keeping the complete12-row objective and all original gates). Save p50/p90, their bracketing wealth states and cumulative weights, occupied upper-asset-bound mass, and age-specific saving. If beta derivatives differ sharply with step size, repeat only the offending center/pair on a locally refined wealth grid. If derivatives are stable, a separate small fixed-beta profile allowing only theta0/theta1 to adjust distinguishes discounting/bequest substitution from a genuine loss penalty. Hold the other six structural coordinates and all targets/weights fixed; retain the separate2.1 normalization. Prestate solve count and cluster time budget before running. This diagnosis alone provides no reason to lower beta's upper bound or remove a target.

## Complete selected target fit

All entries below are copied from the selected table; weights are working minimum-distance weights, not all certified inverse sampling variances.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| initial_normalization | 2.1 | 2.100150605 | 0.0001506046984 | unscored | unscored |
| cps_childlessness | 0.198278751 | 0.2050178343 | 0.006739083314 | 35532.30425 | 1.613708264 |
| cps_exactly_one | 0.2136553252 | 0.2053218834 | -0.008333441786 | 26952.82082 | 1.871772387 |
| nchs_mean_age | 25.97626386 | 26.11316816 | 0.1369042962 | 139.8280678 | 2.620767596 |
| nchs_share30 | 0.249278013 | 0.230558772 | -0.01871924105 | 13866.06543 | 4.858807788 |
| wealth_earnings | 6.145861394 | 5.470594365 | -0.6752670298 | 7.595098473 | 3.463255242 |
| bequest_wealth | 0.0088 | 0.008319317161 | -0.0004806828387 | 5165289.256 | 1.19347103 |
| old_dispersion | 3.515935087 | 4.335051602 | 0.8191165151 | 10.61636153 | 7.123067573 |
| mean_rooms | 5.561097376 | 6.32135461 | 0.7602572339 | 128.0207021 | 73.99482151 |
| ownership_30_55 | 0.6483340343 | 0.5348647079 | -0.1134693264 | 2339.362372 | 30.11996436 |
| first_birth_rooms | 0.7202462624 | 1.113533432 | 0.3932871697 | 137.5652749 | 21.27788109 |
| family_rooms | 0.3470669318 | 0.2130829831 | -0.1339839487 | 280.5280837 | 5.03595558 |
| recent_parent_ownership | 0.1628955092 | 0.1494339007 | -0.01346160849 | 27055.82296 | 4.902918337 |

## Complete selected parameter restrictions

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.9957960668523341 | 0.94 | 0.9995 | False | diagnostic candidate; not a certified estimate |
| kappa_fert | 0.3346343063057974 | 0.02 | 50.0 | True | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.4037305572536571 | 0.02 | 50.0 | True | diagnostic candidate; not a certified estimate |
| chi | 1.045632526446847 | 0.1 | 5.0 | False | diagnostic candidate; not a certified estimate |
| H0 | 8.225499724508786 | 0.2 | 80.0 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.08957524492721414 | 0.0 | 8.0 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.0770953752444599 | 0.02 | 16.0 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.2937155016899614 | 0.0 | 8.0 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.1837953641959 | 0.1 | 2.3 | False | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0.0 | — | — | False | zero restriction |
| psi_child | 0.16525977525321506 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.0463613896121218 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2.0 | — | — | False | externally fixed |

Sources: `initial_calibration_contract/working_contract.json`; `extended_refinement/collected_17378993/{cases.json,jacobian_round_0.json,...,jacobian_round_5.json,selected_target_fit.csv,selected_parameters.csv}`; prior `three_hour_continuation/collected_17376529/cases.json`. All relative paths in this note are beneath the repository or, for result tables, beneath this note's directory.
