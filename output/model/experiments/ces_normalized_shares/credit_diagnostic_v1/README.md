# Mortgage-financing diagnostic in the normalized CES-share experiment

Reference: verified overnight chain 1, price 0.6408361332017416 and housing supply coefficient 7.306962620836552. [Full baseline target fit and all parameter estimates/bounds](../overnight_v1/final_results/README.md). This is an author-requested diagnostic, not an adopted specification or policy result.

Only the uniformly financed mortgage share changes from 0.80 to 0.95. Preferences, all eleven fitted coordinates, the normalized CES-limit denominator, parenthood jump and slope, zero physical housing floor, earnings, entry distributions, taxes, retirement benefits, timing, birth architecture and grids are retained. The house price, rent mapping and housing supply coefficient are held fixed. Down-payment requirements decline from 20% to 5% under the existing joint consumption/saving/purchase budget. Renter unsecured credit remains unchanged.

Two fresh fixed-price lifecycle solutions are planned, one control and one policy case. These are normalized household cross-sections, not a dated transition or a market-clearing policy equilibrium. No recalibration, root search or welfare interpretation is planned. Explicit first-, second- and third-birth flows must be distinguished from topcode-adjusted child units and completed fertility.

The authenticated 3,008-file v5 source remains immutable and read-only. Separate helper scripts are [the driver](../../../../../code/cluster/ces_normalized_shares_calibration/credit_diagnostic.py) and [the launch command](../../../../../code/cluster/ces_normalized_shares_calibration/credit_diagnostic.sh). Reproduce with `bash code/cluster/ces_normalized_shares_calibration/credit_diagnostic.sh` from the repository, with Torch authentication available. The uploaded reference JSON is pinned by SHA-256. Before solving, a zero-lifecycle preflight compares public inputs (only phi may differ), exact grids, and reached consumption-share and material-utility arrays against the normalized utility formula.

Budget: one Torch CPU, 24 GiB, 30 minutes; two lifecycle solves, each subject to the existing native 300-second solve cap, with 600 seconds per case for solving/reporting and a 1,700-second outer timeout. Per-case started/completed/failure receipts and latest-completed summary are saved. Unexpected numerical or acceptance failures terminate without retry or gate changes. Each case retains full numerical arrays, executed parameters, observations and the unchanged 17 standard diagnostic graphs.

Status: both cases solved in job **19161836**; the relaxed case failed the estate gate. Descriptive results, full comparison and limitations follow. Preflight passed exact input/grid and reached utility checks; [receipt](preflight.json).

## Completed descriptive comparison

Verified October 4, New York. Job 19161836 completed both fixed-price solves (26.6 and 15.9 seconds). The 80% control passed the native gate stack; the 95% arm stopped at the no-negative-estates gate. Its saved-array forensic report retains the rejection and makes **zero** additional lifecycle solves. All 17 standard plots exist for both arms. These results are descriptive, not production-accepted policy estimates.

Only financed share changes. In this comparison mortgage relief moves ownership strongly, while births and room demand barely change. This is a finding at this experimental fit and held price; it is not a general claim about financial constraints.

| Quantity | 80% financed | 95% financed | Change |
|---|---:|---:|---:|
| first births per normalized household | 0.0427030936 | 0.0427213491 | 0.0427499% |
| second births per normalized household | 0.0404731232 | 0.040493818 | 0.0511322% |
| third births per normalized household | 0.0289972493 | 0.0290187378 | 0.0741055% |
| Explicit births, all three orders | 0.112173466 | 0.112233905 | 0.0538798% |
| Native topcode-adjusted fertility measure | 2.09999958 | 2.10118829 | 0.0566049% |
| childless_rate_40_44 | 0.310121308 | 0.309823091 | -0.000298216594 (level) |
| exactly_one_among_mothers_40_44 | 0.0744111956 | 0.0742979924 | -0.000113203235 (level) |
| mean_children_ever_born_capped3_age25 | 0.632140381 | 0.633374921 | 0.00123454043 (level) |
| period_mean_age_first_birth | 23.4943372 | 23.4855017 | -0.00883555338 (level) |
| period_share_first_births_age30plus | 0.0764605952 | 0.0762047717 | -0.000255823415 (level) |
| aggregate_mean_occupied_rooms_ahs_uncapped_18_85 | 5.95165759 | 5.9601774 | 0.0085198103 (level) |
| aggregate_mean_occupied_rooms_capped9_18_85 | 5.80923926 | 5.81556133 | 0.00632206959 (level) |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.87299722 | 6.84163183 | -0.0313653962 (level) |
| annual_bequest_flow_to_aggregate_wealth | 0.00731368437 | 0.00719875963 | -0.000114924746 (level) |
| housing_increment_0to1 | 0.683226667 | 0.672406614 | -0.0108200529 (level) |
| old_total_wealth_to_annual_income_median_7684 | 14.0885349 | 13.9311923 | -0.1573426 (level) |
| old_total_wealth_to_annual_income_p90_p50_7684 | 2.82894831 | 2.86089921 | 0.0319508965 (level) |
| own_rate_25_34 | 0.589158295 | 0.7312792 | 0.142120905 (level) |
| own_rate_30_55 | 0.704661732 | 0.77952958 | 0.0748678483 (level) |
| prime30_55_model_dependent_3plus_minus_1to2_rooms_capped9 | 0.405224009 | 0.397019311 | -0.00820469848 (level) |
| recent_parent_minus_no_resident_child_ownership_30_55 | unavailable | unavailable | unavailable (level) |

Birth flows are normalized cross-section flows, not national births. The native topcode-adjusted fertility measure is separate from explicit first-, second- and third-birth events. Ownership differences are proportion changes (multiply by 100 for percentage points). Exactly-one-child fertility is conditional on being a mother in the ages 40–44 observer. Housing rows retain their stated age windows and caps; the uncapped rooms row is the baseline calibration rooms definition. The recent-parent cross-sectional ownership row is unavailable in this helper; it does not replace the separate production flow-response observer.

### Acceptance failure and missing diagnostics

The estate ledger is funded, but its net-negative estate amount is 6.71822679925e-09 model financial units per model period, exceeding the existing 1e-10 no-negative-estate tolerance. Negative-estate death mass is 1.70978995218e-07, or 2.76963 per million deaths. Its cause is not diagnosed here. No tolerance, mortality rule, financing rule, supply coefficient, transfer or utility was changed to obtain acceptance. Checks after the failed estate gate did not execute in that stack; cached reconstruction and descriptive reporting do not confer acceptance.

The first job, 19161617, completed the control and its 17 plots but stopped in JSON reporting. Seven unused permanent-income diagnostics are undefined because this model has no separate permanent-income groups. The report-only correction records these as null with exact paths; requested birth/fertility/housing values are finite. The original failure is preserved. The first cached inspection stopped because saved arrays lack private event-flow fields; the corrected cached reader uses the authoritative forward fertility observer and checks its control flow totals against the saved native totals (error below 1.4e-17). This is a reporting correction, not a new solve.

[Control plots](control_standard_diagnostics/) · [Relaxed plots, estate-gate rejected](relaxed_standard_diagnostics/) · [Control observations](control_summary.json) · [Relaxed cached observations](relaxed_cached_observations.json) · [Estate ledger](estate_ledger.json) · [Native rejection](relaxed_gate_failure.json) · [Zero-solve cached completion](cached_inspection_completed.json). Large numerical arrays and executed inputs remain in the remote job output.

The relaxed arm failed before exporting `executed_P.json`; its 87 stage arrays and pinned primitive inputs remain saved. The cached reader reconstructs those inputs and requires exact equality of every saved shared numerical array before observing them. Original private execution-only fields are not claimed to have been recovered. All 34 downloaded PNG hashes match the remote packets ([hashes](plot_hashes.json)).
