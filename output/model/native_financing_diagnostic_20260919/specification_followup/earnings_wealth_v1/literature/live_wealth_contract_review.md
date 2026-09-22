# Live frozen wealth/income contract review

Date: 2026-09-21
Scope: read-only provenance check of `tmp/income_refinement_local_20260921`; no target, code, or model changes.

## Frozen runtime actually used

The frozen probe imports `intergen_eqscale_seq_optimized.solver` as `model` (`tmp/income_refinement_local_20260921/source/code/model/tools/run_e5f_matched_pf_smoke.py:30`) and the initial probe invokes that stack through `run_e5f_transition_calibration` (`run_e5f_initial_revision_probe.py:22–25, 148`). The relevant wealth observer also imports `intergen_eqscale_seq_optimized` (`e5f_initial_housing_observer.py:115–119`). The active runtime is therefore `INCFAMR`/`EARNINDRRC`-agnostic model code; those variable names occur in the pinned empirical contract, not as runtime data columns.

## Working-age wealth/earnings target

The pinned 13-row contract identifies the scored empirical record as `psid_wealth_gross_earnings_2003_2005`, with data formula
`sum(IW*NETWORTHR)/sum_{age<=65}(IW*EARNINDRRC)` and sample reference persons ages 18–85, selected PSID waves 2003 and 2005. Its pinned model formula is beginning-of-period wealth over all living ages 18–85 divided by four-year-flow income for the 12 working cells covering ages 18–65:
`sum g_asset*(b + owner*p*H) / sum_{j<J_R} g_asset*[P.income*z/((1-tau_pay)*4)]` (`tmp/income_refinement_local_20260921/external_inputs/observer_contract.json`, record `psid_wealth_gross_earnings_2003_2005`; mirrored in `working_contract.json` target row `wealth_earnings`).

The observer requires `period_years=4`, `scale_flows_to_period=True`, `J_R=12`, and complete 18–85 coverage (`e5f_initial_housing_observer.py:120–146, 152–166`). In the live solver, the aggregate diagnostic loops over working cells only and forms annual gross labor earnings from `annual_gross_income_at_state` (`intergen_eqscale_seq_optimized/solver.py:5630–5666`), then divides all-age beginning wealth by that denominator (`solver.py:5681–5691`).

## Property-tax rebate treatment

 `income_at_state` adds `property_tax_lump_sum_transfer` to labor/pension income (`solver.py:228–235`), and `annual_gross_income_at_state` annualizes that total and grosses up working income by `1-tau_pay` (`solver.py:264–272`). Consequently, a nonzero property-tax rebate would enter the runtime denominator named “annual gross labor earnings,” even though the pinned contract formula explicitly shows only `P.income*z/((1-tau_pay)*4)`. The frozen initial contract records the actual rebate as zero, so the discrepancy is inactive in the reviewed baseline. It becomes material if the funded rebate is turned on during an income calibration or wealth-target evaluation.

The property-tax budget code separately computes transfer outlays (`solver.py:5138–5149`); this confirms the rebate is a fiscal transfer, not labor earnings. The current baseline therefore excludes it numerically only because its fixed value is zero, not because the aggregate observer structurally excludes it.

## Retirement-family-income target

The retirement target is a different object. The empirical record `psid_old_wealth_dispersion_2003_2005` uses weighted (Q_{.90}(NETWORTHR/INCFAMR)/Q_{.50}(NETWORTHR/INCFAMR)), living reference persons ages 76–84, `INCFAMR>1000`, completed children observed, and 2003/2005 waves. The pinned contract warns that PSID `INCFAMR` includes transfers, asset income, and other family components, while the model denominator is pension income (plus any configured retirement-income/rebate components); the exact $1,000 cutoff and child-history selection are not reproduced.

The runtime observer uses the model's `annual_gross_income_at_state` for each old-age income state and forms ((b+owner\cdot pH)/income) before weighted quantiles (`e5f_initial_housing_observer.py:177–214`; native solver implementation `intergen_eqscale_seq_optimized/solver.py:5762–5860`). With the frozen baseline `retirement_income_z_scale=0` and `property_tax_lump_sum_transfer=0`, the denominator is the common pension flow annualized over four years. This is not the same economic object as PSID `INCFAMR`; the contract labels it a pension-income proxy and keeps the row non-production-eligible.

## Entrant wealth distribution

The active four-year calibration payload calls `external_entry_wealth_overrides()` (`intergen_eqscale_seq_optimized/calibration.py:1200–1223`), not the separate 18–24 profile. Its five baseline wealth/income nodes are

`[-2.51940697, -0.07907025, 0.10228762, 0.35287169, 3.03955200]`

with weights approximately `[0.2000027, 0.2000966, 0.2000024, 0.1998877, 0.2000107]` (`calibration.py:34–45, 49–56`). The empirical source is PSID young childless renters ages 25–35, using weighted quintile-bin means of `NETWORTH2R/INCFAMR`; reported mean is 0.17922556 and weighted median 0.09996729. The 18–24 nodes are a separate function (`calibration.py:60–96`) used by other profiles, not the active E5F four-year baseline.

At entry, the runtime multiplies these empirical ratios by annual gross model income at the entrant state and linearly distributes the resulting points over the wealth grid (`solver.py:379–399`). For working entrants, annual gross income is after-tax four-year flow divided by four and grossed up by payroll tax; because the entrant state is working-age, a nonzero property-tax rebate would also enter this conversion (`solver.py:264–272`). The entry distribution is therefore externally fixed in ratios but mechanically sensitive to the calibrated income level and any rebate.

## Interpretation for the persistent-plus-iid proposal

The working-age target is already aligned with `EARNINDRRC`-type RP/spouse gross labor earnings in the pinned contract, whereas the entrant and retirement diagnostics retain `INCFAMR`-based empirical ratios. This is a real cross-target denominator mismatch, not an `INCFAMR`/`EARNINDRRC` ambiguity inside the runtime.

For the current frozen baseline, the property-tax issue is numerically dormant because the rebate is zero. Under a proposed persistent-plus-iid income calibration, changing the annual income process changes both the working-age gross-earnings denominator and the mechanically converted entrant wealth levels; the latter can affect initial wealth and therefore wealth moments even when the entrant distribution is treated as external. The retirement row remains a separate pension-versus-family-income proxy and should not be interpreted as an earnings-process target.

An externally fixed entrant marginal can be retained as a diagnostic while keeping the unchanged 13-row objective, provided it is labeled diagnostic and excluded from the scored rows. The pinned scorer already requires all 13 rows and scores 12, with the normalization separate (`external_inputs/score_initial.py` docstring and `working_contract.json`). It does not license silently treating the entrant marginal as identified by the 13-row objective. Any production use would require an explicit contract deciding whether the `INCFAMR` entrant ratios are retained, replaced by a gross-labor denominator, or used only for a diagnostic initialization.

## Bottom-line factual receipt

- Active frozen solver/observer: `intergen_eqscale_seq_optimized`, four-year periods, 12 working cells, full 18–85 wealth stock.
- Working-age scored denominator: model gross labor earnings from four-year after-payroll-tax flows; empirical `EARNINDRRC` aggregate in the pinned target.
- Entrant marginal: five externally fixed `NETWORTH2R/INCFAMR` nodes for young childless renters ages 25–35, not the separate 18–24 profile.
- Retirement target: `NETWORTHR/INCFAMR` ages 76–84; model uses a pension-income proxy and does not reproduce PSID family-income composition.
- Property-tax rebate: enters runtime `annual_gross_income_at_state`, but baseline value is zero; the pinned working-age formula excludes it explicitly. Nonzero rebate would therefore require contract reconciliation.
- Materiality: the denominator and entrant-ratio channels can matter for persistent-plus-iid calibration; the baseline receipt is not already fully corrected for the cross-target income definitions.
