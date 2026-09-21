# Reviewed local calibration panel — September 21

All **28 full objective evaluations passed**, using **168 nested stationary solves**. There are no rejected attempts, missing probes or incomplete solve counts. Two fresh anchor runs match the saved reference and two fresh selected runs match the selected point, exactly in the recorded numeric and checkpoint-array signatures. Smoke ran 03:28:08–04:01:34 UTC and production 04:01:45–05:01:23 UTC, finishing at 01:01 EDT. Raw cluster receipts and plots are retained under [collected](../collected/).

The best tested point has loss **326.983198873**, versus **353.658872914** at the anchor: a **7.54277%** reduction under the identical target/weight system. It raises the existing continuation-fertility taste scale by 5%, from 0.346493693274 to 0.363818377938; all other structural coordinates are unchanged. Each case re-solves prices, stationary population and the maintained fertility normalization. This is a finite local diagnostic, not convergence, a certified calibration, or adoption of the diagnostic earnings candidate.

Most improvement comes from first-birth age (loss contribution falls by 20.3403) and childlessness (by 5.2139). Housing levels barely change: rooms remain 1.11693 above target and ownership 16.4443 percentage points below. Rooms and ownership together account for about 68.2% of the selected loss. There is useful scope within the existing parameters, but the core housing fit remains unresolved.

## Complete selected fit

The first row is an unscored maintained stationary normalization, not a newly measured empirical fertility target. Gap is model minus target; each scored contribution equals weight × gap squared. All 364 target rows and 336 scored contributions across the 28 evaluations were independently checked. [All rows](all_target_fit.csv), [selected CSV](selected_target_fit.csv).

| Moment | Target | Anchor | Selected | Gap | Weight | Contribution |
|---|---:|---:|---:|---:|---:|---:|
| Maintained stationary fertility normalization | 2.1 | 2.10027054 | 2.10027634 | 0.000276343605 | — | — |
| Childless women, ages 40–44 | 0.198278751 | 0.22516949 | 0.222286584 | 0.024007833 | 35532.3042 | 20.4799691 |
| Exactly one child among mothers, ages 40–44 | 0.213655325 | 0.197955732 | 0.199217717 | -0.0144376085 | 26952.8208 | 5.61816832 |
| Period mean first-birth age | 25.9762639 | 26.7000934 | 26.591457 | 0.615193103 | 139.828068 | 52.9196876 |
| First births at age 30+ | 0.249278013 | 0.259828856 | 0.25336475 | 0.00408673663 | 13866.0654 | 0.231582931 |
| Wealth / annual gross labor earnings | 6.14586139 | 6.72243214 | 6.72248587 | 0.57662448 | 7.59509847 | 2.52533827 |
| Annual bequests / aggregate wealth | 0.0088 | 0.00797719829 | 0.0079784584 | -0.000821541603 | 5165289.26 | 3.48621181 |
| Old wealth/income p90 / median, ages 76–84 | 3.51593509 | 4.28268849 | 4.28135774 | 0.765422654 | 10.6163615 | 6.21982726 |
| Mean occupied rooms, capped at 9 | 5.56109738 | 6.67698129 | 6.67802545 | 1.11692808 | 128.020702 | 159.709452 |
| Ownership, heads 30–55 | 0.648334034 | 0.484270888 | 0.48389152 | -0.164442514 | 2339.36237 | 63.2594941 |
| First-birth room response, −1 to +3 | 0.720246262 | 1.01391431 | 1.01653897 | 0.296292708 | 137.565275 | 12.0767687 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347066932 | 0.344338361 | 0.340549789 | -0.00651714268 | 280.528084 | 0.011914911 |
| Recent-parent ownership gap | 0.162895509 | 0.162282188 | 0.158840945 | -0.0040545641 | 27055.823 | 0.444783933 |

## Complete parameter table

The nine structural coordinates use the actual original-search bounds in the hash-pinned `prior_plan.json`. In particular, annual beta is **at its 0.99 upper bound**, despite the raw scorer's generic 0.9995 bound and false near-bound flag. Near-bound means within 1% of the raw level interval; it is a descriptive flag, especially weak for wide log-search intervals, not an economic identification claim. Other rows record their normalization or external restriction. No parameter is adopted. [CSV](selected_parameters.csv).

| Parameter | Selected | Actual lower | Actual upper | Near bound | Restriction/status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.99 | 0.94 | 0.99 | upper | diagnostic candidate; not a certified estimate |
| kappa_fert | 0.197771155 | 0.02 | 50 | lower | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.363818378 | 0.02 | 50 | lower | diagnostic candidate; not a certified estimate |
| chi | 0.953059471 | 0.1 | 5 | no | diagnostic candidate; not a certified estimate |
| H0 | 10.1553845 | 0.2 | 80 | no | diagnostic candidate; not a certified estimate |
| theta0 | 0.170783494 | 0 | 8 | no | diagnostic candidate; not a certified estimate |
| theta1 | 0.0734872299 | 0.02 | 16 | lower | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.167381665 | 0 | 8 | no | diagnostic candidate; not a certified estimate |
| h_P | 2.3 | 0.1 | 2.3 | upper | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | not searched | zero restriction |
| psi_child | 0.191133221 | — | — | not searched | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | not searched | externally fixed |
| pension_period | 2.04636139 | — | — | not searched | budget derived |
| housing_supply_elasticity | 0.63 | — | — | not searched | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | not searched | externally fixed |
| alpha_cons | 0.733 | — | — | not searched | externally fixed |
| sigma | 2 | — | — | not searched | externally fixed |

There are 12 scored moments for nine free structural coordinates plus the separate normalization. This count avoids simple parameter-count underidentification; it does not establish informative independent variation or statistical identification. The saved scores explicitly remain working minimum-distance diagnostics.

## Local response map and step stability

The driver's `finite_difference_jacobian.json` contains 24 **directional per-probe slopes**, not nine paired columns. Here positive and negative probes are paired for seven interior coordinates: $D_j=[m(\theta+h_j)-m(\theta-h_j)]/(2h_j)$. Beta and $h_P$ are at upper bounds and use the inward one-sided slope $D_j=[m(\theta)-m(\theta-h_j)]/h_j$. No missing responses are imputed.

For comparisons and conditioning, columns are $\sqrt W D_j h_j$: one unit of the coordinate is its planned full step, not a log change or one bound span. Column signs follow an increase in the parameter even where the derivative is estimated from the inward move. The relative half-step discrepancy is $\|\sqrt W(D_{j,h/2}-D_{j,h})h_j\|_2/\|\sqrt W D_{j,h}h_j\|_2$. This is descriptive, without an invented acceptance cutoff.

| Parameter | Difference | Full step | Half-step relative discrepancy | Full/half cosine |
|---|---|---:|---:|---:|
| H0 | central | 0.507769224 | 0.0445936409 | 0.999107421 |
| beta_annual | one-sided | 0.001 | 0.906896437 | 0.954880294 |
| chi | central | 0.0476529736 | 0.145318313 | 0.992681949 |
| first_birth_fixed_cost | central | 0.01 | — | — |
| h_P | one-sided | 0.05 | 0.0102314762 | 0.999976718 |
| kappa_fert | central | 0.00988855775 | 0.00037379639 | 0.99999993 |
| kappa_fert_continuation | central | 0.0173246847 | — | — |
| theta0 | central | 0.0085391747 | — | — |
| theta1 | central | 0.00367436149 | — | — |

Beta is step-sensitive: the full inward step lowers loss to 349.5934, whereas the half step gives 353.7171, slightly above the 353.6589 anchor. Its weighted response magnitude changes substantially (90.7% relative discrepancy). Chi has a 14.5% discrepancy; H0 4.46%, child-space requirement 1.02%, and the first fertility taste scale 0.0374%. Four coordinates, including the selected continuation scale, have no half-step test. Exact repetition demonstrates reproducibility, not derivative accuracy; this panel does not distinguish curvature from grid or solver effects.

In the stated scaling the 12-by-9 matrix has condition number **294.44**. Replacing the five tested columns by half-step estimates gives **556.52**, retaining full-step estimates for the other four. Neither number is an identification proof; their sensitivity reinforces caution about using this map as a precise optimization Jacobian. [Paired response data](paired_response_map.csv), [step checks](step_stability.csv), [supplemental heatmap](supplemental_response_map.png).

The local columns show that H0 and the uniform ownership preference move average rooms and ownership in the same direction at this point, while the desired residual correction requires fewer rooms and more ownership. They also move the recent-parent ownership gap. This is a local trade-off, not a proof that joint changes or other regions cannot fit.

## Standard diagnostic figures

The unchanged 17-figure set is retained below, separately from the supplemental response map. The lead visually inspected the three contact sheets. Market clearing is tight (relative residual about 1.04e-7), but that does not certify empirical fit. The plots show concentration of owner housing demand on the ten-room rung, large differences across income states, and continued ownership growth at older ages. They do not contain matched empirical lifecycle overlays, so lifecycle validation remains outstanding. Conditional state policies should not be mistaken for realized population averages. Visible wealth thresholds and some nonmonotone owner-entry segments warrant keeping the existing numerical/fit cautions; no new bug is established by the figure alone.

- [fertility_by_age](../collected/selected/standard_diagnostics/fertility_by_age.png)
- [fertility_policy_by_age_income_state](../collected/selected/standard_diagnostics/fertility_policy_by_age_income_state.png)
- [housing_by_age_income_state](../collected/selected/standard_diagnostics/housing_by_age_income_state.png)
- [housing_market](../collected/selected/standard_diagnostics/housing_market.png)
- [housing_prices](../collected/selected/standard_diagnostics/housing_prices.png)
- [income_state_outcomes](../collected/selected/standard_diagnostics/income_state_outcomes.png)
- [liquid_wealth_by_age_income_state](../collected/selected/standard_diagnostics/liquid_wealth_by_age_income_state.png)
- [market_clearing_by_market](../collected/selected/standard_diagnostics/market_clearing_by_market.png)
- [market_clearing_residuals](../collected/selected/standard_diagnostics/market_clearing_residuals.png)
- [owner_rungs](../collected/selected/standard_diagnostics/owner_rungs.png)
- [ownership_by_age](../collected/selected/standard_diagnostics/ownership_by_age.png)
- [ownership_by_age_income_state](../collected/selected/standard_diagnostics/ownership_by_age_income_state.png)
- [policy_childless_renter_age30](../collected/selected/standard_diagnostics/policy_childless_renter_age30.png)
- [policy_childless_renter_age42](../collected/selected/standard_diagnostics/policy_childless_renter_age42.png)
- [tenure_services](../collected/selected/standard_diagnostics/tenure_services.png)
- [wealth_dist_childless_renter_age30](../collected/selected/standard_diagnostics/wealth_dist_childless_renter_age30.png)
- [wealth_dist_childless_renter_age42](../collected/selected/standard_diagnostics/wealth_dist_childless_renter_age42.png)

[Contact sheet 1](standard_gallery_1.png), [2](standard_gallery_2.png), [3](standard_gallery_3.png). Standard source series and figure metadata: [summary](../collected/selected/standard_diagnostics/summary.json).

## Reproduction and provenance

Run from the repository root: `code/model/.venv/bin/python code/model/tools/analyze_e5f_income_fit_sensitivity.py`. This reads saved receipts only and regenerates checked tables, paired slopes, stability diagnostics and the supplemental heatmap; it does not launch or solve a model. The collected original figure files remain unchanged.

The objective, target/weight fingerprint, actual bounds, selected array signatures and all scorer source fingerprints are preserved in [lead review receipt](review_receipt.json). Remote/source/figure collection evidence is indexed in [collection receipt](../collected/collection_receipt.json). The lead rehashed all 51 selected/repetition PNGs and verified equality of the three 17-figure sets. All 28 scored receipts share the same source fingerprints; the five collected anchor/selected/repetition numerical gate receipts are verified. Source checks here use the certified pinned-run receipts, not a fresh morning rehash of the remote source tree. Raw generic parameter bounds remain visible in the raw files; the reviewed table corrects their annotation without altering evidence.
