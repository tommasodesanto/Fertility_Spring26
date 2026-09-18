# Stationary-state sandbox run: all_switches_psi_fixed

Spec file: `sandbox/specs/all_switches_psi_fixed.yaml`
Overrides applied: {"child_earnings_penalty": [0, 0.2, 0.2, 0.2], "mortgage_origination_only": true, "mortgage_amortization": 0.11, "rental_wedge_intercept": 0.02, "rental_wedge_slope": 0.05, "rental_wedge_knee": 6.0, "hR_max": 11.0, "chi": 1.0, "estate_receiver": "ages_45_65", "bequest_net_of_selling_cost": true, "child_maturation_mode": "parent_age", "mu_young": 0.05, "a_rise": 34.0, "a_full": 62.0, "child_benefit_form": "log"}
Mechanism switches: {"child_earnings_penalty": [0, 0.2, 0.2, 0.2], "child_benefit_form": "log"}

## Timing
- Wall time per Bellman/root evaluation: 7496.93 s (mean over 1 evaluation(s))
- Total wall time: 7497.11 s
- Number of root evaluations (psi-normalization GE solves): 1
- Warm start: price root warm-started from /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/sandbox/baseline/parameters.csv (saved price [0.7931310188535463])
- psi-normalization status: fixed_intercept
- Grid: full production Nb (Nb=120, J=17)

## Loss: 11916.669733

## 13-row target table
| Moment | Target | Model | Gap | Weight | Loss contribution | Source |
|---|---:|---:|---:|---:|---:|---|
| Initial model completed fertility | 2.1 | 0.919937 | -1.18 | - | - | exact_stationary_key |
| Childless women, ages 40–44 | 0.198279 | 0.536989 | 0.3387 | 35532.3042455214 | 4076.440283085953 | approximate_stationary_analogue |
| Exactly one child among mothers, ages 40–44 | 0.213655 | nan | nan | 26952.820824310795 | - | unavailable_in_extract_moments |
| Period mean first-birth age | 25.9763 | 32.4198 | 6.443 | 139.82806784479274 | 5805.455848669935 | approximate_stationary_analogue |
| First births at age 30+ | 0.249278 | 0.609464 | 0.3602 | 13866.065434728798 | 1798.9002402785982 | approximate_stationary_analogue |
| Wealth / annual gross labor earnings | 6.14586 | 4.98241 | -1.163 | 7.595098472533724 | 10.280845481756975 | exact_stationary_key |
| Annual bequests / aggregate wealth | 0.0088 | 0.00839335 | -0.0004066 | 5165289.256198346 | 0.8541533320383308 | exact_stationary_key |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 3.55927 | 0.04334 | 10.616361531535473 | 0.01993971502275099 | exact_stationary_key |
| Mean occupied rooms, capped at 9 | 5.5611 | 5.25321 | -0.3079 | 128.02070205233477 | 12.13571667481115 | approximate_stationary_analogue |
| Ownership, heads 30–55 | 0.648334 | 0.54975 | -0.09858 | 2339.3623724673616 | 22.735993794849428 | exact_stationary_key |
| First-birth room response, −1 to +3 | 0.720246 | 0.516287 | -0.204 | 137.5652749002964 | 5.722650349072365 | approximate_stationary_analogue |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.908769 | 0.5617 | 280.52808370152104 | 88.50909680543218 | approximate_stationary_analogue |
| Recent-parent ownership gap | 0.162896 | 0.222343 | 0.05945 | 27055.822957508266 | 95.6149649829898 | exact_stationary_key |

**Caveat on the 4 fertility-timing rows** (childless women 40-44, exactly-one-child mothers 40-44, mean first-birth age, first births at 30+): the retained target_fit.csv computes these from a specialized dated period/cohort-timing pipeline (`transition_cross_section_moments` / `cohort_timing_moments` in run_e5f_transition_calibration.py), not from the plain stationary `extract_moments()` output this sandbox reuses. The sandbox reports the closest available stationary analogue and flags it 'approximate' below; do not treat these four rows as exact reproductions.

## Parameter table
| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| H0 | 8.11210048056786 | 0.2 | 80.0 | False | structural coordinate (retained, possibly spec-overridden) |
| beta_annual | 0.99 | 0.94 | 0.99 | True | structural coordinate (retained, possibly spec-overridden) |
| chi | 1.0496534423047694 | 0.1 | 5.0 | False | structural coordinate (retained, possibly spec-overridden) |
| first_birth_fixed_cost | 0.26576477618628114 | 0.0 | 8.0 | False | structural coordinate (retained, possibly spec-overridden) |
| h_P | 2.3 | 0.1 | 2.3 | True | structural coordinate (retained, possibly spec-overridden) |
| kappa_fert | 0.33773423411167025 | 0.02 | 50.0 | True | structural coordinate (retained, possibly spec-overridden) |
| kappa_fert_continuation | 0.39785645171756406 | 0.02 | 50.0 | True | structural coordinate (retained, possibly spec-overridden) |
| theta0 | 0.08105103333987912 | 0.0 | 8.0 | True | structural coordinate (retained, possibly spec-overridden) |
| theta1 | 0.08520356663830632 | 0.02 | 16.0 | True | structural coordinate (retained, possibly spec-overridden) |
| psi_child | 0.1489153145785918 | None | None | False | normalized to completed fertility 2.1 unless spec fixes psi |
| payroll_tax | 0.179 | None | None | False | externally fixed (P.tau_pay) |
| housing_supply_elasticity | None | None | None | False | not applicable: this is a dated-transition supply_rule.elasticity value (0.63 in the retained report); the bare stationary solve has no such object |
| tenure_choice_kappa | 0.005 | None | None | False | externally fixed or sandbox switch |
| alpha_cons | 0.733 | None | None | False | externally fixed or sandbox switch |
| sigma | 2.0 | None | None | False | externally fixed or sandbox switch |
| child_benefit_form | log | None | None | False | externally fixed or sandbox switch |
| child_benefit_curvature | None | None | None | False | externally fixed or sandbox switch |
| scale_weighting | None | None | None | False | externally fixed or sandbox switch |
| phi | [0.8, 0.8, 0.8, 0.8] | None | None | False | financed share by parity; solved P.phi (default 0.80, spec-overridable; see apply_spec's phi/n_parity broadcast workaround) |
| lambda_d | 0.0 | None | None | False | unsecured debt-line multiple of mean earnings (default 0.0, spec-overridable) |

