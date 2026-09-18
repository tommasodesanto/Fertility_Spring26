# Stationary-state sandbox run: sw_mortgage_psi_root

Spec file: `sandbox/specs/sw_mortgage_psi_root.yaml`
Overrides applied: {"mortgage_origination_only": true, "mortgage_amortization": 0.11}
Mechanism switches: (none; all at default-off)

## Timing
- Wall time per Bellman/root evaluation: 163.09 s (mean over 10 evaluation(s))
- Total wall time: 1632.49 s
- Number of root evaluations (psi-normalization GE solves): 10
- Warm start: price root warm-started from /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/sandbox/baseline/parameters.csv (saved price [0.7931310188535463])
- psi-normalization status: derived_intercept
- Grid: full production Nb (Nb=120, J=17)

## Loss: 2050.100094

## 13-row target table
| Moment | Target | Model | Gap | Weight | Loss contribution | Source |
|---|---:|---:|---:|---:|---:|---|
| Initial model completed fertility | 2.1 | 2.1 | 1.033e-07 | - | - | exact_stationary_key |
| Childless women, ages 40–44 | 0.198279 | 0.200717 | 0.002438 | 35532.3042455214 | 0.21124146602001106 | approximate_stationary_analogue |
| Exactly one child among mothers, ages 40–44 | 0.213655 | nan | nan | 26952.820824310795 | - | unavailable_in_extract_moments |
| Period mean first-birth age | 25.9763 | 26.1237 | 0.1474 | 139.82806784479274 | 3.039806520378064 | approximate_stationary_analogue |
| First births at age 30+ | 0.249278 | 0.234164 | -0.01511 | 13866.065434728798 | 3.167643753305505 | approximate_stationary_analogue |
| Wealth / annual gross labor earnings | 6.14586 | 5.2106 | -0.9353 | 7.595098472533724 | 6.643480832524264 | exact_stationary_key |
| Annual bequests / aggregate wealth | 0.0088 | 0.00775535 | -0.001045 | 5165289.256198346 | 5.636848645165049 | exact_stationary_key |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 5.04476 | 1.529 | 10.616361531535473 | 24.813598379447065 | exact_stationary_key |
| Mean occupied rooms, capped at 9 | 5.5611 | 5.66007 | 0.09897 | 128.02070205233477 | 1.2540654638115976 | approximate_stationary_analogue |
| Ownership, heads 30–55 | 0.648334 | 0.349344 | -0.299 | 2339.3623724673616 | 209.12728169984143 | exact_stationary_key |
| First-birth room response, −1 to +3 | 0.720246 | 1.0784 | 0.3582 | 137.5652749002964 | 17.64632680436442 | approximate_stationary_analogue |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.220897 | -0.1262 | 280.52808370152104 | 4.465703392969249 | approximate_stationary_analogue |
| Recent-parent ownership gap | 0.162896 | 0.418965 | 0.2561 | 27055.822957508266 | 1774.0940965968305 | exact_stationary_key |

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
| psi_child | 0.18896454850136238 | None | None | False | normalized to completed fertility 2.1 unless spec fixes psi |
| payroll_tax | 0.179 | None | None | False | externally fixed (P.tau_pay) |
| housing_supply_elasticity | None | None | None | False | not applicable: this is a dated-transition supply_rule.elasticity value (0.63 in the retained report); the bare stationary solve has no such object |
| tenure_choice_kappa | 0.005 | None | None | False | externally fixed or sandbox switch |
| alpha_cons | 0.733 | None | None | False | externally fixed or sandbox switch |
| sigma | 2.0 | None | None | False | externally fixed or sandbox switch |
| child_benefit_form | None | None | None | False | externally fixed or sandbox switch |
| child_benefit_curvature | None | None | None | False | externally fixed or sandbox switch |
| scale_weighting | None | None | None | False | externally fixed or sandbox switch |
| phi | [0.8, 0.8, 0.8, 0.8] | None | None | False | financed share by parity; solved P.phi (default 0.80, spec-overridable; see apply_spec's phi/n_parity broadcast workaround) |
| lambda_d | 0.0 | None | None | False | unsecured debt-line multiple of mean earnings (default 0.0, spec-overridable) |

