# Stationary-state sandbox run: cap_eight_psi_root

Spec file: `sandbox/specs/cap_eight_psi_root.yaml`
Overrides applied: {"hR_max": 8.0}
Mechanism switches: (none; all at default-off)

## Timing
- Wall time per Bellman/root evaluation: 91.85 s (mean over 8 evaluation(s))
- Total wall time: 735.61 s
- Number of root evaluations (psi-normalization GE solves): 8
- Warm start: price root warm-started from /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/sandbox/baseline/parameters.csv (saved price [0.7931310188535463])
- psi-normalization status: derived_intercept
- Grid: full production Nb (Nb=120, J=17)

## Loss: 489.334525

## 13-row target table
| Moment | Target | Model | Gap | Weight | Loss contribution | Source |
|---|---:|---:|---:|---:|---:|---|
| Initial model completed fertility | 2.1 | 2.1 | 3.148e-07 | - | - | exact_stationary_key |
| Childless women, ages 40–44 | 0.198279 | 0.20165 | 0.003371 | 35532.3042455214 | 0.40373038122102506 | approximate_stationary_analogue |
| Exactly one child among mothers, ages 40–44 | 0.213655 | nan | nan | 26952.820824310795 | - | unavailable_in_extract_moments |
| Period mean first-birth age | 25.9763 | 26.1256 | 0.1494 | 139.82806784479274 | 3.1202827572069847 | approximate_stationary_analogue |
| First births at age 30+ | 0.249278 | 0.234449 | -0.01483 | 13866.065434728798 | 3.049227579834438 | approximate_stationary_analogue |
| Wealth / annual gross labor earnings | 6.14586 | 5.02749 | -1.118 | 7.595098472533724 | 9.499681424402779 | exact_stationary_key |
| Annual bequests / aggregate wealth | 0.0088 | 0.0081573 | -0.0006427 | 5165289.256198346 | 2.1336072412413802 | exact_stationary_key |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.38238 | 0.8664 | 10.616361531535473 | 7.970028708459002 | exact_stationary_key |
| Mean occupied rooms, capped at 9 | 5.5611 | 5.8049 | 0.2438 | 128.02070205233477 | 7.609546513949581 | approximate_stationary_analogue |
| Ownership, heads 30–55 | 0.648334 | 0.30021 | -0.3481 | 2339.3623724673616 | 283.5078696946895 | exact_stationary_key |
| First-birth room response, −1 to +3 | 0.720246 | 1.35975 | 0.6395 | 137.5652749002964 | 56.26023217446421 | approximate_stationary_analogue |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.257961 | -0.08911 | 280.52808370152104 | 2.2273462535978132 | approximate_stationary_analogue |
| Recent-parent ownership gap | 0.162896 | 0.22768 | 0.06478 | 27055.822957508266 | 113.55297256910919 | exact_stationary_key |

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
| psi_child | 0.18938910648779358 | None | None | False | normalized to completed fertility 2.1 unless spec fixes psi |
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

