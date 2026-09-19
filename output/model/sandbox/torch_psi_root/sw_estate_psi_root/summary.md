# Stationary-state sandbox run: sw_estate_psi_root

Spec file: `sandbox/specs/sw_estate_psi_root.yaml`
Overrides applied: {"estate_receiver": "ages_45_65", "bequest_net_of_selling_cost": true}
Mechanism switches: (none; all at default-off)

## Timing
- Wall time per Bellman/root evaluation: 980.61 s (mean over 8 evaluation(s))
- Total wall time: 7844.97 s
- Number of root evaluations (psi-normalization GE solves): 8
- Warm start: none available (no baseline run found). Bellman/value-function warm start is NOT implemented: run_model_cp_dt (solver.py:1064) has no V-init override hook analogous to p_init_override, so the value function is always cold-started inside the package; only the price root can be warm-started from here.
- psi-normalization status: derived_intercept
- Grid: full production Nb (Nb=120, J=17)

## Loss: 2582.498914

## 13-row target table
| Moment | Target | Model | Gap | Weight | Loss contribution | Source |
|---|---:|---:|---:|---:|---:|---|
| Initial model completed fertility | 2.1 | 2.1 | -4.804e-08 | - | - | exact_stationary_key |
| Childless women, ages 40–44 | 0.198279 | 0.192733 | -0.005546 | 35532.3042455214 | 1.0928759276784945 | approximate_stationary_analogue |
| Exactly one child among mothers, ages 40–44 | 0.213655 | nan | nan | 26952.820824310795 | - | unavailable_in_extract_moments |
| Period mean first-birth age | 25.9763 | 26.4383 | 0.4621 | 139.82806784479274 | 29.852731457392107 | approximate_stationary_analogue |
| First births at age 30+ | 0.249278 | 0.2518 | 0.002521 | 13866.065434728798 | 0.08815938133371161 | approximate_stationary_analogue |
| Wealth / annual gross labor earnings | 6.14586 | 5.08297 | -1.063 | 7.595098472533724 | 8.580423394196375 | exact_stationary_key |
| Annual bequests / aggregate wealth | 0.0088 | 0.00788441 | -0.0009156 | 5165289.256198346 | 4.330080401469346 | exact_stationary_key |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 4.06787 | 0.5519 | 10.616361531535473 | 3.234079621022616 | exact_stationary_key |
| Mean occupied rooms, capped at 9 | 5.5611 | 5.87584 | 0.3147 | 128.02070205233477 | 12.681817954490466 | approximate_stationary_analogue |
| Ownership, heads 30–55 | 0.648334 | 0.446436 | -0.2019 | 2339.3623724673616 | 95.35852739024055 | exact_stationary_key |
| First-birth room response, −1 to +3 | 0.720246 | 1.06121 | 0.341 | 137.5652749002964 | 15.992469506811885 | approximate_stationary_analogue |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.220076 | -0.127 | 280.52808370152104 | 4.523986075760041 | approximate_stationary_analogue |
| Recent-parent ownership gap | 0.162896 | 0.46115 | 0.2983 | 27055.822957508266 | 2406.7637633111067 | exact_stationary_key |

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
| psi_child | 0.18129970325273911 | None | None | False | normalized to completed fertility 2.1 unless spec fixes psi |
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

