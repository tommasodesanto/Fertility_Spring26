# Stationary-state sandbox run: earnings_bgm_psi_fixed

Spec file: `sandbox/specs/earnings_bgm_psi_fixed.yaml`
Overrides applied: {"use_income_types": true, "income_type_transition": "markov", "income_shock_persistence": 0.8635910556159999, "z_grid": [0.2765725263626489, 0.48618547702709247, 0.8546630469076992, 1.5024079456590764, 2.6410755014464775], "z_weights": [0.0625, 0.25, 0.375, 0.25, 0.0625], "Pi_z": [[0.7538457431993805, 0.22071645331912, 0.024233614704291723, 0.0011825490685341885, 2.163970867360049e-05], [0.055179113329780004, 0.7659625505515264, 0.16642425179074063, 0.012138447060819462, 0.0002956372671335472], [0.00403893578404862, 0.11094950119382709, 0.7700231260442486, 0.1109495011938271, 0.00403893578404862], [0.0002956372671335472, 0.012138447060819462, 0.16642425179074066, 0.7659625505515264, 0.055179113329780004], [2.163970867360049e-05, 0.0011825490685341885, 0.024233614704291723, 0.22071645331912, 0.7538457431993805]], "permanent_income_levels_enabled": false, "permanent_income_log_variance": 0.0}
Mechanism switches: (none; all at default-off)

## Timing
- Wall time per Bellman/root evaluation: 79.75 s (mean over 1 evaluation(s))
- Total wall time: 79.80 s
- Number of root evaluations (psi-normalization GE solves): 1
- Warm start: price root warm-started from /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/sandbox/baseline/parameters.csv (saved price [0.7931310188535463])
- psi-normalization status: fixed_intercept
- Grid: full production Nb (Nb=120, J=17)

## Loss: 2096.247478

## 13-row target table
| Moment | Target | Model | Gap | Weight | Loss contribution | Source |
|---|---:|---:|---:|---:|---:|---|
| Initial model completed fertility | 2.1 | 2.04913 | -0.05087 | - | - | exact_stationary_key |
| Childless women, ages 40–44 | 0.198279 | 0.154116 | -0.04416 | 35532.3042455214 | 69.30105749557528 | approximate_stationary_analogue |
| Exactly one child among mothers, ages 40–44 | 0.213655 | nan | nan | 26952.820824310795 | - | unavailable_in_extract_moments |
| Period mean first-birth age | 25.9763 | 27.1689 | 1.193 | 139.82806784479274 | 198.8905560069046 | approximate_stationary_analogue |
| First births at age 30+ | 0.249278 | 0.298041 | 0.04876 | 13866.065434728798 | 32.97122010998203 | approximate_stationary_analogue |
| Wealth / annual gross labor earnings | 6.14586 | 5.48152 | -0.6643 | 7.595098472533724 | 3.352055271558478 | exact_stationary_key |
| Annual bequests / aggregate wealth | 0.0088 | 0.00795241 | -0.0008476 | 5165289.256198346 | 3.7108077296689785 | exact_stationary_key |
| Old wealth/income p90 / median, ages 76–84 | 3.51594 | 2.92005 | -0.5959 | 10.616361531535473 | 3.769670837443433 | exact_stationary_key |
| Mean occupied rooms, capped at 9 | 5.5611 | 6.08682 | 0.5257 | 128.02070205233477 | 35.382786644921396 | approximate_stationary_analogue |
| Ownership, heads 30–55 | 0.648334 | 0.49908 | -0.1493 | 2339.3623724673616 | 52.11346537593633 | exact_stationary_key |
| First-birth room response, −1 to +3 | 0.720246 | 1.2297 | 0.5095 | 137.5652749002964 | 35.70385796856191 | approximate_stationary_analogue |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.347067 | 0.299549 | -0.04752 | 280.52808370152104 | 0.633413538871835 | approximate_stationary_analogue |
| Recent-parent ownership gap | 0.162896 | 0.410625 | 0.2477 | 27055.822957508266 | 1660.4185872536855 | exact_stationary_key |

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
| child_benefit_form | None | None | None | False | externally fixed or sandbox switch |
| child_benefit_curvature | None | None | None | False | externally fixed or sandbox switch |
| scale_weighting | None | None | None | False | externally fixed or sandbox switch |
| phi | [0.8, 0.8, 0.8, 0.8] | None | None | False | financed share by parity; solved P.phi (default 0.80, spec-overridable; see apply_spec's phi/n_parity broadcast workaround) |
| lambda_d | 0.0 | None | None | False | unsecured debt-line multiple of mean earnings (default 0.0, spec-overridable) |

