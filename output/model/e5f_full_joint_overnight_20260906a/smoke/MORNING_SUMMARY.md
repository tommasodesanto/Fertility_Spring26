# Overnight full calibration

Status: smoke_passed. Production remains retained.

Best evaluated loss: 20.695274280. First-birth rooms: 0.515778 against 0.720246.
All eleven parameters were searched; only the first-child upper bound was enlarged to 2.0.
All twelve empirical targets and weights remain unchanged. No global-optimum claim.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.9238178 | 0.0058178405 | 1425.739 | 0.04825737 |
| childless_rate | 0.188 | 0.18910143 | 0.0011014327 | 17180.744 | 0.020842887 |
| mean_age_first_birth | 26.044627 | 26.349204 | 0.30457627 | 44.444444 | 4.1229646 |
| share_first_births_age30plus | 0.2603274 | 0.24314448 | -0.017182918 | 10000 | 2.9525267 |
| housing_increment_0to1 | 0.72024626 | 0.51577822 | -0.20446804 | 137.56527 | 5.751216 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.36769956 | 0.36754102 | -0.00015853784 | 2958.515 | 7.4360047e-05 |
| own_family_gap | 0.16766167 | 0.17572429 | 0.0080626166 | 14229.591 | 0.92500575 |
| own_rate | 0.575472 | 0.54496716 | -0.030504839 | 1207.8461 | 1.1239554 |
| aggregate_mean_occupied_rooms_18_85 | 5.7799705 | 6.3491146 | 0.56914409 | 11.973159 | 3.8784055 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 7.0338073 | 0.16070729 | 6.2876694 | 0.16239058 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.008497223 | -0.00030277696 | 5165289.3 | 0.47352216 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.4481108 | 3.3007964 | -0.14731432 | 56.959772 | 1.2361129 |

| Parameter | Estimate | Lower | Upper | Status | Near bound |
|---|---:|---:|---:|---|---|
| beta_annual | 0.9957387588331625 | 0.94 | 0.9995 | estimated_free_transition_parameter | False |
| kappa_fert | 2.23273003596061 | 0.02 | 50.0 | estimated_free_transition_parameter | False |
| kappa_fert_continuation | 1.7364706586958831 | 0.02 | 50.0 | estimated_free_transition_parameter | False |
| chi | 1.0434717373613915 | 0.1 | 5.0 | estimated_free_transition_parameter | False |
| H0 | 14.562959141565095 | 0.2 | 80.0 | estimated_free_transition_parameter | False |
| theta0 | 0.5544418580619153 | 0.0 | 8.0 | estimated_free_transition_parameter | False |
| theta1 | 0.10635689183849888 | 0.02 | 16.0 | estimated_free_transition_parameter | True |
| hbar_child_rooms | 0.23268814270366286 | 0.1 | 1.8 | estimated_free_transition_parameter | False |
| first_birth_fixed_cost | 4.544052672413152 | 0.0 | 8.0 | estimated_free_transition_parameter | False |
| hbar_first_child_jump | 0.5999999999999999 | 0.0 | 2.0 | estimated_free_transition_parameter | False |
| psi_child_change_2023 | -0.32688089718542696 | -1.5 | 0.2 | estimated_free_transition_parameter | False |
| psi_child_2007 | 0.2997752448535115 | nan | nan | externally_normalized_to_old_completed_fertility | False |
| psi_child_2023 | -0.027105652331915475 | nan | nan | derived_from_old_intercept_and_transition_coordinate | False |
| tenure_choice_kappa | 0.005 | nan | nan | externally_fixed_profile_not_estimated | False |
| housing_supply_elasticity | 0.63 | nan | nan | externally_fixed_profile_not_estimated | False |

Selected case: `/scratch/td2248/projects/Fertility_Spring26_full_joint_overnight_20260906a/output/model/e5f_full_joint_overnight_20260906a/smoke/exact_loop/task_002`.
The standard seventeen diagnostic plots and complete ledgers accompany this result.
