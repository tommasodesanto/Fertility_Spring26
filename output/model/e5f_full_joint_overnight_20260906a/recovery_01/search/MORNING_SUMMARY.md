# Overnight full calibration

Status: failed. Production remains retained.

Best evaluated loss: 19.284439901. First-birth rooms: 0.508447 against 0.720246.
All eleven parameters were searched; only the first-child upper bound was enlarged to 2.0.
All twelve empirical targets and weights remain unchanged. No global-optimum claim.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.9276085 | 0.0096085272 | 1425.739 | 0.13162963 |
| childless_rate | 0.188 | 0.18752607 | -0.0004739325 | 17180.744 | 0.0038590015 |
| mean_age_first_birth | 26.044627 | 26.325894 | 0.28126676 | 44.444444 | 3.516044 |
| share_first_births_age30plus | 0.2603274 | 0.24194081 | -0.018386596 | 10000 | 3.3806692 |
| housing_increment_0to1 | 0.72024626 | 0.50844732 | -0.21179895 | 137.56527 | 6.1710123 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.36769956 | 0.36419295 | -0.0035066042 | 2958.515 | 0.036378708 |
| own_family_gap | 0.16766167 | 0.17267515 | 0.0050134836 | 14229.591 | 0.35766102 |
| own_rate | 0.575472 | 0.57198545 | -0.0034865491 | 1207.8461 | 0.014682607 |
| aggregate_mean_occupied_rooms_18_85 | 5.7799705 | 6.3361231 | 0.55615264 | 11.973159 | 3.703367 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 7.0515164 | 0.17841643 | 6.2876694 | 0.20015176 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.0085325526 | -0.0002674474 | 5165289.3 | 0.36946339 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.4481108 | 3.2913615 | -0.1567493 | 56.959772 | 1.3995212 |

| Parameter | Estimate | Lower | Upper | Status | Near bound |
|---|---:|---:|---:|---|---|
| beta_annual | 0.9956918643790077 | 0.94 | 0.9995 | estimated_free_transition_parameter | False |
| kappa_fert | 2.238195772180536 | 0.02 | 50.0 | estimated_free_transition_parameter | False |
| kappa_fert_continuation | 1.766446315263779 | 0.02 | 50.0 | estimated_free_transition_parameter | False |
| chi | 1.049869529311166 | 0.1 | 5.0 | estimated_free_transition_parameter | False |
| H0 | 14.427261788151577 | 0.2 | 80.0 | estimated_free_transition_parameter | False |
| theta0 | 0.5663517884297847 | 0.0 | 8.0 | estimated_free_transition_parameter | False |
| theta1 | 0.10437601132166745 | 0.02 | 16.0 | estimated_free_transition_parameter | True |
| hbar_child_rooms | 0.2337413868034396 | 0.1 | 1.8 | estimated_free_transition_parameter | False |
| first_birth_fixed_cost | 4.555364635392782 | 0.0 | 8.0 | estimated_free_transition_parameter | False |
| hbar_first_child_jump | 0.606177699084433 | 0.0 | 2.0 | estimated_free_transition_parameter | False |
| psi_child_change_2023 | -0.3236755132052146 | -1.5 | 0.2 | estimated_free_transition_parameter | False |
| psi_child_2007 | 0.29627807420322766 | nan | nan | externally_normalized_to_old_completed_fertility | False |
| psi_child_2023 | -0.027397439001986934 | nan | nan | derived_from_old_intercept_and_transition_coordinate | False |
| tenure_choice_kappa | 0.005 | nan | nan | externally_fixed_profile_not_estimated | False |
| housing_supply_elasticity | 0.63 | nan | nan | externally_fixed_profile_not_estimated | False |

Selected case: `/scratch/td2248/projects/Fertility_Spring26_full_joint_overnight_20260906b/output/model/e5f_full_joint_overnight_20260906a/recovery_01/search/polish_2_joint/task_002`.
The standard seventeen diagnostic plots and complete ledgers accompany this result.
