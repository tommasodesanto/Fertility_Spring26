# Overnight full calibration

Status: smoke_passed. Production remains retained.

Best evaluated loss: 20.358960195. First-birth rooms: 0.513865 against 0.720246.
All eleven parameters were searched; only the first-child upper bound was enlarged to 2.0.
All twelve empirical targets and weights remain unchanged. No global-optimum claim.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.9234125 | 0.0054124997 | 1425.739 | 0.041767243 |
| childless_rate | 0.188 | 0.18919259 | 0.0011925882 | 17180.744 | 0.0244356 |
| mean_age_first_birth | 26.044627 | 26.361366 | 0.31673905 | 44.444444 | 4.4588277 |
| share_first_births_age30plus | 0.2603274 | 0.2438651 | -0.016462299 | 10000 | 2.7100728 |
| housing_increment_0to1 | 0.72024626 | 0.513865 | -0.20638126 | 137.56527 | 5.8593487 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.36769956 | 0.36871816 | 0.0010186051 | 2958.515 | 0.0030696263 |
| own_family_gap | 0.16766167 | 0.17505849 | 0.0073968215 | 14229.591 | 0.77854316 |
| own_rate | 0.575472 | 0.55070765 | -0.024764351 | 1207.8461 | 0.7407395 |
| aggregate_mean_occupied_rooms_18_85 | 5.7799705 | 6.3465852 | 0.56661473 | 11.973159 | 3.8440097 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 7.0367787 | 0.16367875 | 6.2876694 | 0.16845127 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.0085001938 | -0.00029980617 | 5165289.3 | 0.46427552 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.4481108 | 3.2990604 | -0.14905039 | 56.959772 | 1.2654193 |

| Parameter | Estimate | Lower | Upper | Status | Near bound |
|---|---:|---:|---:|---|---|
| beta_annual | 0.9957294031845191 | 0.94 | 0.9995 | estimated_free_transition_parameter | False |
| kappa_fert | 2.238195772180536 | 0.02 | 50.0 | estimated_free_transition_parameter | False |
| kappa_fert_continuation | 1.7322301491336538 | 0.02 | 50.0 | estimated_free_transition_parameter | False |
| chi | 1.0447481691265725 | 0.1 | 5.0 | estimated_free_transition_parameter | False |
| H0 | 14.535717947625159 | 0.2 | 80.0 | estimated_free_transition_parameter | False |
| theta0 | 0.5557589336583452 | 0.0 | 8.0 | estimated_free_transition_parameter | False |
| theta1 | 0.10613495044093449 | 0.02 | 16.0 | estimated_free_transition_parameter | True |
| hbar_child_rooms | 0.23289841116220328 | 0.1 | 1.8 | estimated_free_transition_parameter | False |
| first_birth_fixed_cost | 4.540285143086607 | 0.0 | 8.0 | estimated_free_transition_parameter | False |
| hbar_first_child_jump | 0.6006848485093815 | 0.0 | 2.0 | estimated_free_transition_parameter | False |
| psi_child_change_2023 | -0.3273390564208588 | -1.5 | 0.2 | estimated_free_transition_parameter | False |
| psi_child_2007 | 0.30087533088716095 | nan | nan | externally_normalized_to_old_completed_fertility | False |
| psi_child_2023 | -0.026463725533697857 | nan | nan | derived_from_old_intercept_and_transition_coordinate | False |
| tenure_choice_kappa | 0.005 | nan | nan | externally_fixed_profile_not_estimated | False |
| housing_supply_elasticity | 0.63 | nan | nan | externally_fixed_profile_not_estimated | False |

Selected case: `/scratch/td2248/projects/Fertility_Spring26_full_joint_overnight_20260906b/output/model/e5f_full_joint_overnight_20260906a/recovery_01/smoke/joint_probe_loop/task_001`.
The standard seventeen diagnostic plots and complete ledgers accompany this result.
