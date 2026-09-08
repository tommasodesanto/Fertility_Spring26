# Overnight calibration: morning review, 8 September 2026

Job 17155429 completed at 09:11 EDT after 8h 58m 54s. Selected loss **23.791955**, 9.36% below the previous nested candidate and 21.76% below the sequential control. Both final exact repeats passed.

The accepted simultaneous fertility-nest specification, all twelve targets and weights, and all eleven free-parameter bounds are unchanged. No policy simulations were run; no production promotion occurred.

Ownership is now 57.301% against 57.547%. The first-birth housing response remains 0.443859 rooms against 0.720246: it is slightly above the sequential control (0.439708), but below the previous nested candidate (0.457034). Births after age 30 remain too infrequent (23.563% versus 26.033%). Mean rooms remain too high (6.216 versus 5.780).

There were 162 attempted full histories: 155 valid and seven rejected for market nonconvergence, with no gate relaxation. All 22 final sensitivity probes passed. The best diagnostic probe has loss 23.675506 (final_jacobian_9_plus); it has not received two final exact repeats and does not replace the frozen selection. Further local improvement remains possible.

The weighted residual Jacobian has numerical rank 11 and condition number 415.67 in normalized coordinates. This is a local sensitivity calculation, not proof of strong or global identification.

Comparison limitation: the sequential control shares exhaustive saving, but strict interpolation support and storage conventions differ. These are fixed-code benchmark comparisons, not a decomposition isolating the mathematical choice law.

Verification: independently checked the selected and two repeat summary/fit/parameter hashes, exact equality of all twelve repeat fits, and runtime reference receipts; independently recomputed every comparison loss contribution and both totals. Full checkpoint-array validation was performed by the run adapter, not rerun during this morning readout.

## Complete target comparison

| Moment | Target | Sequential | Nested | Nested gap | Weight | Sequential loss | Nested loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.9229092 | 1.9287128 | 0.010712794 | 1425.739 | 0.034360243 | 0.16362344 |
| childless_rate | 0.188 | 0.18940863 | 0.18905041 | 0.0010504128 | 17180.744 | 0.034090765 | 0.018956665 |
| mean_age_first_birth | 26.044627 | 26.256034 | 26.221748 | 0.17712038 | 44.444444 | 1.9863401 | 1.3942946 |
| share_first_births_age30plus | 0.2603274 | 0.23757923 | 0.23562918 | -0.02469822 | 10000 | 5.1747913 | 6.1000207 |
| housing_increment_0to1 | 0.72024626 | 0.43970795 | 0.44385861 | -0.27638765 | 137.56527 | 10.826627 | 10.50863 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.36769956 | 0.40456653 | 0.3834787 | 0.015779145 | 2958.515 | 4.0211352 | 0.73661526 |
| own_family_gap | 0.16766167 | 0.162116 | 0.16106146 | -0.0066002054 | 14229.591 | 0.43762386 | 0.61987956 |
| own_rate | 0.575472 | 0.54499068 | 0.57300538 | -0.0024666175 | 1207.8461 | 1.122223 | 0.0073487797 |
| aggregate_mean_occupied_rooms_18_85 | 5.7799705 | 6.3185104 | 6.2160316 | 0.43606107 | 11.973159 | 3.4725183 | 2.2766873 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 6.9374156 | 7.0042214 | 0.13112137 | 6.2876694 | 0.026008951 | 0.10810273 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.0084282709 | 0.0085465627 | -0.00025343731 | 5165289.3 | 0.71375256 | 0.33176896 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.4481108 | 3.2361498 | 3.2844302 | -0.16368056 | 56.959772 | 2.5590566 | 1.5260277 |

## Parameters and restrictions

| Parameter | Estimate | Lower | Upper | Free | Status | Near bound flag |
|---|---:|---:|---:|---|---|---|
| beta_annual | 0.9952765791045264 | 0.94 | 0.9995 | True | estimated_free_transition_parameter | False |
| kappa_fert | 2.1681730392479377 | 0.02 | 50.0 | True | estimated_free_transition_parameter | False |
| kappa_fert_continuation | 1.7707705862013274 | 0.02 | 50.0 | True | estimated_free_transition_parameter | False |
| chi | 1.0537270178404272 | 0.1 | 5.0 | True | estimated_free_transition_parameter | False |
| H0 | 13.71604910051014 | 0.2 | 80.0 | True | estimated_free_transition_parameter | False |
| theta0 | 0.570349890219075 | 0.0 | 8.0 | True | estimated_free_transition_parameter | False |
| theta1 | 0.10372395059042462 | 0.02 | 16.0 | True | estimated_free_transition_parameter | True |
| hbar_child_rooms | 0.24708490545163003 | 0.1 | 1.8 | True | estimated_free_transition_parameter | False |
| first_birth_fixed_cost | 4.619731383944025 | 0.0 | 8.0 | True | estimated_free_transition_parameter | False |
| hbar_first_child_jump | 0.4697651033802434 | 0.0 | 0.5 | True | estimated_free_transition_parameter | False |
| psi_child_change_2023 | -0.3213877962199821 | -1.5 | 0.2 | True | estimated_free_transition_parameter | False |
| psi_child_2007 | 0.28938762559349523 | nan | nan | False | externally_normalized_to_old_completed_fertility | False |
| psi_child_2023 | -0.03200017062648686 | nan | nan | False | derived_from_old_intercept_and_transition_coordinate | False |
| tenure_choice_kappa | 0.005 | nan | nan | False | externally_fixed_profile_not_estimated | False |
| housing_supply_elasticity | 0.63 | nan | nan | False | externally_fixed_profile_not_estimated | False |

The first-child housing jump is 0.469765 against an upper bound of 0.5 (6.05% of the bound span remaining). The generic near-bound flag is triggered for theta1. Old-state fertility remains externally normalized to 2.1.

## Search coverage

{
  "smoke": 2,
  "initial_population": 23,
  "de_generation_01": 22,
  "de_generation_02": 23,
  "polish_1_coordinate": 22,
  "polish_1_combined": 12,
  "polish_2_coordinate": 22,
  "polish_2_combined": 12,
  "final_wave_1": 23,
  "final_wave_2": 1
}
