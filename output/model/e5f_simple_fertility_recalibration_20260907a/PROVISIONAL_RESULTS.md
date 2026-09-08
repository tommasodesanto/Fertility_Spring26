# Provisional recalibration readout

Best valid loss **26.249683**, versus starting nested36.371664 and sequential control30.408528. The best is13.68% below sequential; exact final reproductions are pending. Original12targets/11parameters,weights andbounds retained.

36 histories completed with valid receipts (including2 exact anchor smokes); one joint proposal failed the market gate. Job17145615 stopped after1h56m37s before its final repeats. The rejected case is joint005, residual3.588e-4 versus2e-4. It is not the selected candidate. Selected joint012 passes all recorded market/measurement/accounting/population gates and terminal budget/value checks. No policies, global optimum claim, or production promotion.

The reviewed follow-up uses only the2 already-budgeted repeats (39 total attempted histories including the failed case); it does not rerun the failed candidate, change gates or continue search.

## All target fits

| Moment | Target | Weight | Sequential model | Nested model | Sequential gap | Nested gap | Sequential loss | Nested loss |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1425.73899 | 1.92290917 | 1.9308403 | 0.00490917034 | 0.0128403013 | 0.0343602433 | 0.235066347 |
| childless_rate | 0.188 | 17180.7438 | 0.189408631 | 0.187271074 | 0.00140863149 | -0.000728925996 | 0.0340907651 | 0.009128698 |
| mean_age_first_birth | 26.0446273 | 44.4444444 | 26.2560336 | 26.1829092 | 0.211406368 | 0.138281969 | 1.9863401 | 0.849862355 |
| share_first_births_age30plus | 0.260327402 | 10000 | 0.237579234 | 0.233458064 | -0.0227481676 | -0.0268693381 | 5.17479129 | 7.21961331 |
| housing_increment_0to1 | 0.720246262 | 137.565275 | 0.439707954 | 0.457033541 | -0.280538309 | -0.263212722 | 10.8266269 | 9.53065113 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.367699559 | 2958.51499 | 0.404566529 | 0.389985585 | 0.0368669704 | 0.0222860261 | 4.02113519 | 1.46939664 |
| own_family_gap | 0.16766167 | 14229.591 | 0.162115997 | 0.165347902 | -0.00554567338 | -0.00231376798 | 0.437623859 | 0.0761784319 |
| own_rate | 0.575472 | 1207.84609 | 0.544990679 | 0.528988046 | -0.0304813215 | -0.0464839537 | 1.12222303 | 2.60986304 |
| aggregate_mean_occupied_rooms_18_85 | 5.77997048 | 11.973159 | 6.31851039 | 6.21376322 | 0.538539912 | 0.433792733 | 3.47251827 | 2.25306279 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 6.28766943 | 6.93741564 | 6.96360345 | 0.0643156362 | 0.0905034465 | 0.0260089512 | 0.0515015069 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 5165289.26 | 0.00842827094 | 0.00849303408 | -0.000371729063 | -0.000306965921 | 0.713752562 | 0.486715272 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.44811075 | 56.9597722 | 3.23614982 | 3.2880848 | -0.211960932 | -0.16002595 | 2.55905656 | 1.4586432 |

## All parameters and restrictions

| Model | Parameter | Value | Lower bound | Upper bound | Free | Status | Near bound |
|---|---|---:|---:|---:|---|---|---|
| Selected nested | beta_annual | 0.9952765791045264 | 0.94 | 0.9995 | True | estimated_free_transition_parameter | False |
| Selected nested | kappa_fert | 2.1681730392479377 | 0.02 | 50.0 | True | estimated_free_transition_parameter | False |
| Selected nested | kappa_fert_continuation | 1.8057480287692833 | 0.02 | 50.0 | True | estimated_free_transition_parameter | False |
| Selected nested | chi | 1.0434717373613915 | 0.1 | 5.0 | True | estimated_free_transition_parameter | False |
| Selected nested | H0 | 13.71604910051014 | 0.2 | 80.0 | True | estimated_free_transition_parameter | False |
| Selected nested | theta0 | 0.5491891806761955 | 0.0 | 8.0 | True | estimated_free_transition_parameter | False |
| Selected nested | theta1 | 0.10372395059042462 | 0.02 | 16.0 | True | estimated_free_transition_parameter | True |
| Selected nested | hbar_child_rooms | 0.24887678975085575 | 0.1 | 1.8 | True | estimated_free_transition_parameter | False |
| Selected nested | first_birth_fixed_cost | 4.619731383944025 | 0.0 | 8.0 | True | estimated_free_transition_parameter | False |
| Selected nested | hbar_first_child_jump | 0.4649311350610432 | 0.0 | 0.5 | True | estimated_free_transition_parameter | False |
| Selected nested | psi_child_change_2023 | -0.3213877962199821 | -1.5 | 0.2 | True | estimated_free_transition_parameter | False |
| Selected nested | psi_child_2007 | 0.2827242498643133 | nan | nan | False | externally_normalized_to_old_completed_fertility | False |
| Selected nested | psi_child_2023 | -0.038663546355668765 | nan | nan | False | derived_from_old_intercept_and_transition_coordinate | False |
| Selected nested | tenure_choice_kappa | 0.005 | nan | nan | False | externally_fixed_profile_not_estimated | False |
| Selected nested | housing_supply_elasticity | 0.63 | nan | nan | False | externally_fixed_profile_not_estimated | False |
| Sequential control | beta_annual | 0.995116569194981 | 0.94 | 0.9995 | True | estimated_free_transition_parameter | False |
| Sequential control | kappa_fert | 2.1681730392479377 | 0.02 | 50.0 | True | estimated_free_transition_parameter | False |
| Sequential control | kappa_fert_continuation | 1.7364706586958831 | 0.02 | 50.0 | True | estimated_free_transition_parameter | False |
| Sequential control | chi | 1.0434717373613915 | 0.1 | 5.0 | True | estimated_free_transition_parameter | False |
| Sequential control | H0 | 14.562959141565095 | 0.2 | 80.0 | True | estimated_free_transition_parameter | False |
| Sequential control | theta0 | 0.5284284711333161 | 0.0 | 8.0 | True | estimated_free_transition_parameter | False |
| Sequential control | theta1 | 0.10724930821495539 | 0.02 | 16.0 | True | estimated_free_transition_parameter | True |
| Sequential control | hbar_child_rooms | 0.2822101230841891 | 0.1 | 1.8 | True | estimated_free_transition_parameter | False |
| Sequential control | first_birth_fixed_cost | 4.5591384147193255 | 0.0 | 8.0 | True | estimated_free_transition_parameter | False |
| Sequential control | hbar_first_child_jump | 0.3649311350610432 | 0.0 | 0.5 | True | estimated_free_transition_parameter | False |
| Sequential control | psi_child_change_2023 | -0.32871390689556357 | -1.5 | 0.2 | True | estimated_free_transition_parameter | False |
| Sequential control | psi_child_2007 | 0.28801683161268565 | nan | nan | False | externally_normalized_to_old_completed_fertility | False |
| Sequential control | psi_child_2023 | -0.040697075282877915 | nan | nan | False | derived_from_old_intercept_and_transition_coordinate | False |
| Sequential control | tenure_choice_kappa | 0.005 | nan | nan | False | externally_fixed_profile_not_estimated | False |
| Sequential control | housing_supply_elasticity | 0.63 | nan | nan | False | externally_fixed_profile_not_estimated | False |

The first-child jump0.4649 lies0.0351 below the0.5upper bound, although the existing generic near-bound flag is false. Theta1 is flagged near its lower bound. The derived old-state fertility intercept is re-normalized to2.1 in each evaluation.

Collection receipt:search/lead_collection_verification.json. Original sequential control shares exhaustive saving but differs in strict interpolation support/probability storage as disclosed in the specification.
