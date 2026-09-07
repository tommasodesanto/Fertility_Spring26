# Retained-parameter comparison: simple fertility nests

Both full histories completed: new nests in38m10s; sequential exhaustive-saving control in22m06s. These are fixed-parameter objective evaluations, not newly estimated calibrations.

Loss: new **36.371664**; control **30.408528**. All12 targets, weights and11 original parameter values/bounds match exactly. The old-state fertility intercept is separately normalized to2.1 in each model.

All recorded market, measurement, accounting and population gates pass. Terminal budget-violating mass and occupied value drops are zero. The23 focused tests and both lifecycle smokes pass. No policy simulations or exact repeated full histories were run. The control shares exhaustive saving; interpolation-support and probability-storage differences remain, as documented in SPECIFICATION.md.

## Complete target fits

### New fertility nests

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.92285207 | 0.00485207429 | 1425.73899 | 0.0335656383 |
| childless_rate | 0.188 | 0.189413663 | 0.00141366274 | 17180.7438 | 0.0343347258 |
| mean_age_first_birth | 26.0446273 | 26.2573542 | 0.212726908 | 44.4444444 | 2.01123277 |
| share_first_births_age30plus | 0.260327402 | 0.237631466 | -0.0226959357 | 10000 | 5.15105496 |
| housing_increment_0to1 | 0.720246262 | 0.421117724 | -0.299128539 | 137.565275 | 12.3090495 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.367699559 | 0.419762442 | 0.0520628828 | 2958.51499 | 8.01918434 |
| own_family_gap | 0.16766167 | 0.161692962 | -0.00596870801 | 14229.591 | 0.506935941 |
| own_rate | 0.575472 | 0.53833479 | -0.03713721 | 1207.84609 | 1.66582795 |
| aggregate_mean_occupied_rooms_18_85 | 5.77997048 | 6.31039128 | 0.530420795 | 11.973159 | 3.36860303 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 6.93725645 | 0.0641564466 | 6.28766943 | 0.0258803595 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.00842986747 | -0.000370132531 | 5165289.26 | 0.707634765 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.44811075 | 3.2370087 | -0.211102051 | 56.9597722 | 2.53835961 |

### Sequential control

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.92290917 | 0.00490917034 | 1425.73899 | 0.0343602433 |
| childless_rate | 0.188 | 0.189408631 | 0.00140863149 | 17180.7438 | 0.0340907651 |
| mean_age_first_birth | 26.0446273 | 26.2560336 | 0.211406368 | 44.4444444 | 1.9863401 |
| share_first_births_age30plus | 0.260327402 | 0.237579234 | -0.0227481676 | 10000 | 5.17479129 |
| housing_increment_0to1 | 0.720246262 | 0.439707954 | -0.280538309 | 137.565275 | 10.8266269 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.367699559 | 0.404566529 | 0.0368669704 | 2958.51499 | 4.02113519 |
| own_family_gap | 0.16766167 | 0.162115997 | -0.00554567338 | 14229.591 | 0.437623859 |
| own_rate | 0.575472 | 0.544990679 | -0.0304813215 | 1207.84609 | 1.12222303 |
| aggregate_mean_occupied_rooms_18_85 | 5.77997048 | 6.31851039 | 0.538539912 | 11.973159 | 3.47251827 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 6.93741564 | 0.0643156362 | 6.28766943 | 0.0260089512 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.00842827094 | -0.000371729063 | 5165289.26 | 0.713752562 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.44811075 | 3.23614982 | -0.211960932 | 56.9597722 | 2.55905656 |

## Parameters and restrictions

All11 free coordinates below are retained starting values, not estimates obtained from these runs. Near-bound flags use the existing diagnostic.

| Parameter | Retained value | Lower | Upper | Near bound |
|---|---:|---:|---:|---|
| beta_annual | 0.995116569194981 | 0.94 | 0.9995 | False |
| kappa_fert | 2.1681730392479377 | 0.02 | 50.0 | False |
| kappa_fert_continuation | 1.7364706586958831 | 0.02 | 50.0 | False |
| chi | 1.0434717373613915 | 0.1 | 5.0 | False |
| H0 | 14.562959141565095 | 0.2 | 80.0 | False |
| theta0 | 0.5284284711333161 | 0.0 | 8.0 | False |
| theta1 | 0.10724930821495539 | 0.02 | 16.0 | True |
| hbar_child_rooms | 0.2822101230841891 | 0.1 | 1.8 | False |
| first_birth_fixed_cost | 4.5591384147193255 | 0.0 | 8.0 | False |
| hbar_first_child_jump | 0.3649311350610432 | 0.0 | 0.5 | False |
| psi_child_change_2023 | -0.32871390689556357 | -1.5 | 0.2 | False |

Housing choice scale is externally fixed at0.005; housing-supply elasticity at0.63.

| Derived normalization | New nests | Sequential |
|---|---:|---:|
| psi_child_2007 | 0.28733166342792216 | 0.28801683161268565 |
| psi_child_2023 | -0.04138224346764141 | -0.040697075282877915 |

Collection verification: [collection_verification.json](collection_verification.json). Complete source and launch contracts are in the adjacent README and plans.
