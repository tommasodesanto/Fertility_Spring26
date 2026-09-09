# Matched short perfect-foresight equilibria — September 9

Both 12-date paths (2007–2051, terminal state in 2055) now clear markets at the unchanged tolerance and reproduce exactly. These are conditional finite-horizon equilibria at inherited parameters, not new calibrations. Terminal-distance and horizon-stability requirements remain unsatisfied.

The sequential objective is **129.24310527**; the nested objective is **130.25388186**. The latter is about 0.78% higher in this short diagnostic. Both use the same twelve targets, weights and eleven free-coordinate values, with the maintained arm-specific old-state normalization. This does not establish which model will fit better after a proper horizon check and re-estimation.

Maximum absolute market gaps are **0.009987%** and **0.009617%**, respectively, against a **0.02%** tolerance. Every target, parameter, market and measurement file is identical between the selected evaluation and its fresh replay. The complete independent check is [short_equilibria_verification.json](meeting_receipts/short_equilibria_verification.json).

The shared fit problems remain visible: fertility is about 1.773 versus 1.918; childlessness about 23.47% versus 18.8%; the first-birth housing response is 0.455 rooms (sequential) or 0.440 (nested), versus 0.720. These are diagnostic values under the short endpoint. They are not replacements for the retained benchmark.

## Sequential: complete target fit

Unrounded source: [target table](meeting_receipts/historical_root_restart2_sequential_01/sequential/evaluation_003/target_fit.csv).

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.77338617 | -0.144613835 | 1425.73899 | 29.8167094 |
| childless_rate | 0.188 | 0.23471746 | 0.0467174599 | 17180.7438 | 37.4973352 |
| mean_age_first_birth | 26.0446273 | 26.2258076 | 0.181180342 | 44.4444444 | 1.4589474 |
| share_first_births_age30plus | 0.260327402 | 0.243685424 | -0.0166419775 | 10000 | 2.76955416 |
| housing_increment_0to1 | 0.720246262 | 0.455292684 | -0.264953578 | 137.565275 | 9.65713712 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.367699559 | 0.371750591 | 0.00405103204 | 2958.51499 | 0.048551777 |
| own_family_gap | 0.16766167 | 0.111681615 | -0.0559800552 | 14229.591 | 44.5922165 |
| own_rate | 0.575472 | 0.57496865 | -0.000503350432 | 1207.84609 | 0.000306021886 |
| aggregate_mean_occupied_rooms_18_85 | 5.77997048 | 6.21026559 | 0.43029511 | 11.973159 | 2.21687686 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 7.04021008 | 0.167110076 | 6.28766943 | 0.175588058 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.00836309516 | -0.000436904844 | 5165289.26 | 0.985980592 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.44811075 | 3.42762583 | -0.020484926 | 56.9597722 | 0.0239021542 |

### Parameters and restrictions

Unrounded source: [parameter table](meeting_receipts/historical_root_restart2_sequential_01/sequential/evaluation_003/parameters.csv). These are inherited coordinates, not estimates produced by this run.

| Parameter | Value | Lower | Upper | Free coordinate? | Near bound? | Restriction / status |
|---|---:|---:|---:|:---:|:---:|---|
| beta_annual | 0.995276579 | 0.94 | 0.9995 | True | False | fixed or bounded pilot transition parameter |
| kappa_fert | 2.16817304 | 0.02 | 50 | True | False | fixed or bounded pilot transition parameter |
| kappa_fert_continuation | 1.77077059 | 0.02 | 50 | True | False | fixed or bounded pilot transition parameter |
| chi | 1.05372702 | 0.1 | 5 | True | False | fixed or bounded pilot transition parameter |
| H0 | 13.7160491 | 0.2 | 80 | True | False | fixed or bounded pilot transition parameter |
| theta0 | 0.57034989 | 0 | 8 | True | False | fixed or bounded pilot transition parameter |
| theta1 | 0.103723951 | 0.02 | 16 | True | True | fixed or bounded pilot transition parameter |
| hbar_child_rooms | 0.247084905 | 0.1 | 1.8 | True | False | fixed or bounded pilot transition parameter |
| first_birth_fixed_cost | 4.61973138 | 0 | 8 | True | False | fixed or bounded pilot transition parameter |
| hbar_first_child_jump | 0.469765103 | 0 | 0.5 | True | False | fixed or bounded pilot transition parameter |
| psi_child_change_2023 | -0.321387796 | -1.5 | 0.2 | True | False | fixed or bounded pilot transition parameter |
| psi_child_2007 | 0.290051529 | — | — | False | False | externally normalized to old completed fertility |
| psi_child_2023 | -0.0313362669 | — | — | False | False | derived from old intercept and transition coordinate |
| tenure_choice_kappa | 0.005 | — | — | False | False | externally fixed profile not estimated |
| housing_supply_elasticity | 0.63 | — | — | False | False | externally fixed profile not estimated |

The near-bound flag is the inherited generic bound diagnostic; it does not establish weak identification by itself.

## Nested: complete target fit

Unrounded source: [target table](meeting_receipts/historical_root_restart_nested_01/nested/evaluation_006/target_fit.csv).

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.77326857 | -0.144731428 | 1425.73899 | 29.8652201 |
| childless_rate | 0.188 | 0.234741351 | 0.0467413514 | 17180.7438 | 37.5356976 |
| mean_age_first_birth | 26.0446273 | 26.2270614 | 0.182434159 | 44.4444444 | 1.47920989 |
| share_first_births_age30plus | 0.260327402 | 0.243735899 | -0.0165915025 | 10000 | 2.75277957 |
| housing_increment_0to1 | 0.720246262 | 0.440230374 | -0.280015888 | 137.565275 | 10.7863415 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.367699559 | 0.39078966 | 0.0230901011 | 2958.51499 | 1.57734046 |
| own_family_gap | 0.16766167 | 0.112803052 | -0.0548586178 | 14229.591 | 42.8234979 |
| own_rate | 0.575472 | 0.568998283 | -0.00647371699 | 1207.84609 | 0.0506196357 |
| aggregate_mean_occupied_rooms_18_85 | 5.77997048 | 6.20275025 | 0.422779769 | 11.973159 | 2.14011516 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 7.0386783 | 0.165578305 | 6.28766943 | 0.172383846 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.00836500028 | -0.000434999725 | 5165289.26 | 0.977400622 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.44811075 | 3.40764387 | -0.040466886 | 56.9597722 | 0.0932755494 |

### Parameters and restrictions

Unrounded source: [parameter table](meeting_receipts/historical_root_restart_nested_01/nested/evaluation_006/parameters.csv). These are inherited coordinates, not estimates produced by this run.

| Parameter | Value | Lower | Upper | Free coordinate? | Near bound? | Restriction / status |
|---|---:|---:|---:|:---:|:---:|---|
| beta_annual | 0.995276579 | 0.94 | 0.9995 | True | False | fixed or bounded pilot transition parameter |
| kappa_fert | 2.16817304 | 0.02 | 50 | True | False | fixed or bounded pilot transition parameter |
| kappa_fert_continuation | 1.77077059 | 0.02 | 50 | True | False | fixed or bounded pilot transition parameter |
| chi | 1.05372702 | 0.1 | 5 | True | False | fixed or bounded pilot transition parameter |
| H0 | 13.7160491 | 0.2 | 80 | True | False | fixed or bounded pilot transition parameter |
| theta0 | 0.57034989 | 0 | 8 | True | False | fixed or bounded pilot transition parameter |
| theta1 | 0.103723951 | 0.02 | 16 | True | True | fixed or bounded pilot transition parameter |
| hbar_child_rooms | 0.247084905 | 0.1 | 1.8 | True | False | fixed or bounded pilot transition parameter |
| first_birth_fixed_cost | 4.61973138 | 0 | 8 | True | False | fixed or bounded pilot transition parameter |
| hbar_first_child_jump | 0.469765103 | 0 | 0.5 | True | False | fixed or bounded pilot transition parameter |
| psi_child_change_2023 | -0.321387796 | -1.5 | 0.2 | True | False | fixed or bounded pilot transition parameter |
| psi_child_2007 | 0.289387626 | — | — | False | False | externally normalized to old completed fertility |
| psi_child_2023 | -0.0320001706 | — | — | False | False | derived from old intercept and transition coordinate |
| tenure_choice_kappa | 0.005 | — | — | False | False | externally fixed profile not estimated |
| housing_supply_elasticity | 0.63 | — | — | False | False | externally fixed profile not estimated |

The near-bound flag is the inherited generic bound diagnostic; it does not establish weak identification by itself.

## Next numerical step

Longer sequential root 17281784 finished six valid mappings with exact replay, but its maximum gap of 0.07102% remains above tolerance. Verified continuation **17289375** uses the same 28 dates, parameters and targets, with at most four fresh paths including final replay. Expected runtime is about 52 minutes; the internal limit is 70 minutes and the allocation 75 minutes.

The long-run preference continuation and empirical family-group/calendar alignment remain explicit outstanding production decisions. No new policy comparison or production promotion is claimed.
