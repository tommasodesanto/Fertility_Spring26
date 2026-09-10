# Horizon comparison: sequential perfect foresight

The 28-date path (2007–2115; terminal state in 2119) clears markets and reproduces exactly. The maximum market residual is 0.0067368%, below the 0.02% tolerance. Every full target, parameter, market and measurement file agrees exactly with the selected evaluation.

At unchanged parameters, targets and weights, extending from 12 to 28 dates lowers the objective from **129.24310527 to 94.52228097 (26.86%)**. Almost all the net improvement comes from the ownership gap between the model parent groups. The first-birth rooms response deteriorates; average rooms move farther above the target. The empirical/model family-group definitions remain outstanding, so the ownership improvement must not settle that separate measurement issue.

Neither horizon is certified as a sufficiently long approximation. The 28-date final population remains 38.16% from its stationary reference; the final asset-price gap is 12.27%. The historical moments therefore need a substantially longer continuation before re-estimation is treated as final. These are inherited-parameter evaluations, not new calibrations.

## Complete target fit

| Moment | Target | Short model | Long model | Long gap | Weight | Short loss | Long loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| tfr | 1.918 | 1.77338617 | 1.77708613 | -0.14091387 | 1425.73899 | 29.8167094 | 28.3104983 |
| childless_rate | 0.188 | 0.23471746 | 0.233601517 | 0.0456015169 | 17180.7438 | 37.4973352 | 35.7273283 |
| mean_age_first_birth | 26.0446273 | 26.2258076 | 26.2332241 | 0.188596796 | 44.4444444 | 1.4589474 | 1.58083339 |
| share_first_births_age30plus | 0.260327402 | 0.243685424 | 0.24409157 | -0.0162358317 | 10000 | 2.76955416 | 2.6360223 |
| housing_increment_0to1 | 0.720246262 | 0.455292684 | 0.41155194 | -0.308694322 | 137.565275 | 9.65713712 | 13.1088956 |
| prime30_55_parent_3plus_minus_1to2_mean_rooms | 0.367699559 | 0.371750591 | 0.347374623 | -0.0203249356 | 2958.51499 | 0.048551777 | 1.22217144 |
| own_family_gap | 0.16766167 | 0.111681615 | 0.146503259 | -0.0211584106 | 14229.591 | 44.5922165 | 6.37027966 |
| own_rate | 0.575472 | 0.57496865 | 0.587482312 | 0.0120103118 | 1207.84609 | 0.000306021886 | 0.174228885 |
| aggregate_mean_occupied_rooms_18_85 | 5.77997048 | 6.21026559 | 6.32572401 | 0.545753529 | 11.973159 | 2.21687686 | 3.56616847 |
| aggregate_wealth_to_annual_gross_labor_earnings | 6.8731 | 7.04021008 | 7.0573043 | 0.184204299 | 6.28766943 | 0.175588058 | 0.213348319 |
| annual_bequest_flow_to_aggregate_wealth | 0.0088 | 0.00836309516 | 0.00835427886 | -0.000445721136 | 5165289.26 | 0.985980592 | 1.02617423 |
| old_total_wealth_to_annual_income_p90_p50_7684 | 3.44811075 | 3.42762583 | 3.34665242 | -0.101458333 | 56.9597722 | 0.0239021542 | 0.586332124 |

## All parameter values and restrictions

Parameters are exactly identical across the two sequential horizons. No coordinate was re-estimated.

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

## Next bounded diagnostic

Job **17300115** is submitted for one 100-date prescribed-price path (2007–2403; terminal state in 2407), after the startup regression suite. The earlier cached-policy sizing check first reached the terminal distance thresholds at 95 total dates, which motivates this larger horizon; it does not guarantee an equilibrium endpoint by that date.

This first path tests the longer calculation and provides a numerical starting point. Its supplied prices are not yet market-clearing. It preserves the existing seed rule, all parameters, targets, weights, accounting gates and the diagnostic flat preference continuation after 2023. Therefore its fit alone cannot be compared as an equilibrium horizon result.

Run size: 200 Bellman calls, one core, 16 GB, expected about 47–50 minutes from measured 28-date runtimes. The internal limit is 60 minutes; allocation 65 minutes. Progress and complete tables are saved. No price panel, parameter search or policy run has been submitted for this horizon.

The isolated wrapper changes only the maximum date count (40 to 100) and maximum single-probe budget (30 to 60 minutes). Old source D is unchanged. Local tests: 51 passed; four history tests could not run under system Python 3.9 because maintained code uses Python 3.10+ zip(strict=True). The cluster launcher requires all 55 tests to pass under its supported Python before model work.

Unrounded [long fit](meeting_receipts/historical_root_long_restart_01/sequential/evaluation_003/target_fit.csv), [parameters](meeting_receipts/historical_root_long_restart_01/sequential/evaluation_003/parameters.csv), and [verification](meeting_receipts/historical_root_long_restart_01/sequential/verification.json).
