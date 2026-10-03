# Joint birth-count experiment: complete fit comparison

Baseline loss: **13.771131**. Experimental loss: **1672.491299**.
Both cases passed the inherited stationary-GE acceptance and exact-repeat gates.
The ten supplied parameter coordinates are unchanged; this is a specification experiment without recalibration.

Baseline: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/local_solution/cases/20261003T175652812716Z_b1c72f13`
Experiment: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/experiments/birth_count_choice/current_params_v1/cases/20261003T203001912131Z_812b878e`
Target fingerprint: `db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`
Weight fingerprint: `2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`
Baseline provenance check: baseline predates parameter files; recomputed complete saved tables against immutable canonical snapshot.
Experiment provenance check: saved parameter-file provenance plus recomputed complete tables.

The economic change is the joint intended-count menu with existing-age Binomial success probabilities. Birth-order flows count children; recent-parent ownership weights successful households once. Existing target definitions, within-period interpolation, entry, post-interest timing and soft credit are retained.

## All 14 target rows

Gap means model minus target. Blank normalization weights and zero-weight validation rows remain in the table. CSV values retain full precision.

| Moment | Target | Baseline | Experiment | Baseline gap | Experiment gap | Weight | Baseline loss | Experiment loss |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | 2.100 | 2.100 | 2.100 | -1.431e-10 | 3.127e-08 | — | — | — |
| cps_childlessness | 0.198 | 0.199 | 0.263 | 0.001 | 0.065 | 35532.304 | 0.036 | 148.479 |
| cps_exactly_one | 0.214 | 0.216 | 0.185 | 0.002 | -0.029 | 26952.821 | 0.109 | 22.834 |
| nchs_mean_age | 25.976 | 25.968 | 28.420 | -0.009 | 2.443 | 139.828 | 0.010 | 834.739 |
| nchs_share30 | 0.249 | 0.232 | 0.390 | -0.017 | 0.141 | 0.000 | 0.000 | 0.000 |
| wealth_earnings | 6.927 | 6.741 | 6.807 | -0.186 | -0.119 | 7.595 | 0.263 | 0.108 |
| bequest_wealth | 0.007 | 0.007 | 0.007 | -4.558e-04 | -4.955e-04 | 5165289.256 | 1.073 | 1.268 |
| old_dispersion | 3.516 | 2.904 | 2.878 | -0.612 | -0.638 | 0.000 | 0.000 | 0.000 |
| mean_rooms | 5.729 | 5.886 | 4.953 | 0.157 | -0.777 | 128.021 | 3.154 | 77.237 |
| ownership_30_55 | 0.676 | 0.671 | 0.596 | -0.005 | -0.080 | 2339.362 | 0.055 | 15.048 |
| first_birth_rooms | 1.465 | 1.378 | 1.529 | -0.087 | 0.064 | 137.565 | 1.053 | 0.572 |
| family_rooms | 0.385 | 0.299 | 0.170 | -0.086 | -0.215 | 0.000 | 0.000 | 0.000 |
| recent_parent_ownership | 0.128 | 0.124 | -0.016 | -0.004 | -0.143 | 27055.823 | 0.415 | 554.995 |
| early_fertility | 0.810 | 0.534 | 0.395 | -0.276 | -0.415 | 100.000 | 7.602 | 17.210 |

## All 31 parameter records

Restrictions are the inherited reference bounds; they are advisory for this fixed-input GE. Near-bound indicators use the inherited one-percent-of-range screen.

| Parameter | Baseline | Experiment | Lower | Upper | Near bound: baseline | Near bound: experiment |
| --- | --- | --- | --- | --- | --- | --- |
| H0 | 6.406 | 6.406 | 0.200 | 80.000 | False | False |
| beta_annual | 0.966 | 0.966 | 0.940 | 0.990 | False | False |
| chi | 1.050 | 1.050 | 0.100 | 5.000 | False | False |
| first_birth_fixed_cost | 0.305 | 0.305 | 0.000 | 8.000 | False | False |
| kappa_fert | 0.117 | 0.117 | 0.020 | 50.000 | True | True |
| kappa_fert_continuation | 0.401 | 0.401 | 0.020 | 50.000 | True | True |
| theta0 | 0.101 | 0.101 | 0.000 | 8.000 | False | False |
| delta_alpha_jump | 0.000 | 0.000 | — | — | — | — |
| child_benefit_curvature | 0.063 | 0.063 | 0.000 | 0.800 | False | False |
| tenure_choice_kappa | 0.014 | 0.014 | 0.001 | 0.100 | False | False |
| psi_child | 0.179 | 0.179 | 0.010 | 0.500 | False | False |
| child_benefit_CRRA_coefficient | 0.168 | 0.168 | — | — | — | — |
| theta1 | 0.008 | 0.008 | — | — | — | — |
| sigma | 2.000 | 2.000 | — | — | — | — |
| alpha_cons | 0.733 | 0.733 | — | — | — | — |
| delta_alpha | 0.000 | 0.000 | — | — | — | — |
| h_P | 2.594 | 2.594 | 0.100 | 2.600 | True | True |
| utility_reference_rent | 0.110 | 0.110 | — | — | — | — |
| q_annual | 0.020 | 0.020 | — | — | — | — |
| financed_share | 0.800 | 0.800 | — | — | — | — |
| housing_supply_elasticity | 0.630 | 0.630 | — | — | — | — |
| payroll_tax | 0.080 | 0.080 | — | — | — | — |
| pension_period | 0.918 | 0.918 | — | — | — | — |
| annual_depreciation | 0.014 | 0.014 | — | — | — | — |
| period_depreciation | 0.055 | 0.055 | — | — | — | — |
| annual_property_tax | 0.011 | 0.011 | — | — | — | — |
| period_property_tax | 0.042 | 0.042 | — | — | — | — |
| selling_cost | 0.060 | 0.060 | — | — | — | — |
| rental_cap | 6.000 | 6.000 | — | — | — | — |
| wealth_grid_nodes | 120.000 | 120.000 | — | — | — | — |
| income_states | 9.000 | 9.000 | — | — | — | — |

## Convergence and accounting

| Object | Baseline | Experiment |
| --- | --- | --- |
| status | converged | converged |
| price | 0.776 | 0.999 |
| closure_mode | fixed_h0 | fixed_h0 |
| renewal_residual | -6.821e-11 | 1.489e-08 |
| absolute_housing_residual | 0.000 | 0.000 |
| actual_paygo_residual | 2.847e-14 | 1.071e-13 |
| population_scale | 1.000 | 1.394 |
| fixed_h0_population_scale | 1.000 | 1.394 |
| implied_H0_at_population_one | 6.406 | 4.597 |
| adjusted_births_per_normalized_household | 0.130 | 0.130 |
| actual_entry_per_normalized_household | 0.062 | 0.062 |
| normalized_housing_demand | 5.886 | 4.953 |
| physical_housing_supply | 5.886 | 6.902 |
| standard_plot_count | 17.000 | 17.000 |
| target_fit_rows | 14.000 | 14.000 |
| parameter_rows | 31.000 | 31.000 |
| native_population_step.renewal_adjusted_distribution_l1 | 4.314e-15 | 1.727e-14 |
| native_population_step.mass_residual | -5.551e-16 | -6.661e-16 |
| native_population_step.entry_gap | -4.211e-12 | 1.281e-09 |

## Inspection and interpretation

The standard packet is at `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/experiments/birth_count_choice/current_params_v1/cases/20261003T203001912131Z_812b878e/standard_diagnostics`; cached policy and aggregate plots are in that case.
The descriptive arrays attempt_hazard_by_age, first_birth_hazard_by_age, fert_by_age and first_birth_age_distribution use pre-birth exposure rather than the older post-birth exposure. The active target observers are separate and retain their definitions.
The all-zero dead-menu identity for less than 1e-12 occupied household mass preserves the inherited numerical convention; it changes no economic parameter or target. The triggering failed witness was 2.97e-38 household mass.
This first converged experimental steady state does not establish an accepted calibration, identification, grid adequacy, transition readiness or paper adoption.

[Full fit CSV](comparison.csv) · [Full parameter CSV](parameters_comparison.csv)
