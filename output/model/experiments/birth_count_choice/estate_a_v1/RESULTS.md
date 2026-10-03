# Estate-A fixed-parameter comparison

The old single and old multiple baselines isolate Estate-A within each menu; source paths and actual saved estate flags are recorded in comparison.json. Both experimental arms use post-saving net estates W=bp+(1-psi)*P*h, with no extra interest multiplier and no recipient mapping. The single arm caps intended births at one; the multiple arm caps at three. The new wealth target is 4.45838713455674; all other empirical values and weights, including the bequest target, are retained. This is a fixed-parameter GE comparison, not recalibration.

Old native fit tables remain diagnostic and unchanged. The authoritative common new-target comparison is below. Living old-age wealth continues to use beginning b+P*h; death-estate diagnostics use saving bp and chosen housing. The retained SCF bequest target has an outstanding wealth-scope/recipient mapping mismatch and remains provisional. Descriptive count hazards use pre-birth exposure, unlike baseline post-birth descriptions.

| Arm | Birth cap | Estate A | Price | Old diagnostic loss | Common new-target loss |
| --- | --- | --- | --- | --- | --- |
| baseline | 1.000 | False | 0.776 | 13.771 | 53.064 |
| baseline_multiple | 3.000 | False | 0.999 | 1672.491 | 1714.289 |
| single | 1.000 | True | 0.776 | 16.895 | 57.595 |
| multiple | 3.000 | True | 0.999 | 1685.273 | 1728.298 |

## baseline: complete target fit

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
| --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | 2.100 | 2.100 | -1.431e-10 | — | — | normalization |
| cps_childlessness | 0.198 | 0.199 | 0.001 | 35532.304 | 0.036 | scored |
| cps_exactly_one | 0.214 | 0.216 | 0.002 | 26952.821 | 0.109 | scored |
| nchs_mean_age | 25.976 | 25.968 | -0.009 | 139.828 | 0.010 | scored |
| nchs_share30 | 0.249 | 0.232 | -0.017 | 0.000 | 0.000 | validation |
| wealth_earnings | 4.458 | 6.741 | 2.282 | 7.595 | 39.556 | scored |
| bequest_wealth | 0.007 | 0.007 | -4.558e-04 | 5165289.256 | 1.073 | scored |
| old_dispersion | 3.516 | 2.904 | -0.612 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.886 | 0.157 | 128.021 | 3.154 | scored |
| ownership_30_55 | 0.676 | 0.671 | -0.005 | 2339.362 | 0.055 | scored |
| first_birth_rooms | 1.465 | 1.378 | -0.087 | 137.565 | 1.053 | scored |
| family_rooms | 0.385 | 0.299 | -0.086 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.124 | -0.004 | 27055.823 | 0.415 | scored |
| early_fertility | 0.810 | 0.534 | -0.276 | 100.000 | 7.602 | scored |

## baseline: all parameter restrictions

| Parameter | Estimate | Lower | Upper | Near bound | Status |
| --- | --- | --- | --- | --- | --- |
| H0 | 6.406 | 0.200 | 80.000 | False | supplied fixed housing supply coefficient; reference bounds advisory |
| beta_annual | 0.966 | 0.940 | 0.990 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.050 | 0.100 | 5.000 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 0.305 | 0.000 | 8.000 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.117 | 0.020 | 50.000 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.401 | 0.020 | 50.000 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.101 | 0.000 | 8.000 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.000 | — | — | — | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.063 | 0.000 | 0.800 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.014 | 0.001 | 0.100 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.179 | 0.010 | 0.500 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.168 | — | — | — | derived from supplied benefit and curvature |
| theta1 | 0.008 | — | — | — | fixed external restriction |
| sigma | 2.000 | — | — | — | fixed |
| alpha_cons | 0.733 | — | — | — | fixed CEX childless expenditure share |
| delta_alpha | 0.000 | — | — | — | fixed zero under experimental utility contract |
| h_P | 2.594 | 0.100 | 2.600 | True | supplied primitive; reference calibration bounds advisory |
| utility_reference_rent | 0.110 | — | — | — | retained inactive normalization; compensation off |
| q_annual | 0.020 | — | — | — | author-retained 2% annual real rate |
| financed_share | 0.800 | — | — | — | inherited credit contract |
| housing_supply_elasticity | 0.630 | — | — | — | fixed provisional external mapping |
| payroll_tax | 0.080 | — | — | — | derived from adopted pension ratio |
| pension_period | 0.918 | — | — | — | balanced PAYGO |
| annual_depreciation | 0.014 | — | — | — | adopted |
| period_depreciation | 0.055 | — | — | — | compounded |
| annual_property_tax | 0.011 | — | — | — | adopted |
| period_property_tax | 0.042 | — | — | — | linear period convention |
| selling_cost | 0.060 | — | — | — | retained |
| rental_cap | 6.000 | — | — | — | retained provisional |
| wealth_grid_nodes | 120.000 | — | — | — | retained exact grid |
| income_states | 9.000 | — | — | — | retained B15 |

## baseline_multiple: complete target fit

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
| --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | 2.100 | 2.100 | 3.127e-08 | — | — | normalization |
| cps_childlessness | 0.198 | 0.263 | 0.065 | 35532.304 | 148.479 | scored |
| cps_exactly_one | 0.214 | 0.185 | -0.029 | 26952.821 | 22.834 | scored |
| nchs_mean_age | 25.976 | 28.420 | 2.443 | 139.828 | 834.739 | scored |
| nchs_share30 | 0.249 | 0.390 | 0.141 | 0.000 | 0.000 | validation |
| wealth_earnings | 4.458 | 6.807 | 2.349 | 7.595 | 41.905 | scored |
| bequest_wealth | 0.007 | 0.007 | -4.955e-04 | 5165289.256 | 1.268 | scored |
| old_dispersion | 3.516 | 2.878 | -0.638 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 4.953 | -0.777 | 128.021 | 77.237 | scored |
| ownership_30_55 | 0.676 | 0.596 | -0.080 | 2339.362 | 15.048 | scored |
| first_birth_rooms | 1.465 | 1.529 | 0.064 | 137.565 | 0.572 | scored |
| family_rooms | 0.385 | 0.170 | -0.215 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | -0.016 | -0.143 | 27055.823 | 554.995 | scored |
| early_fertility | 0.810 | 0.395 | -0.415 | 100.000 | 17.210 | scored |

## baseline_multiple: all parameter restrictions

| Parameter | Estimate | Lower | Upper | Near bound | Status |
| --- | --- | --- | --- | --- | --- |
| H0 | 6.406 | 0.200 | 80.000 | False | supplied fixed housing supply coefficient; reference bounds advisory |
| beta_annual | 0.966 | 0.940 | 0.990 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.050 | 0.100 | 5.000 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 0.305 | 0.000 | 8.000 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.117 | 0.020 | 50.000 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.401 | 0.020 | 50.000 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.101 | 0.000 | 8.000 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.000 | — | — | — | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.063 | 0.000 | 0.800 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.014 | 0.001 | 0.100 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.179 | 0.010 | 0.500 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.168 | — | — | — | derived from supplied benefit and curvature |
| theta1 | 0.008 | — | — | — | fixed external restriction |
| sigma | 2.000 | — | — | — | fixed |
| alpha_cons | 0.733 | — | — | — | fixed CEX childless expenditure share |
| delta_alpha | 0.000 | — | — | — | fixed zero under experimental utility contract |
| h_P | 2.594 | 0.100 | 2.600 | True | supplied primitive; reference calibration bounds advisory |
| utility_reference_rent | 0.110 | — | — | — | retained inactive normalization; compensation off |
| q_annual | 0.020 | — | — | — | author-retained 2% annual real rate |
| financed_share | 0.800 | — | — | — | inherited credit contract |
| housing_supply_elasticity | 0.630 | — | — | — | fixed provisional external mapping |
| payroll_tax | 0.080 | — | — | — | derived from adopted pension ratio |
| pension_period | 0.918 | — | — | — | balanced PAYGO |
| annual_depreciation | 0.014 | — | — | — | adopted |
| period_depreciation | 0.055 | — | — | — | compounded |
| annual_property_tax | 0.011 | — | — | — | adopted |
| period_property_tax | 0.042 | — | — | — | linear period convention |
| selling_cost | 0.060 | — | — | — | retained |
| rental_cap | 6.000 | — | — | — | retained provisional |
| wealth_grid_nodes | 120.000 | — | — | — | retained exact grid |
| income_states | 9.000 | — | — | — | retained B15 |

## single: complete target fit

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
| --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | 2.100 | 2.100 | 1.759e-11 | — | — | normalization |
| cps_childlessness | 0.198 | 0.199 | 0.001 | 35532.304 | 0.036 | scored |
| cps_exactly_one | 0.214 | 0.216 | 0.002 | 26952.821 | 0.110 | scored |
| nchs_mean_age | 25.976 | 25.968 | -0.008 | 139.828 | 0.009 | scored |
| nchs_share30 | 0.249 | 0.232 | -0.017 | 0.000 | 0.000 | validation |
| wealth_earnings | 4.458 | 6.778 | 2.320 | 7.595 | 40.868 | scored |
| bequest_wealth | 0.007 | 0.006 | -9.297e-04 | 5165289.256 | 4.464 | scored |
| old_dispersion | 3.516 | 2.861 | -0.654 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.881 | 0.152 | 128.021 | 2.942 | scored |
| ownership_30_55 | 0.676 | 0.671 | -0.005 | 2339.362 | 0.053 | scored |
| first_birth_rooms | 1.465 | 1.378 | -0.087 | 137.565 | 1.034 | scored |
| family_rooms | 0.385 | 0.299 | -0.086 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.123 | -0.004 | 27055.823 | 0.473 | scored |
| early_fertility | 0.810 | 0.534 | -0.276 | 100.000 | 7.605 | scored |

## single: all parameter restrictions

| Parameter | Estimate | Lower | Upper | Near bound | Status |
| --- | --- | --- | --- | --- | --- |
| H0 | 6.406 | 0.200 | 80.000 | False | supplied fixed housing supply coefficient; reference bounds advisory |
| beta_annual | 0.966 | 0.940 | 0.990 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.050 | 0.100 | 5.000 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 0.305 | 0.000 | 8.000 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.117 | 0.020 | 50.000 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.401 | 0.020 | 50.000 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.101 | 0.000 | 8.000 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.000 | — | — | — | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.063 | 0.000 | 0.800 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.014 | 0.001 | 0.100 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.179 | 0.010 | 0.500 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.168 | — | — | — | derived from supplied benefit and curvature |
| theta1 | 0.008 | — | — | — | fixed external restriction |
| sigma | 2.000 | — | — | — | fixed |
| alpha_cons | 0.733 | — | — | — | fixed CEX childless expenditure share |
| delta_alpha | 0.000 | — | — | — | fixed zero under experimental utility contract |
| h_P | 2.594 | 0.100 | 2.600 | True | supplied primitive; reference calibration bounds advisory |
| utility_reference_rent | 0.110 | — | — | — | retained inactive normalization; compensation off |
| q_annual | 0.020 | — | — | — | author-retained 2% annual real rate |
| financed_share | 0.800 | — | — | — | inherited credit contract |
| housing_supply_elasticity | 0.630 | — | — | — | fixed provisional external mapping |
| payroll_tax | 0.080 | — | — | — | derived from adopted pension ratio |
| pension_period | 0.918 | — | — | — | balanced PAYGO |
| annual_depreciation | 0.014 | — | — | — | adopted |
| period_depreciation | 0.055 | — | — | — | compounded |
| annual_property_tax | 0.011 | — | — | — | adopted |
| period_property_tax | 0.042 | — | — | — | linear period convention |
| selling_cost | 0.060 | — | — | — | retained |
| rental_cap | 6.000 | — | — | — | retained provisional |
| wealth_grid_nodes | 120.000 | — | — | — | retained exact grid |
| income_states | 9.000 | — | — | — | retained B15 |

## multiple: complete target fit

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
| --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | 2.100 | 2.100 | 2.099e-08 | — | — | normalization |
| cps_childlessness | 0.198 | 0.263 | 0.065 | 35532.304 | 148.442 | scored |
| cps_exactly_one | 0.214 | 0.185 | -0.029 | 26952.821 | 22.782 | scored |
| nchs_mean_age | 25.976 | 28.422 | 2.446 | 139.828 | 836.393 | scored |
| nchs_share30 | 0.249 | 0.390 | 0.141 | 0.000 | 0.000 | validation |
| wealth_earnings | 4.458 | 6.840 | 2.382 | 7.595 | 43.082 | scored |
| bequest_wealth | 0.007 | 0.006 | -0.001 | 5165289.256 | 5.623 | scored |
| old_dispersion | 3.516 | 2.881 | -0.635 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 4.950 | -0.779 | 128.021 | 77.772 | scored |
| ownership_30_55 | 0.676 | 0.596 | -0.080 | 2339.362 | 14.999 | scored |
| first_birth_rooms | 1.465 | 1.531 | 0.066 | 137.565 | 0.597 | scored |
| family_rooms | 0.385 | 0.170 | -0.215 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | -0.016 | -0.144 | 27055.823 | 561.384 | scored |
| early_fertility | 0.810 | 0.395 | -0.415 | 100.000 | 17.224 | scored |

## multiple: all parameter restrictions

| Parameter | Estimate | Lower | Upper | Near bound | Status |
| --- | --- | --- | --- | --- | --- |
| H0 | 6.406 | 0.200 | 80.000 | False | supplied fixed housing supply coefficient; reference bounds advisory |
| beta_annual | 0.966 | 0.940 | 0.990 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.050 | 0.100 | 5.000 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 0.305 | 0.000 | 8.000 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.117 | 0.020 | 50.000 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.401 | 0.020 | 50.000 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.101 | 0.000 | 8.000 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.000 | — | — | — | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.063 | 0.000 | 0.800 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.014 | 0.001 | 0.100 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.179 | 0.010 | 0.500 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.168 | — | — | — | derived from supplied benefit and curvature |
| theta1 | 0.008 | — | — | — | fixed external restriction |
| sigma | 2.000 | — | — | — | fixed |
| alpha_cons | 0.733 | — | — | — | fixed CEX childless expenditure share |
| delta_alpha | 0.000 | — | — | — | fixed zero under experimental utility contract |
| h_P | 2.594 | 0.100 | 2.600 | True | supplied primitive; reference calibration bounds advisory |
| utility_reference_rent | 0.110 | — | — | — | retained inactive normalization; compensation off |
| q_annual | 0.020 | — | — | — | author-retained 2% annual real rate |
| financed_share | 0.800 | — | — | — | inherited credit contract |
| housing_supply_elasticity | 0.630 | — | — | — | fixed provisional external mapping |
| payroll_tax | 0.080 | — | — | — | derived from adopted pension ratio |
| pension_period | 0.918 | — | — | — | balanced PAYGO |
| annual_depreciation | 0.014 | — | — | — | adopted |
| period_depreciation | 0.055 | — | — | — | compounded |
| annual_property_tax | 0.011 | — | — | — | adopted |
| period_property_tax | 0.042 | — | — | — | linear period convention |
| selling_cost | 0.060 | — | — | — | retained |
| rental_cap | 6.000 | — | — | — | retained provisional |
| wealth_grid_nodes | 120.000 | — | — | — | retained exact grid |
| income_states | 9.000 | — | — | — | retained B15 |

## Estate diagnostics

| Arm | Measure | Estate definition | Value |
| --- | --- | --- | --- |
| baseline | annual_positive_death_estate_flow | gross_housing | 0.034 |
| baseline | annual_positive_death_estate_flow_to_living_gross_wealth | gross_housing | 0.007 |
| baseline | old65plus_signed_death_estate_flow | gross_housing | 0.034 |
| baseline | old65plus_positive_death_estate_flow | gross_housing | 0.034 |
| baseline | old65plus_negative_death_estate_flow | gross_housing | 0.000 |
| baseline | old65plus_mean_estate_per_death | gross_housing | 2.211 |
| baseline | old65plus_median_estate_per_death | gross_housing | 1.052 |
| baseline | annual_positive_death_estate_flow | net_selling_cost | 0.030 |
| baseline | annual_positive_death_estate_flow_to_living_gross_wealth | net_selling_cost | 0.006 |
| baseline | old65plus_signed_death_estate_flow | net_selling_cost | 0.030 |
| baseline | old65plus_positive_death_estate_flow | net_selling_cost | 0.030 |
| baseline | old65plus_negative_death_estate_flow | net_selling_cost | 0.000 |
| baseline | old65plus_mean_estate_per_death | net_selling_cost | 1.969 |
| baseline | old65plus_median_estate_per_death | net_selling_cost | 0.869 |
| baseline | ownership_age82 | actual_post_tenure_living_mass | 0.987 |
| baseline | mean_post_saving_bp_age82 | actual_post_tenure_stayer_corrected | -3.286 |
| baseline_multiple | annual_positive_death_estate_flow | gross_housing | 0.034 |
| baseline_multiple | annual_positive_death_estate_flow_to_living_gross_wealth | gross_housing | 0.007 |
| baseline_multiple | old65plus_signed_death_estate_flow | gross_housing | 0.034 |
| baseline_multiple | old65plus_positive_death_estate_flow | gross_housing | 0.034 |
| baseline_multiple | old65plus_negative_death_estate_flow | gross_housing | 0.000 |
| baseline_multiple | old65plus_mean_estate_per_death | gross_housing | 2.220 |
| baseline_multiple | old65plus_median_estate_per_death | gross_housing | 1.199 |
| baseline_multiple | annual_positive_death_estate_flow | net_selling_cost | 0.030 |
| baseline_multiple | annual_positive_death_estate_flow_to_living_gross_wealth | net_selling_cost | 0.006 |
| baseline_multiple | old65plus_signed_death_estate_flow | net_selling_cost | 0.030 |
| baseline_multiple | old65plus_positive_death_estate_flow | net_selling_cost | 0.030 |
| baseline_multiple | old65plus_negative_death_estate_flow | net_selling_cost | 0.000 |
| baseline_multiple | old65plus_mean_estate_per_death | net_selling_cost | 1.967 |
| baseline_multiple | old65plus_median_estate_per_death | net_selling_cost | 0.839 |
| baseline_multiple | ownership_age82 | actual_post_tenure_living_mass | 0.979 |
| baseline_multiple | mean_post_saving_bp_age82 | actual_post_tenure_stayer_corrected | -3.460 |
| single | annual_positive_death_estate_flow | gross_housing | 0.036 |
| single | annual_positive_death_estate_flow_to_living_gross_wealth | gross_housing | 0.007 |
| single | old65plus_signed_death_estate_flow | gross_housing | 0.036 |
| single | old65plus_positive_death_estate_flow | gross_housing | 0.036 |
| single | old65plus_negative_death_estate_flow | gross_housing | 0.000 |
| single | old65plus_mean_estate_per_death | gross_housing | 2.305 |
| single | old65plus_median_estate_per_death | gross_housing | 1.191 |
| single | annual_positive_death_estate_flow | net_selling_cost | 0.032 |
| single | annual_positive_death_estate_flow_to_living_gross_wealth | net_selling_cost | 0.006 |
| single | old65plus_signed_death_estate_flow | net_selling_cost | 0.032 |
| single | old65plus_positive_death_estate_flow | net_selling_cost | 0.032 |
| single | old65plus_negative_death_estate_flow | net_selling_cost | 0.000 |
| single | old65plus_mean_estate_per_death | net_selling_cost | 2.070 |
| single | old65plus_median_estate_per_death | net_selling_cost | 0.869 |
| single | ownership_age82 | actual_post_tenure_living_mass | 0.962 |
| single | mean_post_saving_bp_age82 | actual_post_tenure_stayer_corrected | -3.019 |
| multiple | annual_positive_death_estate_flow | gross_housing | 0.035 |
| multiple | annual_positive_death_estate_flow_to_living_gross_wealth | gross_housing | 0.007 |
| multiple | old65plus_signed_death_estate_flow | gross_housing | 0.035 |
| multiple | old65plus_positive_death_estate_flow | gross_housing | 0.035 |
| multiple | old65plus_negative_death_estate_flow | gross_housing | 0.000 |
| multiple | old65plus_mean_estate_per_death | gross_housing | 2.298 |
| multiple | old65plus_median_estate_per_death | gross_housing | 1.199 |
| multiple | annual_positive_death_estate_flow | net_selling_cost | 0.032 |
| multiple | annual_positive_death_estate_flow_to_living_gross_wealth | net_selling_cost | 0.006 |
| multiple | old65plus_signed_death_estate_flow | net_selling_cost | 0.032 |
| multiple | old65plus_positive_death_estate_flow | net_selling_cost | 0.032 |
| multiple | old65plus_negative_death_estate_flow | net_selling_cost | 0.000 |
| multiple | old65plus_mean_estate_per_death | net_selling_cost | 2.051 |
| multiple | old65plus_median_estate_per_death | net_selling_cost | 0.849 |
| multiple | ownership_age82 | actual_post_tenure_living_mass | 0.950 |
| multiple | mean_post_saving_bp_age82 | actual_post_tenure_stayer_corrected | -3.219 |

[Complete fit CSV](paired_target_fit_new_contract.csv), [complete parameters CSV](paired_parameters.csv), [estate diagnostics CSV](paired_estate_diagnostics.csv).


## Independent Claude-A comparison

Claude-A changed utility only. Its gross bequest-flow report is replaced here by the net flow recomputed from its saved distribution and saving policies; other saved model moments are retained. Both are rescored with the same new wealth target.
Claude price: 0.7761012530937389; Estate-A single price: 0.776101253093739.
[All moment comparisons](claude_a_comparison.csv), [Claude estate and age-82 diagnostics](claude_a_estate_diagnostics.csv).

| Moment | Target | Claude A | Estate A single | Model difference |
| --- | --- | --- | --- | --- |
| initial_normalization | 2.100 | 2.100 | 2.100 | -8.882e-16 |
| cps_childlessness | 0.198 | 0.199 | 0.199 | -5.551e-17 |
| cps_exactly_one | 0.214 | 0.216 | 0.216 | -8.327e-17 |
| nchs_mean_age | 25.976 | 25.968 | 25.968 | 0.000 |
| nchs_share30 | 0.249 | 0.232 | 0.232 | 2.776e-17 |
| wealth_earnings | 4.458 | 6.778 | 6.778 | 1.776e-15 |
| bequest_wealth | 0.007 | 0.006 | 0.006 | 6.939e-18 |
| old_dispersion | 3.516 | 2.861 | 2.861 | 0.000 |
| mean_rooms | 5.729 | 5.881 | 5.881 | -8.882e-16 |
| ownership_30_55 | 0.676 | 0.671 | 0.671 | 1.110e-16 |
| first_birth_rooms | 1.465 | 1.378 | 1.378 | -1.776e-15 |
| family_rooms | 0.385 | 0.299 | 0.299 | -2.665e-15 |
| recent_parent_ownership | 0.128 | 0.123 | 0.123 | 5.551e-16 |
| early_fertility | 0.810 | 0.534 | 0.534 | 0.000 |