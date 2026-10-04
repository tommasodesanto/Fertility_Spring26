# Overnight normalized CES-share calibration results

Verified October 4, 2026 (New York). Array 19133352: all four Slurm tasks COMPLETED, exit 0:0. All four selected candidates passed fresh native GE, exact-repeat arrays and tables, 14 experimental target rows, 31 parameter rows and 17 identical diagnostic PNG hashes. Best chain: 1. No experimental specification or calibration has been adopted.

Experimental differences from post-interest soft chain 13: normalized CES-limit denominator in all family states, free parenthood jump and per-child slope replacing the physical housing floor (h_P=0), and family_rooms promoted to scored. No r*, alpha0 numerator, added birth shock or utility-cost rescaling. Earnings, entry distributions, timing, gross estates, old wealth target and other fixed inputs retained. The national family_rooms target uses inherited 42-metro bootstrap weight; the current-child model observer is a proxy. Eleven scored moments for eleven free coordinates does not certify identification rank.

All searches stopped at the minimum-native-GE start reserve. These are budget-limited local searches, not optimizer convergence or proof of best attainable fit.

| Chain | Completed objective evaluations | Verified loss | Price | Derived H0 |
|---|---:|---:|---:|---:|
| 0 | 82 | 2444.25 | 0.658236 | 7.23323 |
| 1 | 94 | 2336.45 | 0.640836 | 7.30696 |
| 2 | 79 | 2651.2 | 0.628331 | 7.60789 |
| 3 | 68 | 3637.06 | 2.6151 | 1.24491 |

Chain 1 fits aggregate rooms, ownership, wealth and bequests fairly closely, but childlessness is high, exactly-one-child share is low, first births are too early, and the recent-parent ownership response is nearly zero. These four moments account for about 96% of its objective. This establishes a poor fit at the searched points, not that the normalized utility cannot fit.

[Standard 17 diagnostic plots](standard_diagnostics/) are the fresh verified chain-1 packet. The four JSON receipts retain full precision, identity and validation evidence. Near-bound flags below are the native reporting flags; for the two birth-choice shock scales, very wide bounds make the flags sensitive to range normalization.

## Chain 0: full target fit

| Moment | Role | Target | Model | Gap | Weight | Loss contribution |
|---|---|---:|---:|---:|---:|---:|
| initial_normalization | normalization | 2.1 | 2.1 | 8.36089e-09 | — | 0 |
| cps_childlessness | scored | 0.198279 | 0.31197 | 0.113691 | 35532.3 | 459.277 |
| cps_exactly_one | scored | 0.213655 | 0.072524 | -0.141131 | 26952.8 | 536.848 |
| nchs_mean_age | scored | 25.9763 | 23.4193 | -2.55692 | 139.828 | 914.171 |
| nchs_share30 | validation | 0.249278 | 0.0747757 | -0.174502 | 0 | 0 |
| wealth_earnings | scored | 6.92658 | 7.0703 | 0.143712 | 7.5951 | 0.156862 |
| bequest_wealth | scored | 0.00729102 | 0.00729483 | 3.80162e-06 | 5.16529e+06 | 7.46505e-05 |
| old_dispersion | validation | 3.51594 | 2.78378 | -0.732158 | 0 | 0 |
| mean_rooms | scored | 5.72943 | 5.99188 | 0.262448 | 128.021 | 8.81792 |
| ownership_30_55 | scored | 0.67626 | 0.742285 | 0.066025 | 2339.36 | 10.198 |
| first_birth_rooms | scored | 1.465 | 1.00038 | -0.464623 | 137.565 | 29.6969 |
| family_rooms | scored | 0.3851 | 0.516653 | 0.131554 | 280.528 | 4.85493 |
| recent_parent_ownership | scored | 0.127608 | -0.00521582 | -0.132824 | 27055.8 | 477.326 |
| early_fertility | scored | 0.809528 | 0.639097 | -0.170431 | 100 | 2.90466 |

### All eleven free parameters: chain 0

| Parameter | Estimate | Lower | Upper | Near bound (native flag) |
|---|---:|---:|---:|---|
| beta_annual | 0.966407 | 0.94 | 0.99 | False |
| chi | 1.09037 | 0.1 | 5 | False |
| child_benefit_curvature | 0.0692848 | 0 | 0.8 | False |
| delta_alpha | 0.0289665 | 0 | 0.25 | False |
| delta_alpha_jump | 0.0525411 | 0 | 0.25 | False |
| first_birth_fixed_cost | 1.8973 | 0 | 8 | False |
| kappa_fert | 0.16395 | 0.02 | 50 | True |
| kappa_fert_continuation | 0.405409 | 0.02 | 50 | True |
| psi_child | 0.178698 | 0.01 | 0.5 | False |
| tenure_choice_kappa | 0.0135756 | 0.001 | 0.1 | False |
| theta0 | 0.0966017 | 0 | 8 | False |

### Complete parameter accounting: chain 0

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| H0 | 7.23323 | 0.2 | 80 | False | derived housing supply coefficient at N0=1; reference bounds advisory |
| beta_annual | 0.966407 | 0.94 | 0.99 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.09037 | 0.1 | 5 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 1.8973 | 0 | 8 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.16395 | 0.02 | 50 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.405409 | 0.02 | 50 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.0966017 | 0 | 8 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.0525411 | 0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_curvature | 0.0692848 | 0 | 0.8 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.0135756 | 0.001 | 0.1 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.178698 | 0.01 | 0.5 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.166317 | — | — |  | derived from supplied benefit and curvature |
| theta1 | 0.00819308 | — | — |  | fixed external restriction |
| sigma | 2 | — | — |  | fixed |
| alpha_cons | 0.733 | — | — |  | fixed alpha0=.733 in normalized CES-limit share experiment |
| delta_alpha | 0.0289665 | 0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| h_P | 0 | 0 | 0 | True | fixed zero; housing floor removed in this experiment |
| utility_reference_rent | 0.110466 | — | — |  | inactive legacy input; unused by normalized CES-limit utility |
| q_annual | 0.02 | — | — |  | author-retained 2% annual real rate |
| financed_share | 0.8 | — | — |  | inherited credit contract |
| housing_supply_elasticity | 0.63 | — | — |  | fixed provisional external mapping |
| payroll_tax | 0.0802807 | — | — |  | derived from adopted pension ratio |
| pension_period | 0.917784 | — | — |  | balanced PAYGO |
| annual_depreciation | 0.0141614 | — | — |  | adopted |
| period_depreciation | 0.0554538 | — | — |  | compounded |
| annual_property_tax | 0.0105984 | — | — |  | adopted |
| period_property_tax | 0.0423934 | — | — |  | linear period convention |
| selling_cost | 0.06 | — | — |  | retained |
| rental_cap | 6 | — | — |  | retained provisional |
| wealth_grid_nodes | 120 | — | — |  | retained exact grid |
| income_states | 9 | — | — |  | retained B15 |

## Chain 1: full target fit

| Moment | Role | Target | Model | Gap | Weight | Loss contribution |
|---|---|---:|---:|---:|---:|---:|
| initial_normalization | normalization | 2.1 | 2.1 | -4.16014e-07 | — | 0 |
| cps_childlessness | scored | 0.198279 | 0.310121 | 0.111843 | 35532.3 | 444.465 |
| cps_exactly_one | scored | 0.213655 | 0.0744112 | -0.139244 | 26952.8 | 522.586 |
| nchs_mean_age | scored | 25.9763 | 23.4943 | -2.48193 | 139.828 | 861.335 |
| nchs_share30 | validation | 0.249278 | 0.0764606 | -0.172817 | 0 | 0 |
| wealth_earnings | scored | 6.92658 | 6.873 | -0.0535866 | 7.5951 | 0.0218095 |
| bequest_wealth | scored | 0.00729102 | 0.00731368 | 2.26609e-05 | 5.16529e+06 | 0.00265246 |
| old_dispersion | validation | 3.51594 | 2.82895 | -0.686987 | 0 | 0 |
| mean_rooms | scored | 5.72943 | 5.95166 | 0.222223 | 128.021 | 6.32207 |
| ownership_30_55 | scored | 0.67626 | 0.704662 | 0.0284013 | 2339.36 | 1.88701 |
| first_birth_rooms | scored | 1.465 | 0.683227 | -0.781773 | 137.565 | 84.0757 |
| family_rooms | scored | 0.3851 | 0.405224 | 0.0201244 | 280.528 | 0.113611 |
| recent_parent_ownership | scored | 0.127608 | 0.00413398 | -0.123474 | 27055.8 | 412.491 |
| early_fertility | scored | 0.809528 | 0.63214 | -0.177387 | 100 | 3.14662 |

### All eleven free parameters: chain 1

| Parameter | Estimate | Lower | Upper | Near bound (native flag) |
|---|---:|---:|---:|---|
| beta_annual | 0.966101 | 0.94 | 0.99 | False |
| chi | 1.06235 | 0.1 | 5 | False |
| child_benefit_curvature | 0.0569328 | 0 | 0.8 | False |
| delta_alpha | 0.0152666 | 0 | 0.25 | False |
| delta_alpha_jump | 0.0341544 | 0 | 0.25 | False |
| first_birth_fixed_cost | 1.86424 | 0 | 8 | False |
| kappa_fert | 0.15967 | 0.02 | 50 | True |
| kappa_fert_continuation | 0.412938 | 0.02 | 50 | True |
| psi_child | 0.188726 | 0.01 | 0.5 | False |
| tenure_choice_kappa | 0.0127007 | 0.001 | 0.1 | False |
| theta0 | 0.0961321 | 0 | 8 | False |

### Complete parameter accounting: chain 1

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| H0 | 7.30696 | 0.2 | 80 | False | derived housing supply coefficient at N0=1; reference bounds advisory |
| beta_annual | 0.966101 | 0.94 | 0.99 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.06235 | 0.1 | 5 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 1.86424 | 0 | 8 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.15967 | 0.02 | 50 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.412938 | 0.02 | 50 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.0961321 | 0 | 8 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.0341544 | 0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_curvature | 0.0569328 | 0 | 0.8 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.0127007 | 0.001 | 0.1 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.188726 | 0.01 | 0.5 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.177982 | — | — |  | derived from supplied benefit and curvature |
| theta1 | 0.00819308 | — | — |  | fixed external restriction |
| sigma | 2 | — | — |  | fixed |
| alpha_cons | 0.733 | — | — |  | fixed alpha0=.733 in normalized CES-limit share experiment |
| delta_alpha | 0.0152666 | 0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| h_P | 0 | 0 | 0 | True | fixed zero; housing floor removed in this experiment |
| utility_reference_rent | 0.110466 | — | — |  | inactive legacy input; unused by normalized CES-limit utility |
| q_annual | 0.02 | — | — |  | author-retained 2% annual real rate |
| financed_share | 0.8 | — | — |  | inherited credit contract |
| housing_supply_elasticity | 0.63 | — | — |  | fixed provisional external mapping |
| payroll_tax | 0.0802807 | — | — |  | derived from adopted pension ratio |
| pension_period | 0.917784 | — | — |  | balanced PAYGO |
| annual_depreciation | 0.0141614 | — | — |  | adopted |
| period_depreciation | 0.0554538 | — | — |  | compounded |
| annual_property_tax | 0.0105984 | — | — |  | adopted |
| period_property_tax | 0.0423934 | — | — |  | linear period convention |
| selling_cost | 0.06 | — | — |  | retained |
| rental_cap | 6 | — | — |  | retained provisional |
| wealth_grid_nodes | 120 | — | — |  | retained exact grid |
| income_states | 9 | — | — |  | retained B15 |

## Chain 2: full target fit

| Moment | Role | Target | Model | Gap | Weight | Loss contribution |
|---|---|---:|---:|---:|---:|---:|
| initial_normalization | normalization | 2.1 | 2.1 | -2.13869e-07 | — | 0 |
| cps_childlessness | scored | 0.198279 | 0.317333 | 0.119054 | 35532.3 | 503.63 |
| cps_exactly_one | scored | 0.213655 | 0.0672165 | -0.146439 | 26952.8 | 577.985 |
| nchs_mean_age | scored | 25.9763 | 23.3004 | -2.67584 | 139.828 | 1001.19 |
| nchs_share30 | validation | 0.249278 | 0.0708592 | -0.178419 | 0 | 0 |
| wealth_earnings | scored | 6.92658 | 7.15149 | 0.224907 | 7.5951 | 0.384184 |
| bequest_wealth | scored | 0.00729102 | 0.00731707 | 2.60434e-05 | 5.16529e+06 | 0.00350339 |
| old_dispersion | validation | 3.51594 | 2.75093 | -0.76501 | 0 | 0 |
| mean_rooms | scored | 5.72943 | 6.12032 | 0.390881 | 128.021 | 19.56 |
| ownership_30_55 | scored | 0.67626 | 0.772258 | 0.0959977 | 2339.36 | 21.5585 |
| first_birth_rooms | scored | 1.465 | 1.18464 | -0.280363 | 137.565 | 10.8131 |
| family_rooms | scored | 0.3851 | 0.454138 | 0.0690385 | 280.528 | 1.33709 |
| recent_parent_ownership | scored | 0.127608 | -0.00996951 | -0.137578 | 27055.8 | 512.104 |
| early_fertility | scored | 0.809528 | 0.647095 | -0.162432 | 100 | 2.63843 |

### All eleven free parameters: chain 2

| Parameter | Estimate | Lower | Upper | Near bound (native flag) |
|---|---:|---:|---:|---|
| beta_annual | 0.965845 | 0.94 | 0.99 | False |
| chi | 1.11667 | 0.1 | 5 | False |
| child_benefit_curvature | 0.0752102 | 0 | 0.8 | False |
| delta_alpha | 0.0279342 | 0 | 0.25 | False |
| delta_alpha_jump | 0.0851948 | 0 | 0.25 | False |
| first_birth_fixed_cost | 1.89445 | 0 | 8 | False |
| kappa_fert | 0.159266 | 0.02 | 50 | True |
| kappa_fert_continuation | 0.381652 | 0.02 | 50 | True |
| psi_child | 0.173221 | 0.01 | 0.5 | False |
| tenure_choice_kappa | 0.0168232 | 0.001 | 0.1 | False |
| theta0 | 0.101214 | 0 | 8 | False |

### Complete parameter accounting: chain 2

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| H0 | 7.60789 | 0.2 | 80 | False | derived housing supply coefficient at N0=1; reference bounds advisory |
| beta_annual | 0.965845 | 0.94 | 0.99 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.11667 | 0.1 | 5 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 1.89445 | 0 | 8 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.159266 | 0.02 | 50 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.381652 | 0.02 | 50 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.101214 | 0 | 8 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.0851948 | 0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_curvature | 0.0752102 | 0 | 0.8 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.0168232 | 0.001 | 0.1 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.173221 | 0.01 | 0.5 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.160193 | — | — |  | derived from supplied benefit and curvature |
| theta1 | 0.00819308 | — | — |  | fixed external restriction |
| sigma | 2 | — | — |  | fixed |
| alpha_cons | 0.733 | — | — |  | fixed alpha0=.733 in normalized CES-limit share experiment |
| delta_alpha | 0.0279342 | 0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| h_P | 0 | 0 | 0 | True | fixed zero; housing floor removed in this experiment |
| utility_reference_rent | 0.110466 | — | — |  | inactive legacy input; unused by normalized CES-limit utility |
| q_annual | 0.02 | — | — |  | author-retained 2% annual real rate |
| financed_share | 0.8 | — | — |  | inherited credit contract |
| housing_supply_elasticity | 0.63 | — | — |  | fixed provisional external mapping |
| payroll_tax | 0.0802807 | — | — |  | derived from adopted pension ratio |
| pension_period | 0.917784 | — | — |  | balanced PAYGO |
| annual_depreciation | 0.0141614 | — | — |  | adopted |
| period_depreciation | 0.0554538 | — | — |  | compounded |
| annual_property_tax | 0.0105984 | — | — |  | adopted |
| period_property_tax | 0.0423934 | — | — |  | linear period convention |
| selling_cost | 0.06 | — | — |  | retained |
| rental_cap | 6 | — | — |  | retained provisional |
| wealth_grid_nodes | 120 | — | — |  | retained exact grid |
| income_states | 9 | — | — |  | retained B15 |

## Chain 3: full target fit

| Moment | Role | Target | Model | Gap | Weight | Loss contribution |
|---|---|---:|---:|---:|---:|---:|
| initial_normalization | normalization | 2.1 | 2.1 | 2.73799e-07 | — | 0 |
| cps_childlessness | scored | 0.198279 | 0.31323 | 0.114951 | 35532.3 | 469.517 |
| cps_exactly_one | scored | 0.213655 | 0.0690327 | -0.144623 | 26952.8 | 563.737 |
| nchs_mean_age | scored | 25.9763 | 23.1859 | -2.79037 | 139.828 | 1088.73 |
| nchs_share30 | validation | 0.249278 | 0.0614073 | -0.187871 | 0 | 0 |
| wealth_earnings | scored | 6.92658 | 7.28855 | 0.361964 | 7.5951 | 0.995092 |
| bequest_wealth | scored | 0.00729102 | 0.00714521 | -0.000145813 | 5.16529e+06 | 0.109822 |
| old_dispersion | validation | 3.51594 | 2.69972 | -0.816214 | 0 | 0 |
| mean_rooms | scored | 5.72943 | 2.45928 | -3.27016 | 128.021 | 1369.04 |
| ownership_30_55 | scored | 0.67626 | 0.568117 | -0.108143 | 2339.36 | 27.3589 |
| first_birth_rooms | scored | 1.465 | 1.11547 | -0.34953 | 137.565 | 16.8066 |
| family_rooms | scored | 0.3851 | 0.620991 | 0.235891 | 280.528 | 15.6099 |
| recent_parent_ownership | scored | 0.127608 | 0.0723301 | -0.0552783 | 27055.8 | 82.6741 |
| early_fertility | scored | 0.809528 | 0.651865 | -0.157662 | 100 | 2.48574 |

### All eleven free parameters: chain 3

| Parameter | Estimate | Lower | Upper | Near bound (native flag) |
|---|---:|---:|---:|---|
| beta_annual | 0.966598 | 0.94 | 0.99 | False |
| chi | 1.12091 | 0.1 | 5 | False |
| child_benefit_curvature | 0.0592884 | 0 | 0.8 | False |
| delta_alpha | 0.0751034 | 0 | 0.25 | False |
| delta_alpha_jump | 0.164125 | 0 | 0.25 | False |
| first_birth_fixed_cost | 1.86839 | 0 | 8 | False |
| kappa_fert | 0.144058 | 0.02 | 50 | True |
| kappa_fert_continuation | 0.401756 | 0.02 | 50 | True |
| psi_child | 0.176965 | 0.01 | 0.5 | False |
| tenure_choice_kappa | 0.0139955 | 0.001 | 0.1 | False |
| theta0 | 0.0985322 | 0 | 8 | False |

### Complete parameter accounting: chain 3

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| H0 | 1.24491 | 0.2 | 80 | False | derived housing supply coefficient at N0=1; reference bounds advisory |
| beta_annual | 0.966598 | 0.94 | 0.99 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.12091 | 0.1 | 5 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 1.86839 | 0 | 8 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.144058 | 0.02 | 50 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.401756 | 0.02 | 50 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.0985322 | 0 | 8 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.164125 | 0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_curvature | 0.0592884 | 0 | 0.8 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.0139955 | 0.001 | 0.1 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.176965 | 0.01 | 0.5 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.166473 | — | — |  | derived from supplied benefit and curvature |
| theta1 | 0.00819308 | — | — |  | fixed external restriction |
| sigma | 2 | — | — |  | fixed |
| alpha_cons | 0.733 | — | — |  | fixed alpha0=.733 in normalized CES-limit share experiment |
| delta_alpha | 0.0751034 | 0 | 0.25 | False | supplied primitive; reference calibration bounds advisory |
| h_P | 0 | 0 | 0 | True | fixed zero; housing floor removed in this experiment |
| utility_reference_rent | 0.110466 | — | — |  | inactive legacy input; unused by normalized CES-limit utility |
| q_annual | 0.02 | — | — |  | author-retained 2% annual real rate |
| financed_share | 0.8 | — | — |  | inherited credit contract |
| housing_supply_elasticity | 0.63 | — | — |  | fixed provisional external mapping |
| payroll_tax | 0.0802807 | — | — |  | derived from adopted pension ratio |
| pension_period | 0.917784 | — | — |  | balanced PAYGO |
| annual_depreciation | 0.0141614 | — | — |  | adopted |
| period_depreciation | 0.0554538 | — | — |  | compounded |
| annual_property_tax | 0.0105984 | — | — |  | adopted |
| period_property_tax | 0.0423934 | — | — |  | linear period convention |
| selling_cost | 0.06 | — | — |  | retained |
| rental_cap | 6 | — | — |  | retained provisional |
| wealth_grid_nodes | 120 | — | — |  | retained exact grid |
| income_states | 9 | — | — |  | retained B15 |
