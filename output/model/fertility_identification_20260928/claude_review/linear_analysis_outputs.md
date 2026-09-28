# Linear analysis outputs (generated 2026-09-28 16:00 EDT)

Inputs: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/lead_review/jacobian.csv`, `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/lead_review/jacobian_scaled_svd.json`, `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/overnight_calibration_20260928/cluster/final_review/target_fit_primary_rescore.csv`, `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/nchs_natality_timing/first_birth_counts_year_age.csv`.
All derivatives are the saved central differences at the overnight block0506 anchor. "Per log-unit" means the saved derivative multiplied by the anchor parameter value (the contract SVD convention), so entries read as the moment change for a 100 log-percent parameter change; divide by 100 for a 1 percent change.

## Jacobian per log-unit parameter (full step)

| parameter | norm2.1 | childless | one | age1 | share30 | W/E | beq | oldp90 | rooms | own | fbrooms | famrooms | recent | early | psi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | -0.00872 | -0.0216 | +0.0193 | -0.235 | -0.0131 | +0.294 | -0.000122 | -0.811 | +2.97 | +0.201 | -0.0906 | -0.275 | +0.141 | +0.024 | -0.0703 |
| beta | +0.00106 | -0.603 | +0.22 | -17.9 | -1.18 | +94.3 | +0.00193 | -5.87 | +4.14 | +3.13 | +6.38 | +2.69 | -0.307 | +0.853 | -0.546 |
| chi | +9.67e-05 | +0.00217 | +0.00296 | +0.563 | +0.0337 | +1.36 | +0.000211 | -0.0413 | +0.0643 | +2.02 | +0.0312 | +0.0186 | -0.376 | -0.0457 | -0.0305 |
| fbcost | -0.0106 | +0.0977 | -0.101 | -0.46 | -0.0335 | +0.00176 | -1e-06 | +0.013 | -0.0204 | +0.00233 | +0.101 | -0.11 | +0.0485 | -0.0269 | +0.0832 |
| k_first | -1.51e-06 | -0.0386 | +0.0662 | +2.02 | +0.125 | +0.0272 | -2.25e-05 | +0.00105 | -0.00163 | -0.00237 | -0.0136 | +0.0797 | -0.0856 | -0.0827 | +0.00258 |
| k_cont | +0.011 | -0.0845 | +0.0509 | -1.99 | -0.117 | -0.00352 | +1.2e-05 | -0.0103 | +0.0161 | -0.00625 | -0.154 | -0.105 | -0.0565 | +0.152 | -0.0809 |
| theta0 | -1.71e-06 | +0.000218 | -0.000204 | -0.00145 | -9.14e-05 | +0.137 | +0.000828 | -0.0129 | +0.02 | +0.000385 | +0.00253 | +0.0055 | -0.0023 | +2.6e-05 | +0.000588 |
| dalpha | -0.000102 | +0.0113 | -0.00926 | +0.0907 | +0.00471 | +0.17 | -4.64e-05 | -0.118 | +0.309 | +0.0292 | +1.24 | -0.186 | +0.112 | -0.0122 | +0.0288 |
| curv | -1.96e-06 | -0.0035 | +0.00291 | -0.00455 | +1.68e-05 | +0.000208 | -1.05e-07 | -0.000967 | +0.00151 | -8e-06 | -0.00228 | +0.00787 | -0.00123 | +0.00191 | +0.00377 |
| k_ten | +7.55e-05 | -0.00218 | +0.00175 | -0.0226 | -0.00129 | -0.00742 | -0.000176 | +0.122 | -0.21 | +0.00365 | -0.0543 | +0.104 | -0.025 | +0.00271 | -0.002 |

## Jacobian per log-unit parameter (half step)

| parameter | norm2.1 | childless | one | age1 | share30 | W/E | beq | oldp90 | rooms | own | fbrooms | famrooms | recent | early | psi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | +0.000259 | -0.0241 | +0.0181 | -0.25 | -0.0141 | +0.299 | -0.000155 | -0.603 | +2.97 | +0.193 | -0.108 | -0.275 | +0.135 | +0.0265 | -0.0695 |
| beta | +0.00109 | -0.602 | +0.218 | -17.8 | -1.18 | +94.5 | +0.00188 | +30.9 | +4.15 | +3.13 | +6.41 | +2.62 | -0.3 | +0.852 | -0.546 |
| chi | +4.68e-05 | +0.00218 | +0.00303 | +0.567 | +0.034 | +1.37 | +0.000199 | -0.0393 | +0.0612 | +2.03 | +0.0211 | +0.0259 | -0.375 | -0.046 | -0.0306 |
| fbcost | -0.000138 | +0.0948 | -0.102 | -0.484 | -0.035 | +0.000936 | -6.92e-07 | +0.012 | -0.0187 | +0.00235 | +0.0993 | -0.11 | +0.0472 | -0.0232 | +0.0841 |
| k_first | -4.67e-06 | -0.0386 | +0.0662 | +2.02 | +0.125 | +0.0271 | -2.25e-05 | +0.00105 | -0.00163 | -0.00237 | -0.0136 | +0.0797 | -0.0856 | -0.0827 | +0.00258 |
| k_cont | +0.000134 | -0.0815 | +0.0523 | -1.97 | -0.115 | -0.00266 | +1.17e-05 | -0.00925 | +0.0144 | -0.00627 | -0.152 | -0.105 | -0.0552 | +0.148 | -0.0818 |
| theta0 | +1.94e-05 | +0.000212 | -0.000206 | -0.0015 | -9.43e-05 | +0.137 | +0.000816 | -0.0128 | +0.02 | +0.000381 | +0.00251 | +0.0055 | -0.0023 | +3.34e-05 | +0.000589 |
| dalpha | -4.61e-05 | +0.0113 | -0.00923 | +0.0903 | +0.0047 | +0.169 | -3.26e-05 | -0.198 | +0.309 | +0.0302 | +1.23 | -0.185 | +0.114 | -0.0121 | +0.0288 |
| curv | -8.91e-05 | -0.00348 | +0.00292 | -0.00435 | +2.83e-05 | +0.000215 | -1.08e-07 | -0.000959 | +0.00149 | -8.12e-06 | -0.00226 | +0.00787 | -0.00122 | +0.00188 | +0.00377 |
| k_ten | +4.56e-07 | -0.00216 | +0.00175 | -0.0225 | -0.00128 | -0.00747 | -0.000184 | +0.135 | -0.21 | +0.0034 | -0.0537 | +0.105 | -0.0254 | +0.00269 | -0.00201 |

## Step sensitivity: |full - half| / max(|full|, |half|)

| parameter | norm2.1 | childless | one | age1 | share30 | W/E | beq | oldp90 | rooms | own | fbrooms | famrooms | recent | early | psi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H0 | 1.03 | 0.10 | 0.06 | 0.06 | 0.07 | 0.02 | 0.21 | 0.26 | 0.00 | 0.04 | 0.16 | 0.00 | 0.04 | 0.10 | 0.01 |
| beta | 0.03 | 0.00 | 0.01 | 0.00 | 0.00 | 0.00 | 0.02 | 1.19 | 0.00 | 0.00 | 0.00 | 0.03 | 0.02 | 0.00 | 0.00 |
| chi | 0.52 | 0.00 | 0.02 | 0.01 | 0.01 | 0.01 | 0.06 | 0.05 | 0.05 | 0.00 | 0.33 | 0.28 | 0.00 | 0.01 | 0.00 |
| fbcost | 0.99 | 0.03 | 0.01 | 0.05 | 0.04 | 0.47 | 0.31 | 0.08 | 0.08 | 0.01 | 0.02 | 0.00 | 0.03 | 0.14 | 0.01 |
| k_first | 0.68 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| k_cont | 0.99 | 0.04 | 0.03 | 0.01 | 0.01 | 0.24 | 0.03 | 0.10 | 0.11 | 0.00 | 0.01 | 0.00 | 0.02 | 0.02 | 0.01 |
| theta0 | 1.09 | 0.03 | 0.01 | 0.03 | 0.03 | 0.00 | 0.01 | 0.00 | 0.00 | 0.01 | 0.01 | 0.00 | 0.00 | 0.22 | 0.00 |
| dalpha | 0.55 | 0.01 | 0.00 | 0.00 | 0.00 | 0.01 | 0.30 | 0.40 | 0.00 | 0.03 | 0.00 | 0.00 | 0.02 | 0.00 | 0.00 |
| curv | 0.98 | 0.01 | 0.00 | 0.04 | 0.41 | 0.03 | 0.02 | 0.01 | 0.01 | 0.01 | 0.01 | 0.00 | 0.01 | 0.02 | 0.00 |
| k_ten | 0.99 | 0.01 | 0.00 | 0.01 | 0.01 | 0.01 | 0.04 | 0.10 | 0.00 | 0.07 | 0.01 | 0.00 | 0.02 | 0.01 | 0.00 |

Entries above 0.2 mark derivatives that are not stable across the two step sizes (normalization row, old-dispersion row, several bequest entries).

## Primary weights, implied scales and anchor loss contributions

| moment | target | model (block0506) | gap | weight w | 1/sqrt(w) | w*gap^2 |
|---|---:|---:|---:|---:|---:|---:|
| cps_childlessness | 0.198279 | 0.201062 | +0.0027834 | 35532.3 | 0.005305 | 0.2753 |
| cps_exactly_one | 0.213655 | 0.209422 | -0.0042334 | 26952.8 | 0.0060911 | 0.483 |
| nchs_mean_age | 25.9763 | 25.9328 | -0.043483 | 139.828 | 0.084567 | 0.2644 |
| wealth_earnings | 6.92658 | 6.32631 | -0.60027 | 7.5951 | 0.36286 | 2.737 |
| bequest_wealth | 0.00729102 | 0.00705742 | -0.00023361 | 5.16529e+06 | 0.00044 | 0.2819 |
| mean_rooms | 5.72943 | 5.84794 | +0.11851 | 128.021 | 0.088381 | 1.798 |
| ownership_30_55 | 0.67626 | 0.654818 | -0.021443 | 2339.36 | 0.020675 | 1.076 |
| first_birth_rooms | 1.465 | 1.62211 | +0.15711 | 137.565 | 0.08526 | 3.395 |
| recent_parent_ownership | 0.127608 | 0.119548 | -0.0080606 | 27055.8 | 0.0060795 | 1.758 |
| early_fertility | 0.809528 | 0.535426 | -0.2741 | 100 | 0.1 | 7.513 |
| total | | | | | | 19.5813 |
| total excluding early fertility | | | | | | 12.0681 |

## Constrained early-fertility directions (which parameter moves raise early fertility while holding named moments fixed, to first order)

### full step, hold mean first-birth age

Early-fertility gain per unit log-parameter step along the best feasible direction: **+0.0798** (unconstrained gradient norm 0.8724). Log-parameter step for +0.10 early fertility has norm 1.25:

| parameter | log step | multiplier |
|---|---:|---:|
| H0 | +0.200 | 1.221 |
| beta_annual | -0.066 | 0.936 |
| chi | -0.294 | 0.745 |
| first_birth_fixed_cost | -0.770 | 0.463 |
| kappa_fert | +0.224 | 1.250 |
| kappa_fert_continuation | +0.885 | 2.423 |
| theta0 | -0.001 | 0.999 |
| delta_alpha_jump | -0.123 | 0.884 |
| child_benefit_curvature | +0.027 | 1.027 |
| tenure_choice_kappa | +0.026 | 1.026 |

Implied first-order change in every reported moment: norm2.1 +0.0161, childless -0.125, one +0.127, age1 -9.39e-14, share30 +0.0159, W/E -6.6, beq -0.000207, oldp90 +0.237, rooms +0.287, own -0.771, fbrooms -0.82, famrooms -0.203, recent +0.0382, early +0.1, psi -0.107

### full step, hold mean age and childlessness

Early-fertility gain per unit log-parameter step along the best feasible direction: **+0.0618** (unconstrained gradient norm 0.8724). Log-parameter step for +0.10 early fertility has norm 1.62:

| parameter | log step | multiplier |
|---|---:|---:|
| H0 | +0.217 | 1.242 |
| beta_annual | -0.220 | 0.803 |
| chi | -0.627 | 0.534 |
| first_birth_fixed_cost | -0.337 | 0.714 |
| kappa_fert | -0.510 | 0.600 |
| kappa_fert_continuation | +1.319 | 3.739 |
| theta0 | +0.001 | 1.001 |
| delta_alpha_jump | -0.135 | 0.874 |
| child_benefit_curvature | +0.016 | 1.016 |
| tenure_choice_kappa | +0.031 | 1.031 |

Implied first-order change in every reported moment: norm2.1 +0.016, childless +6.13e-15, one +0.0228, age1 -6.86e-14, share30 +0.0292, W/E -21.5, beq -0.000553, oldp90 +1.14, rooms -0.323, own -1.92, fbrooms -1.84, famrooms -0.775, recent +0.271, early +0.1, psi -0.0162

### full step, hold mean age, childlessness, exactly-one

Early-fertility gain per unit log-parameter step along the best feasible direction: **+0.0119** (unconstrained gradient norm 0.8724). Log-parameter step for +0.10 early fertility has norm 8.43:

| parameter | log step | multiplier |
|---|---:|---:|
| H0 | -6.094 | 0.002 |
| beta_annual | -0.156 | 0.855 |
| chi | -1.086 | 0.338 |
| first_birth_fixed_cost | +2.427 | 11.323 |
| kappa_fert | +2.720 | 15.187 |
| kappa_fert_continuation | +3.924 | 50.605 |
| theta0 | -0.032 | 0.968 |
| delta_alpha_jump | -1.739 | 0.176 |
| child_benefit_curvature | +0.975 | 2.650 |
| tenure_choice_kappa | +0.292 | 1.339 |

Implied first-order change in every reported moment: norm2.1 +0.0705, childless +1.07e-13, one +3.22e-14, age1 -1.7e-11, share30 +0.0208, W/E -18.3, beq +0.000197, oldp90 +6.13, rooms -19.4, own -3.98, fbrooms -3.05, famrooms +1.13, recent -0.946, early +0.1, psi +0.391

### half step, hold mean first-birth age

Early-fertility gain per unit log-parameter step along the best feasible direction: **+0.0768** (unconstrained gradient norm 0.8704). Log-parameter step for +0.10 early fertility has norm 1.30:

| parameter | log step | multiplier |
|---|---:|---:|
| H0 | +0.246 | 1.279 |
| beta_annual | -0.066 | 0.936 |
| chi | -0.319 | 0.727 |
| first_birth_fixed_cost | -0.787 | 0.455 |
| kappa_fert | +0.238 | 1.269 |
| kappa_fert_continuation | +0.912 | 2.490 |
| theta0 | -0.001 | 0.999 |
| delta_alpha_jump | -0.132 | 0.876 |
| child_benefit_curvature | +0.028 | 1.029 |
| tenure_choice_kappa | +0.027 | 1.028 |

Implied first-order change in every reported moment: norm2.1 +0.00021, childless -0.126, one +0.134, age1 -6.57e-15, share30 +0.0161, W/E -6.66, beq -0.000222, oldp90 -2.18, rooms +0.416, own -0.82, fbrooms -0.844, famrooms -0.212, recent +0.0491, early +0.1, psi -0.115

### half step, hold mean age and childlessness

Early-fertility gain per unit log-parameter step along the best feasible direction: **+0.0602** (unconstrained gradient norm 0.8704). Log-parameter step for +0.10 early fertility has norm 1.66:

| parameter | log step | multiplier |
|---|---:|---:|
| H0 | +0.267 | 1.306 |
| beta_annual | -0.221 | 0.802 |
| chi | -0.659 | 0.517 |
| first_birth_fixed_cost | -0.347 | 0.706 |
| kappa_fert | -0.499 | 0.607 |
| kappa_fert_continuation | +1.347 | 3.845 |
| theta0 | +0.001 | 1.001 |
| delta_alpha_jump | -0.145 | 0.865 |
| child_benefit_curvature | +0.018 | 1.018 |
| tenure_choice_kappa | +0.033 | 1.033 |

Implied first-order change in every reported moment: norm2.1 +3.45e-05, childless -4.67e-15, one +0.0289, age1 +5.96e-14, share30 +0.0296, W/E -21.7, beq -0.000562, oldp90 -6.94, rooms -0.188, own -1.99, fbrooms -1.87, famrooms -0.782, recent +0.284, early +0.1, psi -0.0225

### half step, hold mean age, childlessness, exactly-one

Early-fertility gain per unit log-parameter step along the best feasible direction: **+0.0132** (unconstrained gradient norm 0.8704). Log-parameter step for +0.10 early fertility has norm 7.56:

| parameter | log step | multiplier |
|---|---:|---:|
| H0 | +4.292 | 73.125 |
| beta_annual | -0.074 | 0.928 |
| chi | -3.165 | 0.042 |
| first_birth_fixed_cost | +3.487 | 32.703 |
| kappa_fert | +3.118 | 22.603 |
| kappa_fert_continuation | +1.461 | 4.310 |
| theta0 | -0.009 | 0.991 |
| delta_alpha_jump | -1.997 | 0.136 |
| child_benefit_curvature | +0.748 | 2.113 |
| tenure_choice_kappa | +0.329 | 1.389 |

Implied first-order change in every reported moment: norm2.1 +0.000611, childless +5.11e-14, one -9.45e-14, age1 +6.96e-12, share30 +0.00995, W/E -10.3, beq -0.00149, oldp90 -4.29, rooms +11.5, own -5.89, fbrooms -3.41, famrooms -1.34, recent +1.37, early +0.1, psi -0.0345

## Linearized Gauss-Newton steps from block0506 (ridge on log steps as a trust region)

| lane | ridge | lane loss at anchor | linear-predicted lane loss | common-primary rescore | other-moment primary loss | max abs log step | predicted early | predicted mean age | predicted W/E | predicted ownership |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| primary | 1000 | 19.58 | 18.15 | 18.15 | 10.66 | 0.01 | 0.536 | 25.917 | 6.59 | 0.659 |
| primary | 100 | 19.58 | 15.58 | 15.58 | 7.97 | 0.07 | 0.534 | 25.942 | 6.72 | 0.659 |
| primary | 10 | 19.58 | 9.52 | 9.52 | 1.80 | 0.30 | 0.532 | 25.958 | 6.84 | 0.668 |
| early10 | 1000 | 87.20 | 84.74 | 19.12 | 11.83 | 0.01 | 0.540 | 25.861 | 6.64 | 0.658 |
| early10 | 100 | 87.20 | 82.77 | 17.35 | 10.08 | 0.07 | 0.540 | 25.856 | 6.64 | 0.658 |
| early10 | 10 | 87.20 | 77.63 | 11.60 | 4.26 | 0.29 | 0.539 | 25.863 | 6.71 | 0.667 |
| early100 | 1000 | 763.38 | 638.83 | 75.91 | 70.22 | 0.14 | 0.571 | 25.373 | 7.07 | 0.647 |
| early100 | 100 | 763.38 | 607.57 | 127.04 | 122.19 | 0.26 | 0.589 | 25.171 | 6.00 | 0.650 |
| early100 | 10 | 763.38 | 604.81 | 142.20 | 137.53 | 0.26 | 0.593 | 25.127 | 5.72 | 0.656 |

Predicted parameter points (ridge 10):

| parameter | anchor | primary GN | early10 GN | early100 GN | bounds |
|---|---:|---:|---:|---:|---|
| H0 | 6.2935 | 6.1527 | 6.1689 | 6.2962 | [0.2, 80] |
| beta_annual | 0.96348 | 0.96902 | 0.96764 | 0.95697 | [0.94, 0.99] |
| chi | 1.0939 | 1.095 | 1.0963 | 1.1063 | [0.1, 5] |
| first_birth_fixed_cost | 0.62094 | 0.45781 | 0.46351 | 0.51037 | [0, 8] |
| kappa_fert | 0.17561 | 0.1359 | 0.13577 | 0.1348 | [0.02, 50] |
| kappa_fert_continuation | 0.33177 | 0.257 | 0.27213 | 0.42462 | [0.02, 50] |
| theta0 | 0.12452 | 0.12955 | 0.13069 | 0.13993 | [0, 8] |
| delta_alpha_jump | 0.13488 | 0.11963 | 0.12168 | 0.1388 | [0, 0.25] |
| child_benefit_curvature | 0.061455 | 0.061996 | 0.061844 | 0.060667 | [0, 0.8] |
| tenure_choice_kappa | 0.012087 | 0.011357 | 0.011412 | 0.011852 | [0.001, 0.1] |

## Linear trade-off frontier: minimum other-moment primary loss for a given early-fertility gain

Minimizes the primary loss of the nine other scored moments (plus ridge 10 on log steps) subject to a first-order early-fertility gain delta. Anchor other-moment loss 12.07.

| early gain | early level | min other-moment loss | mean age | childless | exactly-one | W/E | ownership | rooms | k_first | k_cont | fbcost | beta |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| +0.030 | 0.565 | 44.2 | 25.50 | 0.186 | 0.218 | 6.23 | 0.661 | 5.82 | 0.135 | 0.338 | 0.486 | 0.9624 |
| +0.050 | 0.585 | 105.7 | 25.23 | 0.179 | 0.221 | 5.86 | 0.657 | 5.84 | 0.135 | 0.398 | 0.503 | 0.9585 |
| +0.067 | 0.602 | 179.1 | 25.01 | 0.172 | 0.224 | 5.55 | 0.654 | 5.85 | 0.135 | 0.457 | 0.519 | 0.9552 |
| +0.100 | 0.635 | 377.6 | 24.56 | 0.160 | 0.230 | 4.95 | 0.648 | 5.89 | 0.134 | 0.598 | 0.550 | 0.9489 |
| +0.132 | 0.667 | 640.3 | 24.13 | 0.149 | 0.235 | 4.37 | 0.641 | 5.92 | 0.133 | 0.776 | 0.581 | 0.9427 |
| +0.150 | 0.685 | 818.5 | 23.89 | 0.142 | 0.238 | 4.04 | 0.638 | 5.94 | 0.133 | 0.898 | 0.600 | 0.9393 |
| +0.200 | 0.735 | 1428.4 | 23.22 | 0.123 | 0.246 | 3.13 | 0.628 | 5.99 | 0.132 | 1.349 | 0.655 | 0.9298 |
| +0.274 | 0.809 | 2641.3 | 22.22 | 0.096 | 0.259 | 1.78 | 0.613 | 6.06 | 0.131 | 2.462 | 0.747 | 0.9160 |

Comparison points from actual solves: early_frontier half first-birth scale (+0.067, other-moment loss 491.5), quarter scale (+0.132, other-moment loss 1398.4); overnight identity winner 0526 (+0.121, other-moment loss 550.5 = 552.79 - 100*0.1528^2).

## NCHS 2003-2006 pooled first births: cell shares and mean-age conventions

| model cell | share of first births |
|---|---:|
| <=21 (model cell 18-21, includes ages 12-17) | 0.3417 |
| 22-25 | 0.2156 |
| 26-29 | 0.1935 |
| 30-33 | 0.1426 |
| 34-37 | 0.0751 |
| 38-41 | 0.0263 |
| 42+ | 0.0054 |
| memo: age < 18 | 0.0773 |
| memo: age <= 19 | 0.2070 |
| memo: age >= 30 | 0.2493 |

Mean first-birth age, raw single-year: 25.161. Mean of cell midpoints (contract rule, the actual target): 25.976. Raw mean among mothers 18+: 25.906.

## Time-aggregation bound on model early fertility under data-consistent first-birth timing

Model children ever born at [25,26) = P(birth in cell 18-21) + 0.875 * P(birth in cell 22-25), with at most one birth per four-year cell. With eventual-mother share M = 1 - 0.1983 = 0.8017 and NCHS first-birth cell shares f1 = 0.3417, f2 = 0.2156, and p2 the probability that a cell-1 mother has a second birth in cell 2:

    early = M*f1 + 0.875*(M*f2 + p2*M*f1) = 0.4252 + 0.2397*p2

| p2 | model-consistent early fertility |
|---:|---:|
| 0.0 | 0.425 |
| 0.4 | 0.521 |
| 0.5 | 0.545 |
| 0.6 | 0.569 |
| 0.8 | 0.617 |
| 1.0 | 0.665 |

Target 0.8095; block0506 model 0.5354 (implied p2 about 0.46 under data-consistent timing). The p2 = 1 ceiling is 0.665.

## Near-bound flag arithmetic

- kappa_fert = 0.1756, bounds [0.02, 50]: raw-range rule threshold 0.01*(hi-lo) = 0.500; distance to lower bound 0.156 -> flagged True; in log units the distance to the lower bound is 2.17 of a 7.82 log range (28 percent).
- kappa_fert_continuation = 0.3318, bounds [0.02, 50]: raw-range rule threshold 0.01*(hi-lo) = 0.500; distance to lower bound 0.312 -> flagged True; in log units the distance to the lower bound is 2.81 of a 7.82 log range (36 percent).
