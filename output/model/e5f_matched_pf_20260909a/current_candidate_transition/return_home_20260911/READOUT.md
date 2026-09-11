# Quantitative assessment after the commute — September 11

The utility and pension implementation has produced verified initial and terminal solutions and three cleared short historical paths. Historical preference estimation, continuation-horizon verification, and policies under the revised calibration remain unfinished. No further search was submitted during this assessment.

## Main conclusions

- Social Security: the actual-budget repair is verified in the initial economy, all three trial terminal equilibria, and all three short historical paths. This does not certify the still-solving longer paths.
- Initial calibration: the unrestricted reference loss is 158.0764. Fixing annual beta at0.98 and reoptimizing gives222.5244, with two exact repetitions. The0.99 profile reached168.2011 before a population-mass check stopped it; it has no final exact repetitions or joint refinement and is not an optimum.
- Lower beta has a real trade-off. At0.98 the first-birth rooms response is0.7575 against0.7202, but wealth/earnings falls to3.8823 against6.1459 and ownership at30–55 to49.60% against64.83%. All rows below must be considered together.
- Preference sensitivity: all three short roots clear and their reported rate profiles reproduce exactly. These are fixed amplitudes, not adaptive estimation. Correct female/maternal-age measurement and the outer preference update remain missing.
- Longer paths: the central28-date run finished its eight-mapping budget without convergence (maximum housing residual1.1056%, pension residual0.028509%). Exact replay/checkpoint checks pass. Both alternatives were still solving at collection. None supplies horizon certification.
- Presentation: the revised specification and empirical evidence can be developed now. New numerical historical-fit and policy claims are not ready. Do not combine old policy results with the revised utility/pension calibration.

## Reproduced short-path fertility diagnostics

**These columns do not measure exactly the same object.** Data are equal-weight averages of published annual female period TFR. Model columns sum four-year birth flows divided by model-household mass in each age cell. Their numerical proximity is suggestive, not a certified data fit. The paths also use only six dates, so horizon sensitivity remains open.

| Birth-data window | Model decision | Female TFR, data | Smaller decline −0.025 | Middle decline −0.05 | Larger decline −0.10 |
|---|---:|---:|---:|---:|---:|
| 2008–2011 | 2007 | 1.974875 | 2.007444 | 1.924405 | 1.745358 |
| 2012–2015 | 2011 | 1.861 | 1.942204 | 1.822482 | 1.572013 |
| 2016–2019 | 2015 | 1.755375 | 1.883424 | 1.732933 | 1.430065 |
| 2020–2023 | 2019 | 1.64575 | 1.842701 | 1.669684 | 1.333686 |

The2023 decision generates2024–2027 births; it is not compared with the2020–2023 data block. Different2007 outcomes reflect anticipation of different announced future paths despite the same inherited state.

## Complete initial-calibration fits and parameters

All three profiles use the same12 scored moments/weights and separate2.1 normalization. These are the approved working minimum-distance diagnostics, not a claim of fully certified empirical SMM or global identification. Baseline has nine free structural coordinates; fixed-beta profiles have eight. Near-bound flags reproduce the saved convention and need not indicate a binding constraint.

### Unrestricted reference: reproduced r5_joint_09

Loss: **158.076390747**.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.100151 | 0.0001506047 | — | — |
| Childless women, ages 40–44 | 0.1982788 | 0.2050178 | 0.006739083 | 35532.3 | 1.613708 |
| Exactly one child among mothers, ages 40–44 | 0.2136553 | 0.2053219 | -0.008333442 | 26952.82 | 1.871772 |
| Period mean first-birth age | 25.97626 | 26.11317 | 0.1369043 | 139.8281 | 2.620768 |
| First births at age 30+ | 0.249278 | 0.2305588 | -0.01871924 | 13866.07 | 4.858808 |
| Wealth / annual gross labor earnings | 6.145861 | 5.470594 | -0.675267 | 7.595098 | 3.463255 |
| Annual bequests / aggregate wealth | 0.0088 | 0.008319317 | -0.0004806828 | 5165289 | 1.193471 |
| Old wealth/income p90 / median, ages 76–84 | 3.515935 | 4.335052 | 0.8191165 | 10.61636 | 7.123068 |
| Mean occupied rooms, capped at 9 | 5.561097 | 6.321355 | 0.7602572 | 128.0207 | 73.99482 |
| Ownership, heads 30–55 | 0.648334 | 0.5348647 | -0.1134693 | 2339.362 | 30.11996 |
| First-birth room response, −1 to +3 | 0.7202463 | 1.113533 | 0.3932872 | 137.5653 | 21.27788 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.3470669 | 0.213083 | -0.1339839 | 280.5281 | 5.035956 |
| Recent-parent ownership gap | 0.1628955 | 0.1494339 | -0.01346161 | 27055.82 | 4.902918 |

| Parameter | Value | Lower | Upper | Near bound | Restriction/status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.9957961 | 0.94 | 0.9995 | False | diagnostic candidate; not a certified estimate |
| kappa_fert | 0.3346343 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.4037306 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| chi | 1.045633 | 0.1 | 5 | False | diagnostic candidate; not a certified estimate |
| H0 | 8.2255 | 0.2 | 80 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.08957524 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.07709538 | 0.02 | 16 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.2937155 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.183795 | 0.1 | 2.3 | False | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | False | zero restriction |
| psi_child | 0.1652598 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.046361 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2 | — | — | False | externally fixed |

Source tables: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_target_fit.csv` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/selected_parameters.csv`.

### Annual beta0.98: completed bounded search and two exact repetitions

Loss: **222.524449191**.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.1 | 3.665736e-07 | — | — |
| Childless women, ages 40–44 | 0.1982788 | 0.2124979 | 0.01421914 | 35532.3 | 7.184064 |
| Exactly one child among mothers, ages 40–44 | 0.2136553 | 0.2046513 | -0.009004038 | 26952.82 | 2.185138 |
| Period mean first-birth age | 25.97626 | 26.09289 | 0.1166271 | 139.8281 | 1.901923 |
| First births at age 30+ | 0.249278 | 0.2346194 | -0.01465864 | 13866.07 | 2.979482 |
| Wealth / annual gross labor earnings | 6.145861 | 3.882343 | -2.263518 | 7.595098 | 38.9136 |
| Annual bequests / aggregate wealth | 0.0088 | 0.008732953 | -6.704722e-05 | 5165289 | 0.02321968 |
| Old wealth/income p90 / median, ages 76–84 | 3.515935 | 4.926317 | 1.410382 | 10.61636 | 21.11781 |
| Mean occupied rooms, capped at 9 | 5.561097 | 6.322453 | 0.7613551 | 128.0207 | 74.20869 |
| Ownership, heads 30–55 | 0.648334 | 0.4959938 | -0.1523402 | 2339.362 | 54.29085 |
| First-birth room response, −1 to +3 | 0.7202463 | 0.7574763 | 0.03723005 | 137.5653 | 0.190676 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.3470669 | 0.1711666 | -0.1759003 | 280.5281 | 8.679799 |
| Recent-parent ownership gap | 0.1628955 | 0.1428707 | -0.0200248 | 27055.82 | 10.84919 |

| Parameter | Value | Lower | Upper | Near bound | Restriction/status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.98 | 0.98 | 0.98 | False | Externally fixed in this profile |
| kappa_fert | 0.2848706 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.3262545 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| chi | 1.05093 | 0.1 | 5 | False | diagnostic candidate; not a certified estimate |
| H0 | 8.47204 | 0.2 | 80 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.06311569 | 0 | 8 | True | diagnostic candidate; not a certified estimate |
| theta1 | 0.0283618 | 0.02 | 16 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.1553662 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.3 | 0.1 | 2.3 | True | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | False | zero restriction |
| psi_child | 0.1577083 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.046361 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2 | — | — | False | externally fixed |

Source tables: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/current_candidate_transition/return_home_20260911/beta_098/selected_target_fit.csv` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/current_candidate_transition/return_home_20260911/beta_098/selected_parameters.csv`.

### Annual beta0.99: stopped first-round profile; unrepeated best candidate

Loss: **168.201074488**.

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.100222 | 0.0002222855 | — | — |
| Childless women, ages 40–44 | 0.1982788 | 0.2070217 | 0.008742968 | 35532.3 | 2.716071 |
| Exactly one child among mothers, ages 40–44 | 0.2136553 | 0.2042096 | -0.009445759 | 26952.82 | 2.404794 |
| Period mean first-birth age | 25.97626 | 26.01837 | 0.04210529 | 139.8281 | 0.2478949 |
| First births at age 30+ | 0.249278 | 0.2269837 | -0.02229432 | 13866.07 | 6.891946 |
| Wealth / annual gross labor earnings | 6.145861 | 4.853467 | -1.292395 | 7.595098 | 12.68597 |
| Annual bequests / aggregate wealth | 0.0088 | 0.00857263 | -0.0002273696 | 5165289 | 0.2670296 |
| Old wealth/income p90 / median, ages 76–84 | 3.515935 | 4.51595 | 1.000015 | 10.61636 | 10.61668 |
| Mean occupied rooms, capped at 9 | 5.561097 | 6.298926 | 0.7378286 | 128.0207 | 69.69332 |
| Ownership, heads 30–55 | 0.648334 | 0.5181074 | -0.1302266 | 2339.362 | 39.6732 |
| First-birth room response, −1 to +3 | 0.7202463 | 0.9811284 | 0.2608822 | 137.5653 | 9.362626 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.3470669 | 0.1853837 | -0.1616832 | 280.5281 | 7.333415 |
| Recent-parent ownership gap | 0.1628955 | 0.1476262 | -0.01526933 | 27055.82 | 6.30813 |

| Parameter | Value | Lower | Upper | Near bound | Restriction/status |
|---|---:|---:|---:|---|---|
| beta_annual | 0.99 | 0.99 | 0.99 | False | Externally fixed in this profile |
| kappa_fert | 0.3346343 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| kappa_fert_continuation | 0.4037306 | 0.02 | 50 | True | diagnostic candidate; not a certified estimate |
| chi | 1.045633 | 0.1 | 5 | False | diagnostic candidate; not a certified estimate |
| H0 | 8.2255 | 0.2 | 80 | False | diagnostic candidate; not a certified estimate |
| theta0 | 0.08957524 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| theta1 | 0.07709538 | 0.02 | 16 | True | diagnostic candidate; not a certified estimate |
| first_birth_fixed_cost | 0.2937155 | 0 | 8 | False | diagnostic candidate; not a certified estimate |
| h_P | 2.227911 | 0.1 | 2.3 | False | diagnostic candidate; not a certified estimate |
| hbar_child_rooms | 0 | — | — | False | zero restriction |
| psi_child | 0.1648933 | — | — | False | normalized to 2.1 |
| payroll_tax | 0.179 | — | — | False | externally fixed |
| pension_period | 2.046361 | — | — | False | budget derived |
| housing_supply_elasticity | 0.63 | — | — | False | externally fixed |
| tenure_choice_kappa | 0.005 | — | — | False | externally fixed |
| alpha_cons | 0.733 | — | — | False | externally fixed |
| sigma | 2 | — | — | False | externally fixed |

Source tables: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/current_candidate_transition/return_home_20260911/beta_099/selected_target_fit.csv` and `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/current_candidate_transition/return_home_20260911/beta_099/selected_parameters.csv`.

## Next decision

Keep the announced-path experiment and finish its measurement objective and adaptive preference fit. Use the completed trial results as initial evaluations. Decide between the high-beta reference and a completed0.99 profile using the full fit, rather than demanding another initial-calibration optimum before historical fitting. Establish the horizon accuracy of the selected history before promoting policy results.

The slides task can independently revise the exposition, utility-floor equation and pension description, while scientific claims and numerical replacements are obtained from the technical tasks.
