# 2023 snapshot: model and data

These cross-sectional comparisons are validation moments, not fitted targets. CPS is June 2024; first-birth timing is NCHS 2023; resources and room distributions use PSID 2019.

| Validation moment | Data | Model | Model − data |
|---|---:|---:|---:|
| Children ever born, ages 22–25 (model 3+ adjusted) | 0.3089 | 0.3052 | -0.0037 |
| Childlessness, ages 22–25 | 0.7891 | 0.7280 | -0.0611 |
| Children ever born, ages 40–44 (model 3+ adjusted) | 1.9184 | 1.5968 | -0.3216 |
| Childlessness, ages 40–44 | 0.1882 | 0.2993 | +0.1111 |
| First-birth mean age (band midpoints) | 28.0902 | 26.9731 | -1.1171 |
| First-birth share age 30+ | 0.3766 | 0.2872 | -0.0894 |
| Mean rooms, all households (PSID) | 5.2279 | 5.7537 | +0.5258 |
| Mean rooms, owners (PSID) | 6.2821 | 6.4698 | +0.1877 |
| Mean rooms, renters (PSID) | 3.7893 | 4.0155 | +0.2261 |
| Mean total net wealth / own mean earnings | 5.4616 | 5.0353 | -0.4262 |
| Negative financial position, all ages | 0.3606 | 0.3192 | -0.0414 |

## Complete target system for the retained transition

Only the 2020–2023 window was fitted. These are the original household-rate fertility statistics, distinct from children-ever-born stocks above.

| Birth window | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| 2008–2011 | 1.974875 | 1.626619 | -0.348256 | 0 | 0 |
| 2012–2015 | 1.861000 | 1.641307 | -0.219693 | 0 | 0 |
| 2016–2019 | 1.755375 | 1.641700 | -0.113675 | 0 | 0 |
| 2020–2023 | 1.645750 | 1.643134 | -0.002616 | 1 | 6.8433246e-06 |

Estimated shock: psi_child=0.1199969464, bounds [0.0017892072, 0.3578414413], away from either bound. [All 31 supplied baseline parameters and reference bounds](../transition/solution/baseline_parameters.csv).

## CPS data check

The weighted five-year age groups reproduce Census Table 1 population totals and children-count shares to its published rounding. Ages 35–45 have about 600–700 respondents per single age. The annual-age line is a cross-section of different cohorts, not a trajectory for the same women. No monotonicity is imposed.

Mean graphs use the retained model 3+ weight, 3.602359422009, and CPS public counts (five or more coded five). Distribution bars still group 3+. Applying the completed-fertility tail mean at younger ages and within income groups is a reporting approximation: fourth and later birth dates are not separately modeled. The previous capped-at-three view was comparable on its own terms, but omitted this existing measurement adjustment.

[Official Census Table 1](https://www2.census.gov/programs-surveys/demo/tables/fertility/2024/am-women-fertility/t1.xlsx). Census documentation identifies PRTAGE as the masked public-use age variable; masking is a possible contributor to single-age irregularity, not an established explanation of these specific jumps. Mean graph definitions now include the existing top-bin adjustment; model behavior and calibration targets are unchanged.
