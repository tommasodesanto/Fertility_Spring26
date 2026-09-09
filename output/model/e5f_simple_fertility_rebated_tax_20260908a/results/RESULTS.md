# Verified rebated property-tax impact

Annual property tax rises from 1% to 2%, with equal household rebates in both equilibria. The same calibrated 2023 population is inherited.

| Outcome | Rebated 1% | Rebated 2% | Change |
|---|---:|---:|---:|
| Births per household | 0.08767873 | 0.08812525 | +0.509263 % |
| Rooms per household | 6.32524837 | 6.13207918 | -3.053938 % |
| All-age ownership | 0.65550158 | 0.65571787 | +0.021629 pp |
| Young ownership | 0.33365784 | 0.34601743 | +1.235959 pp |
| Young rooms | 5.42091395 | 5.18879543 | -4.281907 % |
| Asset price | 0.73052835 | 0.69545103 | -4.801637 % |

Births increase only modestly. Young ownership improves, while young housing services fall. No causal channel attribution is established before the eight-cell decomposition.

| Group | Ownership change (pp) | Rooms change (%) |
|---|---:|---:|
| all_ages / dependent_children | -0.177713 | -2.906482 |
| all_ages / no_dependent_children | +0.138984 | -3.148325 |
| all_ages / childless_parity_zero | +1.246385 | -4.002175 |
| young_whole_nodes_25_34 / dependent_children | +0.106470 | -3.476742 |
| young_whole_nodes_25_34 / no_dependent_children | +2.249957 | -5.240267 |
| young_whole_nodes_25_34 / childless_parity_zero | +2.476625 | -5.289859 |

Young model nodes are ages 26, 30 and 34; annual-age ACS alignment remains unresolved. Family groups are defined by realized post-choice status, so changes combine behavior and composition. Households without dependent children include empty nesters; parity-zero households are also reported separately.

Jobs 17222674, 17222675 and 17222676 completed successfully. Both coupled price/rebate roots and fresh full-policy replays pass. Root residual tolerance 1e-4, fiscal absolute tolerance 2.5e-5, market tolerance 2e-4; mass, feasibility, budget, probability and occupied-value gates pass. The lead verified all remote receipt artifact hashes, including checkpoints, before collecting these tables.

No long-run effect, population forecast, recalibration or production promotion is implied. Full active target-fit and free-parameter tables: [morning review](../../e5f_simple_fertility_overnight_20260907a/morning_review/MORNING_REVIEW.md).
