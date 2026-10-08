reference for baseline: saved packet
## baseline

Legacy screen vs saved packet: 96 solution arrays compared, 0 differ [].
Fix vs saved packet: 62 of 96 arrays differ.

- saved residuals: {"renewal_residual": -1.3829918321661694e-09, "absolute_housing_residual": 0.0, "rebate_T": 0.1921730548652237, "rebate_residual": -7.341640333171904e-07}
- legacy residuals: {"renewal_residual": -1.3829918321661694e-09, "absolute_housing_residual": 0.0, "rebate_T": 0.1921730548652237, "rebate_residual": -7.341640333171904e-07}
- fix residuals: {"renewal_residual": 1.6965993201223384e-05, "absolute_housing_residual": 0.0, "rebate_T": 0.19214659162389147, "rebate_residual": -1.380585231350422e-10}

keys: moment, target, model, gap, weight, loss_contribution, role

| moment | target | model saved | model legacy | model fix | fix - saved | weight | loss saved | loss fix |
|---|---|---|---|---|---|---|---|---|
| initial_normalization | 2.1 | 2.1 | 2.1 | 2.10004 | +3.56e-05 |  | 0 | 0 |
| cps_childlessness | 0.19827875100684264 | 0.202193 | 0.202193 | 0.202161 | -3.2e-05 | 35532.3042455214 | 0.5444 | 0.5356 |
| cps_exactly_one | 0.21365532522014702 | 0.217103 | 0.217103 | 0.217125 | +2.18e-05 | 26952.820824310795 | 0.3203 | 0.3244 |
| nchs_mean_age | 25.976263860992496 | 25.9628 | 25.9628 | 25.9631 | +0.00028 | 139.82806784479274 | 0.02523 | 0.02419 |
| nchs_share30 | 0.2492780130410667 | 0.237511 | 0.237511 | 0.237502 | -9.27e-06 | 0.0 | 0 | 0 |
| wealth_earnings | 4.45838713455674 | 4.58349 | 4.58349 | 4.58314 | -0.00035 | 7.595098472533724 | 0.1189 | 0.1182 |
| bequest_wealth | 0.007291023472616158 | 0.00712342 | 0.00712342 | 0.0071237 | +2.84e-07 | 5165289.256198346 | 0.1451 | 0.1446 |
| old_dispersion | 3.51593508651872 | 3.30769 | 3.30769 | 3.30769 | +0 | 0.0 | 0 | 0 |
| mean_rooms | 5.729434240102641 | 5.81601 | 5.81601 | 5.81522 | -0.000797 | 128.02070205233477 | 0.9596 | 0.942 |
| ownership_30_55 | 0.6762604168538028 | 0.669864 | 0.669864 | 0.670244 | +0.00038 | 2339.3623724673616 | 0.09572 | 0.08467 |
| first_birth_rooms | 1.465 | 1.26385 | 1.26385 | 1.2647 | +0.000849 | 137.5652749002964 | 5.566 | 5.519 |
| family_rooms | 0.38509964969278165 | 0.293947 | 0.293947 | 0.294026 | +7.9e-05 | 0.0 | 0 | 0 |
| recent_parent_ownership | 0.12760836356692162 | 0.126658 | 0.126658 | 0.12445 | -0.00221 | 27055.822957508266 | 0.02444 | 0.2699 |
| early_fertility | 0.8095276384290021 | 0.552574 | 0.552574 | 0.552556 | -1.86e-05 | 100.0 | 6.603 | 6.603 |

Total loss saved 14.402424, fix 14.566397.

reference for phi095: legacy rerun (no saved packet)
## phi095

Legacy screen vs saved packet: 96 solution arrays compared, 0 differ [].
Fix vs saved packet: 62 of 96 arrays differ.

- saved residuals: {"renewal_residual": -0.014170585107748157, "absolute_housing_residual": 0.0, "rebate_T": 0.19066292745880686, "rebate_residual": -1.22576082960658e-07}
- legacy residuals: {"renewal_residual": -0.014170585107748157, "absolute_housing_residual": 0.0, "rebate_T": 0.19066292745880686, "rebate_residual": -1.22576082960658e-07}
- fix residuals: {"renewal_residual": -0.011140609073998164, "absolute_housing_residual": 0.0, "rebate_T": 0.19098093000006644, "rebate_residual": -6.975011122745857e-09}

keys: moment, target, model, gap, weight, loss_contribution, role

| moment | target | model saved | model legacy | model fix | fix - saved | weight | loss saved | loss fix |
|---|---|---|---|---|---|---|---|---|
| initial_normalization | 2.1 | 2.07024 | 2.07024 | 2.0766 | +0.00636 |  | 0 | 0 |
| cps_childlessness | 0.19827875100684264 | 0.214986 | 0.214986 | 0.211815 | -0.00317 | 35532.3042455214 | 9.919 | 6.511 |
| cps_exactly_one | 0.21365532522014702 | 0.215768 | 0.215768 | 0.216728 | +0.000959 | 26952.820824310795 | 0.1204 | 0.2544 |
| nchs_mean_age | 25.976263860992496 | 25.9334 | 25.9334 | 25.9469 | +0.0135 | 139.82806784479274 | 0.2569 | 0.1203 |
| nchs_share30 | 0.2492780130410667 | 0.236824 | 0.236824 | 0.237484 | +0.000659 | 0.0 | 0 | 0 |
| wealth_earnings | 4.45838713455674 | 4.56445 | 4.56445 | 4.5399 | -0.0246 | 7.595098472533724 | 0.08545 | 0.05046 |
| bequest_wealth | 0.007291023472616158 | 0.00708214 | 0.00708214 | 0.00711187 | +2.97e-05 | 5165289.256198346 | 0.2254 | 0.1658 |
| old_dispersion | 3.51593508651872 | 3.30769 | 3.30769 | 3.30769 | +0 | 0.0 | 0 | 0 |
| mean_rooms | 5.729434240102641 | 5.77031 | 5.77031 | 5.77994 | +0.00962 | 128.02070205233477 | 0.2139 | 0.3265 |
| ownership_30_55 | 0.6762604168538028 | 0.843979 | 0.843979 | 0.806141 | -0.0378 | 2339.3623724673616 | 65.81 | 39.46 |
| first_birth_rooms | 1.465 | 1.22457 | 1.22457 | 1.22231 | -0.00226 | 137.5652749002964 | 7.952 | 8.102 |
| family_rooms | 0.38509964969278165 | 0.29129 | 0.29129 | 0.291427 | +0.000137 | 0.0 | 0 | 0 |
| recent_parent_ownership | 0.12760836356692162 | -0.00528238 | -0.00528238 | -0.0344049 | -0.0291 | 27055.822957508266 | 477.8 | 710.2 |
| early_fertility | 0.8095276384290021 | 0.548806 | 0.548806 | 0.54979 | +0.000984 | 100.0 | 6.798 | 6.746 |

Total loss saved 569.180067, fix 771.908811.

## baseline, price re-solved at the same parameters (fix, fixed_h0 closure)

price 0.7794387507545272, population scale 1.0001768654431398, renewal residual -1.2214204092586556e-09, housing residual 0.0

| moment | target | saved | fix (GE price) | diff | weight | loss saved | loss fix |
|---|---|---|---|---|---|---|---|
| initial_normalization | 2.1 | 2.1 | 2.1 | +3.39e-10 |  | 0 | 0 |
| cps_childlessness | 0.198279 | 0.202193 | 0.202172 | -2.13e-05 | 35532.3042455214 | 0.5444 | 0.5385 |
| cps_exactly_one | 0.213655 | 0.217103 | 0.217128 | +2.56e-05 | 26952.820824310795 | 0.3203 | 0.3251 |
| nchs_mean_age | 25.9763 | 25.9628 | 25.9632 | +0.00038 | 139.82806784479274 | 0.02523 | 0.02383 |
| nchs_share30 | 0.249278 | 0.237511 | 0.237507 | -3.69e-06 | 0.0 | 0 | 0 |
| wealth_earnings | 4.45839 | 4.58349 | 4.58314 | -0.00035 | 7.595098472533724 | 0.1189 | 0.1182 |
| bequest_wealth | 0.00729102 | 0.00712342 | 0.00712367 | +2.52e-07 | 5165289.256198346 | 0.1451 | 0.1447 |
| old_dispersion | 3.51594 | 3.30769 | 3.30762 | -6.89e-05 | 0.0 | 0 | 0 |
| mean_rooms | 5.72943 | 5.81601 | 5.8151 | -0.000912 | 128.02070205233477 | 0.9596 | 0.9395 |
| ownership_30_55 | 0.67626 | 0.669864 | 0.67024 | +0.000376 | 2339.3623724673616 | 0.09572 | 0.0848 |
| first_birth_rooms | 1.465 | 1.26385 | 1.26471 | +0.00086 | 137.5652749002964 | 5.566 | 5.519 |
| family_rooms | 0.3851 | 0.293947 | 0.294032 | +8.54e-05 | 0.0 | 0 | 0 |
| recent_parent_ownership | 0.127608 | 0.126658 | 0.124425 | -0.00223 | 27055.822957508266 | 0.02444 | 0.2742 |
| early_fertility | 0.809528 | 0.552574 | 0.552541 | -3.35e-05 | 100.0 | 6.603 | 6.604 |

Total loss saved 14.402424, fix at re-solved price 14.571723.
