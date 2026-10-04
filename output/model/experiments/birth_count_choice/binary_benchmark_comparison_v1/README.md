# One-birth Estate-A recovery versus working revised-timing anchor

Reproduce from saved reports only: `python3 build_comparison.py`. This performs no model solve. Both selected points passed fresh native and exact-repeat checks; neither has an optimizer-convergence, grid-convergence, or identification-rank certificate. The Estate-A point is experimental and is not adopted as the paper baseline.

## Native and cross-contract scores

| Saved model moments | New wealth target 4.458387 | Old wealth target 6.926584 |
|---|---:|---:|
| One-birth Estate-A chain 1 | 21.275413361 (native) | 46.146943506 (arithmetic rescore) |
| Working revised-timing chain 13 | 53.064444038 (arithmetic rescore) | 13.771131463 (native) |

An arithmetic rescore replaces the target in `weight × (saved model moment − target)²`. It does not re-estimate parameters or solve either economy. The two rows differ in the estate utility and death-flow definition, bequest measurement, wealth target, beta search bound, and all ten estimated parameter values. Thus neither the native losses across columns nor cross-scored losses identify a causal Estate-A effect or establish which economic specification is better.

The working chain-13 anchor uses the author-adopted post-interest transaction timing and the old PSID wealth/earnings target. Both points retain the same consumption/housing flow-utility form: constant consumption share `α=0.733`, physical parent housing floor `h_P>0`, and equivalence scale `e(m)=((2+0.7m)/2)^0.7`. The separate normalized-CES, parent-dependent-share, `h_P=0` experiment is not either point here. The experimental chain-1 point retains the one-birth menu, soft financing with financed share 0.8, and the other fixed economic inputs. The old estate values its housing leg at gross `Ph′`; Estate A instead uses terminal estate wealth `W = b′ + (1−ψ)Ph′` with selling cost `ψ = 0.06` and no extra interest on `b′`. The net housing leg enters utility, native death-flow accounting, and the empirical bequest observer. The new aggregate wealth/earnings target is 4.45838713455674 versus 6.92658379107299, with unchanged numerical weight 7.595098472533724. The lower search bound for annual beta changes from 0.94 to 0.93. The recovered point re-estimates all ten free coordinates. Entry wealth/income distributions, earnings, transfers/floors, housing/financing rules, other targets and weights are reported unchanged in the recovery status and source receipt; they were not individually replayed for this comparison. SCF wealth scope and estate recipient/creditor mapping remain provisional.

## Complete moment comparison

`new` denotes Estate A; `old` denotes the working anchor. Empty weight/loss is the normalization row; zero-weight rows are validation moments. Full precision and both cross-contract gap/loss columns are in [moments_full.csv](moments_full.csv).

| moment | role | new_target | old_target | weight | new_model | new_gap_new_target | new_loss_new_target | old_model | old_gap_old_target | old_loss_old_target |
|---|---|---|---|---|---|---|---|---|---|---|
| initial_normalization | normalization | 2.1 | 2.1 | — | 2.10000013 | 1.32416349e-07 | — | 2.1 | -1.4313839e-10 | — |
| cps_childlessness | scored | 0.198278751 | 0.198278751 | 35532.3042 | 0.19937083 | 0.00109207865 | 0.0423770971 | 0.199280303 | 0.00100155206 | 0.0356426866 |
| cps_exactly_one | scored | 0.213655325 | 0.213655325 | 26952.8208 | 0.219833417 | 0.00617809196 | 1.02875737 | 0.215670678 | 0.00201535235 | 0.109472793 |
| nchs_mean_age | scored | 25.9762639 | 25.9762639 | 139.828068 | 25.973636 | -0.00262789247 | 0.000965627305 | 25.9676516 | -0.00861227749 | 0.0103712329 |
| nchs_share30 | validation | 0.249278013 | 0.249278013 | 0 | 0.238300649 | -0.0109773639 | 0 | 0.232470048 | -0.0168079646 | 0 |
| wealth_earnings | scored | 4.45838713 | 6.92658379 | 7.59509847 | 5.0291101 | 0.570722962 | 2.47391117 | 6.74051972 | -0.186064071 | 0.262941083 |
| bequest_wealth | scored | 0.00729102347 | 0.00729102347 | 5165289.26 | 0.00646095316 | -0.000830070308 | 3.55897064 | 0.00683518502 | -0.000455838452 | 1.07328871 |
| old_dispersion | validation | 3.51593509 | 3.51593509 | 0 | 3.24958085 | -0.266354236 | 0 | 2.9036903 | -0.612244788 | 0 |
| mean_rooms | scored | 5.72943424 | 5.72943424 | 128.020702 | 5.65136914 | -0.0780650985 | 0.780178591 | 5.88639743 | 0.15696319 | 3.15410273 |
| ownership_30_55 | scored | 0.676260417 | 0.676260417 | 2339.36237 | 0.679406098 | 0.00314568128 | 0.0231487176 | 0.671391172 | -0.00486924511 | 0.0554652244 |
| first_birth_rooms | scored | 1.465 | 1.465 | 137.565275 | 1.24629123 | -0.218708772 | 6.58023227 | 1.3775252 | -0.0874748032 | 1.05262764 |
| family_rooms | validation | 0.38509965 | 0.38509965 | 0 | 0.302526362 | -0.082573288 | 0 | 0.299060671 | -0.0860389784 | 0 |
| recent_parent_ownership | scored | 0.127608364 | 0.127608364 | 27055.823 | 0.128658659 | 0.00105029543 | 0.0298458329 | 0.123692354 | -0.00391600917 | 0.414904503 |
| early_fertility | scored | 0.809527638 | 0.809527638 | 100 | 0.549584836 | -0.259942802 | 6.75702605 | 0.533804682 | -0.275722956 | 7.60231486 |

## Complete parameter comparison

The ten `estimated/free` bounds come from the Estate-A start contract and old native parameter table. `H0` is derived for household scale `N0=1`; its displayed interval is advisory, not an optimization bound. `Near` is each native report's flag, not a new threshold calculation. Full native status wording is in [parameters_full.csv](parameters_full.csv).

| parameter | role | new_estimate | new_lower | new_upper | new_native_near_bound | old_estimate | old_lower | old_upper | old_native_near_bound |
|---|---|---|---|---|---|---|---|---|---|
| H0 | derived | 6.1245552591467405 | 0.2 | 80.0 | False | 6.40569359569417 | 0.2 | 80.0 | False |
| beta_annual | estimated/free | 0.9444437995936287 | .93 | .99 | False | 0.9663191380998087 | 0.94 | 0.99 | False |
| chi | estimated/free | 1.088196976146629 | 0.1 | 5.0 | False | 1.0500762402240174 | 0.1 | 5.0 | False |
| first_birth_fixed_cost | estimated/free | 0.2827067667825205 | 0.0 | 8.0 | False | 0.3045994545418478 | 0.0 | 8.0 | False |
| kappa_fert | estimated/free | 0.11612031168396098 | 0.02 | 50.0 | True | 0.11652185618155607 | 0.02 | 50.0 | True |
| kappa_fert_continuation | estimated/free | 0.4164292250768279 | 0.02 | 50.0 | True | 0.40070847699255835 | 0.02 | 50.0 | True |
| theta0 | estimated/free | 0.14389227339225366 | 0.0 | 8.0 | False | 0.10097400014250629 | 0.0 | 8.0 | False |
| delta_alpha_jump | fixed/retained | 0.0 | — | — | — | 0.0 | — | — | — |
| child_benefit_curvature | estimated/free | 0.10313883854657178 | 0.0 | 0.8 | False | 0.0629608477756522 | 0.0 | 0.8 | False |
| tenure_choice_kappa | estimated/free | 0.015769047253080707 | 0.001 | 0.1 | False | 0.014123854930856623 | 0.001 | 0.1 | False |
| psi_child | estimated/free | 0.17788162876169034 | 0.01 | 0.5 | False | 0.17892072066041628 | 0.01 | 0.5 | False |
| child_benefit_CRRA_coefficient | derived | 0.15953512417243712 | — | — | — | 0.1676557204030058 | — | — | — |
| theta1 | fixed/retained | 0.008193084126995582 | — | — | — | 0.008193084126995582 | — | — | — |
| sigma | fixed/retained | 2.0 | — | — | — | 2.0 | — | — | — |
| alpha_cons | fixed/retained | 0.733 | — | — | — | 0.733 | — | — | — |
| delta_alpha | fixed/retained | 0.0 | — | — | — | 0.0 | — | — | — |
| h_P | estimated/free | 2.5426884388222315 | 0.1 | 2.6 | False | 2.593759507364224 | 0.1 | 2.6 | True |
| utility_reference_rent | fixed/retained | 0.11046592704873838 | — | — | — | 0.11046592704873838 | — | — | — |
| q_annual | fixed/retained | 0.020000000000000018 | — | — | — | 0.020000000000000018 | — | — | — |
| financed_share | fixed/retained | 0.8 | — | — | — | 0.8 | — | — | — |
| housing_supply_elasticity | fixed/retained | 0.63 | — | — | — | 0.63 | — | — | — |
| payroll_tax | derived | 0.08028070961950022 | — | — | — | 0.08028070961950022 | — | — | — |
| pension_period | derived | 0.917784047463731 | — | — | — | 0.917784047463731 | — | — | — |
| annual_depreciation | fixed/retained | 0.01416143718381309 | — | — | — | 0.01416143718381309 | — | — | — |
| period_depreciation | derived | 0.05545379079326218 | — | — | — | 0.05545379079326218 | — | — | — |
| annual_property_tax | fixed/retained | 0.010598360773872594 | — | — | — | 0.010598360773872594 | — | — | — |
| period_property_tax | derived | 0.042393443095490375 | — | — | — | 0.042393443095490375 | — | — | — |
| selling_cost | fixed/retained | 0.06 | — | — | — | 0.06 | — | — | — |
| rental_cap | fixed/retained | 6.0 | — | — | — | 6.0 | — | — | — |
| wealth_grid_nodes | fixed/retained | 120.0 | — | — | — | 120.0 | — | — | — |
| income_states | fixed/retained | 9.0 | — | — | — | 9.0 | — | — | — |

## Provenance and limits

New selected chain 1: array 19141024, inventory `974a14e243da6a2ad0572bb9825b47ab349828f9144cbbad69e74b40a4408b22`, target fingerprint `c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70`, weight fingerprint `f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4`. Old selected chain 13: 48-chain collection, target fingerprint `db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`, weight fingerprint `2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`. [Input hashes](source_hashes.json) pin the exact six saved files used here.

Sources: `CALIBRATION_STATUS.md` (working anchor and Estate-A recovery sections); `output/model/experiments/birth_count_choice/estate_a_recovery_20261004_v1/collection/RESULTS.md` and `binary/provenance/start_contract.json`; `output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/collection.json`. These are local saved reports, not a new numerical verification.
