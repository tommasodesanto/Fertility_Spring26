# Estate-A cap two at the verified binary winner

One fixed-parameter diagnostic, not a recalibration or adopted specification. The sole menu change is the intended-birth cap from one to two per four-year period. The ten estimated coordinates and every other fixed economic input are held at the verified binary chain-1 winner. Its derived population-one housing coefficient, `H0=6.1245552591467405`, is held fixed while the native stationary GE solves the birth-renewal price and allows population scale to adjust. Estate A, post-interest timing, soft credit, target 4.45838713455674, all weights and numerical gates are retained. The bequest target has the already noted scope/recipient mismatch.

**Result:** native exit `0` in 94.98 seconds; nine lifecycle solves. New-contract weighted loss **1386.601504576907**, compared with **21.275413361071312** for the cap-one winner under the same target and weight contract. The native old-wealth-target diagnostic loss is **1407.6485378721213**; it is a distinct scoring contract and is not used for that comparison. This comparison changes only the birth menu with fixed model inputs and H0; it does not compare separate recalibrations.

Price 0.999335757238; household population scale 1.378613699598; implied H0 at population one 4.442546349956. Birth-renewal residual -1.841e-07; housing residual 0.000e+00; PAYGO residual 8.634e-14. Native selected/repeat target and parameter tables match exactly, as do SHA-256 hashes of all 17 standard figures; `verification_receipt.json` records the hashes. The packet also contains eight policy and seven aggregate plots.

## Complete target fit

| Moment | Target | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| initial_normalization | 2.1 | 2.099999613459218 | -3.8654078204913844e-07 |  |  |
| cps_childlessness | 0.19827875100684264 | 0.26240143302399127 | 0.06412268201714863 | 35532.3042455214 | 146.09882735113194 |
| cps_exactly_one | 0.21365532522014702 | 0.18373304093069032 | -0.029922284289456702 | 26952.820824310795 | 24.132022072394665 |
| nchs_mean_age | 25.976263860992496 | 28.009497178375813 | 2.033233317383317 | 139.82806784479274 | 578.0545071930503 |
| nchs_share30 | 0.2492780130410667 | 0.3733799858578359 | 0.12410197281676918 | 0.0 | 0.0 |
| wealth_earnings | 4.45838713455674 | 5.131117370885647 | 0.6727302363289072 | 7.595098472533724 | 3.437283114084193 |
| bequest_wealth | 0.007291023472616158 | 0.006264556980557837 | -0.0010264664920583205 | 5165289.256198346 | 5.442321587389018 |
| old_dispersion | 3.51593508651872 | 3.2256506225378447 | -0.29028446398087526 | 0.0 | 0.0 |
| mean_rooms | 5.729434240102641 | 4.78740990362855 | -0.9420243364740912 | 128.02070205233477 | 113.60683207037748 |
| ownership_30_55 | 0.6762604168538028 | 0.6219412632857847 | -0.054319153568018086 | 2339.3623724673616 | 6.902453474817177 |
| first_birth_rooms | 1.465 | 1.3327239721634658 | -0.13227602783653425 | 137.5652749002964 | 2.406972398285271 |
| family_rooms | 0.38509964969278165 | 0.21759770868009465 | -0.167501941012687 | 0.0 | 0.0 |
| recent_parent_ownership | 0.12760836356692162 | -0.007470528511200003 | -0.13507889207812163 | 27055.822957508266 | 493.66885412151805 |
| early_fertility | 0.8095276384290021 | 0.45103870981531524 | -0.35848892861368686 | 100.0 | 12.85143119385891 |

## All 31 parameter records

The native parameter table reports advisory original beta bounds 0.94–0.99. The retained Estate-A search restriction was 0.93–0.99; beta is held at 0.9444437995936287 here. Every numeric estimate, including H0, equals the binary selected packet exactly.

| Parameter | Estimate | Lower | Upper | Near bound | Status |
|---|---:|---:|---:|---|---|
| H0 | 6.1245552591467405 | 0.2 | 80.0 | False | supplied fixed housing supply coefficient; reference bounds advisory |
| beta_annual | 0.9444437995936287 | 0.94 | 0.99 | False | supplied primitive; reference calibration bounds advisory |
| chi | 1.088196976146629 | 0.1 | 5.0 | False | supplied primitive; reference calibration bounds advisory |
| first_birth_fixed_cost | 0.2827067667825205 | 0.0 | 8.0 | False | supplied primitive; reference calibration bounds advisory |
| kappa_fert | 0.11612031168396098 | 0.02 | 50.0 | True | supplied primitive; reference calibration bounds advisory |
| kappa_fert_continuation | 0.4164292250768279 | 0.02 | 50.0 | True | supplied primitive; reference calibration bounds advisory |
| theta0 | 0.14389227339225366 | 0.0 | 8.0 | False | supplied primitive; reference calibration bounds advisory |
| delta_alpha_jump | 0.0 |  |  |  | fixed zero under experimental utility contract |
| child_benefit_curvature | 0.10313883854657178 | 0.0 | 0.8 | False | supplied primitive; reference calibration bounds advisory |
| tenure_choice_kappa | 0.015769047253080707 | 0.001 | 0.1 | False | supplied primitive; reference calibration bounds advisory |
| psi_child | 0.17788162876169034 | 0.01 | 0.5 | False | supplied primitive; reference calibration bounds advisory |
| child_benefit_CRRA_coefficient | 0.15953512417243712 |  |  |  | derived from supplied benefit and curvature |
| theta1 | 0.008193084126995582 |  |  |  | fixed external restriction |
| sigma | 2.0 |  |  |  | fixed |
| alpha_cons | 0.733 |  |  |  | fixed CEX childless expenditure share |
| delta_alpha | 0.0 |  |  |  | fixed zero under experimental utility contract |
| h_P | 2.5426884388222315 | 0.1 | 2.6 | False | supplied primitive; reference calibration bounds advisory |
| utility_reference_rent | 0.11046592704873838 |  |  |  | retained inactive normalization; compensation off |
| q_annual | 0.020000000000000018 |  |  |  | author-retained 2% annual real rate |
| financed_share | 0.8 |  |  |  | inherited credit contract |
| housing_supply_elasticity | 0.63 |  |  |  | fixed provisional external mapping |
| payroll_tax | 0.08028070961950022 |  |  |  | derived from adopted pension ratio |
| pension_period | 0.917784047463731 |  |  |  | balanced PAYGO |
| annual_depreciation | 0.01416143718381309 |  |  |  | adopted |
| period_depreciation | 0.05545379079326218 |  |  |  | compounded |
| annual_property_tax | 0.010598360773872594 |  |  |  | adopted |
| period_property_tax | 0.042393443095490375 |  |  |  | linear period convention |
| selling_cost | 0.06 |  |  |  | retained |
| rental_cap | 6.0 |  |  |  | retained provisional |
| wealth_grid_nodes | 120.0 |  |  |  | retained exact grid |
| income_states | 9.0 |  |  |  | retained B15 |

## Evidence and reproduction

- `launch_receipt.json`, `terminal_receipt.json`, and `run.log` record the single bounded local run.
- `cases/20261004T150401752556Z_701c4377/` is the complete saved case with `input_contract.json`, full CSVs, native exact-repeat packet, and 32 plots.
- Binary source: `../estate_a_recovery_20261004_v1/collection/binary/provenance/search_completed.json`, SHA-256 `2e104f260050b802370de3d8bdd0d0792c8a7af9e13be849fadfdd6171106f71`; verified selected and repeat receipts are adjacent.
- Zero-solve preflight: `output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/experiments/birth_count_choice/run_cap2_at_binary_winner.py --preflight`.
- Single-run command: `NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/experiments/birth_count_choice/run_cap2_at_binary_winner.py` (external process limit 1,800 seconds; native limit 1,500 seconds).
