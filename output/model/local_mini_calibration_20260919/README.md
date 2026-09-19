# Local mini-calibration — September 19, 2026

Completed: baseline replay, twelve nearby candidates in two waves of six processes, and a verified independent repeat of the selected candidate. Search plus repeat: 13.90 minutes; sampled combined process-tree peak RSS: 11.74 GiB. Baseline stationary solve: 217.01 seconds. Setup and the earlier baseline run are excluded from the search timer.

Baseline loss 179.298424248025; selected diagnostic loss 177.566967759517 (0.966% lower). Selected case: `joint_2`. This is a bounded neighborhood test of the frozen paper-baseline initial objective, not a converged production calibration or a hardware-only speed comparison. No matched serial run of all twelve candidates was performed.

The contract retains twelve scored moments, nine structural coordinates, and the separate completed-fertility normalization. Moment count alone does not establish local rank identification. The objective and source pins are unchanged; beta uses the author cap 0.99 and child housing requirement uses the cap 2.3. No policy closure or live production calibration was changed.

## Complete target comparison

| Moment | Target | Baseline | Selected | Selected gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|---:|
| Initial model completed fertility | 2.1 | 2.1000008 | 2.1000012 | 1.2185078e-06 | — | — |
| Childless women, ages 40–44 | 0.19827875 | 0.19625838 | 0.19592793 | -0.0023508222 | 35532.304 | 0.19636449 |
| Exactly one child among mothers, ages 40–44 | 0.21365533 | 0.21643897 | 0.21652078 | 0.0028654525 | 26952.821 | 0.22130471 |
| Period mean first-birth age | 25.976264 | 26.095246 | 26.063249 | 0.086984666 | 139.82807 | 1.0579856 |
| First births at age 30+ | 0.24927801 | 0.23215522 | 0.23043785 | -0.018840161 | 13866.065 | 4.9217829 |
| Wealth / annual gross labor earnings | 6.1458614 | 4.893427 | 4.83202 | -1.3138414 | 7.5950985 | 13.110502 |
| Annual bequests / aggregate wealth | 0.0088 | 0.0085824095 | 0.0086073081 | -0.00019269192 | 5165289.3 | 0.1917881 |
| Old wealth/income p90 / median, ages 76–84 | 3.5159351 | 4.4998003 | 4.384293 | 0.86835795 | 10.616362 | 8.0052199 |
| Mean occupied rooms, capped at 9 | 5.5610974 | 6.4246524 | 6.4054046 | 0.84430719 | 128.0207 | 91.260151 |
| Ownership, heads 30–55 | 0.64833403 | 0.53989617 | 0.52677358 | -0.12156046 | 2339.3624 | 34.568627 |
| First-birth room response, −1 to +3 | 0.72024626 | 0.98080602 | 0.98307248 | 0.26282622 | 137.56527 | 9.5026819 |
| Rooms: 3+ versus 1–2 resident children (model dependent proxy) | 0.34706693 | 0.17204278 | 0.17377595 | -0.17329098 | 280.52808 | 8.4241919 |
| Recent-parent ownership gap | 0.16289551 | 0.14410266 | 0.14787236 | -0.015023154 | 27055.823 | 6.1063682 |

Full provenance, sample definitions and uncertainty fields: [target comparison](search_v2/comparison_target_fit.csv). The normalization row is separate and unscored.

## Parameters and restrictions

| Parameter | Baseline | Selected | Lower | Upper | Near a bound |
|---|---:|---:|---:|---:|---|
| beta_annual | 0.99 | 0.9895 | 0.94 | 0.99 | True |
| kappa_fert | 0.33773423411167025 | 0.3360455629411119 | 0.02 | 50.0 | True |
| kappa_fert_continuation | 0.39785645171756406 | 0.39984573397615186 | 0.02 | 50.0 | True |
| chi | 1.0496534423047694 | 1.0444051750932455 | 0.1 | 5.0 | False |
| H0 | 8.11210048056786 | 8.07153997816502 | 0.2 | 80.0 | False |
| theta0 | 0.08105103333987912 | 0.08064577817317972 | 0.0 | 8.0 | False |
| theta1 | 0.08520356663830632 | 0.08477754880511479 | 0.02 | 16.0 | True |
| first_birth_fixed_cost | 0.26576477618628114 | 0.2644359523053497 | 0.0 | 8.0 | False |
| h_P | 2.3 | 2.2885 | 0.1 | 2.3 | True |
| hbar_child_rooms | 0.0 | 0.0 |  |  | False |
| psi_child | 0.1489153145785917 | 0.14824063921279804 |  |  | False |
| payroll_tax | 0.179 | 0.179 |  |  | False |
| pension_period | 2.0463613896121218 | 2.0463613896121218 |  |  | False |
| housing_supply_elasticity | 0.63 | 0.63 |  |  | False |
| tenure_choice_kappa | 0.005 | 0.005 |  |  | False |
| alpha_cons | 0.733 | 0.733 |  |  | False |
| sigma | 2.0 | 2.0 |  |  | False |

Blank bounds identify normalized, derived or fixed objects rather than search coordinates. Full classifications and transformations: [parameter comparison](search_v2/comparison_parameters.csv). Near-bound flags follow the frozen scorer, except beta is evaluated against the actual 0.99 author cap.

## Verification and performance evidence

- Every one of the twelve proposals passed the frozen scored-run checks. Baseline full score agrees with the saved replay within absolute 1e-8 and relative 1e-10 tolerances (loss difference about 2e-13). The selected full score passed the same repeated-run comparison.
- All 17 standard diagnostic graphs are present for the baseline and selected candidate. Both packets were visually inspected. Selected housing-market relative residual is about 9.6e-7. Policy plots retain the existing discrete-grid irregularities; this test is not a new audit of model economics.
- Single-thread NumPy/BLAS/Numba per process; six concurrent processes; 420-second case ceiling; 25-minute controller ceiling; 24 GiB combined descendant RSS ceiling. No limit fired in the successful search.
- Timings are observed on the native M5 Pro installation. No old-Mac or cluster mini-calibration was run for this test.

Selected fiscal/normalization checks:

- Completed fertility: 2.10000121851, target 2.1.
- Pension relative gap: 3.19223e-12.
- Property-tax rebate relative gap: 8.19602e-07.

## Artifacts and reproduction

- [Controller summary](search_v2/summary.json), [selected receipt](search_v2/selected.json), [all target rows](search_v2/all_target_fits.csv), [all parameter rows](search_v2/all_parameters.csv).
- Full case checkpoints, logs and the original 17 graphs: `search_v2/cases/joint_2/case/evaluation/`. Baseline evidence: `run_v3/cases/baseline_smoke/`.
- Staging: `code/model/tools/stage_local_paper_calibration.py`; launcher: `code/model/tools/run_local_mini_calibration.py`; launcher tests: `code/model/tools/test_run_local_mini_calibration.py`.
- To regenerate the selected diagnostic packet and score, rerun the frozen scored driver with `search_v2/cases/joint_2.json`, the staged template/helper/joint paths in `staged/README.json`, and a fresh output directory. The driver generates all 17 graphs automatically. Install the temporary pickle aliases as described in `docs/workflow/local_laptop_setup.md` before replay; remove them afterward.
- Frozen objective SHA-256: `4440ea07f4de957740ca6c04961d2806d9b9ef782c7a0e7dad4ce73e1db651b1`. The exact source inventory and normalized seed are verified by the scored wrapper before model execution.

Earlier launch attempts failed on archived relative paths and NumPy pickle import compatibility. They are retained for audit in `run`, `run_v2`, `run_v3`, and `search`; only the baseline from `run_v3` and the completed `search_v2` candidates are reported above. No failed numerical candidate was silently promoted or scored.

After installing the temporary compatibility module, this command regenerates
all selected-case graphs and the full score (use a fresh output path):

```sh
code/model/.venv/bin/python -B output/model/local_mini_calibration_20260919/staged/recipe/run_e5f_joint_rebated_initial_scored.py --helper "$PWD/output/model/local_mini_calibration_20260919/staged/recipe/run_e5f_rebated_initial_overnight.py" --helper-sha256 d9aa97b890442d45971ec622b4b41687da10ffa26eedaf0198f352c7e6ecb790 --joint "$PWD/output/model/local_mini_calibration_20260919/staged/recipe/run_e5f_joint_rebated_initial_probe.py" --joint-sha256 9eee3bca39f2a98f4a58cf18196d695b6a9db3e93ef18f8eaa2bf4dbb1243bbb --template "$PWD/output/model/local_mini_calibration_20260919/staged/template" --proposal "$PWD/output/model/local_mini_calibration_20260919/search_v2/cases/joint_2.json" --output "$PWD/output/model/local_mini_calibration_20260919/selected_replay"
```
