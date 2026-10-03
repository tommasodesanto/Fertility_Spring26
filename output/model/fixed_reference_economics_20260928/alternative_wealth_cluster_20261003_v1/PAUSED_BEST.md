# Paused experimental new-wealth search: best retained candidate

Best saved search loss: **22.141841386410267**, chain **6**, label **0060_nm**.
This is a provisional search checkpoint under the experimental new wealth target. Fresh native selected-point verification, exact repeat, and optimizer convergence are pending; this is not an adopted calibration or certified best point. No model solves were performed for this readout.

Source: [pause_20261003.json](pause_20261003.json), object `chains[chain=6].files.best_so_far.best`.
Source SHA-256: `f5eea8bab420330e173afdad4b519c400e9c0138dd78a543539f0005e8cf5657`.
Pause snapshot verified epoch: `1791061133.5716393` (2026-10-03T20:58:53.571639+00:00).
Target fingerprint from chain-6 start contract: `c7a3d185668122e508a6c322bc5ef0715ebb0ecb23948c8d9b184ee25d1cde70`.
Weight fingerprint from chain-6 start contract: `f762ebb5684ab30487b3b8b64fc10977fda396b520035d91c0c5c803255f88e4`.
Saved candidate report identity: `/work/results/run/0060_nm/phase_b_ge/selected_root`.
Price: 0.7928926715538719; derived H0: 6.04347890577259; population normalization: 1.0.

The empirical aggregate wealth/earnings target is 4.45838713455674. Bounds below are the actual chain-6 start-contract bounds. Near-bound means the estimate is within 1% of the total permitted interval from either endpoint; it is computed here, not inferred from an earlier calibration.

## Complete target fit

| Moment | Target | Model | Gap | Weight | Loss contribution | Role |
| --- | --- | --- | --- | --- | --- | --- |
| initial_normalization | 2.100 | 2.100 | 4.267e-09 | — | — | normalization |
| cps_childlessness | 0.198 | 0.204 | 0.006 | 35532.304 | 1.353 | scored |
| cps_exactly_one | 0.214 | 0.215 | 0.002 | 26952.821 | 0.080 | scored |
| nchs_mean_age | 25.976 | 26.036 | 0.060 | 139.828 | 0.504 | scored |
| nchs_share30 | 0.249 | 0.241 | -0.008 | 0.000 | 0.000 | validation |
| wealth_earnings | 4.458 | 5.029 | 0.571 | 7.595 | 2.477 | scored |
| bequest_wealth | 0.007 | 0.007 | -4.558e-04 | 5165289.256 | 1.073 | scored |
| old_dispersion | 3.516 | 3.173 | -0.343 | 0.000 | 0.000 | validation |
| mean_rooms | 5.729 | 5.629 | -0.100 | 128.021 | 1.288 | scored |
| ownership_30_55 | 0.676 | 0.677 | 9.062e-04 | 2339.362 | 0.002 | scored |
| first_birth_rooms | 1.465 | 1.234 | -0.231 | 137.565 | 7.358 | scored |
| family_rooms | 0.385 | 0.295 | -0.090 | 0.000 | 0.000 | validation |
| recent_parent_ownership | 0.128 | 0.122 | -0.006 | 27055.823 | 0.904 | scored |
| early_fertility | 0.810 | 0.543 | -0.267 | 100.000 | 7.104 | scored |

## All ten estimated parameters

| Parameter | Estimate | Actual lower bound | Actual upper bound | Within 1% of bound |
| --- | --- | --- | --- | --- |
| beta_annual | 0.945 | 0.930 | 0.990 | False |
| chi | 1.087 | 0.100 | 5.000 | False |
| first_birth_fixed_cost | 0.323 | 0.000 | 8.000 | False |
| kappa_fert | 0.122 | 0.020 | 50.000 | True |
| kappa_fert_continuation | 0.406 | 0.020 | 50.000 | True |
| theta0 | 0.106 | 0.000 | 8.000 | False |
| h_P | 2.497 | 0.100 | 2.600 | False |
| child_benefit_curvature | 0.107 | 0.000 | 0.800 | False |
| tenure_choice_kappa | 0.014 | 0.001 | 0.100 | False |
| psi_child | 0.187 | 0.010 | 0.500 | False |

Checks: 14 complete fit rows, ten scored moments and ten free coordinates; saved loss equals the summed scored contributions, and each scored gap/contribution matches its saved target, model and weight. These arithmetic checks do not replace native verification.
