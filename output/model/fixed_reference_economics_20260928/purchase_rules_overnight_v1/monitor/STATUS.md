# Overnight calibration search monitor

Final snapshot 2026-10-02T08:55:58.683501+00:00. All submitted Torch and local jobs have terminal receipts. No Slurm launcher exit was nonzero. The two best points passed fresh native selected-point postchecks with exactly matching losses, 14 target rows, and effective parameters. They remain experimental and **optimizer convergence is not certified**.

| Purchase rule | Best loss | Source | Search stop reason |
|---|---:|---|---|
| Hard | **97.011220** | Torch continuation chain 11, case 0047_nm | native or global evaluation budget exhausted |
| Quarter saving | **51.556036** | Local chain 54, case 0094_nm | native evaluation budget exhausted |

Full tables: [hard 14-row target fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/monitor/hard_verified_best_chain11_target_fit.csv), [hard 31-parameter bounds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/monitor/hard_verified_best_chain11_parameters.csv), [quarter 14-row target fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/monitor/quarter_verified_best_chain54_target_fit.csv), [quarter 31-parameter bounds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/monitor/quarter_verified_best_chain54_parameters.csv).

| Source | Terminal / selected postchecked | Hard valid / budget-uncomputed | Quarter valid / budget-uncomputed | Other scientific status |
|---|---:|---:|---:|---|
| Original Torch 19009131 | 48 / 48 | 891 / 24 | 1533 / 24 | — |
| Five Torch continuation arrays | 25 / 25 | 813 / 40 | 301 / 8 | — |
| Broader regions Torch 19019947 | 16 / 15 | 121 / 12 | 105 / 9 | 13 quarter numerically inadmissible; quarter slot 15 no selected candidate |
| Local original and resume | 10 / 9 | 295 / 5 | 563 / 5 | Hard chain 50 no selected candidate |
| Local restart_v2 | 5 / 5 | 199 / 18 | 79 / 3 | — |

Across all saved case receipts: hard 2,319 valid full GE evaluations and 99 budget-exhausted uncomputed attempts; quarter 2,581 valid, 49 budget-exhausted uncomputed, and 13 numerically inadmissible attempts. These sums include repeated parameter evaluations across origins, so they are not counts of unique parameter vectors. No active calibration search job or checkpoint remains; this monitor does not report policy-job status.

Broader quarter slot 15 illustrates the count distinction: its `completed_full_ge=14` records evaluator attempts, but all 13 completed cases were numerically inadmissible because the price root was unbracketed, and the final case stopped at the budget. None had a computed base loss, so no candidate was selected.
