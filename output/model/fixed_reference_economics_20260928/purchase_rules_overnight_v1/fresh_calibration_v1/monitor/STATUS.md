# Fresh 80% calibration monitor

Snapshot 2026-10-02T19:29:10.676413+00:00. Remote root: `/scratch/td2248/projects/purchase_fresh_calibration_v1`.
Local design SHA-256: `463fce700a4d3279afa3a0e837b6bcc271589206435e3b117ad713442e206b20`. Remote design SHA-256: `463fce700a4d3279afa3a0e837b6bcc271589206435e3b117ad713442e206b20`.

Slurm job: `19040483`; active rows: 0.

| Arm | Valid scored cases | Budget-uncomputed | Numerically inadmissible | Other uncomputed | Terminal / 12 | Fresh postchecks / 12 |
|---|---:|---:|---:|---:|---:|---:|
| hard | 909 | 27 | 0 | 0 | 12 | 12 |
| quarter | 995 | 21 | 0 | 0 | 12 | 12 |

A completed model attempt is counted as valid only when its case status is `passed` and it has a computed loss. The other statuses are listed separately.

## Hard

Verified candidate: loss **88.588403**, slot 8; [14-row target fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor/hard_verified_slot8_target_fit.csv), [31-parameter bounds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor/hard_verified_slot8_parameters.csv).
Provisional candidate: loss **88.588403**, slot 8; [14-row target fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor/hard_provisional_slot8_target_fit.csv), [31-parameter bounds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor/hard_provisional_slot8_parameters.csv).

## Quarter

Verified candidate: loss **48.319938**, slot 19; [14-row target fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor/quarter_verified_slot19_target_fit.csv), [31-parameter bounds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor/quarter_verified_slot19_parameters.csv).
Provisional candidate: loss **48.319938**, slot 19; [14-row target fit](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor/quarter_provisional_slot19_target_fit.csv), [31-parameter bounds](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor/quarter_provisional_slot19_parameters.csv).

Fresh selected-point verification confirms the quoted numerical point. It does not certify optimizer convergence or adopt the experimental economics.
