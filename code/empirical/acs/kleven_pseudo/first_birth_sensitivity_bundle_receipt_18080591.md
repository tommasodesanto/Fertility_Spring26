# ACS sensitivity bundle receipt: job 18080591

Slurm job `18080591` completed successfully on 2026-09-20 in 2:34 with 4 CPUs and 64 GB. Both independent stages returned exit code 0. The prior staging failure, job `18080547`, and its outputs remain preserved.

The short-window housing stage used verified true-ACS source years 2005--2019, fixed implied support cohorts 2007--2016, event times (-2,-1,0,1,2,3), and (-2) as reference. It reused the deterministic verified source-key join and did not rematch observations. The prefit gate passed all 360 outcome-by-state-by-gender-by-cohort groups, with positive outcome-valid support in all six event cells. The joined analysis had 571,718 rows, and the six fits wrote level curves, +3 minus (-1) full-covariance contrasts, fit checkpoints, weighted event-(-2) baselines, and Kish effective sample sizes.

The transformed second-birth stage found exact raw-field concordance for 2,190,987 packet rows: sex, age, race, education, marital status, and state had zero mismatches and zero one-sided missing values. Its builder receipt reports 106,212 strict eligible rows, 5,478 event-zero anchors, 166,435 one-child donors, 1,453 full-pre-window anchors, 5,047 reference anchors, and 17,552 donor targets. Exact matching with replacement produced 31,312 matched links. This remains a support diagnostic; it does not establish causal identification or a housing estimate.

Compact local copies of the receipts are under:

`code/empirical/acs/kleven_pseudo/output/first_birth_sensitivity_bundle_18080591/`

The remote output directories are recorded in `short_window/bundle_stage_status.txt`. The short-window output includes `short_window_support_gate.csv`, `short_window_baseline_ess.csv`, `short_window_contrasts.csv`, `short_window_curves.csv`, and `short_window_result_receipt.json`; the transformed-support output includes `diagnostic_manifest.json`, `raw_key_concordance.csv`, `builder_audit.csv`, `anchor_support.csv`, `target_support.csv`, and `donor_weights_by_event_cell.csv`.
