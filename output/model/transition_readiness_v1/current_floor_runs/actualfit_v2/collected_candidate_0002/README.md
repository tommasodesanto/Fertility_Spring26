Candidate 0002: diagnostic lower preference, not fitted
======================================================

This completed candidate uses permanent 2007 preference
0.16985485829436683, equal to the initial level times exp(-0.01). It is a
derivative proposal, not a fitted shock. The complete four-window target table
is `target_fit.csv`; `target_fit.json` retains the original source-rich schema.
Targets are 1.974875, 1.861, 1.755375 and 1.64575; model moments are
2.084491400612583, 2.0852768411819356, 2.085378206973357 and
2.0855981163693076. Weights are 0, 0, 0, 1, giving weighted loss
0.19346636547362792. Earlier windows are validation. The retained household-rate
model analogue differs from published female-exposure TFR.

`shock_parameter.csv` records the proposal, absolute search bounds and
near-bound flag; `baseline_parameters_all31.csv` preserves all baseline
parameter estimates, bounds, status and original near-bound flags. The exact
executed V2 plan was authenticated against SHA256
`6b6b28afbf798562f2f28fd8c38eb33d844a36363bc53c45d4b35b87aeb34945`.

The 12/16-date path roots converge in four/five mappings. Housing residuals are
8.776309779227423e-6 and 4.312083167398036e-6; fiscal residuals are
1.14143317433806e-5 and 2.9547684132548232e-6. These pass the authorized dated
fiscal tolerance 2e-5. Fresh root replay gaps are zero, all 28 dated accounting
audits pass, and policy/projection errors are zero. Macro/fertility horizon
comparison passes, but **both terminal checks fail**. The physical 2023 state
comparison also fails: normalized distribution L1 is 0.0015694665555649686,
above its unchanged 0.001 tolerance. `state_2023_horizon_comparison.json` is
preserved in the collected candidate directory; value integrity remains true.

The actual 2023 checkpoint was downloaded without opening or changing it:
`candidate_0002/state_2023_checkpoint/actual_2023.pkl.gz`, 32,869,383 bytes,
SHA256 `9794f7ec276208d6564879876235febb10060edda31ec59ca1100f2a436d4dab`.
Checkpoint receipt SHA256 is
`60ca1a724af64883b6ea6b6f313a8fbd6b556a6fe2f7f1bf1bd71d0eaabf4b82`.
It is the actual period-4 state with both native queues, forecasts and saved
continuation; no rescaling or reconstruction occurred. Its original receipt
keeps `shock_fit_complete=false`, `state_experiment_ready=false`, and
`full_path_certified=false`. No fitted experiment-ready state is claimed.

All 60 compact collected files match remote SHA256 pins in
`remote_inventory.json`. `verification.json` records the result gates, input and
table pins, state comparison and outstanding limitations. `derive.py` executes
only the archived deterministic `fit_rows` AST after verifying its source hash;
it does not import model code. No model call, plot, source edit or job submission
occurred during collection. Standard graph review and full 104/128 production
certification remain outstanding, and the ongoing fit remains the lead's work.
