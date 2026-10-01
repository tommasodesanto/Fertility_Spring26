Candidate 0001: verified baseline, fit still in progress
======================================================

This is the unchanged-preference candidate, not the fitted one-shock transition.
The permanent preference level is 0.17156192800028292; its four-window model
fertility measure is 2.100001297131149 in every window. The complete target table
is `target_fit.csv` (source-rich schema in `target_fit.json`); the weighted loss
is 0.2063442409453312. Only 2020–23 is scored; the preceding windows validate.
`shock_parameter.csv` gives the estimate, absolute search bounds and near-bound
flag. `baseline_parameters_all31.csv` preserves every baseline parameter and
its original bounds/status. The target measurement is the retained household
rate analogue, whereas NCHS published TFR uses female exposure.

Both diagnostic horizons, 12 and 16 dates, converge and reproduce exactly.
Their maximum housing residuals are respectively 3.077128737338753e-7 and
5.92966108775315e-7; maximum fiscal residuals are 4.0475104971193814e-7 and
7.818728362592892e-7. All dated accounting gates pass; policy reproduction and
projection errors are zero. Initial endpoint renewal residual is
6.176814453251467e-7. Both terminal checks pass and the early macro/fertility
comparison is exactly stable. These diagnostic horizons do not establish the
required 104/128 production path certification.

The authentic 2023 checkpoint is downloaded under
`candidate_0001/state_2023_checkpoint/actual_2023.pkl.gz` (32,871,448 bytes), SHA256
`7e1f22a14477d82fab618845b6e9f7aa18b45cf7589f642c61cfbd587b8a1c4d`.
The checkpoint receipt SHA256 is
`4e702d5c96c0754965ca532b857fdda8cbc9d14b5858d3270a686b52f3bfb0b1`.
It contains the actual period-4 state, both native queues and saved forecasts
and continuation, without reconstruction or rescaling. Its receipt explicitly
keeps `state_experiment_ready=false`, `shock_fit_complete=false` and
`full_path_certified=false`. It is a usable serialized baseline artifact, not
an approved fitted state for experiments.

`remote_inventory.json` authenticates all 40 collected files. `verification.json`
records the gates, tables, identity, input pins and outstanding limitations.
`collect.py` downloads only completed candidate receipts and the checkpoint;
`derive.py` checks copied hashes, reproduces the original archived deterministic
`fit_rows` function and builds the readout. Neither imports model code, opens the
native pickle, invokes a model solve, submits a job, or plots diagnostics.
The ongoing fit is owned by the lead; this collection does not certify its
later candidates. Standard plots remain uncollected and visually unreviewed.
