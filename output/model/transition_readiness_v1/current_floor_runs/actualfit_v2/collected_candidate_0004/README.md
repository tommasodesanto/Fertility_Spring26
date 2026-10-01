Candidate 0004: rejected horizon comparison
==========================================

The permanent preference proposal is 0.14766471988070498. This is not a fitted
shock. Both 12/16-date equilibrium roots converge and replay exactly, but the
unchanged horizon measurement criterion rejects the candidate. No fabricated
`complete.json` was created. All fit, state-experiment and production readiness
flags remain false; no tolerance changed during collection.

The 2020–23 target is 1.64575. Its model moments are 1.885774640307272 for 12
dates and 1.8892769682539958 for 16 dates, with gaps 0.240024640307272 and
0.2435269682539958. Weighted losses are 0.05761182795463525 and
0.05930538426698267. `target_fit_horizon_012.csv` and
`target_fit_horizon_016.csv` contain every target, model moment, gap, weight and
loss contribution; JSON versions preserve original measurement and provenance.
Weights are 0, 0, 0, 1, so earlier windows validate. `shock_parameter.csv`
contains the proposal, absolute bounds and near-bound flag;
`baseline_parameters_all31.csv` preserves every baseline parameter.

Both roots use six mappings. Maximum housing residuals are
1.5837490347934126e-5 and 7.489571387665621e-5; fiscal residuals are
1.4205472898428462e-6 and 7.0234059016671135e-6. All 28 dated accounting audits
pass; replay, policy reproduction and projection gaps are zero. These results
do not override the rejected horizon comparison: last-window fertility differs
by 0.003502327946723671 across horizons, both terminal checks fail, and the
physical 2023 distribution comparison fails. Original comparison receipts are
retained in `candidate_0004/` and embedded in `verification.json`.

The authentic long-horizon 2023 checkpoint is downloaded at
`candidate_0004/state_2023_checkpoint/actual_2023.pkl.gz` (32,809,849 bytes),
SHA256 `1292ff812f62d2ef75e8c8aee380891e2d5768a3bf594c5f7fc18441819780aa`.
Its receipt SHA256 is
`a3e996e74c0de2d53b1e66e7f9de4901cb2cb6950a99ed2a9d723a9539db1de1`.
The state was not opened, reconstructed, rescaled or approved for experiments.

All 71 collected files match their remote SHA256 inventory. The pinned executed
V2 plan, archived deterministic `fit_rows` source and empirical input hashes
were verified before table extraction. Collection made zero model calls and
submitted no jobs; it changed only this output folder. The retained household
fertility measure remains an analogue to female-exposure NCHS TFR, not the same
empirical denominator. Standard graph inspection and 104/128 production
certification remain outstanding.
