Normalized restart: baseline candidate preserved, shock unfitted
===============================================================

Candidate 1 holds the new normalized baseline preference at
0.17198899419542374. Its four model fertility moments remain
2.099999999935852; this is not a fitted shock. `target_fit.csv` reports all four
targets, model moments, gaps, weights and contributions. The 2020–23 target is
1.64575 and its gap is 0.4542499999358518. Earlier-window weights are zero and
the final window weight is one. `target_fit.json` preserves the original
source-rich measurement schema. `shock_parameter.csv` gives the proposal,
absolute bounds and near-bound flag.

The restart uses the newly constructed baseline with fixed housing supply
normalization and initial population N0=1. The fresh reference reconstruction
receipt and actual new 14-row target and 31-row parameter reports are retained
under `native_reference/`; `baseline_parameters_all31.csv` copies the complete
new parameter report. No old baseline certificate was reused. The original
candidate and comparison receipts preserve the strict 12/16-date checks, but
`shock_fit_complete=false`, `state_experiment_ready=false`, and
`full_path_certified=false` remain explicit. The deadline ended the run before
a fitted shock was obtained.

The actual 2023 checkpoint is preserved under
`candidate_0001/state_2023_checkpoint/actual_2023.pkl.gz` (32,878,977 bytes),
SHA256 `6e6ecf9f6375c7513c2a4ea7a6a80cda8883f1891239cb3c9b8267e4f8e7322d`.
Checkpoint receipt SHA256 is
`90618da0bc02349e2c5ce3fbddd357acd2ba7defd8fa083f5711a4941d1d3010`.
The authentic export was copied without opening, rebuilding or rescaling it.

Existing 17 standard plots belong to the freshly reconstructed reference report,
not an invented dated transition report. No 2011 policy graph was produced or
claimed during this collection. Graphs still require the lead's visual review.
`remote_inventory.json` records copied hashes and `verification.json` contains
the audited tables, map gates, comparisons and preserved readiness limits.
Collection performed zero model solves and submitted no jobs; only this output
folder was written, with previous results retained.

Collection stopped at a report identity mismatch: the fresh reconstruction's
31-row parameter CSV has SHA256 `ac6f961dd939026e785021ebd3f88b850b7f5a6fe8053a38e8e86de35685d94c`,
whereas the canonical new normalized report has SHA256
`c94aed6f3a2acbc12d79b6101022d3a99f36d1257b8ca80c6b704e22df14036d`.
A narrow comparison found all 31 parameter estimates identical, but the full
CSV identity difference remains for the lead to review. The authenticated
canonical table is `baseline_parameters_all31.csv`; the original fresh table
is preserved separately. `verification.json` explicitly reports the mismatch,
and `derive.py` stops on it instead of silently certifying the report.
