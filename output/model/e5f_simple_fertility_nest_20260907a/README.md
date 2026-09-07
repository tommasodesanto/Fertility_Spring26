# Simple simultaneous fertility-nest experiment

**September 7, return from commute: both retained-parameter objectives complete.**
Simple fertility nests17142457 completed in38m10s, loss36.37166360862253;
sequential exhaustive-saving control17142458 completed in22m06s,
loss30.408527701170645. New loss is19.6101% higher at the identical11 retained
coordinates. All12 targets and weights, source bundle and plans verified;
all24 fit-row losses independently recomputed. The old-state fertility intercept
is separately normalized to2.1 as maintained. No recalibration or policy run.

Both five-date histories pass recorded market/measurement/mass/population gates;
terminal budget-violating mass and occupied negative value steps are zero.
New first-birth housing response0.421118 versus control0.439708, target0.720246.
The3+- versus1–2-child rooms gap rises0.404567 to0.419762, target0.367700;
this contributes4.00 of the5.96 additional loss. TFR and childlessness barely
change. Interpretation: functional, tractable experiment with worse housing
fit at retained parameters; neither optimal attainable fit nor policy validity
has been established. No exact repeated full history yet. Controls share
exhaustive saving but interpolation support/storage differs. Large checkpoints
remain remote; six collected summary/fit/parameter artifact hashes verified.
All jobs ended; no automatic monitor. Full tables and verification:
`output/model/e5f_simple_fertility_nest_20260907a/RESULTS.md`.

[Complete result tables](RESULTS.md). Earlier launch notes follow.

Source: isolated branch `codex/fertility-nest-computation`, commit `88b21439`.
Specification: [SPECIFICATION.md](SPECIFICATION.md). Production is unchanged.

Cluster verification job 17142456 is running; all23 tests passed there in0.85s.
New objective17142457 and sequential exhaustive-saving control17142458 depend
on successful verification. Each objective evaluates only the retained task010
parameter vector with all12 targets and all11 original estimated coordinates.
No search, policies, illustrations or monitoring automation.

The verification has a30-minute internal limit and35-minute allocation.
Each full-history objective has a2-hour Slurm cap. Plans and submission scripts
here pin source, targets, input hashes, scope, solve forecasts and stop gates.
Latest observed state is recorded in launch_receipt.json, not a live monitor.

Local exact default-off reproduction passed all10 arrays. The first lifecycle
harness omitted sequential calendar operator wiring; its failed receipt is
preserved. Corrected local attempt hit600-second cap during matched-sequential
verification, without a completed lifecycle receipt. Cluster harness now records
individual phases and periodic Python stacks; no numerical runtime setting was
changed. Cluster checks new nested choice first, then sequential control.
No complete objective or new fit is available at submission.

Remote root: `/scratch/td2248/projects/Fertility_Spring26_simple_fertility_20260907a`.
Outputs are under its `output/model/e5f_simple_fertility_nest_20260907a/`.
After collection require cluster_probe/summary.json complete, both historical
validation receipts, full target_fit_long.csv and parameter_table.csv before
interpreting fit. Any failed dependency blocks both objective jobs.


**Commute launch update:** verification17142456 COMPLETED successfully. All23
cluster tests and both full-grid fixed-price lifecycle/accounting/budget/value
checks pass. New solve28.14s, sequential exhaustive21.24s; current-population
L1 errors1.55e-14 and1.85e-14, birth gaps below1.1e-15. Neither fixed price
clears the housing market; these are code checks, not calibration fits. Both
objective dependencies released, waiting for scheduling at last check. At the
observed cluster timing,50–100price evaluations are roughly24–47minutes for
the new model, plus history/audit overhead; the2-hour caps remain unchanged.
Local slowness was not reproduced on the cluster; its specific cause remains
unresolved. Collected receipts: same output folder, `cluster_probe/`.

