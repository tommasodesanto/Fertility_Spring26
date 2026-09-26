# Saved floor-linear inspection — stopped duplicate diagnostic

The coordinating lead explicitly stopped this diagnostic after completing and
verifying the requested calculation centrally (reported main commit
`91965de9`, Torch job `18602589`). This folder preserves only the inspection
already completed here. It contains no reconstructed incentive-gap results,
population-weighted quantiles, new plots or new model solves.

Our Torch inspection job `18603092` completed with exit `0:0` in 15 seconds.
It authenticated the original comparison contract and selected checkpoint,
loaded the frozen source and saved objects, and recorded their configuration
and schemas. The checkpoint hash matched the original receipt exactly.
The saved policy has sequential births and independent child counts enabled;
joint nesting, fertility nesting and readiness gating are disabled. Matching
saved `evaluation.g_pre` and policy probabilities are present. These findings
agree with the central diagnostic's stated extraction method; no material
discrepancy was found in this inspection.

`inspection.json` and `inspection_status.json` are the saved-object and job
receipts; `inspection_submission.json` records what ran. `inspect_saved.py`
is the Torch-only inspector. The two source-evidence JSON files retain the
authenticated Bellman mapping and the distribution timing excerpts. An
independent read-only source review confirmed that `evaluation.g_pre` is the
pre-attempt distribution, while `solution.g_beginning_distribution` already
includes current-period births. No new forward pass was necessary or run.

The original `../supplemental_birth_wait/availability.json` remains unchanged:
its earlier stored-action-values-only authorization differs from the later
explicit permission to reconstruct gaps from probabilities. There is no
competing diagnostic or soft-floor experiment in this folder.

Central result location, as supplied by the coordinating lead:
`output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/utility_fertility_rationale_review/saved_floor_linear_20260926/`.
The central lead owns its mathematical review, interpretation and subsequent
authorized experiment.
