# Independent overnight verification, 4–5 September 2026

The author authorized up to 12 hours of further review and exploratory quantification on 4 September. Work began at 00:48 UTC on 5 September; the ceiling is 12:48 UTC. This is a ceiling, not a requirement to consume the full allowance. The production selection remains the repaired `task_010` until a separate author decision.

## Scope and sequence

1. Retrieve the completed 23-case, eleven-parameter young-ownership coordinate panel. Reconstruct the original twelve-row and experimental thirteen-row objectives separately. Trace the renewal-identity rejection before changing any validator.
2. Reconstruct the selected historical path and 2023 solution on Torch. Save the standard diagnostic packet and audit occupied-state budget floors and global versus local saving choices. Do not substitute a stationary solution for the dated state.
3. Check ownership measurement definitions and the transition from Census-raked historical ages to endogenous household entry. Keep alternative measurements explicitly diagnostic.
4. Conditional on the preceding checks, evaluate a bounded set of diagnostic calibration proposals or quantification exercises. Preserve the model, original target system, numerical gates, and production files. New target weights or closure assumptions must be labeled experimental.
5. Produce a short morning PDF with supporting complete fit/parameter tables, numerical receipts, and a list of decisions for the author. Update the longer audit where new evidence changes its conclusions.

## Budgets and stop rules

- Read-only measurement worker: one `explorer_fast` worker, full context, ten-minute limit; return the available receipts at the limit, without model runs or edits.
- Panel analysis: no solves. Stop proposals if provenance differs, the scientific contract fails, or effective rank/side consistency does not justify the proposed direction.
- Dated-state reconstruction and numerical audit: one exact-loop smoke on Torch before expanding; estimate later stages using its observed time. Every case writes a checkpoint and latest/best-so-far summaries. Maximum initial stage: two hours, including queue time, before reassessment.
- Exploratory quantification: at most two bounded rounds and at most 24 new full evaluations in total. At the observed 7–11 minutes per full calibration case this is 3–4.5 serial CPU-hours, plus queueing and reconstruction; parallel scheduling can reduce elapsed time. Every submitted job has an explicit wall limit and isolated output directory.
- Investigate any active computation with no heartbeat/checkpoint for 30 minutes. Stop on accounting, feasibility, mass, target-fingerprint, or source mismatch. Never convert a failed run into a pass by weakening a gate.
- Reserve the final hour for collection, source verification, PDF rendering, and backup. Finish earlier if the remaining decisions require author judgment or additional computation adds little information.

## Provenance and files

Canonical baseline: `CALIBRATION_STATUS.md`; production candidate `output/model/e5f_transition_calibration_eta063_kappa005_rooms_repair_gauss_newton_coord_20260904a/task_010`.

Inherited diagnostic panel: `e5f_transition_calibration_eta063_kappa005_rooms_repair_youngown_overid_coord_20260904a`, job 16963855. Collector 16964275 completed all 23 cases; ridge-plan job 16964564 rejected a renewal identity. Original artifacts remain immutable.

Torch snapshot: `/scratch/td2248/projects/Fertility_Spring26_pf_20260823`; account `torch_pr_570_general`. Use existing launch conventions and one numerical thread per allocated CPU.

New results: `output/model/e5f_overnight_independent_verification_20260905a/`. Supporting temporary prompts: `tmp/e5f_overnight_20260905a/`. The author-owned `latex/JMP_DS_draft/` subtree is read-only. Existing unrelated working-tree changes are outside this task.

This file records the authorization and initial design. Completed work, job IDs, estimates, and any departures will be recorded in the output folder's progress record and morning report.

## Designs after the first receipts

The full inherited panel has relative numerical rank ten at the declared
\(10^{-3}\) threshold in both target systems. No ridge proposal or calibration
search is justified by that diagnostic. The two bounded numerical questions
below replace further optimization.

1. **Saving maximization.** The 50,000-draw oracle identifies rare non-global
   choices. A separate `reviewer_strong` worker checks the analytic oracle and
   full kernels, read-only, with the wrapper's twenty-minute limit. The lead
   checks its conclusions and reruns the independent smooth-optimizer tests.
   The exact-loop smoke comprises two complete backward solutions at a common
   price and one globally optimized market-clearing loop. The subsequent
   comparison, conditional on smoke passage, comprises six separately cleared
   2023 cases (two numerical methods times baseline, supply +20%, and family
   credit) and three global-method evaluations at the corresponding local-method
   prices. Each case records reported-budget excess and retains the standard
   visual diagnostic set. The production reported-consumption floors are
   inherited and measured; their nonzero excess is not declared feasible.
   Job `16981492` was cancelled while still pending, with zero solve time,
   to add these receipts and the market-loop smoke. Revised smoke job:
   `16983266`, thirty-minute wall limit (failed and superseded as recorded
   below). Intermediate pending job `16983006`
   was also cancelled with zero compute time to correct the dated diagnostic
   unit conversion and add compressed case checkpoints after visual review.
   Full-stage runtime and submission must
   be based on its observed cost.

2. **Finite-step stability.** Seven complete historical evaluations use the
   unchanged thirteen-row diagnostic profile: anchor, annual patience plus/minus,
   owner premium plus/minus, and bequest shifter plus/minus. The step is 0.005
   in the original transformed coordinate, one quarter of the inherited panel's
   step. The original 23-slot generator is used only for task IDs
   `1,2,3,8,9,14,15`; this is explicitly a seven-case subset, not a full panel or
   an eleven-column Jacobian. Job `16982612`, task 2, is the exact-loop smoke.
   Only after all gates and artifacts pass may the other six tasks be submitted,
   with at most six concurrent cases and twenty-five minutes per case. The
   observed seven-to-eleven-minute historical evaluation cost implies
   49–77 serial CPU-minutes. Reassess after two hours including queue time.

Both designs pin their inputs and commands in the output directory's JSON
submission receipts. The source snapshot for the saving experiment is
`/scratch/td2248/projects/Fertility_Spring26_independent_audit_20260905`,
constructed from the production-matching Git revision
`ac676c2a2c6b27f92319625f77affde2fd6b4b74` (`4fa4c2e^`); the small-step
panel uses the unchanged diagnostic source in the original Torch snapshot.
Together the planned expansions remain below the twenty-four-evaluation budget.

The final seventeen-graph packet for the reconstructed production state is
`output/model/e5f_overnight_independent_verification_20260905a/diagnostics_dated_units/`.
It preserves the established graph set and layouts, while dividing housing
quantities by actual household mass and correcting the old per-adult label.
The output folder's README provides the one-command `--graphs-only`
regeneration from the hash-pinned checkpoint. The source file is
`code/model/tools/run_e5f_independent_numerical_audit.py` and must be run from
the frozen production snapshot. It rejects a mismatched checkpoint or live
code bundle before any diagnostic computation.

## Completed execution and disposition

Both bounded questions are complete. No further numerical panel or calibration
search was launched. The two PDF sources incorporate the completed results.

| Stage | Torch job | Outcome and observed cost |
|---|---|---|
| Dated reconstruction and exact-loop smoke | 16981272 | Complete, 6:43; all 75 historical entries exactly reproduced |
| 50,000-draw interval oracle | 16981404 | Complete, approximately 49 seconds from the checkpoint |
| First checkpoint-based global smoke | 16983266 | Failed exact local reproduction after two backward solutions; approximately 52 seconds |
| Corrected global smoke | 16984897 | Complete, approximately 54 seconds in the driver; exact local reproduction and global dominance pass |
| Seven-case smaller-step smoke | 16982612, task 2 | Complete, 9:47; every scientific and numerical gate and artifact check passes |
| Remaining six smaller-step cases | 16985059 | All complete, 9:59–10:51 per task, within 25-minute limits |
| Nine-case impact saving comparison | 16985060 | Complete, 14:06, within the 60-minute limit; all six cleared cases meet requested residual `1e-5` |

The first global smoke failed because the new audit adapter reloaded the saved
state without reinstalling the three sequential calendar operators. This was
an audit-driver error, not a production-kernel failure. A zero-solve replay with
the correct operators reproduces all saved distributions and aggregate
quantities exactly. The corrected driver installs and asserts those operators
before every case. The successful reconstruction and saved-array audits used
the correct production initialization and remain valid. The failed output is
retained under `global_saving_smoke_v3/`; the corrected result is under
`global_saving_smoke_v4/`.

One bounded local fallback smoke was attempted during the priority delay, in an
isolated temporary environment matching Torch NumPy/Numba versions. It completed
one backward solution and was stopped when the Torch smoke exposed the adapter
error. It was not retried locally and was not used as numerical evidence. The
environment, elapsed time, and termination receipt are retained. The two pending
cancelled versions consumed no compute. Counting reconstruction, smoke attempts,
the stopped local evaluation, and both completed expansions conservatively gives
23 scenario evaluations, within the cap of 24. Price-loop backward solves are
counted separately in each case receipt and are not hidden in that scenario count.

The saving comparison uses a common inherited 2023 distribution and tighter
current market precision, with no history re-estimation. Global saving changes
the supply birth response from `1.290982%` to `1.290873%`, and the family-credit
response from `0.017845%` to `0.017850%`. The smaller-step panel confirms the
owner-premium versus first-birth-rooms trade-off and material derivative step
dependence. Its seven tasks are a diagnostic subset, not a complete new Jacobian.
No incumbent, target, weight, external restriction, production code, or gate was
changed. Remaining work requires a stated measurement/identification strategy
or a separately designed discretization comparison, so another panel would not
be an appropriate use of the remaining twelve-hour allowance.

The authoritative completed readouts are `global_saving_receipts/summary.json`,
`smallstep_receipts/summary.json`, and `panel_receipts/summary.json` in the
overnight output directory. They pin inputs, complete fit/parameter tables,
operator identities, nine compressed case checkpoints, and standard graph
counts. The output README gives the regeneration commands; the morning brief
and full audit explain the economic conclusions and their limits.

## Authorized continuation after the first PDF delivery

At approximately 04:00 UTC the author explicitly invited further safe work while
he sleeps. The completed 23-evaluation stage above remains closed. The following
separate discretization experiment has a new, bounded allowance of **seven
additional scenario evaluations**, within the original 12:48 UTC overall ceiling.
It does not authorize a calibration search or production changes.

The lead will embed the exact inherited 2023 pre-choice measure on nested liquid
wealth grids. Subdividing each original interval into two pieces gives 239 nodes;
four pieces gives 477. Every old node and its mass remain exactly in place, all
new nodes initially have zero mass, and the endpoints remain -12 and 30. This
changes backward-solution and current transaction interpolation resolution,
without reconstructing history or changing the economic housing choice set.
The independently checked global saving kernels are used throughout, so their
known local-optimizer difference is not confounded with grid resolution.

1. **Exact-loop smoke:** one 120-node global baseline evaluation at its previously
   cleared price, reproducing the completed comparison; one 239-node baseline
   through the full existing market loop. Both write a checkpoint, all seventeen
   standard graphs, accounting and probability checks, and numerical exposures.
2. **First comparison:** only after the smoke passes, two 239-node cleared cases:
   housing supply +20% and dependent-child 95% LTV. Baseline comes from the smoke.
3. **Optional second comparison:** at most three 477-node cases (the same baseline
   and policies), only if the first comparison materially changes the interpretation
   or leaves appreciable resolution dependence. Its own first baseline is the
   size-specific exact-loop smoke; the policy cases stop if it fails. This is not
   permission to keep doubling grids until a preferred answer appears.

Before the 239-node policy results are available, the second-round trigger is
specified as any of: a change above 2% of the coarse supply birth effect; above
5% of the coarse credit birth effect; above 1% in baseline births per household;
or above 0.1 percentage point in baseline ownership. These are disclosed
diagnostic materiality thresholds, not sampling tests or acceptance standards.
If none is triggered and all first-round gates pass, stop at 239 nodes.

Observed 120-node global backward solves take roughly 11–14 seconds. Doubling
states and interpolation intervals suggests roughly four times that cost at 239
nodes: 45–60 seconds per backward solution, or about 10–20 minutes per market
case at 14–20 iterations. The three fine cases are provisionally 30–60 serial
CPU-minutes; a possible 477-node round is provisionally two to four CPU-hours.
Reassess these estimates after the smoke. Allocate one CPU and 16 GB for the
first round, with a 90-minute smoke limit and 90-minute policy-stage limit;
the optional round cannot exceed four hours or encroach on the final report hour.

Stop on any frozen-source, helper, inherited-checkpoint or smoke-receipt hash
mismatch; nonfinite/negative mass; changed initial measure or economic parameters;
failed population/feasibility gate; or failure to meet the requested market
residual of 1e-5. The production acceptance thresholds remain documented separately.
Retain failed-case evidence. Write runtime heartbeat every minute and completed
case summaries/checkpoints immediately; investigate active silence at 30 minutes.
There is no best calibration candidate: `best_so_far.json` explicitly records
the latest admissible resolution comparison and no production selection.

One independent read-only `explorer_fast` worker has a ten-minute limit to trace
the ownership age conventions. The lead will separately quantify the occupied
mass exposed to the existing saving/first-birth irregularities from saved arrays.
Neither read-only task requires a model solve. Useful findings will be added to
both PDFs with the conditional nature of the grid experiment stated prominently.

The numerical driver is `code/model/tools/run_e5f_nested_wealth_grid_audit.py`;
the isolated launcher is `code/cluster/submit_e5f_nested_wealth_grid_audit.sh`.
Two tests verify preservation of the joint inherited measure, including
nonlinear state payoffs on an unequal grid. The independent read-only collector
is `code/model/tools/build_e5f_continuation_receipts.py`. It pins the saved
source and state receipts and independently checks first-birth attribution
against the post-fertility childless stock difference.

Smoke `16986180` completed in 16 minutes with a 15.3-minute fine-grid baseline,
15 backward solves, residual `2.289275e-6`, zero feasibility projection, and no
occupied negative wealth-value steps. Both checkpoint hashes and 34 standard
graphs were verified; the seventeen fine-grid graphs were visually inspected.
The remaining two 239-node cases completed as job `16986566` in about 29 minutes.
All gates and standard graph packets pass. None of the predeclared refinement
triggers is crossed; `nested_grid/stop_decision.json` records the calculations.
This branch stops after four scenarios, without a 477-node launch.

## Bounded rental-cap mechanism sensitivity

At 04:53 UTC, after the fine-grid baseline and supply cases passed, the lead
specified one additional model-sensitivity question within the author's broad
overnight authorization: does the six-room rental limit materially govern the
contrast between supply and family credit? Around 40% of dependent-child renters
reach that limit under the supply policy. This is a concrete, occupied margin.

The diagnostic raises the continuous rental upper bound from six to **eight
rooms**, the penultimate existing owner size. Eight is an explicit diagnostic
value, not an empirical estimate or production default. The owner menu and its
option-specific taste distribution do not change. All estimated parameters,
external restrictions, wealth endpoints, fiscal objects, and the inherited 2023
measure remain fixed. The experiment uses the already tested **239-node** grid
and global saving. It changes an economic housing opportunity, so it must never
be described as numerical refinement, re-estimated SMM, or a new production result.

The additional allowance is **three cases**: one cap-eight baseline through the
complete market loop as the exact-loop smoke; then, only if it passes, supply
+20% and dependent-child LTV 95%. Stop on the same source, checkpoint, population,
probability, feasibility, requested market precision, and occupied value-monotonicity
gates. Retain reported-budget excess separately. Each case must save the full
seventeen-graph packet, checkpoint, latest case, and explicit no-selection record.
Heartbeat is once per minute. Use one CPU and 8 GB, 90 minutes per stage. The
observed 239-node cost gives an initial 45-60 serial CPU-minute estimate for all
three cases. No additional cap values, owner-menu variants, calibration, or
historical re-estimation are included. The first-birth rooms target requires its
matched 2019-2023 branch and is not remeasured by this current-impact exercise.

This brings the maximum continuation allowance to ten scenarios (seven grid,
three rental-cap), without extending the original 12:48 UTC ceiling or the
reserved final reporting hour. The cap baseline may run alongside the remaining
grid-credit case in the isolated snapshot; they read common immutable inputs and
write separate output folders. Economic conclusions will be reconciled only
after both comparisons are collected. A failed smoke stops its own branch.

Cap-eight smoke `16987040` completed in 16 minutes. Its checkpoint and source
hashes, sole public parameter override, numerical gates, and all seventeen
standard graphs were verified; the contact sheet was visually inspected.
Policy job `16987657` was submitted at 05:27 UTC for the final two cases,
with a 90-minute limit and an expected 30-40 serial minutes. The smoke receipt
SHA256 is `49c235b071b45b1bb5e552dca0857588f9b4c51cc42a1b41985467375386d7a0`.
No further numerical or economic branches will be launched in this continuation.


## Continuation completed

Both final cap policies completed as job `16987657` in 33 minutes at 06:00 UTC.
All original gates pass; maximum cap-eight market residual is `5.2349e-6`.
The continuation stops after **seven** scenarios: four nested-grid evaluations
and three rental-cap cases. No 477-node resolution or further economic variant
was launched. Every checkpoint hash and standard seventeen-graph packet is
verified. The initial source and helper hashes remain pinned.

The read-only collector verifies exact inherited input distributions and grids.
An additional demand for exact identity of processed pre-choice arrays rejected
the cap-eight supply case: existing feasibility preparation moves `8.3098e-31`
units of mass, yielding an L1 difference of `1.6620e-30`. That rejection remains
recorded. The collector now separately verifies the intended exact input identity
and exactly replays the existing preparation operator from saved policies in all
three cases. No production gate was altered. See
`rental_cap/processed_measure_identity_investigation.json` and
`continuation_receipts/rental_inherited_measure_checks.json` in the output packet.

The same read-only collector independently aggregates the existing ACS/MMS
housing-stock cells: renter seven-plus-room share `9.1661%`, with source hash and
sample differences recorded. This is descriptive saved-cell evidence, not a raw
microdata rebuild, an empirical standard error, or a new calibration target.

The revised PDFs include age-window alignment, zero current exposure of the
visible first-birth spike, the wealth-grid stopping decision, and the complete
rental-cap sensitivity. Production remains unchanged. Remaining decisions concern
measurement, identifying evidence, first-birth housing fit, housing opportunities,
and demographic/fiscal closure, rather than another diagnostic search tonight.
