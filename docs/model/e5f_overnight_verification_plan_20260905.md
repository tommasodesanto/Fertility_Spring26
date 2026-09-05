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
