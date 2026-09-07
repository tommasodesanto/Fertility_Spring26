# Codex worker task

## Goal

Assess what the saved evidence establishes about three market-nonconvergence rejections in the active experimental simultaneous-choice search. Distinguish established facts from possible numerical/economic explanations. This is a bounded independent diagnosis for the lead's final PDF, not a request to change code or rerun calibration.

## Scope

Read the active immutable experimental source through the local worktree `tmp/e5f_joint_nested_full_20260906a/code/model/`, especially tools/run_dynamic_population_transition.py (clear_scalar_housing_market), tools/run_e5f_open_population_transition.py (retry), and the joint_nested.py choice support. Logs, receipts and outputs were copied to `output/model/e5f_joint_nested_full_20260906a/parallel_search_o/output/model/joint_nested_overnight/search/initial_population/`. Inspect case_023.log, case_024.log, case_027.log and task_023, task_024, task_027; read the sibling rejects_ledger.csv. Cases5 and10 hit the unchanged one-hour cap; summarize their last actual progress only if helpful. Inspect the existing experiment memo docs/model/e5f_joint_nested_experiment.md and readout only for model restrictions. Main production files may differ from the experimental worktree: do not confuse them.

## Context

Full project startup required, using root memory/AGENT_MEMORY.md, latest daily, CALIBRATION_STATUS.md in order. Latest state at11:45UTC supersedes older progress counts: Torch job17106283 remains running on32CPUs352GB; first population ended27 valid new histories+3 market failures(cases23,24,27)+2 timeouts(cases5,10), plus4 imported smoke histories. None improved loss450.7052931460765. No DE/polish fit the budget. The controller has moved to24 final histories (22Jacobian probes+2exactrepeats) concurrently with4 full policy paths. No source or contract changes are authorized for this worker. Numeric residuals reported at failure: case24 about.001286 after an earlier.003226 attempt; case23 .004361; case27 .002696. All exceed unchanged.0002 gate. Selected-source smoke previously passed all4histories and8policy dates. The root lead retains all economics/identification/calibration judgment.

## Do not touch

Read-only worker: no source/document/output edits, model solves, job submissions/cancellations, git mutations, thread messages, or production changes. Never touch protected latex/JMP_DS_draft. No remote tool or SSH use is needed; local copied evidence is enough. Do not assert market nonexistence, a coding bug, or unreachable targets from a failure receipt alone. Do not propose relaxing gates or silently replacing targets. Do not audit the entire codebase.

## Required output

A concise final report (wrapper saves it to the assigned file) with: (1) three-case table giving last completed date, failing period if established, exact residual/iteration evidence and retry behavior; (2) code-line references explaining the solver algorithm and what failure certifies; (3) whether the saved files distinguish a demand discontinuity, failed bracket, max-iteration exhaustion, or another mechanism; (4) smallest future diagnostic that would resolve the uncertainty, without running it; (5) explicit separation of findings from conjectures. The lead will inspect the cited lines and reconcile with source/target tables before any conclusion is accepted. Do not include a broad refactor plan or new calibration recommendation.

## Verification

Read-only logs, JSON/CSV and source checks only. No model run. State what cannot be checked from the saved files. Avoid dumping full files or unfiltered tracebacks.

## Stop and report if

Stop with the best supported report by the wrapper's default30-minute limit. Stop earlier if evidence is insufficient rather than expanding into simulations or asking the sleeping user. No automatic retries or new workers.
