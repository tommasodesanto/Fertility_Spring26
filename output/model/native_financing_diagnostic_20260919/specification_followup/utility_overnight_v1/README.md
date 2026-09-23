# Overnight utility comparison — September 22–23

**Delivery update:** The author returned before10AM and explicitly requested the report immediately. The reviewed53-page PDF is released on that instruction; all experiment jobs are terminal. Artifact/source backup: `b6032d94`.

## September 23, 09:40 EDT — terminal results collected and reviewed

Torch access returned at09:22EDT. All four smokes passed;20 of40 production workers completed and20 failed after retaining earlier cases. The reviewed inventory has **444 scored/20 failed/272 unrun objectives** and **2,590 started/2,571 completed/19 unfinished stationary solves**. Nineteen failures were the strict initial housing-equilibrium gate; one was the unchanged exact-zero budget-mass gate. All experiment jobs are terminal. No new run, retry or relaxed gate.

The selected losses under the identical objective are B_floor464.8522, B_shares851.2713, D_floor416.4814 and D_shares1524.0382. B means one persistent earnings process; D adds iid risk. All retain inherited heterogeneous wealth. The main floor/shares cases predict ownership39.49%/51.33% versus64.83%, rooms6.48/7.09 versus5.56, and recent-parent ownership gaps12.50/5.51pp versus16.29pp. Neither is an adopted fit. All four selected entry marginals are preserved without censoring.

Only D_floor completed both selected repetitions. Independent original-selected-versus-both checks agree exactly in price, finite V/g bytes, all finite moments, psi and full tables. Seven matching NaNs are unused fixed-type diagnostic defaults; they caused the raw finalizer JSON serialization to fail. That raw failure remains preserved alongside the separate reviewed comparison. The other six repetitions remain unrun/unavailable. No grid-convergence or harmless-continuation claim.

Working recommendation: retain persistent earnings, heterogeneous entry wealth and the floor as the comparison baseline, while fixing the chosen birth-timing/model-data observer before a further long refit. Do not choose utility solely by weighted loss. All six Opus5.5 phases and empirical checks are complete; target changes remain unadopted. The53-page PDF includes all13 target rows,17/19 parameter rows and17 original diagnostics per selected cell. Delivery is held for10:00EDT.

[Final lead review](final_lead_review.json), [complete selected target table](final_readout/selected_target_fits.csv), [complete selected parameter table](final_readout/selected_parameters.csv), [failure evidence](execution_failure_review.json), [checkpoint comparisons](checkpoint_review/independent_lead_support.json).

Author authorizes **40 production workers**, with fixed current targets while a separate Opus5.5/lead target review proceeds. Deadline: **September23 10:00 America/New_York**. Laptop remains open for local review and reporting; Torch computation runs independently.

## Submitted experiment

- 18 workers: persistent income, parenthood-only housing floor, retained power equivalence scale.
- 18 workers: same income, zero housing floor, housing expenditure share changes with **current dependent children**, retained power equivalence scale.
- 2 workers per utility: persistent plus iid income, supplementary comparison.
- All arms use the same inherited heterogeneous wealth marginal from the passing September22 B/D battery, with disclosed diagnostic rank coupling; **no zero-entry fallback**.
- Main income approximation15 persistent nodes; supplementary7x3=21 nodes. Main15 improves on the unvalidated seven-node approximation; supplementary21 remains exploratory because prior45-state D failed. Wealth160 nodes with upper3000 unchanged. This is not convergence evidence.
- Same deterministic age profile, stationary entrant income, four-year timing, current-income purchase eligibility, transfers, full target/weight fingerprint and household/market gates. No fixed income types. Utilities are experimental alternatives, not adoption.

The floor specification has nine searched structural coordinates including h_P; the share specification has ten, replacing h_P with two bounded share tilts. Fertility utility scale normalizes separately in each evaluation. Equal worker counts do not establish equally optimized fits or identification.

## Submission

Submitted **23:15EDT September22**, immutable bundle `tmp/utility_overnight_20260923_v1`, remote `/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_overnight_20260923_v1`, results under `results/`. Complete pins and job IDs: `manifest.json`, `submission.json`.

| Arm | Smoke | Production array | Two repeats | Final comparison |
|---|---:|---:|---:|---:|
| Persistent + floor |18304038|18304073 (18)|18304074|18304077|
| Persistent + shares |18304078|18304079 (18)|18304080|18304081|
| Persistent+iid + floor |18304082|18304083 (2)|18304085|18304086|
| Persistent+iid + shares |18304087|18304088 (2)|18304089|18304090|

Each cell has its own successful-smoke dependency. Production includes up to18 adaptive proposals per worker, with6 local/6 medium/6 broad primary chains and one local/one medium supplementary chain. Each selected point is repeated twice in parallel after its array; the finalizer compares the original selected checkpoint against BOTH repeats in prices, V, g, all saved moments, normalized psi and full tables. Only checked checkpoint/path provenance and wall time are excluded from numerical equality. Equal random innovations do not imply matched parameters after centers diverge.

Cap:8 smoke+720 production+8 repeat objectives =736 objectives and at most5,888 stationary solves. Time budgets can leave many proposals unrun. Observed objective times imply roughly20 minutes per primary point and28 per supplementary point; no fixed realized count promised. Each worker uses1CPU/16GiB, maximum40 production CPUs. Previous D21 peak memory5.1–5.5GiB supports this allocation. Smoke allocations3h, production10h but absolute stop08:00, repeat allocation1h and absolute stop09:00. No laptop-dependent job chaining.

Validation:10 utility tests,23 frozen accounting tests in separate processes,3 new driver tests including8 local zero-solve template/dynamic preflights,4 remote zero-solve preflights, Bash syntax and source review passed. Combined older tests initially collided through Python module caching; separate intended-source runs pass, with no model-source correction. One staging check ran before transfer finished and correctly returned failure with zero jobs; successful post-transfer check preserved separately.

## Preparation and stop rules

Sol owns the new preparation/runner/submission/test files; Luna's independent read-only formula/runtime/entry audit is `utility_audit.json`. Lead reviews core formula and contract wiring before frozen zero-solve preflights, exact-loop smokes, then dependent production. Existing untracked preference files are read-only; any adapter compatibility patch must be vendored, hashed and disclosed. Targets and gates may not be relaxed. Preserve failed evidence, no automatic retries or duplicate jobs.

Search stops by08:00EDT. Reserve selected-versus-two-exact-repetition verification and collection through09:00, report QA through09:50. Raw native progress/latest/best must be observable; investigate30minutes without progress. Full13 targets, complete parameter tables/actual bounds (17 or explicitly expanded rows),17 unchanged diagnostic plots per selected arm, entry and continuation diagnostics, fingerprint checks and honest attempted/completed/rejected/incomplete counts are required. Selection includes smoke results. A finite search is not a converged calibration.

Latest runtime evidence: inherited persistent15 objective1187 seconds and inherited21-state objective1699 seconds, each six stationary solves. Actual scope/solve-count cap will be written in the frozen launch manifest before submission. No additional policy runs or author draft/mock/slide edits.

## Review and delivery

Target research plan and live Opus receipt: sibling `target_review_v1/overnight/README.md`. Follow-up `overnight-utility-and-target-review` checks at :00/:20/:40 and delivers at10:00EDT. Stay quiet during healthy unchanged execution; notify meaningful failure/action. Commit/push only this task's verified sources and artifacts, preserving unrelated work.

## Collection

Reviewed collector: `code/model/tools/collect_e5f_utility_overnight.py`; six focused checks pass after lead provenance/count corrections. Four-template local empty readout validates pins only, not scored-case collection. Remote collector and `audit_e5f_earnings_entry_checkpoint.py` are staged under sibling `utility_overnight_collection_v1/`. Use the frozen remote manifest with `--manifest`, remote `results/` with `--results-root`, and a **new time-labelled directory** with `--output`; wait for completion before rsync of compact outputs. Never download giant checkpoints. Verify first actual smoke collection before relying on final readout. `collector_lead_review.json` records remaining checks. Source launch commit ca90cc77 is pushed.

**23:28EDT operational interruption:** Torch SSH master closed and fresh authentication is required (standard wrapper confirms permission denied). All four smoke jobs and all dependent production/verification jobs were successfully submitted earlier and run independently. Author was asked asynchronously to refresh `ssh torch`; no credential workaround or job retry. Local Opus/data review continues. First collector version is remote; latest timeout-classification-only patch failed to upload and must be rsynced after login returns. Preserve the live bundle. After refresh, maintain a detached read-only connection so the configured600-second idle control expiry does not strand20-minute follow-ups. A10AM report must explicitly state uncollected results if access remains unavailable.
