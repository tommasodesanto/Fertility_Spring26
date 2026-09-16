# Stage 3 profiling assessment

Status: **BLOCKED for native measurement; Kimi completed a bounded final synthesis after the initial exploratory run was stopped.**

The local baseline inventory at `output/model/e5f_refactor_baseline_20260916/baseline/git_status_full.txt` contains only the captured status list. The recovery manifest and receipts point to `/scratch/td2248/...`, which is not available in this workspace. The exact frozen source/input bundle needed to run the current `original_queue` / long successive refit path is therefore unavailable. No solve, synthetic benchmark, or wrong-default benchmark was run.

Kimi K3 was dispatched through OpenCode Go with read-only permissions. Repository access worked: it read the requested source, README, manifests, receipts, and status files. The initial run was stopped after context/tool expansion, then continued once with all tools denied and a final-only prompt. The final synthesis is in `kimi_final.jsonl`; the exploratory stream is in `kimi_raw.jsonl`.

## Evidence and recommendations

Historical receipts record 104-date mapping evaluations of 2,568.8–3,412.2 seconds (about 25–33 seconds per date), 24-evaluation stationary roots around 659–691 seconds, and warm continuation steps around 19 seconds. These are historical original-queue/recovery values, not current local measurements; old September 5 timings are a different source. The average seconds/date does not establish a Bellman hotspot. The active exact policy cache is additive and exact-call keyed (`code/model/tools/e5f_exact_policy_cache.py`), but existing receipts contain no cache counters, phase splits, RSS, or thread counts.

Ranked next steps, pending native inputs:

1. Add an observer wrapper around the current driver that records fresh-process versus warmed state (with JIT status explicitly recorded), per-date Bellman/remaining-evaluator/measurement/serialization timings, call counts, peak RSS, thread count, and exact-cache hit/miss/eviction counters. The first mapping is not automatically an all-JIT or cold observation. Subtracting Bellman time from `evaluate_period` gives remaining evaluator time, not a pure forward timing, unless the forward seam is timed directly.
2. Run one representative native mapping in one process and one thread, then repeat warm; compare rows, accounting residuals, feasibility, first-period replay, and serialized hashes under the existing gates.
3. Only after that pilot, target the largest measured phase. Do not optimize based on old September 5 timings or change algorithm, grids, precision, gates, diagnostics, or cache-key semantics.

The pilot must fail closed: the certified staged source hash `ffcb2334…` and native inputs under `/scratch/td2248/...` are absent locally; the local helper is a later variant. The project `.venv` had a startup failure, while the bundled runtime remains available for pure harness checks. Actual flags, errors, and tolerances must be restored from the pinned contract rather than hard-coded from worker recollection. Keep raw samples and hardware/runtime metadata in this directory.

## Lead review note

The raw Kimi final contains stale assertions about baseline availability and interpreter status. Those assertions are preserved as worker evidence in `kimi_final.jsonl` but are not adopted here. This report adopts only the receipt-backed historical timings and the qualified profiling recipe above; no phase hotspot, cache benefit, or numerical tolerance is established until the pinned native contract and inputs are restored.
