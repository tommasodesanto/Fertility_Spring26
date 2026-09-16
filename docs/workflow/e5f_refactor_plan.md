# E5F behavior-preserving refactor plan

## Stage 1 — freeze baseline inventory (inventory done; freeze pending)

Use the last certified recovery runtime as the reference. Its manifest, source hashes, local availability, dirty-tree status, and unresolved scaled-step discrepancy are recorded in `output/model/e5f_refactor_baseline_20260916/baseline/manifest.json` and `REPORT.md`. Stage 1 remains outstanding until the exact `ffcb…de3c4` helper and required scratch inputs are restored or explicitly re-certified.

## Stage 2 — comparison tests (focused self-checks pass; native pending)

Use the no-solve harness around separately pinned reference and candidate sources. The reviewed dependency-light focused harness passes 12 tests, including distinct source pins, refreshed-hash mutations, input/source mismatch, missing/escape paths, self-check rejection, object-array rejection, and malformed receipts. Preserve function signatures, array ordering, calendar, tolerances, sentinels, projection, final certification, and serialization. Native cross-version replay remains pending restoration of missing inputs.

## Stage 3 — measured costs (Kimi review complete; native benchmarks deferred)

Measure matched local workloads when available; defer native benchmarks until the exact pinned runtime and inputs are restored. Kimi’s bounded profiling review completed with provider-reported completed-step cost 0.415824; this is an analysis receipt, not a native benchmark. Do not launch a long run from the dirty working tree.

## Stage 4 — one bounded Kimi implementation pilot (review gate)

Run one bounded Kimi implementation pilot on the measured hotspot in isolated files after review; stop on access, context, or time failure.

## Stage 5 — selected structural/performance patches (review gate)

Patch only selected hotspots after review. Separately pin reference and candidate source hashes while requiring the same scientific and input contract. Target acceptance is at least 10% repeated runtime reduction or 20% memory reduction on a matched workload, with no material numerical or behavioral regression.

## Stage 6 — integrated verification (author/review gate)

Promote only a reviewed diff with reproducible receipts, comparison tests, and the required Torch gate. Missing inputs or any hash mismatch fails closed and remains outstanding.

Current invariant/test contract: source identities may differ between reference and candidate, but input fingerprints and scientific metadata must match; outputs, ordering, tolerances, sentinels, gates, and serialized receipts must remain unchanged. Local checks validate harness behavior and saved artifacts. Native equivalence and performance claims require Torch execution against the certified snapshot. Today’s no-model-run result reflects missing baseline inputs and is an execution restriction for this package, not a timeless Stage 3 rule.

Selective commit plan (after lead release): include `code/model/tools/verify_e5f_refactor.py`, its focused test, `code/model/README.md`, this plan, and compact evidence reports/JSON only. Exclude raw Kimi streams, the full dirty-tree status log, and large/generated assets; output artifacts remain ignored.

Future-worker rule: use small explicit file bundles with bounded reading and tool calls; if context expands, stop and synthesize existing evidence rather than repeating broad lookups.
