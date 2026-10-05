# Fertility storage cleanup — October 5, 2026

Project-only Torch cleanup has removed **210,731 regular files and one symlink**, accounting for **685,764,225,024 allocated regular-file bytes (685.764 GB)**. Only final-link removals count as reclaimed blocks. No current account/quota total is inferred.

The latest author-approved retirement removed **64 entire old experiment roots and 29 newer obsolete output trees**, releasing **328.566 GB**. All 93 exact targets are absent. This includes the September 6–7 joint-choice deployment family shown in the author's screenshot, June–July revisions, old full-joint/simple-fertility/matched-history runs, and obsolete transition v6/v7/v9 results. Later parent source/code/runtime remain. The separate earlier checkpoint and bested-chain passes released 357.198 GB.

[The legacy catalogue](legacy_catalogue.md) records the experiment families and exact paths. Its size tables are historical inventory attribution, not current presence or reclaim guarantees. [The cumulative receipt](cleanup_summary.json) pins all completion evidence; the latest detailed receipt is local at `legacy_retirement_v1/93_completion_receipt.json`.

**CSI and unrelated projects were untouched.** Current October 4 calibration/search work, current H24/H32 packets/maps/full tracebacks, immutable v5, frozen September 14 reference, selected presentation results and the author manuscript remain protected. Explicit retained input identities matched before and after every retirement family. No model/pipeline/economics/gates changed; no model job or changed-budget run was launched.

Lead verification independently reconciled all 93 exact targets, 153,491 unique regular-file removals, per-stream compressed/raw hashes and final block sums. The pass unlinked 5.484 GB of shared-link allocations that were **not** counted as freed because another link survived. Directory/symlink blocks and tiny temporary execution staging are excluded. Detailed per-file receipts were losslessly compressed to 3.52 MB locally; no new remote archive was created. Every command completed once without failure, cap stop or retry. Never rerun deletion commands.

Remaining mixed boundaries: unlisted candidate-path batches, 43 native-financing directories from the saved layout, three initial-population output boundaries, and empirical room-target results. Required input/source boundaries are retained; current presence/size of these remaining mixed trees was not rescanned. Their exact known paths are recorded locally in `legacy_retirement_v1/known_retained_mixed_boundaries.json`.

The recent Estate calibration writer evaluates each guess and saves selected/repeat solution/shared arrays before comparing calibration losses (`cluster_calibrate.py:158`, improvement check at164; `model/native_price.py:54–60`). Existing case samples show roughly73-MB packets. This explains non-improving guess retention; it is not a measured terabyte total. A future pipeline change could keep compact guess/loss/moments/status/provenance and save full solutions only on improvements, with atomic winner publication and required final/reference states. That code change has not been implemented.

The old quota receipt displayed 4.83 TB/5 TB before later deletions. It does not establish present usage or how much account storage belongs to Fertility; CSI was excluded from this work.
