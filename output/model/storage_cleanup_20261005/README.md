# Fertility storage cleanup — October 5, 2026

The author-authorized cleanup removed **1,694 nonbest generated checkpoint files**, accounting for **309,594,630,656 allocated bytes (309.595 GB)**, from four ended historical Torch experiments.

| Historical experiment | Files removed | Allocated GB removed | Protected checkpoint entries |
|---|---:|---:|---:|
| September 25 nightpair | 680 | 110.037 | 40 |
| September 23 utility search | 379 | 64.213 | 56 |
| September 28 gated search | 520 | 111.533 | 23 |
| September 7 search | 115 | 23.811 | 38 |

Worker/stage bests (including ties), latest completed checkpoints, selected/reference/final cases, pinned dependencies, repeats/probes, and ambiguous cases remain. Each trial's compact parameters, loss, moments, status and provenance remain, along with plots, logs, source, inputs and runtime. **CSI and unrelated projects were not touched.** Current transition and calibration evidence, the frozen September 14 reference and immutable packages remain protected. No model or pipeline changed; no model job was launched.

Every successful unlink matches the exact manifest. Regular single-link files passed metadata and path guards; every protected checkpoint metadata record matched before and after. Nightpair's measured directory-size decrease equals its removal total and a scratch write/fsync passed. Later narrow directory-size commands exceeded their 60-second bounds; those removal totals come from successful per-file unlink receipts, without a separate directory-size claim. This is metadata verification, not a new scientific or large-array content audit.

The last account quota report displayed **4.83 TB / 5.00 TB (96.52%)** before the final September 7 deletion. It does not reconcile exactly to all four file-level passes; no inferred current account total is asserted.

[The concise cumulative receipt](cleanup_summary.json) pins each detailed completion receipt and the quota record. Detailed manifests, guard scripts and execution logs remain local under this folder. **All four passes have executed once; do not rerun their execution commands.**

Future retention proposal, not implemented: keep compact parameters/loss/moments/status/provenance for every guess; retain full current and previous best per independent chain plus selected final/reference solutions and lightweight optimizer restart state. Publish the new best atomically before retiring the old version. Broader data engineering and pipeline work remain a separate discussion.
