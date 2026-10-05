# Fertility storage cleanup — October 5, 2026

The author-authorized cleanup removed **2,054 generated checkpoint files**, accounting for **356,459,412,480 allocated bytes (356.459 GB)**, from five ended historical Torch experiment families in six passes. The latest round removed **360 files / 46.865 GB**.

| Historical experiment | Total files removed | Allocated GB removed | Current protected checkpoint entries |
|---|---:|---:|---:|
| September 25 nightpair | 680 | 110.037 | 40 |
| September 23 utility search | 430 | 72.795 | 5 |
| September 28 gated search | 537 | 115.179 | 6 |
| September 7 search | 152 | 31.515 | 1 |
| September 11–12 recent parent cache | 255 | 26.933 | 10 |

The second round tightened retention for completed searches: keep selected full solutions for each distinct contract, explicit pinned inputs, and ambiguous dependencies. Utility, gated and September 7 checkpoint entries fell from **117 to 12**, with retained allocation falling from **21.655 GB to 1.723 GB**. These 12 entries comprise nine physical states and three selected-export symlinks. Nightpair's original 40 protected states remain; it was not revisited in the tighter round.

Each trial's compact parameters, loss, moments, status and provenance remain, along with plots, logs, source, inputs and runtime. **CSI and unrelated projects were not touched.** Current transition and calibration evidence, the frozen September 14 reference and immutable packages remain protected. No model or pipeline changed; no model job was launched.

Every successful unlink matches its exact manifest. Candidates were regular single-link files with unchanged metadata and guarded paths. Protected checkpoint identities matched before and after; the recent-parent postcheck confirms all 255 candidates absent and all 10 selected states intact. Nightpair's measured directory-size decrease equals its removal total. Later reclaimed-byte totals come from successful per-file unlink receipts, without a separate directory-size claim. This is metadata verification, not a scientific or large-array content audit.

**Remaining limits:** 262 multiply-linked nonselected states in the recent-parent cache were left untouched. Their named-path allocation may count shared blocks repeatedly and is not a reclaimable-space estimate. The tighter round did not prune plots. This receipt does not claim that every historical output folder is now lean.

The last account quota report displayed **4.83 TB / 5.00 TB (96.52%)** before the final first-round deletion. It predates the latest round; no inferred current account total is asserted.

[The concise cumulative receipt](cleanup_summary.json) pins each detailed completion receipt and the old quota record. Detailed manifests, guard scripts and execution logs remain local under this folder. **All six passes have executed once; do not rerun their execution commands.**

Future retention proposal, not implemented: keep compact parameters/loss/moments/status/provenance for every guess; save full solutions when loss improves, plus selected final/reference solutions and lightweight optimizer restart state. Keep a previous valid best until the replacement is published atomically. This reduces serialization and disk I/O; it does not reduce the equilibrium solve itself. Broader data engineering and pipeline work remain a separate discussion.
