# fable_analysis: soft financing, interest timing, and early fertility

Independent read-only analysis (Fable 5.1, October 3, 2026). No model solve, no
production, calibration, Slurm, Git, or manuscript change.

## Verified facts (source-linked)

| Fact | Value | Source |
|---|---|---|
| Original-timing winner | chain 15, loss 18.4453051941, price 0.718172 | `input_snapshot/original_fit.csv`, `original_summary.json` |
| Alternative-timing winner | chain 13, loss 13.7711314635, price 0.776057 | `input_snapshot/alternative_fit.csv`, `alternative_summary.json` |
| Early fertility, both arms | 0.5330788 / 0.5338047 vs target 0.8095 | reproduced from saved arrays in `analysis/out/native_grid_analysis.json` |
| Age-25 share with any birth | model 0.442, CPS 0.457 | same |
| Age-25 children per mother (capped 3) | model 1.21, CPS 1.77 | same |
| Observer bound conditional on current 18–21 first-birth flow and 22–25 hazard (not a global ceiling) | 0.676 | same |
| Held-coordinate timing swap | ownership 30–55 0.665→0.777, early fertility 0.5276→0.5249 | `../smoke_collection/readout.json` |
| Saved distributions | all three are post-fertility; no pre-fertility array saved | identity check in JSON |
| Young childless renters passing implemented closing screen | Original timing: 95–99% at ages 18/22/26; alternative timing: 93–97% | JSON, Section 3 of memo; `revision1/POST_REVIEW_CORRECTIONS.md` |
| New wealth-target arm (entrant wealth unchanged) | loss 48.1707 under a different contract; wealth row contributes 33.60 | `input_snapshot/new_wealth/` |

## Contents

- `ECONOMIC_MEMO.md`: central claim, evidence, counterevidence, early-fertility decomposition, limitations, falsifiable checks.
- `analysis/native_grid_analysis.py`: reproducible computations (run with the project venv, `MPLBACKEND=Agg`).
- `analysis/out/native_grid_analysis.json`: all numbers for both arms.
- `analysis/out/F1..F4_*.png`: four supplemental figures. The standard 17-plot packets per arm are untouched in their native roots.
- `analysis/out/F5_native_asset_policy_mass.png`: native-node first-birth, conditional renter consumption/saving policies, and the corresponding pre-fertility wealth mass for childless inherited renters at ages 22–25 in income states 4 and 6. `analysis/plot_native_asset_states.py` rebuilds it from the hash-pinned saved arrays without a solve. Mass is conditional on arm and income state; policy markers are shown only where mass exceeds \(10^{-14}\).
- `input_snapshot/`: pinned inputs with SHA-256 in `manifest.json`.
- `revision1/`: preserved pre-review files, same-session revision receipts, independent saved-array hashes, and the audited arithmetic correction note.

Run:

```sh
PROJECT=/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
MPLBACKEND=Agg "$PROJECT/code/model/.venv/bin/python" \
  "$PROJECT/output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/fable_analysis/analysis/native_grid_analysis.py"
```
