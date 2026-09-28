# Codex worker task

## Goal
Read-only compatibility review: identify required changes for the existing dated transition runtime to consume the exact September 28 block0506 saved checkpoint without preference or measurement changes.

## Scope
Load mandatory project context in order with bounded text reads. Then code/model/tools/e5f_current_transition_runtime.py, run_e5f_current_transition.py, run_e5f_current_transition_smoke.py and their direct wrappers/helpers; output/model/daytime_calibration_20260927/credit_transition/README.md; fixed_reference_manifest.json and measurement_audit_v1/README.md under output/model/fertility_identification_20260928/. Focus on checkpoint loading, frozen module overlays, dated household/KFE maps, earnings/pension/estate/entry/housing closures, stationary vs transition validation, source pinning.

## Context
Full startup required. Label: 2007 stationary reference — block0506, September 28 verified export. Never replace saved 260 parameters with constructor defaults. Current checkpoint is on Torch, no local full checkpoint. Author explicitly chose fixed physical housing stock/supply with prices/rents clearing and removal of artificial borrowing/down-payment limits retaining lifetime solvency and repayment. Lead owns economic definition and model-critical implementation; you locate compatibility evidence and hazards.

## Do not touch
Read-only. No edits, model imports, hashes, tests, solves, rendering, SSH, jobs, downloads, git actions or additional workers. Do not inspect conversation archives or overlap the historical specification review.

## Required output
At most 1800 words with concrete file/line evidence: reusable code; required adapter hooks; defaults that would silently alter economics; existing tests/gates; isolated overlay design and smallest zero-shock smoke. Distinguish no constraint from phi=1 and grid clipping. Lead verifies all model-critical evidence.

## Verification
Trace actual call sites and parameter mutations, not just docstrings. No computation.

## Stop and report if
Hard 20-minute limit within the lead's 45-minute preparation budget. Return best evidence and gaps; no retries, broader audit or implementation.
