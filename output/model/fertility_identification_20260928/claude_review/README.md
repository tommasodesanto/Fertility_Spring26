# claude_review — independent review of the September 26–28 calibration

Read-only review written September 28, 2026 (16:10–16:45 EDT) against the snapshot recorded in `review_report.md` section 0. Nothing outside this folder was changed; no model code was run.

- `review_report.md` — the report: executive assessment, answers to the eight review questions, three recommended checks, missing evidence, direct answer.
- `appendix_tables.md` — full contract, target-fit (block0506 and lane winners), 31-parameter, lane-snapshot, resume-manifest, early-frontier, failure-diagnosis and lifecycle tables; the September 26–28 decision ledger with status-note line references; the code line references behind the specification claims.
- `linear_analysis_outputs.md` — every Jacobian-derived number in the report (pivots, step sensitivity, implied weights, constrained directions, Gauss–Newton benchmark, linear trade-off frontier, NCHS cell shares, the time-aggregation bound, near-bound arithmetic).
- `linear_analysis.py` — stdlib-only script that regenerates `linear_analysis_outputs.md` from the saved Jacobian, SVD, block0506 rescore table and NCHS counts file.

A PDF rendering of the report should be produced on Torch, not on the Mac.
