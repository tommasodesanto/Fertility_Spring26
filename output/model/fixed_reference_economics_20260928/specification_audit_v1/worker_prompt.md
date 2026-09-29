# Codex worker task

## Goal

Urgent low-cost audit for undisclosed or unresolved economic assumptions like
the inherited renter debt taper. Cover the full live frozen model, not every
obsolete repository file. Return an auditable inventory and short findings.

## Scope

Start from `output/model/fertility_identification_20260928/fixed_reference_manifest.json`,
its contract and authoritative `resume_v1/selected_export/primary/` tables,
then the pinned `intergen_eqscale_seq_optimized` package and production runtime
sources named by that contract. Inspect all economic mechanisms: preferences
and utility clipping; tenure/housing menus and costs; liquid wealth and credit;
earnings and entrant wealth/debt; conception/birth/timing and dependent exits;
survival/estates/recipient/creditor accounting; fiscal/pension/tax rules;
population renewal/entry queues/geographic closure; supply/rents/equilibrium;
measurement and normalization. Separate numerical approximations with economic
effects from economic primitives. Inspect current transition code only as a
separate module: do not treat old/experimental transitions as frozen results.

Compare implementation to actual `latex/september_14_presentation.tex`, current
`latex/JMP_DS_draft/sections/03_model.tex`, current slides indexed in
`latex/README.md`, `docs/model/ACTIVE_DECISION_LEDGER.md`, Sept15 consolidated
review and Sept26 `borrowing_negative_equity_resolution.md`. Advisor Google Doc
read this session does NOT disclose taper: ID1hxESCRA89O028-Kx4LmBjdM19R_CobIkn_GkJdMgbbo.
It says adopted purchase/incumbent rules and repayment at sale/death need
explanation; do not infer approval from generic labels or a Git author name.
Frozen reference checkout's deck is an earlier source, not actual Sep14 PDF.
If evidence of human approval cannot be found narrowly, report unknown.

## Context

Full startup required: AGENTS.md, memory/AGENT_MEMORY.md, latestdaily2026-09-29,
CALIBRATION_STATUS.md, code/model/README.md. Reference remains **2007 stationary
reference — block0506, September 28 verified export**. Reference solver SHA
b637a655a9344b63f4461ee0fa4796c04bd98188477c4e6ace2c48ae0fc8aec1;
parameters66f86697c2c58ca3864305bf13dd2be71a008905b2beb573f1a4ebafabef5464;
kernels639c9a21797dbc9f2a0e9a891f283c115353c2edfcb89c959a7fe9f32b86ca27.
These were freshly verified on Torch by lead. Mark other source identity gaps.

## Do not touch

Read-only. No shared source/document edits, Google writes/messages, code imports,
numeric tests/solves/rendering, checkpoint reads/downloads, Git operations,
recalibration, target changes, automatic cleanup. No subworker fanout. Use
bounded rg/text reads; no whole session archive scans. Twenty-minute hard cap.

## Required output

Your final message, saved by wrapper to `audit_findings.md`, must include:
1. At most one-page prioritized findings table: actual rule, source/line,
   disclosure/decision evidence, consequence, uncertainty, smallest next check.
2. Compact mechanism coverage matrix: examined / partly examined / not examined.
3. Distinguish undocumented, documented-but-unresolved, approved-but-stale, and
   numerical approximation. Do not call all oddities bugs or claim global clean bill.
4. At most three highest-priority new concerns beyond the known renter taper.
Evidence appendix may be a concise rule inventory; no long report or prose dump.

## Verification

Read exact decisive source lines and current values/flags from exported tables
or contract, not defaults alone. Defaults for inactive branches are not baseline
features. Lead independently checks flagged mathematics before reporting.

## Stop and report if

Time limit reached, missing source identity, conflict, economic choice required.
Report coverage gaps and ambiguity without inventing facts or approval. No retry.
