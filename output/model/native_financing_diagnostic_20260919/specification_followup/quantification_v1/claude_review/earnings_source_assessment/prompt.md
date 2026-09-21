# Goal
Provide an independent assessment of the earnings-source decision for Tommaso's fertility/housing paper. The author specifically requests Fable. Read-only research assessment; do not implement or run models.

# Context and scope
Read mandatory startup context in order: memory/AGENT_MEMORY.md (bounded relevant excerpts are sufficient), latest memory/daily/2026-09-21.md, CALIBRATION_STATUS.md. Use SESSION_DIARY only for a concrete historical gap. Then inspect:
- docs/model/ACTIVE_DECISION_LEDGER.md, earnings E1-E7 and September 21 additions.
- docs/model/m5_recalibration_contract_20260716.md, especially explicit July 16 literature-process approval.
- docs/model/eqscale_calibration_reconciliation_20260722.md, especially section 3.6 author-confirmed July 23 Floden-Linde plus HSV decision.
- docs/prompts/HANDOFF_fable_eseries_20260723.md.
- code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/README.md and referenced parameter tables as needed.
- docs/model/structural_model_specification_fable_20260920.md and ../claude_review/lead_review.md (actual path relative to this output directory's parent).
- output/model/native_financing_diagnostic_20260919/specification_followup/decision_packet.md and quantification_v1/morning_view.md as needed.
Read narrowly relevant implementation to establish current income concept/aggregation, not a broad solver audit.

The ledger E2 currently overstates own PSID as maintained source. Lead has corrected this verbally: historical explicit decisions favored literature estimates; current author is happy to reconsider and doubts OUR PSID pipeline. Assess afresh, don't treat E2 as binding. Proposed architecture is age profile + persistent AR(1) + iid transitory, without permanent types; this is a diagnostic candidate, not final adoption. Legacy Sept 14 reference remains frozen. The local 20-minute refinement failed its 600-second two-repetition preliminary timeout before any proposals; this is not model failure and not evidence for source selection.

# Questions
1. What source choices were explicitly approved, when, and why? Separate historical approvals from later implementations or inferred approvals. Identify why own estimates entered the current candidate if evidence permits; report gaps instead of inventing reasons.
2. Assess Floden-Linde plus HSV, Sommer/Sullivan/Verbrugge or Sommer fertility, and Boar-Gorea-Midrigan against the CURRENT economic income object. Which estimate measures wages, individual/household earnings, equivalized disposable income, fixed heterogeneity, persistent risk, iid risk or measurement error? Distinguish published PSID estimates from our PSID pipeline and truly independent datasets.
3. Is importing the recommended published process defensible for this model? What needs adjustment (labor supply, taxes/transfers, household composition/selection, age profile), and what would double-count taxes or risk? Do not shrink variance just to fix fit.
4. What four-year mapping is supported by the chosen reference? Distinguish endpoint AR(1) rho^4 and summed innovation variance, four-year average/sum income, and iid averaging. Do not claim a transformation is standard without a primary source. De Nardi 2004 directly forms five-year PSID cells; Doepke-Kindermann 2019 three-year model does not establish a generic AR1 conversion.
5. Give one recommended baseline source/concept, one meaningful robustness alternative, and the minimal remaining author decisions and numerical validation before calibration. Respect the author's request to stop endlessly reopening choices. Explain tradeoffs if the recommendation conflicts with no permanent types. Data and model observers must match when regression moments are compared.

# Required output
Return a self-contained memo <=1800 words with: recommendation first; compact evidence table (source, income concept/sample, stochastic components, model compatibility); historical decision timeline; explicit verdict on July choice; exact minimal next steps. Cite local paths/line numbers and primary-paper URLs/page/table numbers. Distinguish verified facts, inference and unresolved items. No need to repeat calibration loss tables because this task does not assess fit outcomes. Do not adopt architectural-impossibility or family ownership-premium claims from prior reviews.

# Verification
Use Read/Grep/Glob for local evidence and WebSearch/WebFetch for primary literature verification. Verify any quantitative parameter you quote; otherwise omit or qualify it. No code executions, no numerical runs, no git actions, no file edits. Your final response is captured by the launcher as final.md.

# Stop condition
Maximum 20 minutes and 40 turns. Return best bounded assessment with missing evidence identified. No retries, extra workers, policy experiments, or expansion to unrelated specification choices.
