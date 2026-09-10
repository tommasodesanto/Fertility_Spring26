# Calibration architecture: decision research

Completed 10 September 2026. Research and empirical diagnostics only: no model,
target/weight contract, parameter estimate, numerical gate or cluster job changed.

Read [DECISION_REPORT.md](DECISION_REPORT.md) or the eight-page PDF at
`output/pdf/calibration_design_decision_20260910.pdf`. The first two pages give the
recommendation and main economic issue; the remaining pages document targets,
measurement, literature, computation and the complete inherited numerical tables.

## Recommendation, not adoption

Calibrate an approximate pre-announcement economy using NCHS 2003-2006, CPS
2004/2006, ACS 2005/2006 and PSID 2003/2005 observations; retain the reviewed pooled
birth-housing response and explicit external restrictions. Assess initial fit and
identification, then freeze common parameters and fit one announced 2007-2023
preference-shock amplitude to a correctly dated late fertility-rate window. Use
the resulting 2023 distribution for matched policy scenarios. The proposed early
system has 13 restrictions for 10 common parameters plus the initial preference
level; the local weighted Jacobian is still required.

The recommendation replaces the earlier *joint dated fit* proposal as the lead's
preferred design, following newly established early-data feasibility. It does not
silently change the active production objective. Both designs remain coherent;
joint estimation with the same stationary initial state does not reconstruct old
cohort histories. The old completed-fertility mean is a visible diagnostic in the
proposed period-oriented initial benchmark, never a production target removed by
this review. New parameter/target assignments and unresolved observers are explicit.

## Evidence map

- [Early housing review](early_housing.md), [all target candidates](housing/early_housing_target_candidates.csv),
  [saved metro components](housing/early_housing_metro_components.csv),
  [source receipt](housing/early_housing_source_receipt.json),
  [original-current-sample reproduction](housing/current_2023_overlap_verification.json).
  The lead added the strictly pre-announcement 2005/2006 window to the same saved
  components; it requires no second ACS raw pass. Bootstrap covariance and draws
  for that window are included. All city/sample percentages refer to their stated
  denominators: 48.9% in 2012 and 4.54% in 2023 are the old filter's excluded shares
  of the full 42-city head weight, not additions relative to admitted weight.
- [Early wealth review](initial_wealth.md), [aggregate results](wealth/aggregate_wealth_results.csv),
  [old-wealth results](wealth/old_wealth_results.csv), [verification](wealth/verification.json).
  All six windows are retained, including strictly pre-announcement 2003/2005 and
  2005 alone. Person-cluster bootstrap uncertainty preserves the authoritative
  builder's definitions; this is not a complete complex-survey variance estimator.
- [Early fertility extraction](../parameter_target_audit/fertility/README.md),
  [cohort/period coherence review](../parameter_target_audit/design_research/early_fertility.md),
  [independent design critique](../parameter_target_audit/design_research/fertility_design_check.md).
  The mean first-birth age uses the established boundary-collapsed model midpoint
  bins; the age-30+ share uses actual maternal age. These are period birth-flow
  observations, not the timing histories of the old CPS cohorts.
- [Computational assessment](computation.md), including the additional renewal
  accounting obligations, empirical initial-state alternative and four-year
  calendar mapping. Early-data feasibility in the decision report supersedes the
  assessment's earlier conditional paragraph about missing initial housing rows.
- [Final existing numerical receipts](computation/final_replay/collection_receipt.json),
  [root summary](computation/final_replay/summary.json),
  [full fit](computation/final_replay/evaluation_003/target_fit.csv),
  [all parameters](computation/final_replay/evaluation_003/parameters.csv).
  The finite 100-date price path converged at inherited parameters. Horizon
  certification, new calibration and matched policy results remain outstanding.
- [Lead verification](lead_verification.json) records independent component,
  bootstrap, hash and numerical arithmetic checks. [PDF QA](pdf_qa.json) records
  full-page visual inspection and complete numerical-cell checks.
- [Public-only Claude research](claude_public_literature_review.md),
  [safe prompt](claude_public_prompt.txt), [execution receipt](claude_public_execution.json).
  This is a lead sheet, not canonical research truth. The decision report uses only
  lead-verified primary sources. The attempt to send a project-specific research
  prompt was rejected before execution; the public-only replacement had no project
  access or details and used two research agents. It did not attach to the user's
  separately opened GUI instance.

## Reproduction

Run `build_decision_report.py` with the bundled Python runtime to regenerate the
PDF and Markdown from the saved receipts. No model code is executed. The PDF
builder requires ReportLab and the system Arial fonts. The original empirical
scripts and their commands are documented in the housing/wealth/fertility notes.
Use `housing/summarize_early_housing.py` to regenerate the ACS candidate table and
covariance from saved components without rereading raw data.

Current unresolved prerequisites: certified female fertility exposure and maternal
age mapping; calendar allocation; completed-cohort and family-group observers;
new CPS/timing uncertainty and objective weights; the initial 2.1-linked renewal
law, queues and entry gates; source/role of housing supply restrictions; local
identification; horizon certification; and the population/fiscal/geographic policy
contract. National fertility/PSID and metro housing are an explicit maintained
geographic approximation, not a new empirical equivalence claim.
