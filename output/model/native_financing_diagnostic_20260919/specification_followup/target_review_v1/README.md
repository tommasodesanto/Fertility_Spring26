# Target validation before the utility comparison

## September 23 — national depreciation and property tax adopted

Author adopts annual depreciation **1.416%** (2007 national BEA/Fed)
and property tax **1.060%** (2007–2011 ACS national owner-tax/value ratio).
Definitions, period/aggregation differences from DUE, retained official data,
full-precision values and reproducible calculation are in
[national housing inputs](national_housing_inputs/README.md). Implementation and paper/slide reconciliation pending. No model run, frozen
target or objective changes.

## September 23 — curvature retained with further review

Author retains sigma = 2. Review its interpretation and literature support for
the consumption–housing composite, floor and child-dependent shares later.
No numerical change. Bequest source/input workers are complete, but the 2007
estate-flow calculation has not been performed; do not mark remeasurement done.

## September 23 — provisional first-birth rooms entry

Author records 0.770 rooms provisionally and takes the deferred timing and
reference-period issue to Claude. This labels the original-date common-sample
+3 coefficient, not the corrected-date estimate or a finalized target. Frozen
objectives and weights remain unchanged. The handoff is
`docs/prompts/first_birth_timing_claude_handoff.md` at the repository root.

## September 23 — housing supply review remains open

[Source and implementation review](housing_supply_review.json) verifies that
elasticity 0.63 is present in the checklist, while its empirical provenance is
unresolved. H0 is the supply level at the reference user cost. DUE uses observed
initial rents to anchor its supply level; KMV uses construction employment.
The proposed quantity/cost anchor is for discussion, not an adopted target or
specification. Capped mean rooms and physical supplied stock are distinct.
No new model runs, target/weight changes or paper edits.

## September 23 — provisional estate-flow target

Author requests showing 0.880% annual bequests / aggregate wealth as a
provisional literature benchmark while 2007 recomputation proceeds. Do not
label it a remeasured 2007 estimate or completed empirical work. Theta1 stays
externally fixed as recorded below; continue reviewing the other targets.

## September 23 — externally fixed theta1 chosen

Author fixes the bequest wealth shift using DUE's 1%-of-median-annual-earnings
restriction: theta1 = 0.008193084126995582 (display 0.008), mapped from the
current 15-state persistent-income model's gross working-income median
0.8193084126995582 and mean 1. Source: the saved actual-parameter extraction
in `overnight/empirical_entry/common_scale_candidate/actual_frozen_parameters.json`.
Use the age/state weighted inverse-CDF median; do not multiply by four.
This replaces the proposed internal theta1 calibration for the next run,
not the frozen experiments. Theta0 stays internally calibrated; its target
remeasurement remains pending. Final target/weight bookkeeping is still open.

Low-priority follow-up: compare against a run estimating theta1 internally,
with a justified estate-shape target and other parameters allowed to adjust.
No new run or paper edit.

## September 23 — beta data-year decision

Author chooses pooled PSID 2005/2007: aggregate net wealth / annual gross
labor earnings 6.927, full precision 6.92658379107299, existing bootstrap
SE 0.417310142072186. Source: `initial_2005_2007` in the September9 design
research `wealth/aggregate_wealth_results.csv`. This supersedes the pending
vintage and 6.146 discussion below for the next target set. Weights remain
open; no frozen experiment, model run or rescoring changes.

## September 23 — beta target approach accepted

Retain aggregate net wealth / annual gross labor earnings as the main beta
moment, including housing net of debt. Existing value 6.146 is unchanged;
the data-year choice remains open. A supplementary literature review is
delegated, separately from the discussion document, for lead review before
its findings are accepted. References belong in
`docs/literature/economic_decisions.md`. No model runs, weights or rescoring.

## September 23 — saved-search identification check

[Descriptive evidence](fertility_saved_search_identification_check.json) covers
all193 B_floor evaluations. Mean first-birth age and the age30+ share have
correlation0.999 in this sample. This raises a concrete question about their
distinct information but is not a formal identification test. Saved isolated
coordinate comparisons exist for H0, beta and chi, not the three fertility
parameters. No new solves or rescoring; retain chosen targets and open
identification status.

## September 23 — age-grouping decision confirmed

The author accepts the current common four-year grouping and bunching births
before age18 into the youngest model age group as an explicit approximation.
Keep first-birth mean age 25.976 and age30+ share 24.928%. No numerical target
changes. Separate parameter identification remains under review; cost removal
and a common fertility taste scale remain two distinct deferred checks.

## September 23 — author requires measurement and identification review

Agreement on proposed fertility outcomes does not close their age/population
mapping or prove separate parameter identification. Review both before closing
the block. [Bounded arithmetic and code review](fertility_measurement_identification_review.json)
reproduces the existing grouped mean (25.976) and distinguishes raw completed-year
ages (25.161) from an annual midpoint proxy (25.661). First births below age18
are 7.731% of the data and are collapsed into the youngest model cell; this is
a substantive support approximation, separate from common four-year bins.

Main parameter–target associations are hypotheses about informative variation.
Neither reproducing scalars nor having more moments than parameters proves
identification. The cost-removal and common-taste-scale checks remain distinct.
No targets, weights, frozen sources or model runs changed in this review.

## September 23 — decision: retain birth timing; two deferred checks

Retain mean first-birth age 25.976 years and first births at ages 30+ of
24.928%, from NCHS 2003–2006 with the existing model age grouping. The author
accepts these targets alongside the two CPS shares and 2.10 normalization;
full-precision values and frozen experiments are unchanged.

Keep two distinct follow-ups after recalibration: test whether the first-child
cost can be removed, and assess whether the first-child and additional-child
taste scales (kappa_1 and kappa_C) can be replaced by one common scale. Let
remaining parameters adjust when evaluating simpler specifications. Retain
both scales and the cost now; no new run or parameter removal.

The author also requested a bounded Luna check of the continuing slides'
fertility formulas against implemented code. This is read-only and separate
from the deferred recalibration tests.

## September 23 — decision: retain both CPS fertility shares

The author retains 19.828% childlessness among women aged 40–44 and
21.366% exactly one child among mothers aged 40–44, with the existing
pooled June 2004/2006 CPS definitions and sampling weights. Objective weights
and birth-timing targets remain separate decisions.

Deferred check: after recalibration, compare fit and behavior with a version
that sets the first-child cost to zero and permits the remaining parameters to
adjust. Consider removing the cost if it proves unnecessary. Keep it for now;
no new calculation, target-value change or frozen-bundle edit.

## September 23 — decision: retain the 2.10 normalization

The author retains 2.10 mean children ever born at ages 46+ as the imposed
completed-fertility normalization for the reference steady state, determining
the fertility utility scale separately. It is not a measured 2007 fertility
statistic or period TFR. This closes the normalization decision only; other
fertility targets remain open. The existing decision PDF already recommends
retaining this value.

Add alternative completed-fertility normalizations and transitions from different
initial conditions to the long-term, low-priority list. No calculation now and
no change to frozen target values or experiments.

## September 23, 14:31 EDT — decision: use the national housing sample

The author chooses a nationally weighted U.S. housing sample instead of the
42 selected metros. Keep the 2007 reference. Other measurement choices remain
open; do not overwrite the frozen experimental targets or present the existing
fits as national calibrations.

The four national point estimates with the existing definitions are already in
[`../housing_profiles_v1/full/target_recomputed.json`](../housing_profiles_v1/full/target_recomputed.json).
The author requests these four values in the discussion document and mock draft
now. Recalculating uncertainty and choosing weights remain open; losses and
rescoring are deferred until the next calibration run. Continue the remaining
target decisions without starting those calculations. The decision PDF now
shows the national values; the historical selected-fit CSV and frozen
numerical contracts retain their original 42-metro targets.

## September 23, 12:08 EDT — reviewed lunch packet ready

Read the [five-page decision PDF](../../../../pdf/calibration_lunch_decisions.pdf)
first. It gives six recommended choices, all 13 current targets and their
parameter connections/confidence flags, followed by a supplemental native
first-birth timing diagnostic. These are recommendations for author decisions,
not adopted changes. The July economic definitions are the starting point.

- [Complete current targets and both selected fits](lunch_targets_and_selected_moments.csv)
  includes model values, gaps, actual weights and loss contributions; the
  selected fits and objective are unchanged.
- [Decision source](lunch_decision_review.json) separates the choice, evidence,
  consequence and operation required after the choice.
- [First-birth comparison](lunch_first_birth_comparison.csv),
  [supplemental figure](lunch_first_birth_branch_comparison.pdf),
  [scientific review](lunch_first_birth_final_review.json), and
  [PDF/table QA](lunch_pdf_final_qa.json) retain exact evidence.

The active scored destination response reproduces exactly at selected B_floor
and B_shares. Immediate responses are 0.902596503 / 0.997422993 rooms; destination
responses are 0.987250816 / 1.043654597. Moving the native reading to the origin
therefore reduces it by 0.084654313 / 0.046231604, but neither native reading
automatically equals the PSID -1/+3 regression. A common simulated-panel
estimator remains unimplemented. Do not describe this as having shown that the
empirical miss survives (or disappears under) identical measurement.

The helper uses actual selected checkpoints, g_pre, and the active
begin/finish observer with destination continuation births. It does not use
the different legacy native statistic. A one-cohort smoke and two full
fixed-policy measurements all completed 0:0; no household/GE solves or searches
were performed. All source, checkpoint, helper and loader identities are
recorded. `lunch_first_birth_helper_review.json` preserves rejected draft
approaches and the independent review. No jobs remain active.

Rebuild the PDF and supplemental native-date plot (the original 17-plot sets
remain unchanged):

```sh
code/model/.venv/bin/python code/model/tools/build_e5f_target_review_pdf.py \
  --packet output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1 \
  --decision-review --output output/pdf/calibration_lunch_decisions.pdf
```

## September 23, 11:29 EDT — bounded first-birth measurement calculation

**Author priority, clarified before departure:** The post-lunch meeting must
start with decisions, not another request for two hours of review. Overnight
research is complete. An additional 20-minute agent, `lunch_decision_evidence`,
assembles all 13 rows and cross-cutting choices from existing evidence into
`lunch_decision_evidence.json`; the lead supplies recommendations. For every
choice state what can be decided immediately, the consequence, and the exact
calculation (if any) required after that choice. An unfinished first-birth
measurement calculation must not delay decisions on the remaining targets.

Grounding is complete; receipts are `lunch_first_birth_empirical_inventory.json`
and `lunch_first_birth_model_inventory.json` (agent write in progress).
Joint PSID histories exist locally; selected model checkpoints remain remote,
and Torch access was freshly confirmed. The first bounded model calculation
will reproduce the saved destination birth/control contrast and extract the
immediate-origin contrast from the same birth-weighted cohorts. This is a
measurement diagnostic, not yet an identical empirical regression. Agent
`first_birth_model_match` owns the new measurement helper and focused test;
lead must review its core and verify saved-moment reproduction before running
it remotely or interpreting any new contrast. No automatic annual projection
or causal claim. Empirical four-year pairs may contain intervening interviews;
do not unnecessarily select only four-year gaps between adjacent records.

The author is away for approximately 90 minutes and requests a useful calculation
before resuming calibration decisions after lunch. The lead selects one question:
how much of the first-birth housing-fit discrepancy reflects unlike model/data
measurements? Deadline for a reviewed result or precise limitation: about 13:00 EDT.

Two ten-minute read-only agents inventory the empirical sample/estimator and
saved-model propagation respectively (`first_birth_empirical_match` and
`first_birth_model_match`). The lead then specifies the common statistic before
inspecting its model fit, reviews measurement code, and authorizes only bounded
candidate measurement using existing data and saved policies. Preserve the
original regression/target and objective. No new household/equilibrium solve,
calibration search, target adoption, document rewrite, or broad literature review.

A model comparison must state its birth/interview timing and interpolation;
there is no automatic exact eight-year bridge. Education cannot silently be
removed only from the model estimator. Marginal event-time support cannot stand
in for joint panel histories. Report any reduced empirical estimator as a new
candidate alongside the unchanged original. If checkpoints/access or valid
observation histories are unavailable, stop that dependent calculation and
report the evidence available; do not substitute the inherited seed for selected
floor/share checkpoints or improvise a convenient target. Aim for one comparison
table, one supplemental event-time plot, and a recommendation, with uncertainty
and unresolved ingredients stated. Existing six-phase review is reused.

## September 23 — author-led calibration decisions

The author requests a roughly two-hour first block reviewing the calibration
strategy and every target, with targeted recomputation and rescoring where
justified, followed by a decision tree for earnings and utility by the evening.
The mock quantification section is now explicitly authorized as provisional.
This does not adopt any experimental specification or authorize another search.

Reuse the completed six-phase Opus review and its independently checked
empirical calculations in [the overnight packet](overnight/README.md). Preserve
the rationale of signed-off July choices; reopen a row only for a specific
unresolved definition, changed reference sample, or model-data discrepancy.

The advisor-facing [working checklist](https://docs.google.com/document/d/1hxESCRA89O028-Kx4LmBjdM19R_CobIkn_GkJdMgbbo/edit?tab=t.0) now links to a [Parameters and targets details tab](https://docs.google.com/document/d/1hxESCRA89O028-Kx4LmBjdM19R_CobIkn_GkJdMgbbo/edit?tab=t.umcf44buj7j3) in the same document. It separates completed evidence from pending decisions, lists all 13 working values, and records four next actions. Native readback verified both links and exact preservation of every other main-list paragraph and its formatting/list metadata. Sharing was unchanged.

Working order:

1. State the reference economy, population and units, distinguishing cohort
   fertility, period timing and cross-sectional household moments. Different
   data sources need justified mappings, not mechanically identical samples.
2. Review all five fertility rows, including the separate mean-fertility
   normalization; all three wealth/bequest rows; and all five housing rows.
   Each row needs its empirical object, model measurement, uncertainty,
   economic role and explicit retain/repair/decision status.
3. Review external restrictions, normalization, free parameters and informative
   variation jointly. A count of moments is not an identification test. Explain
   the mixed working weights and any proposed replacement independently of
   which specification fits better.
4. Recompute only the quantities needed by agreed changes. Distinguish a
   target/weight rescore of saved moments, a new measurement from saved model
   outputs, and a change requiring a fresh solution or normalization. Preserve
   all old scores and label candidate contracts separately.
5. Compare the retained reference and experimental fits under common definitions
   where recoverable. Diagnose remaining fit and numerical issues before choosing
   another search. Do not attribute a changed ranking to utility alone when
   earnings, initial distributions, timing or numerical grids also differ.

The two-hour checkpoint should deliver row-level decisions and named blockers,
not an unsupported promise that a new simulated-panel estimator or full refit is
already complete. The evening decision tree should distinguish measurement
repairs, numerical validation, and genuine economic specification choices.

Parallel preparation is bounded: a target decision sheet from existing evidence
(10-minute collection task) and a provisional draft of
`latex/JMP_DS_mock/sections/04_quantification.tex` (20-minute drafting task).
The lead owns identification, target choices, and review. Author draft and
slides are unchanged; no empirical candidate becomes an active target without
an explicit recorded decision.

Preparation completed: [13-row decision sheet](target_decisions_20260923.csv),
with current values, unresolved choices, required calculations and evidence;
the lead checked every row and corrected three overstatements. The mock uses the
reviewed July prose and organization while describing the current 2007 stationary
calibration. The author explicitly corrected the earlier literal restoration:
obsolete July economics and the historical 14/15 tables are removed. The mock
now contains all 13 current working targets and the current parameter roles,
with model-fit values left blank pending specification decisions. The decision
sheet remains current. Both preparation
and restoration checks are recorded in `target_decisions_20260923_review.json`.
These artifacts organize the discussion; they do not mark the target decisions
complete.

**Later author clarification, September 22:** complete target reconciliation is the first task tomorrow morning, September 23. Tonight's experimental comparisons may retain the frozen targets and weights while income, initial wealth and numerical checks take priority. The earlier launch hold in the dated report is superseded; its substantive findings are unchanged. See the latest `CALIBRATION_STATUS.md`.

September 22, 2026. The author made target validation the first step before any further overnight calibration. No jobs, target changes, weight changes or paper edits were made in this review.

The complete 13-row frozen target system (12 scored rows plus the separate fertility normalization) is numerically consistent with the actual scorer. This does not settle whether every empirical and model measurement is comparable. Main remaining decisions are the national versus 42-metro population; the first-birth room-response estimator and matching model measurement; child-group definitions; and wealth/income concepts. Existing room caps, fractional age masks and repaired recent-parent observer are implemented and were checked in the actual scored source. Working scales mix sampling uncertainty, temporal variation and an external tolerance, so they are not uniformly sampling standard errors.

Primary deliverable: [six-page PDF](../../../../pdf/calibration_target_review.pdf). [Complete target catalogue](target_catalogue.csv) gives all target values, weights, scales, empirical definitions, model measurements and assessments. [Reviewed narrative](review.json) drives the report; [integrity check](lead_integrity_review.json) verifies the frozen contract and recorded runtime. Detailed source receipts are fertility_receipt.json, housing_receipt.json and wealth_receipt.json; the latter two include lead corrections to obsolete implementation warnings. [PDF verification](pdf_qa.json) records final artifact identity and checks.

The selected simple-process heterogeneous-entry case took 483.8 seconds for one complete normalized evaluation: six stationary solves totaling 448.4 seconds, averaging 74.7 seconds, plus overhead. This timing is observed, not a guarantee for every parameter point.

Next: choose the quantitative population, settle the event-study comparison, reconcile child and wealth definitions, explicitly record accepted approximations and working weights, then freeze a source/target/observer contract before reconsidering the two-utility search. No target is dropped or replaced without preserving the identifying information for affected parameters.

Rebuild from the repository root:

```sh
python3 code/model/tools/build_e5f_target_review_pdf.py --packet output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1 --output output/pdf/calibration_target_review.pdf
```

This review relies on saved empirical receipts and implemented source; it does not certify a fresh raw-data replication or structural identification. Geography remains pending the author's response.

## July decisions retained

The [bounded July decision trace](july_decision_review.json) confirms explicit July24 saving/bequest sign-off and a deliberately constructed18–24 entrant-wealth proxy with1835 family-years. The lead withdraws the broad before-launch reconstruction recommendation. Compare new income/timing implementation against established definitions; the family-income conversion issue is a compatibility question whose quantitative materiality has not yet been established. No target or entry change follows from this review.
