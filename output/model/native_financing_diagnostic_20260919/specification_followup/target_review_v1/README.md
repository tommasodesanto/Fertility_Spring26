# Target validation before the utility comparison

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
