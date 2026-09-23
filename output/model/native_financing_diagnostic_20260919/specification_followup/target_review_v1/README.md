# Target validation before the utility comparison

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
