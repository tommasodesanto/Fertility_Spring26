# Post-presentation issues ledger

## Author decision — September 12, 2026

For the September presentation, use the historical May empirical plots and
original regression specification. Tommaso has not had time to review the
subsequent measurement and design revisions. This is a presentation-version
decision, not certification that the historical measurement is correct.
Suggested disclosure: **Original May specification; measurement revisions
under review.** Do not describe that specification as timing-corrected.

The May rooms graph is reproducible: all 18 saved points match the recognized
original specification to numerical precision. The rooms survey-year assignment
problem is verified against the assembled source panel, rather than merely
conjectured. Preserve both findings when discussing the historical results.

The new annual and binned comparisons remain diagnostics. Do not replace the
calibration target or its weight, relabel the historical graph, or promote a
new specification without reconciling the empirical and model definitions.
Empirical follow-up is deferred; the completed data jobs and their collection
automation are stopped/paused. No new empirical regressions are launched by
this decision. Resume the items below after the presentation review.

## Numbers that must not be conflated

| Object | Rooms | Interpretation |
|---|---:|---|
| May graph coefficient at +3 | 0.796858555 | Original normalization omits both −2 and −6. |
| May graph +3 minus estimated −1 | 0.740737457 | Four-year contrast calculated from that same historical curve. |
| Current pinned calibration target | 0.720246262 | Later August specification; +3 minus −1, different sample/weighting/control construction. |
| Historical slide scalar | 0.66 / 0.664 | Does not equal the saved May graph's +3 coefficient; source remains to be reconciled. |

The comparable-horizon scalars 0.740737457 and 0.720246262 differ by
0.020491195 rooms (about 2.8% of the May contrast). This numerical proximity
does not establish equivalent estimands or validate either specification.
Tommaso's preferred pre-birth reference remains −2; subtracting −1 answers a
different question. A −2-to-+3 comparison spans five calendar years, whereas
the currently pinned model mapping spans four years.

## Deferred empirical work

All rows remain **open / deferred**, except the reproduction result recorded
above. Close each with a reproducible result and an explicit specification
decision; smoothness or statistical significance alone is not a criterion.

| Priority | Issue and established evidence | Required resolution |
|---|---|---|
| 1 | Rooms answers are stored one interview early in the merged panel. Every assigned value in the 352,250-row prepared common sample matches its survey-year source after alignment. | Reconstruct rooms directly from the source-year crosswalk; recover the 52,317 source-observed values not restored by shifting; validate year-specific response codes. Preserve May reproduction separately. |
| 1 | Annual event-time support alternates after PSID becomes biennial. Some cohorts lack the intended −2 reference; May also omits −6. | Choose exact years versus two-year windows, retain a reference before the final pre-birth year, and document cohort support and aggregation weights. The tested baseline −3/−2 is a changed object, not exactly −2. |
| 1 | The first binned diagnostic is much smaller: 149,402 fitted rows versus 345,751 in the prior annual common fit. | Review the restrictions individually: exclude 2019+ dates, exclude no-recorded-birth-year observations, require all six displayed windows, then estimator exclusions. Do not attribute the full difference to binning or data transfer. |
| 1 | The binned driver excludes 101,848 rows without a recorded first-birth year after its date restriction. These may include childless people and unknown histories. | Reconstruct history status before deciding which observations can be controls. Compare admissible last-treated and confirmed-childless designs separately; report their populations and identifying assumptions. |
| 1 | A last-treated 2019 cohort cannot remain an untreated comparison group in 2019 and later. | Either restrict admissible dates/cohort comparisons or choose and justify another control group. Preserve this change separately from outcome timing. |
| 2 | Ludovica code constructs rounded frequency weights and uses them in csdid blocks; its separate Sun–Abraham command reproducing May has no weight argument. August uses direct IW probability weights. | Compare unweighted and correctly defined survey-weighted estimates on identical samples; state the population represented. Do not claim the original analysis generally forgot weights. |
| 2 | August restricts to women who are reference persons/spouses, selects one per household-year, and excludes multi-family-unit dwellings. | Decide the observational unit and population; quantify each restriction, repeated household outcomes, and appropriate weighting/clustering. These are not interchangeable with mechanical corrections. |
| 2 | August defines first biological birth across 20 child records and all available histories, unlike the original first-child field. | Compare birth dates and exclusions person by person; distinguish measurement corrections from a changed biological-parent population. |
| 2 | HOMEOWN has a recovered 41-wave source-year mapping; no evidence justifies applying the rooms shift to it. Original regressions exclude tenure 'neither owns nor rents'. | Numerically validate the mapping, decide the ownership denominator, then rerun ownership with the chosen event design. For ownership transitions, use the previous observed interview rather than calendar-year L. |
| 2 | Moving variables lack a recovered upstream year crosswalk. In 2019 the question covers moves since January 2017 and dates the most recent move; recall conventions vary by vintage. | Recover/validate source mapping, move dates and recall intervals; distinguish most recent move from any move, first-mentioned reason from all reasons, non-movers from missing responses. Bin only after establishing the outcome clock. |
| 2 | Reason code 3 includes expansion/better housing; code 6 includes neighborhood, schools and proximity to friends/relatives. | Check vintage-specific codes and use accurate outcome names. Revisit move gating, missing-to-zero errors and reason-response denominators. |
| 3 | Presentation numbers, empirical contrasts and calibration mapping differ. | Reconcile 0.66/0.664, 0.797, 0.741 and 0.720; choose the intended horizon and population; regenerate the target estimate, covariance-based uncertainty, weight and provenance together before recalibration. |

## Evidence and reproduction

- Recognized original code: `/Users/tommasodesanto/Desktop/Projects/Fertility/Codes/code_per tommi_addingcontrolsandfixingthings.do`.
- May rooms graph: `/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Graphs/rooms_f_c_y_all.png`.
- May saved coefficient table: `/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Tables/rooms_f_c_y_all_estimates.dta`.
- Consolidated audit: [first-birth correction review](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/README.md).
- Binned diagnostic: [results and sample flow](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/binned_rooms/summary.csv), [sample exclusions](../../code/data/psid_followup_mar2026/output/first_birth_correction_review/binned_rooms/sample_flow.json).
- Source-year ownership construction: `/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/Construction_Files/Code/01 Collect housing variables.do`, lines 95–118.
- Moving definitions: [PSID 2019 family codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/FAM2019ER_codebook.pdf), pp. 54–56; [1984 family codebook](https://psidonline.isr.umich.edu/documents/psid/codebook/FAM1984_codebook.pdf), p. 161.

Presentation asset changes are owned by the existing September presentation
task. This data task communicated the author's choice and the distinction
between the historical graph and the pinned calibration scalar. The ledger
does not certify restoration of other historical figures whose source assets
have not yet been verified.
