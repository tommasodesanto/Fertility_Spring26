# Handoff: resolve the deferred PSID first-birth housing review

Work in `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

I need clarity on an empirical issue we postponed just before the September
presentation. We used Sun–Abraham throughout. Correcting the interview dates
assigned to PSID rooms reduced the response sharply. We chose to present the
original figures and revisit this afterward. That review is now the task.

Start by explaining, plainly, what was wrong, what is established, and what
remains undecided. Do not start with a broad literature review, model/data
matching, or a new estimator. Do not equate replication with validity. Push
back if our existing interpretation is wrong.

Current author instruction: record **0.770 rooms provisionally** while we
discuss this. This is not approval of its timing or a finalized target. Do not
change model targets, weights, frozen results, paper text or Google documents.

Verified evidence to inspect:

- Original-date and timing-corrected fits use the same 345,751 observations,
  36,026 person clusters and Sun–Abraham command. Their +3 coefficients are
  0.770 and 0.403. The original code omits BOTH -2 and -6, so do not simply
  label this a clean -2 comparison.
- Their -1-to-+3 differences are 0.758 and 0.238: different statistics from
  the +3 coefficients. The author had expressly preferred a -2 reference.
- The rooms shelf assigned interview outcomes to the preceding interview;
  source validation supports moving them forward. Inspect the evidence and
  its limits, including missing coverage.
- The later household-based regression already corrects room dates and
  reports 0.720 for -1-to-+3. It differs in reporters, controls, weights,
  birth-history construction and event indicators. No completed controlled
  decomposition explains the gap from the original corrected regression.
- An overnight calculation reused the later regression and averaged changes
  within the same 27 birth cohorts using common endpoint weights, yielding
  0.798–0.811. These are candidate aggregations, not a reversal of the timing
  correction and not adopted targets. Distinguish Sun–Abraham cohort-share
  aggregation from PSID survey weights. Changing cohort shares is not itself
  an estimator error.
- Some cohorts have no observed -2 interview. A saved normalization test of
  the later regression changed the reported target from 0.720 to 0.730 while
  preserving fitted outcomes. Assess what this establishes and what it does
  not establish; do not interpret it as an estimate of total bias.

Read the mandatory repository startup, then these sources (relative paths):

1. `code/data/psid_followup_mar2026/output/first_birth_correction_review/README.md`
   — read the later completed-review sections as well as the initial audit.
2. In that folder: `reference_m2_summary.csv`,
   `timing_only_comparison_summary.csv`, `reference_m2_support.json`,
   `all_wave_validation_receipt.json`, and the `binned_rooms/` results.
3. `code/data/psid_followup_mar2026/audit_original_rooms_timing.do`
   and `sa_rooms_first_birth_household_aligned_v1.do` in the same code folder.
4. `output/model/e5f_first_birth_measurement_review_20260905a/reference_cluster/full/reference_comparison.json`.
5. `output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/README.md`
   and its `covariance_replay/analysis/README.md` for completed covariance work.
6. `CALIBRATION_STATUS.md`, September 12 entries, and
   `docs/model/POST_PRESENTATION_ISSUES.md` for the author's deferral.

Give me one clear comparison table: specification, actual date assignment,
reference support, comparison group, sample, survey weighting, reported
contrast, estimate and uncertainty. Then explain the verified timing result
with one concrete interview example. Separate demonstrated errors from
defensible empirical choices. Recommend how to handle biennial interviews
and missing reference support while retaining Sun–Abraham where appropriate.
Explain any proposed sample restriction and the population it would exclude.

If further regressions are necessary, specify the smallest controlled sequence
that isolates the unresolved differences and what each run would tell us.
Discuss that plan with me before expensive work. Do not rerun completed checks,
choose a specification because its effect is larger, or silently substitute
a different pre-birth reference. Keep model matching for a subsequent discussion.
