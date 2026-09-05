# First-birth housing measurement review

**Completed decision: retain the 0.720246-room target.** Paired Torch job
`17007732` completed successfully in 8m03s (exit 0:0). Its original fit reproduces
the authoritative point to 4.4e-16. Changing the 1986 cohort's reference gives
0.730459 rooms, a +0.010213 difference, with identical 49,457 estimation rows,
4,112 clusters and fitted outcomes equal within 4.41e-12. The three predeclared
common-cohort point summaries are 0.811582, 0.797759 and 0.804369 rooms and are
invariant across these normalizations to machine precision. No new standard
error or replacement target is claimed. The concern does not explain the
model's shortfall. No model solve or policy path was run in this follow-up.

The final aggregate result is `reference_aggregation_check.json`; the paired
coefficients and sample receipts are in `reference_cluster/full/`. The original
`manifest.json` belongs to the earlier audit saved in commit `39fa668`; it does
not pin the revised follow-up driver. `reference_followup_verification.json`
records the completed follow-up, frozen-source hashes and temporary-sample
cleanup. The chronology below preserves the failed local attempts explicitly.

Read `docs/model/e5f_first_birth_measurement_review.md` or
`output/pdf/e5f_first_birth_measurement_review.pdf` from the repository root.
This is a diagnostic review; the primary model, target, weights, and calibration
are unchanged. The review identifies a normalization concern, not a corrected
empirical target or an estimate of its bias.

## Completed and incomplete work

- The primary PSID source/code/output provenance receipt verified successfully.
- `smoke/` reproduces 49,872 input rows, 4,527 women, and 2,486 treated women in
  28.1 seconds and exports only cohort/event aggregates.
- `full/` records one identical regression replay, expected at 512 seconds from
  the prior receipt, stopped at its declared 900-second limit (902.5 seconds
  including termination). No cohort coefficient or covariance objects were
  exported. It is a FAILED regression replay, not a verified re-estimation.
- `design_check/smoke/` is a subsequent sample-only check of the newly identified
  normalization hypothesis, not a second regression. It verifies nonmissing
  regression inputs, person-constant birth cohorts, and the exact event-indicator
  partition; all checks pass in 84.4 seconds. Its input support matches the first
  smoke. The independent read-only algebra review is retained separately.
- `post_run_verification.json` confirms the primary source and every primary
  output remained unchanged after the timed-out regression.
- `measurement_summary.json` contains the aggregate findings, checked null
  direction and saved-model comparison. `normalization_null_directions.csv`
  retains every cohort's endpoint shares and the relevant null loading.
- `all_target_fits.csv` and `all_parameters.csv` retain all rows for the 35
  coordinate/joint cases. No new equilibrium, policy path or calibration ran.

The three byte-identical input-support exports are stored once in Git as
`input_cohort_event_support.csv`; the originals remain available locally. The
post-run verification records their equality and shared hash.

The original plan was one full unchanged regression, hard stop 900 seconds,
with no automatic follow-on specification or calibration. Progress was written
at least every two seconds. The second sample check tests a different hypothesis
using exact design identities and has a 180-second limit. No individual records
were exported. Source hashes and generated do-files identify each actual replay;
the reusable driver subsequently gained the design assertions and failure-path
hash checks, so its current bytes are not the original full replay source.

## Reproduction

Regenerate all reported aggregates and supplemental figures from saved receipts:

```sh
python3 code/model/tools/analyze_e5f_first_birth_measurement.py
```

Run a new sample-only check in a fresh folder:

```sh
python3 code/data/psid_followup_mar2026/audit_first_birth_rooms_measurement.py smoke --output OUTPUT_DIRECTORY --timeout 180
```

The original full replay command is retained for diagnosis, not recommended as
an automatic retry:

```sh
python3 code/data/psid_followup_mar2026/audit_first_birth_rooms_measurement.py full --output OUTPUT_DIRECTORY --timeout 900
```

Rebuild the PDF after checking the source:

```sh
python3 code/model/tools/build_e5f_independent_audit_pdf.py --source docs/model/e5f_first_birth_measurement_review.md --output output/pdf/e5f_first_birth_measurement_review.pdf --date '5 September 2026' --heading 'FERTILITY / HOUSING MEASUREMENT' --no-source-index
```

The two new plots are supplemental. The standard seventeen model figures remain
in the earlier complete calibration packet. The follow-up recovers cohort point
coefficients and verifies the actual reference change. Alternative-contrast
uncertainty and a common empirical-to-model population/time contract remain
outstanding. A new target contract requires an explicit author decision; none
was made here.

## Decisive reference-choice check (September 5, after author challenge)

Scope: exactly two interacted least-squares fits using the primary sample,
weights, controls, fixed effects and cohort/event columns. First reproduce the
existing point estimate. Then explicitly omit F1event for the 1986 first-birth
cohort, whose missing K=-2 observation creates the verified null direction.
This changes a normalization within the same regression column space. Compare
every fitted housing value and require identical estimation membership. No
alternative empirical target or standard error will be proposed from this test.

The expensive cohort-share covariance calculation is skipped; cohort shares are
computed from the same marked input sample and interacted regression is the
same reghdfe call. This is a materially changed method from repeating the failed
full eventstudyinteract covariance pipeline. Baseline point reproduction must
pass within 2e-6 rooms and maximum fitted-value difference within 2e-6 rooms.
The exact two-fit loop passed on a synthetic panel in 8.1 seconds with maximum
fitted-value discrepancy 5.83e-16. Full input remains 49,872 rows; primary fitted
sample must be 49,457 rows and 4,112 clusters. Estimate: two fits, conservatively
up to about 17 minutes using twice the prior 512-second full-estimator runtime;
the skipped covariance should reduce this. Hard total budget 1,200 seconds.
Local licensed Stata and local raw data avoid a cluster staging operation for
this under-30-minute check. No follow-on search is launched.

Commands:

```sh
python3 code/data/psid_followup_mar2026/audit_first_birth_rooms_measurement.py reference_smoke --timeout 120
python3 code/data/psid_followup_mar2026/audit_first_birth_rooms_measurement.py reference_full --timeout 1200
```

The local `reference_full/` run wrote a heartbeat every two seconds but timed
out before its first coefficient export. The completed-fit receipts, coefficient
exports and final comparison were subsequently recovered in
`reference_cluster/full/`.
The shared helper is
`code/data/psid_followup_mar2026/audit_first_birth_rooms_reference_invariance.do`.
No individual residuals or other microdata are exported.

Before either full coefficient export is available, predeclare one supporting
arithmetic comparison: average within-cohort (+3 minus -1) changes among cohorts
observed at both endpoints in the fitted sample, using common weights at both
endpoints. Report three explicit descriptive weighting rules: original -1 shares,
original +3 shares, and their mean, each normalized over the common cohorts.
Check that these point summaries agree under both fitted normalizations. This
requires no extra regression and does not supply a new empirical target, its
standard error or proof of a causal interpretation. It helps separate a formal
normalization failure from a large numerical change under consistent aggregation.

The local `reference_full/` run timed out at 1,200 seconds without a first-fit
checkpoint. `reference_cluster/` retains the completed paired test on Torch, job
17007732, with a one-hour hard limit, 8 CPUs and 32 GB. It used Stata 19 under
`version 17.0` with a pinned copy of the installed ado sources. Synthetic
exact-loop job 17007687 completed both fits (maximum fitted difference
7.49e-16); its shell marker check failed because rg is unavailable. That check
alone was replaced with grep and verified against the successful log before
launching the full job. A local rsync option incompatibility was resolved before
any full job launched; see `cluster_preflight.json` for the receipt.

A temporary 2.19 MB minimal analysis sample was built by the unchanged primary
sample code, verified at 49,872 rows and 4,527 women, and transferred to the
author-owned private Torch directory with a matching SHA-256. It contains the
person identifiers needed for fixed effects and clustering and therefore stays
outside Git; only aggregate outputs are retained in this report. This explicit
cluster staging supersedes the earlier local-only no-microdata-export plan for
this follow-up, without changing the estimation sample. The private working
copies are temporary and the original PSID source is untouched.

After collecting the remote `full/` folder into `reference_cluster/full/`, run:

```sh
python3 code/model/tools/analyze_e5f_first_birth_measurement.py --reference-check-only --reference-directory output/model/e5f_first_birth_measurement_review_20260905a/reference_cluster/full
```

That command requires a completed Stata marker, the original row/cluster counts,
fitted-value invariance, and reproduction of the retained 0.720246 point before
it reports the paired estimates and descriptive common-cohort comparisons.
