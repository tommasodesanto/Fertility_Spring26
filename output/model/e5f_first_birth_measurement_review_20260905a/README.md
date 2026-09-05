# First-birth housing measurement review

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
in the earlier complete calibration packet. Full-coefficient recovery, actual
omitted-column inspection, alternative-contrast uncertainty and a common
empirical-to-model population/time contract remain outstanding. A new target
contract requires an explicit author decision; none was made here.
