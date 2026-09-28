# Overnight two-page memo renderer

`build_report.py` adapts the evening report's ReportLab layout. It reads saved
artifacts only; it never imports or solves the model. **Torch Slurm only.**
The PDF and two page PNGs still need visual review before delivery. Authoring
this builder is not evidence that the report has rendered or the run passed.

## Inputs and authentication

Default run: `output/model/overnight_calibration_20260928/gated_v1/search`.
It must have `complete.json` with `status=bounded_search_complete`, `selected`,
`common_primary_key`, and the real ledger. Three or four selected candidates
are supported. It recognizes repeat cases by the driver's `repeat_` prefix,
not the `design` field (which holds the selected candidate key).

Pass an independently prepared `authentication.json` with this exact schema:

```json
{
  "status": "authenticated_final_export",
  "complete_sha256": "SHA256 of search/complete.json",
  "contract_sha256": "same SHA256 as complete.json",
  "common_primary_best": {"case": "same selected case at common_primary_key"},
  "primary_fit_sha256": "SHA256 of --fit CSV",
  "lanes": {
    "primary": {"all14_31_byte_exact": true,
      "all17_repeat_export_plot_hashes_exact": true,
      "plots_visually_reviewed": true},
    "identity": {"all14_31_byte_exact": true,
      "all17_repeat_export_plot_hashes_exact": true,
      "plots_visually_reviewed": true},
    "block": {"all14_31_byte_exact": true,
      "all17_repeat_export_plot_hashes_exact": true,
      "plots_visually_reviewed": true}
  },
  "selected_files_sha256": {
    "selected_export/primary/target_fit.csv": "SHA256",
    "selected_export/primary/parameters.csv": "SHA256",
    "selected_export/primary/receipt.json": "SHA256",
    "selected_export/primary/export_receipt.json": "SHA256",
    "selected_export/primary/standard_diagnostics/ACTUAL_NAME.png": "SHA256"
  }
}
```

Fill **every** selected lane, including `common_primary` if present, and every
one of its 17 actual PNGs plus the four named text receipts/tables. Keys are
relative to the search directory. No example placeholder is accepted.
Only set review flags after checking the two original/fresh-repeat table pairs,
plot hashes and visual panels. The builder verifies these pins but does not
replace the lead's independent scientific/source review.

`--fit` is the common winner's full 14-row table rescored under primary weights:
`moment,target,model,gap,weight,loss_contribution,role`. Physical values must be
byte-identical fields to the selected export's table. All 31 parameter rows
are retained; 10 bounded searched rows plus normalized `psi_child` appear in
the PDF. Fit gaps, contributions, totals and parameter bounds are rechecked.

`--findings` is concise lead-reviewed prose (plain text, not markup):

```json
{
  "reviewed_by_lead": true,
  "economic_fit": "About 45 words: principal raw moment misses and tradeoffs.",
  "numerical_review": "About 55 words: search failures, repeats, gates, memory and recovered old timeouts; distinguish runtime failure from equilibrium nonexistence.",
  "price_diagnostic": {
    "status": "pending",
    "text": "About 40 words; pending is honest until results are reviewed."
  },
  "early_diagnostic": {
    "status": "pending",
    "text": "About 65 words: early target/model, other-moment primary loss, and sacrifices."
  },
  "next_steps": "About 40 words: concrete unresolved decisions, no automatic new run."
}
```

Allowed diagnostic statuses: `reviewed`, `pending`, `failed`, `unverified`.
For `reviewed`, additionally supply `evidence_path` and `evidence_sha256`; the
path must exist on the rendering host. Do not call an incomplete probe verified.

Early diagnostic source: `<stage>/early_frontier_v1/run_v1/complete.json` has
`accepted`, `anchor_checks`, `records`, `summary.ranking`, and `summary.nondominated_ids`.
`early_gap_ranking.csv` contains `early_model`, `early_gap`, and
`other_moments_primary_loss`; full target/parameter comparisons are separate.
Only `accepted=true` supports verified exploratory interpretation. It does not
supply fresh repeats for promotion as a calibrated candidate.

Price diagnostic source: `<stage>/price_diagnostic_v1/run_v1/complete.json`
has `cases` with process status; each arm's `result.json` reports actual gates.
`child_finished` alone is not numerical acceptance. These fixed-benefit-level
price-start arms do not normalize completed fertility or certify full policies.

## Render and inspect on Torch

Use the already available report dependencies/approved container. If the
runtime requires the repository's original logical path, use the existing bind
mount as for previous reports. In a Slurm allocation:

```sh
python build_report.py \
  --project-root /scratch/td2248/projects/fertility_night_calibration_20260928_v1/project \
  --authentication /ABS/PATH/final_review/authentication.json \
  --fit /ABS/PATH/final_review/target_fit_primary_rescore.csv \
  --findings /ABS/PATH/final_report/findings.json
```

`--run` and `--pdf` optionally override the search/output paths. Defaults write
`output/pdf/overnight_calibration_20260928.pdf`, plus `page_1.png`, `page_2.png`,
`summary.json`, `report_qa.json`, and supporting CSVs under this report directory.
The builder aborts on overflow rather than silently dropping table rows.
Before the first actual PDF authoring operation, run the PDF skill's
`container_tools/mark_artifact_operation_started.mjs --operation-kind create
--expected-output-count 1 --output-format pdf` once using its available runtime.

No provisional switch is provided: incomplete or unauthenticated final results
must not accidentally receive the passing-candidate language. For a layout test,
use a separate scratch project with an explicitly synthetic, fully consistent
fixture, and label its output clearly; never publish it as the final report.
If the main run fails, write an honest failure memo separately rather than
weakening this renderer's authentication checks.

Final delivery: `../../../pdf/overnight_calibration_20260928.pdf` (two pages; both rendered pages checked). Full target/parameter tables are in `../cluster/final_review/`; all 17 selected standard plots are in `selected_standard_diagnostics/`. Price and early-fertility diagnostic reviews are complete. No new model solves or specification changes.
