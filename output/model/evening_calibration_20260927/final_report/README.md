# Final evening calibration memo

The author-facing report is `output/pdf/evening_calibration_20260927.pdf`, exactly
two pages. `page_1.png` and `page_2.png` are its Torch-rendered visual QA pages.
The memo reports the common-primary best across all accepted search cases,
`initial_0347_block`, rather than comparing incompatible lane objective scalars.

`summary.json` retains all 14 targets and all 31 parameter rows at full precision.
`common_best_target_fit.csv` contains the primary-weight comparison. The three
lane-prefixed target, parameter, and export-receipt files preserve alternatives.
`report_qa.json` records page/count checks and the independently recomputed loss.
The separate standard 17 diagnostic figures remain under the authenticated final
review packet; this memo does not replace them.

## Source and verification

Source: the Torch project under
`/scratch/td2248/projects/fertility_evening_calibration_20260927_v1/project`,
`output/model/evening_calibration_20260927/gated_v2/{search,final_review}`.
Contract v4 SHA256:
`4453e92f712b1b6b6d3b9a11a4b16c0231a05ce314d33b207ec50d3a4c5475be`.

`build_report.py` checks completed status, 360 search cases, six successful fresh
repeats, the common-primary winner, complete target/parameter counts, all gaps,
all primary loss contributions, and the final authentication flags. It then
renders with ReportLab and PyMuPDF on Torch. No model is imported or solved.
The underlying checkpoint and source hashing was performed by the independent
final collector, not duplicated in this report builder.

The first rendering-only job, 18683896, failed immediately on a missing closing
bracket in the comparison-table list comprehension. It generated no report and
did not touch model outputs. Corrected rendering job 18683921 produced the first
valid two-page PDF. Final rendering job 18683951 applies only wording precision:
validation comparisons are not labelled formal failures; memory is rounded to
102.1 GiB. Rendering logs are retained here. Both final pages require visual QA
after collection; the lead receives their paths.

Reproduce on Torch with ReportLab and PyMuPDF available:

```bash
python build_report.py --project-root /scratch/td2248/projects/fertility_evening_calibration_20260927_v1/project
```

No optimum, identification-rank, transition, or policy claim follows from this
bounded search. Further work requires the next author decision.
