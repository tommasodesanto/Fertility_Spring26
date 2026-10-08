# Use the October 7 tables in JMP Slides

The verified generator is `code/model/tools/build_jmp_result_tables_v6.py`. The verified fragments are in this v6 folder. V1–v5 outputs are retained intermediate artifacts; use v6.

The continuing JMP deck compiles from `latex/`. For that working directory, include a fragment with the `JMP_slides/` prefix:

```tex
\input{JMP_slides/tables/results_20261007_v6/calibration_cap6}
\input{JMP_slides/tables/results_20261007_v6/calibration_cap5}
\input{JMP_slides/tables/results_20261007_v6/calibration_cap5_hP2_illustration}
\input{JMP_slides/tables/results_20261007_v6/calibration_benchmark}
```

Place each desired fragment inside its own frame. Each fit fragment contains all 14 moments in `Moment | Target | Model` form. The parent-space-floor-two case is an illustration. The README's shorter `tables/...` include syntax assumes compilation from `latex/JMP_slides/`; prefix `JMP_slides/` for the established build directory.

The existing JMP preamble already loads `booktabs` and `tabularx`. The fragments scope their font and spacing locally. The deck itself was not edited.

For every other source table, find its original heading in `manifest.json` and select the corresponding `rendered_fragments`. Horizontal parts repeat row identifiers; row parts retain the original row order. Carry the relevant source definitions, denominators and treatment notes into the frame: the manifest retains the full source text. In particular, `fp_price110` rebalances the rebate, whereas `A6_P110` holds it fixed.

Run the generator from the project root into a new, absent destination:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 code/model/tools/build_jmp_result_tables_v6.py --output-dir /absolute/new/destination
```

Repeat `--project-table PATH` for explicitly identified Markdown or CSV tables. The additional project-file tables remain pending identification.

Verification: eight authenticated source files, 41 tables, zero uncovered source cells; all 344 fragments are byte-identical to the v5 payload compiled with zero overfull boxes. The lead independently checked 4,211 numeric occurrences, literal price multipliers, repeated identifiers and six representative layouts. See `output/model/codex_review_tables_20261007_v1/FINAL_VERIFICATION.json` and its linked QA receipts.
