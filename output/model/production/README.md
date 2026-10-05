# Standard saved model outputs

- `2007/`: working 2007 calibration snapshot and model–data PDF.
- `2023/`: saved 2023 snapshot from the overnight one-shock transition and the same model–data PDF.
- `transition/`: retained transition results and fertility-path PDF.

One plotting file reads these saved solutions. From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/plot_saved_outputs.py
```

No model solve. No separate page images. `solution` links point to the original retained results, which stay in place; replacing a link selects another stored run. These folders are the standard place to inspect outputs; experiments retain their own directories.

The current 2007 and 2023 snapshots use different parameter points. The one-shock transition remains experimental with the limitations recorded in its source README. Folder naming does not change certification or calibration targets. Data vintages are printed on each PDF.

`2023/model_data_fit.md` contains the concise model–data comparison, all transition target rows, parameter-table link and CPS source check. `2023/target_fit.csv` preserves the original complete transition target system. The plotting file refreshes these with the PDFs.

Mean fertility plots use the existing model 3+ weight (3.602 children); data retain public counts through five or more. Count-distribution bars continue to group 3+. This is a reporting adjustment, with no change to household decisions or population accounting. The completed-fertility tail weight at younger ages and within income groups is an explicit approximation.

## Fast saved-results readout

`readout_index.json` selects explicit saved sources. These commands read only that
small index and the requested CSV; they do not inspect solution arrays, run the
model, search folders, or contact a cluster:

```sh
python3 code/model/tools/read_results.py list
python3 code/model/tools/read_results.py show base_2007
python3 code/model/tools/read_results.py fit base_2007
python3 code/model/tools/read_results.py parameters base_2007
python3 code/model/tools/read_results.py plots base_transition
```

All responses are compact JSON. `fit` returns every saved target row, including
target, model, gap, weight, and contribution fields; `parameters` returns every
saved parameter row. `show` and `list` state the source identity, contract,
classification, pointer readiness, missing sources and limitations. A missing
CSV or plot fails visibly. `base_2007` is the October 3 local working snapshot;
`base_2023` and `base_transition` use the retained one-shock exercise, with a
different parameter point. For those one-shock slots, the `parameters` command
also returns the fitted `psi_child` row as `supplemental_rows`; the 31-row base
parameter table alone does not contain that shock value. `best_estate_a` and
`best_original` are provisional
saved search points; their losses are not native accepted results. The
`transition_latest` candidate is separate from the retained one-shock path and
has three saved first-shock diagnostic panels. Read its status before interpreting it.

Routine questions about these selected results can use the readout immediately.
The index does not test whether a newer run exists or whether source contents
changed after selection. Update a slot only after reviewing its source and
contract. To update it explicitly, copy the current slot object into a small
JSON manifest, add `"slot": "base_2007"` (or the relevant slot), edit its
source and artifact paths, then run:

```sh
python3 code/model/tools/read_results.py select base_2007 /path/to/reviewed_manifest.json
```

Selection checks named artifact paths, writes the index atomically, and retains
the previous slot in `selection_history`. It does not promote a candidate or
change calibration status automatically. Recheck underlying evidence only when
the source identity changes, evidence conflicts, or a requested fact is new.
