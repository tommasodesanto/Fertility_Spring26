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
