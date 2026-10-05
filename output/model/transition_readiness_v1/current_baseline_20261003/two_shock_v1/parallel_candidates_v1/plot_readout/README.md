# Saved c03_h24 transition plots

The [usual fertility plot](fertility_2007_2063.png) shows the saved first-shock path from 2007 through 2063 alongside the four US fertility targets. The [fertility and demography panels](fertility_demography.png) and [housing panels](housing.png) show all 24 saved four-year dates through 2099. Matching PDFs are in this folder. [Plotted values](series.csv) and [source hashes and gates](provenance.json) make the figures auditable.

This is the fixed-preference c03_h24 diagnostic at first-shock preference 0.14736308634876963, from Torch job 19195645. The H24 root and final replay pass. The horizon comparison, terminal certification, scalar estimation and second shock remain pending; the continuation after 2015 is not a fitted two-shock history. The plots perform zero model solves. The saved inputs were copied read-only from `/scratch/td2248/projects/current_estate_parallel_candidates_20261004_v1/results/c03_h24/run/`: `candidate/latest_completed_full.json` to `c03_h24_native_dated.json`, and `latest_completed.json` to `c03_h24_latest_completed.json`. The plotted US values match the four-row target contract recorded in the parent `README.md` and are saved in `historical_fertility_targets.json`.

Regenerate in seconds with the project plotting environment:

```sh
output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python \
  code/model/tools/build_current_estate_saved_transition_panels.py \
  --native output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/parallel_candidates_v1/plot_readout/c03_h24_native_dated.json \
  --receipt output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/parallel_candidates_v1/plot_readout/c03_h24_latest_completed.json \
  --targets output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/parallel_candidates_v1/plot_readout/historical_fertility_targets.json \
  --output output/model/transition_readiness_v1/current_baseline_20261003/two_shock_v1/parallel_candidates_v1/plot_readout
```

The plotter accepts the current-estate saved native dated JSON and its matching completion receipt through `--native`, `--receipt` and `--output`. It verifies row identity, four-year dates, finite plotted series, positive adult population and root/replay gates. It requires a constant dated preference for a first-shock receipt; a receipt explicitly marked `two_shock_result: true` may have a dated preference path. This is the supported receipt contract, not a claim that arbitrary native JSON can be plotted. `--allow-incomplete` permits an explicitly labeled diagnostic from an incomplete saved path. `--targets` adds observed values when a matching contract is available.
