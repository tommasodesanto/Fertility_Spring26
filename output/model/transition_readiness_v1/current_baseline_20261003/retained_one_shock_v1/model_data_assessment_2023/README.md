# Routine model–data comparison: saved 2023 transition

Uses exactly the standard assessment's 19 panels and three-page layout.
Run from the repository root:

```sh
output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python code/model/tools/plot_transition_model_data.py
```

The command reads the saved full 2023 household solution and survey caches. It does not estimate a transition or clear an equilibrium. The first preparation recovered one date's policies from the original frozen engine and saved continuation values; reproduction is checked at floating-point precision. Future runs load `assessment_solution.pkl.gz`.

- `assessment_solution.pkl.gz`: 2023 parameters, values, consumption, housing, savings, tenure/location/fertility policies, current and beginning distributions, prebirth distribution and first-birth flows.
- `../state_2023/actual_2023.pkl.gz`: original inherited state, both entry queues, continuation/terminal values and price/pension/preference forecasts.
- `state_receipt.json`: original checkpoint hash and policy reproduction check.
- `model_data_assessment.pdf`, `page_1.png`–`page_3.png`: standard plots; `plotted_long.csv` and `metadata.json`: plotted numbers and source details.

Data dates are explicit: CPS June 2024 (nearest local fertility supplement), NCHS 2023, ACS 2023, and the latest PSID wave with measured age, weight and wealth fields (2019). The local shelf has 2021 records but those required fields are missing. This is descriptive validation, not a change to calibration targets.

This routine uses only the 2023 cross-section. Original experimental certification limitations remain in the parent README.

Snapshot check: distributions come from `actual_2023.initial_state.g_pre`, not the stationary observer arrays embedded in parameters. First births by age, population, ownership and housing demand reproduce the original dated record. The plot loader requires this receipt. An earlier plotting cache used stale stationary distributions; it has been replaced.
