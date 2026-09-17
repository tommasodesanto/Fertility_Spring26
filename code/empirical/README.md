## Empirical Layout

Updated: `2026-04-21`

This folder holds the non-core empirical and data-side scripts that do not need
to live in the top-level calibration surface.

### Subfolders

- `housing/`
  ACS room-distribution audit and model-vs-data comparison scripts used by the
  live housing-fit pipeline.
- `acs/`
  ACS supporting regressions and plots that feed empirical notes and side facts.
- `roundup/`
  PSID / empirical-roundup scripts and compilers.

### ACS fertility identification feasibility

The September 17 design review is in
[`acs_fertility_pseudopanel_feasibility.md`](../../docs/model/acs_fertility_pseudopanel_feasibility.md),
with a verified [sample and variable inventory](../../docs/model/acs_fertility_pseudopanel_inventory.md).
It separates pooled fertility instruments from cohort validation; no new estimation
or calibration target is approved by the memo.

### Historical Calibration Plotting

The old MATLAB plotting path is archived, not active model code.

- archive:
  `calibration_archive/model_history_2026-05-07/legacy_matlab_2026-05-07/`
