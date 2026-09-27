# Final overnight readout - September 27, 2026

The two-page author memo is `output/pdf/fertility_overnight_memo_20260927.pdf`
from the repository root. The identical build copy is `memo.pdf` here.
This is a reproducible provisional calibration point, not an established
optimum or a completed economic calibration.

## Evidence and selection

All numerical work finished by 07:02 EDT. There were 436 successful search
evaluations: 240 original main, 96 main continuation, 60 unit-weight diagnostic,
and 40 early-fertility-weight diagnostic evaluations. Sixteen smoke/repeat
checks appear in `candidates.csv`; two additional cold-start checks are kept
outside that comparison. All four searches passed two final exact repetitions
of all 14 target rows and 31 full parameter rows. No accepted-run numerical
failure, timeout, inadmissible result or resource intervention was recorded.

Actual capacity was up to ten local evaluators and zero cluster search workers.
Expired SSH access blocked the requested 24 cluster workers. Two earlier local
preparation smokes were cancelled after pin drift; a zero-solve cluster
submission failed. These are preserved and excluded from accepted-run counts.
The legacy Torch replay and cross-host checks remain uncollected.

The main selection is continuation `de_0093`, primary loss 42.281937042645964.
The main model has no housing floor and delivers a first-birth rooms response
of 1.610 against 1.465, but ownership, wealth and children born by age 25 remain
below target. No economic specification, target, grid or scientific gate changed
during search. Alternative weights are diagnostics, not adopted specifications.

## Files

- `target_fit.csv`: all 14 restrictions, including the unscored fertility
  normalization, with targets, model values, gaps, primary weights and losses.
- `parameters.csv`: all ten fitted parameters, bounds and bound flags; nine
  are searched and the positive child-benefit level is normalized internally.
- `weighting_comparison.csv`, `weighting_parameters.csv`, `weighting_scores.csv`:
  each weight system's own selected point, rescored with common primary weights.
  The early-weight winner scores 57.608; re-ranking all its saved candidates
  would instead select another point at 57.444. These are different selections.
- `target_provenance.json`: provenance retained from the frozen target contract.
- `candidates.csv`, `summary.json`, `monitor.json`: full local search audit.
- `report_qa.json`: final PDF, table and layout checks.
- `../primary_final/`: selected receipt, estate ledger, exact-repeat verification,
  full 31-row parameter table, and the unchanged 17 standard diagnostic PNGs.
- `../early_weight_final/`, `../identity_final/`, `../primary_initial_final/`:
  independent final verification and plot evidence for the other searches.
- `../fertility_overnight_support_20260927.zip`: portable main diagnostic packet,
  all fit/parameter/weight comparison tables, provenance and verification.

The weighted-gap SVG is supplemental. It does not replace any standard graph.
Rendered `page-*.png` files are layout checks, not additional report pages.
Two complete PDF pages and the main and alternative final plot packets were
visually inspected. Wide wealth axes limit detailed boundary inspection.

## Regeneration

From the repository root, rebuild the core memo and tables into a fresh output
directory using the reportlab-enabled bundled Python:

```sh
/Users/tommasodesanto/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  code/model/tools/build_e5f_overnight_memo.py \
  --index output/model/supervised_calibration_20260927/report_index.json \
  --output tmp/pdfs/fertility_overnight_rebuild
```

The following command regenerates all 17 standard plots from authenticated
saved arrays without a model solve. The output directory must not exist:

```sh
EXPECTED_UTILITY_OVERNIGHT_SHA256=3b770d8c8c22d2b0449b34a575d6353b063bc015d74ce11016dad7e22ed7ca5e \
E5F_LOCAL_EXECUTION_AUTHORIZATION=tommaso_authorized_20260927_local_primary_continuation_v1 \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 MPLBACKEND=Agg \
code/model/.venv/bin/python code/model/tools/plot_e5f_calibration_case.py \
  --contract tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/production_contract.json \
  --case tmp/e5f_overnight_local_20260927/portable/night_launch_v4/primary_continuation/search/de_0093/case \
  --output tmp/pdfs/fertility_overnight_diagnostics_rebuild
```

Preserve the immutable portable source, contracts and checkpoints. Neither
regeneration command authorizes another calibration search.

## Next decisions

Check age-25 fertility measurement and child-benefit curvature jointly with the
two fertility taste scales, retaining every fertility target. Review the credit
constraint, 2% interest rate, tenure scale, rental choices and income/wealth
inputs before a wider search. Keep the finer housing grid and conception
schedule as separate tests. Estate/SCF recipient scope, older-wealth income
denominators and the housing event-study observer remain explicit approximations.
