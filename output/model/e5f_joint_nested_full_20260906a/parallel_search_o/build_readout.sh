#!/bin/bash
set -euo pipefail
# Optional first argument is a NEW output path; the builder never overwrites.
export PYTHONPATH="/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/pdfs/plot_dependencies${PYTHONPATH:+:$PYTHONPATH}"
"/Users/tommasodesanto/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3" "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/build_e5f_joint_nested_review.py" \
  --selected-dir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/support_repair_m/smoke/smoke_histories/task_004" \
  --reference-fit "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_experiment_20260906a/reference_target_fits.csv" \
  --reference-parameters "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_experiment_20260906a/reference_parameters.csv" \
  --policy-results "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/parallel_search_o/output/model/joint_nested_overnight/equilibrium_path" \
  --allow-partial-policies \
  --search-verification "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/parallel_search_o/output/model/joint_nested_overnight/search/final_verification.json" \
  --narrative "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/parallel_search_o/morning_narrative_final.json" \
  --market-diagnostic-figure "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/parallel_search_o/market_trace/results/price_trace.png" \
  --display-diagnostics "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_joint_nested_full_20260906a/parallel_search_o/display_graphs_readout_final" \
  --fixture-label 'EXPERIMENTAL - RETAINED STARTING POINT' \
  --output "${1:-/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/joint_nested_review_20260907_final.pdf}"
