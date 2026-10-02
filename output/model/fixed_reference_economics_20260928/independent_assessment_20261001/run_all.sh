#!/bin/sh
# Regenerates every tabulation in this folder. Read-only on saved outputs; zero model solves. Run from the repository root.
set -e
D=output/model/fixed_reference_economics_20260928/independent_assessment_20261001
PY=code/model/.venv/bin/python
for s in 01_search_cloud_jacobian 02_purchase_tests_saved_arrays 03_early_fertility_bound_and_selection 04_financing_arms_target_table; do
  $PY $D/$s.py > $D/$s.out.txt 2>&1
done
