# Per-child housing need probe (experimental, mechanism test only; NOT a calibration)

Repo: /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
Write ONLY inside this folder: output/model/fixed_reference_economics_20260928/per_child_need_probe_v1/
Do not edit any file outside it. Python: code/model/.venv/bin/python. Set NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1. Local, one core.
Time limit 75 minutes total. Stop condition: four solves done and table written, OR any blocker -> stop and write RECEIPT.md explaining exactly what blocked.

## Reference point (do not change any other parameter)
Quarter-saving purchase rule, selected overnight point (chain 54, loss 51.556036, price 0.6744838540900874, H0 7.288573389887633).
Arrays/receipts: output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/local_runtime/runs/local10_v1/chain54/postcheck/selected_postcheck/phase_b_ge/selected_repeat/stage/ (summary.json, normalization.json, solution_arrays.npz) and chain54/search, chain54/postcheck JSONs.
Parameter table: latex/corina_progress_20260930/evidence/overnight_quarter_parameters.csv.
Engine: code/model/experiments/quarter_saving_solvency/source/small_credit_lab/engine/.
Existing fixed-price drivers to copy patterns from: output/model/fixed_reference_economics_20260928/run_fixed_price.py, quarter_saving_solvency_v2/, purchase_rule_comparison_v1/ (each case packet has a run.py that sets P.phi before the solve).

## Step 0 (gate)
Reproduce the reference point at FIXED price 0.6744838540900874 (no GE root) with financed share 0.8. Confirm first_birth_flow and ownership_30_55 match the selected point's saved values to ~1e-6 (or explain the difference). If you cannot build the exact parameter object, STOP.

## Experiment: 2x2 at fixed price, otherwise identical
Child floor in engine/shared.py:286-290 is h_bar = hbar_first_child_jump + hbar_child_rooms * nk. Report what nk is (children ever born or at home) and the current values of both parameters (expected: jump ~2.4917, child_rooms 0).
Arms:
- need0: current values.
- need1: hbar_child_rooms = 1.0 and hbar_first_child_jump = (current floor at one child) - 1.0, so the ONE-child floor is unchanged, two children +1 room, three +2.
Crossed with financed share (P.phi / financed_share) in {0.8, 1.0}. Price, H0 and every other parameter fixed. No recalibration.

## Report (RECEIPT.md, concise, plus results.csv)
For each of the 4 cases: first_birth_flow; second- and third-birth flows (from parity_birth_flows_by_age or equivalent observer); completed fertility (initial_normalization / TFR); childlessness 40-44; ownership_30_55; own_rate_25_34; first_birth_rooms; mean_rooms; share of renter parents at the 6-room cap by number of children.
Matched-entrant cell age 18-21 (same entrants in every case): first-birth probability and ownership.
Then the financing effect (phi 1.0 minus 0.8) within need0 and within need1, for every row above, and the difference-in-differences.
List the exact parameter objects used, solve times, and any warnings/convergence flags. Do not interpret economics.
