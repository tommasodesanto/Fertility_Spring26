# Tenure-barrier probe (experimental mechanism test; NOT a calibration; nothing adopted)

Repo: /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26
Write ONLY inside output/model/fixed_reference_economics_20260928/tenure_barrier_probe_v1/ locally, and in a NEW dated folder under the Torch scratch project tree. Edit no repo source file. Python locally: code/model/.venv/bin/python, one thread (NUMBA/OMP/OPENBLAS/MKL_NUM_THREADS=1).
Cluster: read docs/workflow/delegation_and_cluster_playbook.md (section on Torch; `code/cluster/torch.sh`, ssh alias `torch`, user td2248). Check for existing jobs before submitting. Smoke-test first. If SSH fails, run the solves locally instead (each solve is ~7 s) and say so.
Budget: 2.5 hours wall total; Torch jobs <= 1 hour walltime, one small job or array (<= 16 tasks, 1 core each). Stop condition: Part A + Part B tables written, or a blocker -> stop and write RECEIPT.md stating exactly what blocked.

## Starting point (reuse, do not rebuild from scratch)
Previous probe, verified bit-identical to the saved overnight quarter point: output/model/fixed_reference_economics_20260928/per_child_need_probe_v1/ (run_probe.py, params_used.json, RECEIPT.md, NOTE.md, case folders need*_phi*/ with stage/solution_arrays.npz). Copy and adapt run_probe.py. Fixed price 0.6744838540900874, fixed H0 7.288573389887633, no GE root, no recalibration.
"need1" = hbar_child_rooms 1.0, hbar_first_child_jump 1.49169824815624 (one-child floor unchanged). "need0" = reference (2.49169824815624, 0.0).

## Part A — zero-solve decomposition on saved arrays (local)
Use need1_phi10 saved arrays (and need1_phi08 for comparison). Population: renter parents (children at home m>=1) at the 6-room cap, weighted by mass (pre-tenure weights as in the previous probe: g_post_fertility x tenure_probs). For this population report, by m and age cell:
1. Share for whom buying each owner rung (H_own values) is FEASIBLE under the purchase screen and borrowing floors (code: engine kernels.py ~380-420, 880-920; household.py ~890-915). Cite lines.
2. If the arrays allow: tenure choice probabilities over {rent, each rung}, and the choice-specific values (rent vs best owner rung) — mass-weighted mean gap and its distribution. If choice-specific values are not saved, say so and skip.
3. Mean net financial position b, income state, and whether they are near the wealth-grid bounds.

## Part B — counterfactual solves (Torch preferred), 12 solves
All at fixed price/H0 above; each arm at financed share 0.8 and 1.0:
- B1 need1 control (must reproduce previous need1_phi08 / need1_phi10 first_birth_flow to ~1e-10: identity gate).
- B2 need1 + selling_cost = 0.
- B3 need1 + finer owner grid H_own = (2,3,4,5,6,7,8,9,10). FIRST verify the engine supports an owner grid of any length (search for hard-coded rung counts, n_tenure, array shapes). If it does not, skip B3/B4/B6 and report why.
- B4 need1 + finer grid + selling_cost = 0.
- B5 need1 + tenure_choice_kappa = 0.001 (lower bound of its search range).
- B6 need0 + finer owner grid.
Report for every solve: first/second/third birth flows; completed fertility; childlessness 40-44; ownership_30_55; own_rate_25_34; first_birth_rooms; mean_rooms; renter-parent cap share by m=1,2,3; ownership of parents by m; matched entrant (age 18-21) first-birth probability and ownership; fp.gates; residuals; solve time.
Then the financing effect (phi 1.0 minus 0.8) for every row within each arm, and each arm's financing effect minus B1's.

## Deliverables
results.csv, partA_tables (csv/md), RECEIPT.md (Outcome / Verification / Artifacts-job IDs / Unresolved / Reported cost). Do not interpret economics.
