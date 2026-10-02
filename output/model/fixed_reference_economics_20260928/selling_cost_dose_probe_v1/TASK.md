# Selling-cost dose-response + lock-in check (experimental; NOT a calibration; nothing adopted)

Repo: /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26. Write ONLY inside output/model/fixed_reference_economics_20260928/selling_cost_dose_probe_v1/. Edit no repo source. Local, one thread (NUMBA/OMP/OPENBLAS/MKL_NUM_THREADS=1), code/model/.venv/bin/python. (Torch SSH is failing; do not try it.) Budget 75 min. Stop when tables written or on blocker (write RECEIPT.md saying what blocked).

Start from output/model/fixed_reference_economics_20260928/tenure_barrier_probe_v1/ (run_barrier_probe.py, params_used.json; arms B1 and B2 verified). Same fixed price/H0, per-child need "need1" (hbar_child_rooms 1.0, hbar_first_child_jump 1.49169824815624).

## Part 1: dose response, 8 solves
selling cost P.psi in {0.00, 0.02, 0.04, 0.06} x financed share {0.8, 1.0}. psi=0.06/0.0 arms must reproduce B1/B2 first_birth_flow exactly (identity gate).
Rows: first/second/third birth flows; completed fertility; childlessness 40-44; ownership_30_55; own_rate_25_34; first_birth_rooms; renter-parent cap share m=1,2,3; matched entrant (age 18-21) first-birth prob and ownership. Financing effect (1.0 minus 0.8) per psi.

## Part 2: lock-in diagnostic from saved arrays of psi 0.06 and 0.00, both phis (zero extra solves)
Among CHILDLESS households (n=0, m=0) by age cell 18-37 (age cells 0-4): mass by beginning tenure state (renter, owner of each rung), and the first-birth attempt/birth probability by beginning tenure state (mass-weighted). Also for childless owners who have a first birth: share that keep the same rung, upsize, downsize, switch to renting (from saved tenure policies after the birth branch if available; if not saved say so).
Report as tables; do not interpret.

Deliverables: results.csv, part2 tables, RECEIPT.md (Outcome / Verification / Artifacts / Unresolved / Reported cost).
