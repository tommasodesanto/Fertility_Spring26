# Owner-ladder floor + split transaction cost probe (experimental; NOT a calibration; nothing adopted)

Repo: /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26. Write ONLY inside output/model/fixed_reference_economics_20260928/owner_floor_split_cost_probe_v1/ (locally) and a new dated folder on Torch scratch if you use Torch. Do NOT edit any repo source file; any engine modification happens in a COPY of the engine inside this folder.
Compute: Torch is reachable again (code/cluster/torch.sh; docs/workflow/delegation_and_cluster_playbook.md; check existing jobs first; array 19040483 is someone else's production run — do not touch). Solves take 7-14 s, so a single small Torch job (<= 1 h, <= 16 tasks) or local one-thread runs (code/model/.venv/bin/python; NUMBA/OMP/OPENBLAS/MKL_NUM_THREADS=1) are both acceptable; say which you used. Budget 2.5 h wall. Stop when deliverables are written or on blocker (RECEIPT.md says what blocked).

Start from output/model/fixed_reference_economics_20260928/tenure_barrier_probe_v1/run_barrier_probe.py (verified). Fixed price 0.6744838540900874, fixed H0 7.288573389887633, no GE root, no recalibration. need0 = reference floor (hbar_first_child_jump 2.49169824815624, hbar_child_rooms 0); need1 = (1.49169824815624, 1.0).

## Part 0 (lookup, 10 min max)
Find the DUE paper PDF (Greaney, Parkhomenko, Van Nieuwerburgh 2025, "Dynamic Urban Economics"; search ~/Desktop, ~/Downloads, repo latex/ and docs/ for the PDF). Quote (<= 2 short sentences) how its transaction cost psi = 0.06 is charged: on sale, on purchase, or split; on the whole house value or otherwise; with page/equation number. If not found, say so. Do not guess.

## Part A: owner-ladder floor, 12 solves
Owner grid H_own variants x financed share {0.8, 1.0} x need {need0, need1}:
- A1 (2,4,6,8,10) control. need0 arms must reproduce per_child_need_probe_v1 need0_phi08/need0_phi10 first_birth_flow exactly; need1 arms must reproduce B1.
- A2 (4,6,8,10)
- A3 (6,8,10)
Selling cost unchanged (0.06, seller pays).

## Part B: split transaction cost (engine copy; 4 solves + 2 identity checks)
Copy the executed quarter engine into ./engine_split/. Add a buyer cost psi_buy charged as psi_buy * p * H_new at every purchase (renter->owner and owner->different rung), alongside the existing seller cost psi * p * H_old on every sale. It MUST be applied identically in (i) the household Bellman budget/feasibility (incl. the purchase screen and borrowing floors that reference cash at closing) and (ii) the KFE/distribution wealth transition and any wealth statistics. List every edited line (file:line, before/after) in RECEIPT.md.
Identity checks: psi_buy = 0 must reproduce A1 bit-for-bit (need0, both phis).
Runs: psi = 0.03, psi_buy = 0.03, need0, phi {0.8, 1.0}; and the same with need1.

## Rows for every solve
first/second/third birth flows; completed fertility; childlessness 40-44; ownership_30_55; own_rate_25_34; first_birth_rooms; mean_rooms; renter-parent cap share m=1,2,3; matched entrant (18-21) first-birth prob and ownership; childless owners aged 22-29 by rung (mass shares); target_fit loss if computed; fp.gates; residuals; solve time.
Financing effect (phi 1.0 minus 0.8) for each arm, and each arm's effect minus its control.

Deliverables: results.csv, RECEIPT.md (Outcome / Verification / Artifacts-job IDs / Unresolved / Reported cost), the engine diff as engine_split.diff. Do not interpret economics.
