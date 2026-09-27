# Current production transition integration

The author requested all accepted recent changes in the production transition and a same-parameter credit experiment starting from the exact calibrated steady state. Economic changes remain explicit: baseline uses current calibrated credit; the experimental option replaces artificial borrowing limits with grid-supported feasible continuation and nonnegative net estates. All preferences, earnings, entry wealth, grids and targets remain fixed. No fertility renormalization or initial-population reset.

## Verified native integration

The current native solver includes conditional entrant wealth/income coupling, the exact transaction grid, income-aware purchase screens, exact allocation reporting, and the optional experimental solvency-credit mode. Current preferences, B15 earnings, survival, fiscal and supply inputs are read from the authenticated calibrated checkpoint. PF population accounting uses the adopted split 16/20-year entry queue and a dated estate ledger funding the actual next entrant cohort. See `preparation/PRODUCTION_INTEGRATION.md` for the complete 14-object reconciliation with September 14.

- `native_baseline_v1/`: free-price baseline replay; all 14 target gaps below 2.331e-6 and all 31 parameters exact. Its price moved 4.336e-7 inside the root tolerance; retained as evidence.
- `native_baseline_fixed_price_v1/`: fixed original price isolates implementation changes. ALL 107 saved arrays exact, including values, policies, distributions, maps and shared objects; all targets/parameters exact. Solve 20.647 seconds.
- `native_endpoint_v1/`: native experimental endpoint reproduces ALL 107 arrays and all 14/31 rows exactly. Solve 22.392 seconds.
- `native_baseline_fixed_price_v2/` and `native_endpoint_v2/`: final source receipts after adding actual dated rents to the existing plot routine. Again all 107 arrays exact; all 17 standard PNGs byte-identical to the visually inspected v1 packets. Solves 21.072 and 21.966 seconds. No economic or numerical gate changes.
- `baseline_approval_v2.json`, `terminal_approval_v2.json`: lead-reviewed numerical/source/plot evidence. These approve native reproduction, not a full transition.

Complete target/weight/loss tables and all 31 parameters with restrictions/bounds live in each `case/target_fit.csv` and `case/parameters.csv`. Baseline loss 42.282; endpoint loss 57.193. The endpoint is a counterfactual, not a better calibration. The calibrated initial distribution and the terminal normalized distribution are different objects: multiply only the latter by its recorded population scale for terminal-distance comparisons.

## Dated verification and root entrypoint

`native_smoke_v1/` is the bounded two-date native mapping test. Reference prices are constant; experimental prices genuinely change between the calibrated and closed-endpoint prices. It checks actual next-entry estate funding, household budgets, feasibility, population accounting and exact cache equality. Only the reference no-shock arm is market/fiscal-gated. Experimental residuals are reported because a prescribed price/pension path is not an equilibrium solution. The plan has one worker, at most 20 minutes total and 10 minutes per arm, a 2 GiB cache and immutable source snapshots. Both arms PASSED. Every cache-on/cache-off output is exact. Reference market residual is 7.445e-7 and fiscal residual 7.512e-13. Changing-price credit arm takes 51.828 seconds without cache and 27.749 with cache (two actual solves, two hits); housing residual 0.030187 is explicitly not equilibrium. No negative estates or funding/budget violations. The whole prescribed-path test is a verification of the mapping, not a solved price path.

The prior `smoke_v3/` is a completed frozen-runtime prototype, not native certification. Earlier `smoke_v1/` and `smoke_v2/` failures remain preserved: tiny saved frontier projection and an entry-prehistory mismatch, respectively. Neither was repaired by changing scientific gates.

`code/model/tools/run_e5f_current_transition.py` implements the joint price/pension root at fixed payroll tax. It requires a pinned plan, approved native replays and exact-loop smoke, explicit horizon/evaluation/time budgets and numerical starting paths. It saves latest/best/final mappings and separate root/terminal gates. A converged final mapping exports the standard 17 plots at first/middle/final dates using actual dated rents. Long-horizon and horizon-extension runs remain outstanding; the root driver has not yet been numerically certified. No full equilibrium transition has launched.

The initial stationary fertility normalization has a tiny birth/entry mismatch. Queue prehistory reproduces actual entry; future actual births enter the queue without adjustment. Longer no-shock paths can therefore drift slightly after the four/five-date maturation lag. Do not silently recalibrate preferences, reset population or relabel that drift as exact stationarity.

## Reproduction

Use a fresh interpreter, the current project virtual environment, and one thread for each numerical library. Baseline command is `e5f_current_transition_runtime.py --output NEW_DIRECTORY --replay --fixed-reference-price --seconds 300`; endpoint is `run_e5f_current_endpoint_replay.py --output NEW_DIRECTORY --seconds 300`. The latter is a fixed-price replication of the authenticated endpoint, not another root search. Never overwrite completed directories. Source and checkpoint pins are validated by the drivers; source snapshots are retained beside each run.

These drivers still use the authenticated saved contract and observer provenance under the local portable reference bundle. Cluster packaging must preserve that dependency closure and adapt paths with explicit authentication; copying only the new driver is insufficient. No cluster portability claim has been made for this new runtime.
