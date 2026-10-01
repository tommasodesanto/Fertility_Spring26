# Local floor Nelder–Mead continuation

Experimental searches, with no production adoption. Reuse the verified round2
native evaluator and its read-only local source overlay. The inherited economy
uses 120 wealth nodes and 9 income states, nonnegative mean-preserving five-bin
entrant wealth, zero unsecured credit, 2% annual interest, fixed H0 and psi,
and the floor utility already documented in round2. Target weights are unchanged:
14 reported moments, including 10 scored moments, for 9 free coordinates.

`run_nm.py` applies SciPy bounded adaptive Nelder–Mead to coordinates scaled to
[0,1]. The simplex uses beta steps of .002 and floor steps of .1; other steps
are 10% of the starting parameter or explicit physical steps recorded in each
`search_contract.json`. Six starts reuse round2 indexes 0,1,5,6,3,7: winner,
historical/curvature, low floor, high floor, lower cost, and zero birth cost.

Each chain has a four-hour hard deadline, 500 objective-call maximum, two threads,
four GiB RSS cap, exact-vector caching and a final fifteen-minute reserve. The
reserve accommodates the original native path's unchanged reporting guards.
Six chains
have at most 3000 search GEs; this is a cap, not a promise of completed work.
The old native evaluator initially retains its intrinsic selected-price repeat
and full report for each objective. There is no new baseline/repeat gate.
An optional reviewed `fast_objective.make_evaluator` can omit intermediate plots
and the intrinsic selected repeat while retaining root/accounting/parameter gates.
Final selected verification always uses the original native reporting path in
a separate fresh supervised process, after the search process exits. This
avoids repeating frozen-runtime authentication in one process.
Optimizer success alone never certifies convergence or economic identification.

Every computed objective writes `latest_completed.json`, `best_so_far.json`,
`cases.json`, all 14 target rows and all 31 effective parameter rows. Explicit
uncomputed numerical-root rejections receive a 1e12 penalty and a receipt;
accounting/source/feasibility failures stop the chain. A numerical rejection is
never presented as a valid model loss. New task-owned intermediate array/cache
payloads are removed after parsing; only current-best standard plots are kept.
All old packets are preserved. Free disk below 20 GiB stops the new chains.
Supervisors refresh progress every 60 seconds and enforce caps every 3 seconds.

Preparation check: run the real SciPy loop with the toy evaluator (zero model
solves). Explicit author steering skips additional model smoke gates. Launch
only after lead review:

```sh
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/utility_floor_nm_v1/launch.py
```

The launcher returns immediately after starting six detached supervisors. Their
launch/terminal/watchdog receipts and native logs live under `chain_0` through
`chain_5`. Any selected postcheck that cannot fit its remaining budget is reported
as provisional; no automatic restart or extension is performed.
