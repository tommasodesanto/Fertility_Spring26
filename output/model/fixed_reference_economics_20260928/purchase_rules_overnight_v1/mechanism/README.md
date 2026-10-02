# Dated financed-share mechanism, two fitted purchase rules

This isolated driver starts each experiment from that rule's **fresh postchecked**
80% financing equilibrium. It requires the read-only collector's
`collection/readout/selected_hard.json` or `selected_quarter.json`, the matching
`results/chain_N/postcheck/completed.json`, both native selected reports, and
the selected repeat state arrays. It refuses a search best-so-far point. The
source and target pins are checked before reconstructing and exactly repeating
the chosen 2007 baseline. The 120×9 grid, mean-preserving nonnegative entrant
distribution, four-year timing, earnings, housing floors, transfers, taxes and
all ten fitted parameters stay within the chosen arm's contract.

The sole dated policy change is the owner financed share from 0.80 to 1.00.
`control` holds 0.80; `temporary` uses 1.00 at the first date and 0.80
thereafter; `permanent` holds 1.00 at every date. `dated_phi.py` applies the
date's share to both the backward household solve and the forward policy
reproduction. The original exact policy cache hashes the complete parameter
object, including the date's financed share. `integration.py` installs the
isolated dated function only during the original audited mapping, retaining
its household budget, purchase, estate, policy-array, population and queue
checks. The child preference remains at its fitted 80% value.

The calibrated housing-supply coefficient $H_0$ stays fixed. Population and
both birth-to-entry queues evolve endogenously. The permanent case uses the
existing closed-demography terminal stationarity condition: price clears
birth renewal and the terminal population scale equals physical supply divided
by per-household demand. This pins **long-run replacement fertility**, while
dated first-birth flows and hazards remain free to respond. A failed terminal
root, dated price/pension root, exact replay or terminal state check is a failed
case, never a calibrated result.

Run one arm, one case and one horizon per process; use both 12 and 16 dates:

Before submitting the dated array, run one authenticated, native one-date
control smoke per arm with the same arguments and engine, adding
`--kind control --horizon 1 --smoke-one-date`. Its `completed.json` must say
`status: passed`, `kind: smoke_one_date`, and `horizon: 1`. This runs a fresh
baseline reconstruction and the actual dated backward/forward callbacks;
it is a launch gate, not a reported policy result.

```sh
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/mechanism/run_case.py \
  --arm hard --kind temporary --horizon 12 \
  --selected-json output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/collection/readout/selected_hard.json \
  --selected-completed output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/results/chain_N/postcheck/completed.json \
  --out output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/mechanism/results/hard_temporary_h12 \
  --deadline-epoch DEADLINE_EPOCH --maximum-policy-calls CALL_CAP
```

`DEADLINE_EPOCH` and `CALL_CAP` must be finite launch-plan values. Use one
thread for Numba, BLAS and OpenMP, and the same source tree on Torch. The
driver writes a run contract, per-mapping heartbeat, latest completed case,
best-so-far case, native dated audits and a terminal or failure receipt. A
failed case must retain its folder for diagnosis. Standard 17 plots remain in
the selected baseline ROOT/REPEAT reports; date 0, midpoint and final-date
diagnostic packets are saved by the native observer.

The first-birth response uses the native dated `birth_flow_first`, at-risk
childless mass and `first_birth_hazard` by age. For each arm and horizon,
compare temporary and permanent date-0 and later values with that arm's
control. Young renter-to-owner closing finance diagnostics are prepared in
`../buyer_diagnostics/`; do not infer financing inability from remaining a
renter. The quarter-rule inability measure needs branch-specific saving
feasibility and is not produced by this driver.

Zero-solve verification:

```sh
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/mechanism/test_zero_solve.py
```

This packet is experimental; no baseline or policy result is adopted.
