# Fixed-credit runtime validation v1 — prepared, not submitted

**Incomplete and unsubmitted; do not run.** The lead rejected this partial
controller adaptation. A smaller harness is being prepared in sibling
`runtime_validation_v2/`; only its reviewed launch receipt can authorize a run.
This skeleton and its source evidence are preserved, with zero model solves.
Its affordability JSON is also unaccepted: `income.z1` has a transcription error
(1.209 instead of 0.209), and purchase cash omitted division of income by
`R_gross`. Use the actual-function checkpoint diagnostic in v2, not this JSON.

Reference: **2007 stationary reference — block0506, September 28 verified export**.

The only proposed economic change is the explicitly selected scalar
`unsecured_credit_limit=0.0`. Prices, rents, fiscal objects, entry mass,
preferences, fertility normalization, and the 160-node grid remain frozen.
The source overlay is pinned module-by-module: `parameters.py` `d7b4d23c…`,
`solver.py` `cefc1627…`, and `kernels.py` `379d179a…`.

The two age-18 renter cells at wealth index 44 are infeasible at zero credit:
cash before any rent, consumption, or saving is -0.134050 and -0.067770 for
income states z0 and z1. Purchase cash before a downpayment is also negative.
With zero transfers and no estate/birth grant, no attainable tenure or
fertility choice restores feasibility. Therefore the strict case is a
zero-lifecycle diagnostic, not a model solution.

The only conditional runtime action left is one fresh-process control with the
overlay loaded but `unsecured_credit_limit=None`. It must reproduce retained
native policy, birth, and distribution arrays exactly and pass all gates before
any future strict run could be considered. No strict solve, GE, price sweep,
positive credit limit, entry mutation, or retry is authorized here.

Before any submit, pin `driver_sha256` in `plan.json`, stage only this directory
and the existing `overlay/` directory to the stated paths, then run:

```sh
sbatch output/model/fixed_reference_economics_20260928/credit_no_taper_v1/fixed_credit_contract_v1/runtime_validation_v1/launch_runtime_validation.sh --self-test-only
```

Production submission is deliberately not prepared as an executable command
until lead review pins the staged source and resolves the controller's
single-case adaptation.
