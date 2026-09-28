## Conclusion

The existing runtime is **not drop-in compatible** with the September 28 verified export. It hard-pins the September 27 `de_0093` contract/checkpoint, old frozen observers, static-elastic supply, and an experimental credit mode that conflicts with block0506’s saved DUE repayment flag.

The new reference must be the actual September 28 verified-export checkpoint, SHA `b15ba92d…`, not the original search-case checkpoint SHA `1c4345b…`. The manifest distinguishes them explicitly. [fixed_reference_manifest.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/fixed_reference_manifest.json:1432) [fixed_reference_manifest.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/fixed_reference_manifest.json:2004)

## Required adapter surface

1. Replace hard-coded reference identity with an explicit immutable `ReferenceSpec`.

   Current constants point only to the old portable contract and checkpoint (`3b770d…`, `090c9e…`). [e5f_current_transition_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_current_transition_runtime.py:27)

   The new spec must pin, before unpickling:

   - September 28 contract `68323a…`;
   - source manifest `07d843…`;
   - verified-export checkpoint `b15ba…`, whose physical Torch path is recorded in the manifest;
   - selected-export parameter, target, objective, and observer artifacts.

   These identities are documented in the measurement audit. [measurement_audit_v1/README.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/measurement_audit_v1/README.md:18)

2. Preserve the serialized parameter object wholesale.

   The existing setup does the right basic operation—unpickle first, then `copy.deepcopy(selected['parameters'])`—and checks the serialized wealth and entry grids. [e5f_current_transition_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_current_transition_runtime.py:242)

   Retain that pattern. The adapter may add only a logged, allow-listed transition overlay. It must fingerprint/check `Pi_z`, `z_grid`, survival, pension/tax, entry arrays, supply objects, credit flags, and grids after every overlay. The new manifest explicitly authenticates every serialized instance field and says no constructor defaults were substituted. [fixed_reference_manifest.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/fixed_reference_manifest.json:1705)

3. Replace the old frozen measurement bundle as one overlay closure.

   The runtime currently loads fertility/housing/accounting/estate tools from the September 27 portable source root, while allowing only narrow byte-verified changes to the current reporter and recent-parent observer. [e5f_current_transition_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_current_transition_runtime.py:231) [e5f_current_transition_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_current_transition_runtime.py:249)

   That is not sufficient for a September 28 measurement-preserving claim: the runtime’s frozen ancestry is still `3b770d…`, whereas block0506’s verified source manifest is `07d843…`. Use a new, isolated import namespace/source snapshot from the September 28 manifest; reject mixed local/frozen modules. Do not “update” observers individually or inherit September 27 reporting exceptions.

4. Make the housing closure explicit.

   The current path passes `packet['supply_rule']` unchanged. [run_e5f_current_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_current_transition.py:193) That saved rule is static-elastic in the old transition workflow. The dated PF framework already supports `HousingSupplyRule("fixed-stock", …)`, returning the fixed initial stock independently of price. [run_dynamic_population_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_dynamic_population_transition.py:152)

   Build a transition-only fixed-stock rule from the authenticated date-0 stock and retain it in every packet/receipt. Prices still clear demand against that constant quantity. Do not modify `H0` or `xi_supply`; the existing reporter’s population scaling is display-only, not an economic fixed-stock closure. [e5f_current_transition_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_current_transition_runtime.py:139)

## Dated household/KFE closure hazards

- The existing initial state uses saved pre-choice mass but reconstructs the entry queue from current entrant mass, saved births, and hard-coded \(1/2.1\). [run_e5f_current_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_current_transition.py:161) The PF implementation deliberately seeds adjusted prehistory with observed entry even if the stationary birth/entry gap is nonzero. [run_e5f_perfect_foresight_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_perfect_foresight_transition.py:100)

  This is compatible only if the adapter records the saved block0506 entry gate and queue convention. The current README warns that longer no-shock paths can drift because future births enter without adjustment. [credit_transition/README.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/daytime_calibration_20260927/credit_transition/README.md:23)

- Fixed tax is structurally supported: the root solves the dated pension path while holding payroll tax fixed. [run_e5f_current_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_current_transition.py:307) But all pension starting values/bounds must be pinned in the new plan; no stationarity constructor or analytic pension default may enter. The saved block0506 pension is `0.917784…`; the manifest also records the actual balanced tax and accounts. [fixed_reference_manifest.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/fixed_reference_manifest.json:1227) [fixed_reference_manifest.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/fixed_reference_manifest.json:1842)

- Estate funding is only provisional in the current runtime and uses the actual next entrant cohort in dated audits. [e5f_current_transition_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_current_transition_runtime.py:190) The new overlay must pin an origin-specific estate auditor compatible with block0506’s DUE flag; it cannot reuse a general stationary ledger.

- A genuine transition also requires an independently authenticated closed terminal checkpoint. Current production validates exact grid, 31 reported parameters, many primitives, terminal receipt, and identical model pins. [run_e5f_current_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_current_transition.py:257) No such block0506 2023 endpoint is identified in the supplied September 28 packet, so production transition work remains blocked after the reference adapter.

## Credit: critical incompatibility

Block0506 saves `native_due_stayer_credit=True` and \(\phi=(0.8,0.8,0.8,0.8)\). [fixed_reference_manifest.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/fixed_reference_manifest.json:1175) [fixed_reference_manifest.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fertility_identification_20260928/fixed_reference_manifest.json:1279)

The current root unconditionally sets `native_solvency_credit=True`. [run_e5f_current_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_current_transition.py:158) The solver explicitly rejects simultaneous natural-solvency and DUE stayer modes. [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:3217) Thus the existing natural-credit root will fail or requires a new lead-reviewed credit overlay; flipping flags silently would change repayment economics.

Also:

- “No down-payment constraint” is **not** \(\phi=1\). \(\phi=1\) makes \((1-\phi)pH=0\), but still permits only debt down to \(-\phi pH=-pH\). [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:3705)
- Current natural credit removes the purchase threshold by setting it to \(-\infty\), but still clips purchase support at `b_grid[0]` and retains strict continuation/death solvency. [solver.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:3418)
- Grid support is not an economic borrowing rule: it remains a numerical feasibility boundary and must be audited separately. The existing audit reports both lower/upper-grid mass and outside-grid mass. [e5f_solvency_credit_benchmark.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/e5f_solvency_credit_benchmark.py:146)

## Existing gates and smallest smoke

Existing useful gates include pin/approval checks, exact cache-on/off equality, budget/estate/purchase audits, mass and policy reproduction, market/fiscal checks, and terminal distribution/queue/price/rent checks. [run_e5f_current_transition_smoke.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_current_transition_smoke.py:235) [run_e5f_perfect_foresight_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_perfect_foresight_transition.py:437)

Tests are mostly structural/synthetic: runtime tests perform no model imports or solves; root tests use fabricated artifacts. [test_e5f_current_transition_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/test_e5f_current_transition_runtime.py:1) [test_run_e5f_current_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/test_run_e5f_current_transition.py:12)

Smallest safe smoke: a two-date, fixed-stock, saved-credit **reference-only** mapping with identical price, pension, \(\psi\), distribution, and terminal checkpoint; require exact cache equality and all dated gates. Do this before the credit experiment. Removing borrowing/down-payment limits is itself an economic shock, so it cannot be called a zero-shock invariance test; test it separately against the same authenticated reference with explicit DUE/solvency semantics.