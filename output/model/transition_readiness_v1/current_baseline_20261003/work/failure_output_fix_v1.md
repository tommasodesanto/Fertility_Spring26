# Stationary inherited-state diagnostic routing: unapplied fix

**Status:** isolated patch and five zero-solve mock cases passed on October 4, 2026. This patch is for a future separately pinned continuation only. No active source, immutable deployment or running local/cluster job was changed. It changes diagnostic output routing and records that routing; it changes no economic primitive, state, solve equation, projection, acceptance gate or tolerance.

## Evidence and cause

[The monitored failure snapshot](../monitor/remote_health.json) records v6 tasks 10 and 11 failing with `OSError: [Errno 30] Read-only file system` at the saved baseline's `native/phase_b_ge/root_06/stage/inherited_state_failures`. Their last completed endpoint evaluations were valid, with scores `0.1397790114598909` and `0.2537838919257691`, respectively. The next endpoint evaluation failed before a completed record was returned. Each runtime counted 33 actual policy calls; the controller had accounted for 31. These failed tasks do not establish a fitted result.

A bounded read of task 10's remote driver log confirmed the call chain:

```
/scratch/td2248/projects/current_estate_transition_20261003_v6/results/panel_guess_10/driver.log
one_shock_floor.NativeAdapter._endpoint.evaluate
  -> CurrentEstateARuntime.stationary
  -> FloorRuntime.stationary, floor_runtime.py:437
  -> current birth_count_evaluate_period.evaluate_period
  -> gate_pre_fertility_distribution
  -> _require_exact_inherited_distribution, run_dynamic_population_transition.py:449
  -> destination.parent.mkdir -> errno 30
```

The current saved case carries `native_inherited_distribution_evidence_dir` assigned when `model/native_price.py:32` created its stationary stage. `FloorRuntime.stationary` deep-copies those parameters at line 416, changes `psi_child`, and calls the native solve and period evaluation without replacing that obsolete output path. The diagnostic write therefore addresses the saved case, which is read-only in the cluster container. The dated scaffold already redirects this field to its active mapping output at `pinned_tools/run_e5f_preference_transition.py:315–319`; the missing redirect is specific to the stationary path.

The feasibility helper writes evidence for **both** rejected mass and a retained tail below the existing tolerance. Its unchanged criteria are total occupied dead mass above `DEAD_MASS_TOL=1e-12`, or any occupied nonfinite value, with the native value cutoff `-1e9`. It writes the evidence before raising `InheritedDistributionInfeasible`. Consequently, errno 30 masks the classification: task 10/11 logs alone cannot establish whether the attempted diagnostic described a true feasibility failure or a permitted tail. Their dead mass, affected support and census remain unknown because the write did not complete. The patch cannot make a genuinely infeasible endpoint pass; it permits the existing classification and evidence to complete.

The two-call controller/runtime discrepancy is retained as evidence. The exception occurs before `_endpoint.evaluate` calls `_account(record)`, so actual native calls can be missing from completed controller records; this patch does not change counters or certify their equality after failure.

## Minimal proposed change

[The unapplied patch](failure_output_fix_v1.patch) adds four lines to the copied parameter object in `FloorRuntime.stationary`, before the native callbacks:

- preserve the previous diagnostic path for an audit record;
- route evidence to the current endpoint evaluation's `folder/inherited_state_evidence`;
- write `folder/output_override.json` with the saved/effective paths and `economic_change=false`.

The saved parameter object remains unchanged. Each endpoint evaluation has its own directory. No fallback, retry, gate relaxation or population repair is introduced. The inherited controller and its native solve sequence remain unchanged. A future deployment must regenerate source pins and use an isolated output/stage; do not apply this patch to sources underlying current jobs.

## Verification

[The mock script](failure_output_fix_v1_mock.py) applies the proposed diff only in memory, executes the actual stationary setup prefix through its copied parameters, and executes the unchanged calendar feasibility helper extracted from its syntax tree. It neither imports the model runtime nor calls any numerical solve. A mocked read-only saved folder reproduces errno 30 in the old routing even for below-tolerance mass. The proposed routing passes these five cases:

1. Below-tolerance dead mass: writable diagnostic, retained, zero projection.
2. Above-tolerance dead mass: writable diagnostic, unchanged `InheritedDistributionInfeasible` rejection.
3. Occupied nonfinite value: writable diagnostic, unchanged rejection.
4. No bad values: no diagnostic file, unchanged mass.
5. Old below-tolerance route: reproduces errno 30, demonstrating that errno alone does not identify a rejection.

Tests assert independent active output paths, unchanged saved parameters and input distribution, exact original cutoff/tolerance, and unchanged active-source SHA-256 values. [The mock receipt](failure_output_fix_v1_mock_receipt.json) records all five passes, source/patch hashes and zero native calls. The first mock execution exposed a missing `ctx` stub in the test harness; adding that stub resolved it without source changes.

Run from the repository root:

```sh
output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python output/model/transition_readiness_v1/current_baseline_20261003/work/failure_output_fix_v1_mock.py
git apply --check output/model/transition_readiness_v1/current_baseline_20261003/work/failure_output_fix_v1.patch
```

This verifies routing and preservation of the existing feasibility classification with synthetic arrays. It is not a fresh endpoint replay, a proof that either failed task is feasible, or numerical validation of a future continuation.
