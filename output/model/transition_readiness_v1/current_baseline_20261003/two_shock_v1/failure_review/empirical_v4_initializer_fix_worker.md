Implemented the scoped numerical initialization correction.

- Cold paths now initialize `prices` and `pensions` as constant \(q_T\) and \(b_T\) arrays; no state, queue, gate, bounds, or solver controls changed.
- `warm_start.json` records cold method `stationary_endpoint_flat_initialization`; warm labels remain unchanged.
- Added cold and warm seam regressions verifying endpoint-flat vs. copied warm paths/Jacobian, inherited-state mapping, and key root controls.

Changed files:

- [two_shock_runtime.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/experiments/birth_count_choice/two_shock_runtime.py:95)
- [test_two_shock.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/experiments/birth_count_choice/test_two_shock.py:469)

Verification passed:

```sh
output/model/publication_refactor_20260929/local_env_v1/venv313/bin/python \
  code/model/experiments/birth_count_choice/test_two_shock.py
```

Result: `Ran 27 tests ... OK`.

No native model/transition smoke was run. Native feasibility and acceptance remain unvalidated until the separately authorized fresh tiny package smoke.