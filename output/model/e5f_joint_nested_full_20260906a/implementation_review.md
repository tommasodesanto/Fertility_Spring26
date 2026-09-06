Implemented dated joint-choice plumbing in the two owned scripts, default-off.

Changed functions:

- `PolicyBundle`, `policy_from_solution`, `solve_policy`, `evaluate_period`, `reconstruct_stationary_pre_fertility` in [run_dynamic_population_transition.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_dynamic_population_transition.py).
- `first_birth_housing_response`, `begin_dated_first_birth_housing_branch`, `finish_dated_first_birth_housing_branch` in [run_e5f_transition_calibration.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_transition_calibration.py).
- Added [test_e5f_joint_nested_integration.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/test_e5f_joint_nested_integration.py).

Key behavior:

- Policies own `joint_choice`; joint evaluation calls `factor_distribution(..., mode="natural")` and advances using the resulting distribution-specific tenure kernel.
- Dated first-birth branches use treated/control selection kernels at origin, then natural/wait kernels at destination.
- The older same-policy diagnostic now routes through the same dated joint branch logic.
- Non-joint paths and target/gate definitions remain unchanged.

Verification:

- `python -m py_compile` passed for both scripts and the new test.
- `git diff --check` passed.
- The synthetic test could not run: this worktree has no model virtual environment and the system Python lacks NumPy. No full model run was attempted.

Caveat: this is plumbing verification only, not numerical certification. Remaining outside-scope constructor needing joint attachment if used with `joint_nested_choice=True`: `policy_from_objects` in `run_e5f_perfect_foresight_transition.py`. Existing calendar/open-population advancement already consumes `evaluation.policy`, so it receives the effective kernel from dated evaluation.