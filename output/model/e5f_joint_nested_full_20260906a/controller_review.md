Implemented the bounded controller and pure tests:

- [run_e5f_joint_nested_long_search.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py)
- [test_e5f_joint_nested_long_search.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/test_e5f_joint_nested_long_search.py)

Validation passed: syntax compilation, 3 pure controller tests, and `git diff --check`. No model run, job submission, commit, or push.

Lead contract requirements: schema `e5f_joint_nested_long_v1`; inline `base_plan` plus canonical `base_plan_sha256`; `controller_sha256`, `adapter_sha256`, `planner_sha256`; seed/reference paths and SHA-256s; source/bundle/target/domain values matching `base_plan`; fixed budget fields; `absolute_finish_epoch`; optional `reused_helper_sha256` and finalizer path/hash.