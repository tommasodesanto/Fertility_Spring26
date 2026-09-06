Implemented the scoped controller and pure-test files:

- [run_e5f_joint_overnight_search.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_e5f_joint_overnight_search.py)
- [test_run_e5f_joint_overnight_search.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/test_run_e5f_joint_overnight_search.py)

Verified static compilation, whitespace checks, and pure tests:

```text
NUMBA_DISABLE_JIT=1 PYTHONPATH=code/model/tools python code/model/tools/test_run_e5f_joint_overnight_search.py
Ran 4 tests — OK
```

No model solves, cluster actions, commits, or unrelated files were changed. The controller expects the lead-created `run_e5f_joint_overnight_case` adapter.