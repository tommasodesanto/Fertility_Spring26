Implemented `parallel32_overlap` only in the three authorized experimental-clone files.

- New 32-worker profile with `5400 = max(3600, 4200) + 1200` final reserve.
- Exactly reuses `parallel32_fixed` initialization and frozen-selection behavior.
- Starts one pinned 4-process policy finalizer before the 24-history final batch; polls it during batch progress; awaits and validates that same process only after both exact repeats pass.
- Early finalizer failure stops/reaps active process groups. Finalizer output remains unverified until repeats pass.
- Contract builder pins the overlap/frozen-selection fields, 4-worker policy forecast requirement, and 4200-second policy cap.
- Added pure-controller tests plus short real subprocess overlap and orphan-process-group cleanup tests; bad/missing repeats block final success without changing selection.

Verification passed:

```text
code/model/.venv/bin/python …/test_e5f_joint_nested_long_search.py
Ran 28 tests in 1.439s — OK
```

Also passed `py_compile` and `git diff --check`. No model solves, cluster work, commits, or pushes. Remaining: lead’s requested Torch compiled tests and pinned full preflight.