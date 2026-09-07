Implemented the scoped parallel finalizer change.

- Default remains serial when `policy_workers` is absent; only exact integer `1` or `4` is accepted.
- `4` uses four spawned processes, one per fixed policy case. Each process retains the original per-case monkeypatch/audit/save body unchanged and writes only its own folder.
- Parent collects results in canonical `CASES` order, preserving receipt schema/order. Unexpected child failure triggers parent termination, timed reap, then kill/reap if needed.
- `prepare()` and existing integration behavior remain unchanged.

Added focused pure test:
`tmp/e5f_joint_nested_full_20260906a/code/model/tools/test_e5f_joint_nested_finalize_parallel.py`

Verified:

- New pure orchestration tests: 3 passed.
- Existing integration plumbing tests: 4 passed.
- `py_compile` passed.

No scientific/full compiled policy smoke was run. The two task files are intentionally ignored by the repository’s `/tmp/e5f_joint_nested_full_20260906a/` exclude rule, so they do not appear in normal Git diff/status; lead should inspect them directly before pinning the updated finalizer hash and adding `policy_workers: 4` to the contract.