Blocking finding:

- [`run_e5f_joint_nested_long_search.py:528`](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:528) crashes if the initial batch drains with no valid receipt: `available` becomes `[None]`, then lines 531–532 subscript it. Thus three initial timeouts (or any all-rejected initial drain) ends as `failed` rather than reaching `final_assessment()`’s intended `no_valid_completed_case` path at line 564. This does not affect the valid-incumbent path, but it is an unhandled timeout-stop outcome.

Otherwise, the intended behavior is implemented correctly:

- Triple timeout sets a sticky reason and writes `search_stop.json` ([339–343](.../run_e5f_joint_nested_long_search.py:339)); `can_fit` and `submit` then refuse all non-final work ([243–248](.../run_e5f_joint_nested_long_search.py:243), [375–381](.../run_e5f_joint_nested_long_search.py:375)).
- Existing futures are drained rather than killed; every returned complete receipt is recorded before the batch exits ([383–399](.../run_e5f_joint_nested_long_search.py:383)).
- Unstarted cases are explicitly saved with plan/SHA metadata and are not counted as completed or rejected ([400–403](.../run_e5f_joint_nested_long_search.py:400)).
- The called collector accepts a partial plan, validates every completed receipt, records missing IDs, and does not fabricate rows or losses ([build_e5f_bounded_refinement_plan.py:29](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/build_e5f_bounded_refinement_plan.py:29)–89).
- Final Jacobian/repeats still use the unchanged absolute budget and are permitted despite the sticky stop ([571–579](.../run_e5f_joint_nested_long_search.py:571)). Final-repeat rejection remains fatal because it is a smoke batch ([394–397](.../run_e5f_joint_nested_long_search.py:394)).
- Unexpected scientific failures remain fatal at line 395; they do not become ordinary rejected proposals.

Checks: read the required current context excerpts; inspected the current scoped controller and its 14-test file, the invoked `planner.collect`, and re-ran the scoped Git diff (no tracked diff was emitted for these snapshot paths). I attempted the pure controller suite: it ran 14 tests, but this read-only sandbox cannot create `TemporaryDirectory`, yielding 8 environment-only errors; six non-temp tests passed.