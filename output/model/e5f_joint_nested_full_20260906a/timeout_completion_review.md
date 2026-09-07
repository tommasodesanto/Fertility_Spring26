Blocked: the prepared `search_stop_reason` revision is not yet visible in either scoped file.

The current controller still raises at [run_e5f_joint_nested_long_search.py:336](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:336). That enters `batch`’s broad exception handler, which calls `stop_active()` ([lines 391–394](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:391)), so it still kills active cases and prevents `final_assessment()`.

Thus, as currently visible, it does not stop only new search proposals, does not drain asynchronous completions, and cannot preserve the desired receipts/final repeats. Budget reservation and scientific-failure fatality remain as before, but the required operational behavior is absent.

Checks: mandatory context files and git status; reread the scoped controller and test file; inspected `_reject`, `can_fit`, `batch`, `search`, and `final_assessment`. No tests or model jobs run; no files modified.