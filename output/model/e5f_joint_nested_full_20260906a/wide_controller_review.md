Review found two launch blockers.

1. Policy receipts will always be rejected.

   [run_e5f_joint_nested_long_search.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:112) requires:
   `family-ltv-95`, `property-tax-double`.

   The unchanged finalizer emits:
   `dependent-child-ltv95`, `property-tax-2pct-no-rebate` in [run_e5f_joint_nested_finalize.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_finalize.py:21).

   This rejects both the imported policy-smoke receipt and the final four-path receipt despite valid gates. Repair: use the finalizer’s exact four canonical names in `validate_policy_receipt`, and add a pure receipt-validation test.

2. The advertised 640-history ceiling is not schedulable in the stated 12-hour allocation.

   [build_e5f_joint_nested_long_contract.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/build_e5f_joint_nested_long_contract.py:68) records the implied wall time. With measured \(2364\) seconds/history,
   \[
   640 \times 2364 / 32 / 3600 = 13.13\text{ hours},
   \]
   before Jacobian repeats, policy paths, and overhead—above the 12-hour Slurm limit in [submit_e5f_joint_nested_long.sh](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/cluster/submit_e5f_joint_nested_long.sh:6).

   Repair: treat 640 as a hard ceiling only, record a computed time-feasible attempt cap (with final-policy reserve), and fail contract construction if the declared planned cap exceeds the deadline.

Additional finding: [run_e5f_joint_nested_long_search.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:551) can label the run `complete_verified` after exact repeats even when the finalizer was skipped for insufficient time ([line 603](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py:603)). Use a distinct “calibration repeats verified; policy finalization incomplete” status unless a complete final policy receipt exists.

No blocker found in the proposed 64-to-two-32-case chunking itself: global case IDs flow through the ledger and DE lookup, each plan stays under the adapter’s 40-case cap, and the shared executor remains capped at 32 workers. Add a mocked pure test asserting plans `{1..32}` and `{33..64}`, merged results `{1..64}`, and two planner collections.

Imported-smoke provenance is materially stronger now: original contract/plan controller–adapter–planner pins, four receipt-backed cases, artifact validation, exact anchor replay, and selected-summary-bound policy receipt checks are all present. I would additionally compare imported plan `helper_sha256`, `target_set`, and `first_child_jump_upper` with the original contract’s base plan to close the remaining stale-plan gap.

Verification: controller test execution was attempted; five tests ran, but two tempfile-based tests could not create a directory because this review environment is read-only. No model, submission, or writes were performed.