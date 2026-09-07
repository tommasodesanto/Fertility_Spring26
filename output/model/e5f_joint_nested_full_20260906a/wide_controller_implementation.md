Implemented the `wide32` controller profile; no contract was built, no SHA invented, and no job submitted.

Changed:

- [run_e5f_joint_nested_long_search.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/run_e5f_joint_nested_long_search.py)
  - Adds `wide32`: 32 workers, 64 population, 640 histories, 8 DE generations, 2 polish rounds, 9h/12h ceilings.
  - Implements deterministic anchor + 47 grid points (omits `(2, .8)`) + 16 radius-`.12` perturbations.
  - Derives population/DE loop sizes from the selected profile.
  - Reserves 22 Jacobian cases and 2 byte-identical final repeats.
  - `wide32` requires hash-pinned imported smoke; v1’s existing local-smoke path remains supported.

- [build_e5f_joint_nested_long_contract.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/build_e5f_joint_nested_long_contract.py)
  - Requires an explicit profile, runtime estimate status (`measured`/`provisional`), and imported-smoke root/proof/original-contract hashes.
  - Labels the runtime estimate without asserting that 640 histories fit nine hours.

- [submit_e5f_joint_nested_long.sh](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/cluster/submit_e5f_joint_nested_long.sh)
  - Requests 32 CPUs / 192GB / 12h.
  - Requires contract path and SHA via environment; removes stale hardcoded hash.

- [test_e5f_joint_nested_long_search.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_joint_nested_full_20260906a/code/model/tools/test_e5f_joint_nested_long_search.py)
  - Adds coverage for 64 deterministic bounded proposals, unique grid coverage, profile/budget accounting, stale imported-contract rejection, and byte-exact repeat-plan copying.

Imported smoke verification: the controller hashes the original contract and smoke proof, checks source/bundle/target/domain against the active contract, then validates every original plan/receipt/artifact before accepting it. Policy smoke proof and receipt hashes are also rechecked. Non-final stages reserve 24 histories, while timeout checks use full per-wave case caps; missing final repeats cannot yield `complete_verified`.

Tests passed:

- `python3 code/model/tools/test_e5f_joint_nested_long_search.py` — 7 passed
- Python compilation, shell syntax check, and `git diff --check` — passed

Only the four authorized files are modified.