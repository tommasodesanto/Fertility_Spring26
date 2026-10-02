# Fixed-coordinate parenthood housing floor diagnostic v2

This isolated restart packet preserves the v1 diagnostic contract and fixes its
replay-gate process handling. V1's gate called the old-bound and extended-bound
native solves sequentially in one Python interpreter; the first solve succeeded
and the second raised `RuntimeError: Fresh interpreter required: model already
imported`. In v2, a parent gate starts one child interpreter for each arm. Both
arms receive the same absolute deadline, set once at gate start, so process
startup and the first solve count against the shared 1,800-second budget.

The gate still requires the exact incumbent old-bound target table, 31 parameter
rows, selected price, derived H0, and native ROOT/REPEAT report with 17 plots;
then it compares the extended-bound arm's target table, all estimates, price,
H0, loss, and every non-h_P parameter field. The gate receipt remains
`status=incumbent_replay_passed`, with `old` and `extended` result objects and
`source_manifest_sha256`. The four points remain h_P = 2.3, 2.4, 2.5, and 2.6.
Only the local h_P upper bound changes from 2.3 to 2.6. The same 14-target,
31-parameter, ROOT/REPEAT, and 17-plot checks remain active. No target, weight,
parameter, or economic input was changed.

`incumbent.json` is byte-identical to v1. The manifest pins this v2 runner and
the same normalized-calibration-v2 code, plan, source pins, target fingerprint,
weight fingerprint, and authenticated incumbent. The stage and submit scripts
use the unique remote directory `/scratch/td2248/projects/normalized_floor_extension_v2`;
all authenticated v2 source and input mounts are reused. Staging and submission
have not been run.

The zero-solve mock uses the production parent/child command path for both gate
arms, confirms distinct fresh-process receipts with the same deadline, and
retains negative target-drift and non-floor-parameter-drift rejection checks.
Run from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/normalized_floor_extension_v2/run.py --mode mock --out /tmp/normalized_floor_extension_v2_mock
```

For lead review only: after independently reviewing the source and mock receipt,
the lead can stage and submit with the v2 scripts. This agent did not stage,
submit, cancel, or run model code.
