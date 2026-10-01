Saved native dated diagnostics
=============================

`floor_report.py` authenticates the selected calibration handoff and a pinned
`mapping.json` or `native_record.json`, verifies every saved native diagnostic
packet hash, and calls the original renderer under actual native bindings. It
never reconstructs the stationary reference or evaluates a transition. Native
`solve_*` entry points are temporarily forbidden while reading/rendering packets.
The actual dated pickle supplies parameters, including the shocked preference.

Use the original calibration deployment mounts and a fresh snapshot of this
helper; do not modify the source inventory of a submitted job. After a mapping
has completed, invoke the following inside that authenticated environment:

```
python code/model/experiments/transition_readiness/floor_report.py \
  --handoff HANDOFF.json --handoff-sha256 HANDOFF_SHA256 \
  --mapping-record ACTUAL_MAPPING/native_record.json \
  --mapping-record-sha256 ACTUAL_MAPPING_SHA256 \
  --run-identity PINNED_RUN_PLAN.json --run-identity-sha256 RUN_PLAN_SHA256 \
  --output FRESH_REPORT_DIRECTORY
```

The optional run-plan arguments validate its numerical identity and source-file
pins; callers should supply them when the run plan is available. The report
preserves the first, middle and last dated packets and verifies the unchanged
17 PNG names and hashes at each date. `report.json` records all packet/input
pins, numerical identity, original mapping gates and residuals, calendar years,
PNG paths and hashes, and zero native calls. Failed map gates are preserved:
reporting is diagnostic only and does not certify a fitted or production path.
Visual review remains pending until the caller inspects the generated figures.

Five zero-model-call tests passed. Actual native rendering must be checked after
the readiness job supplies its dated packets; no native run or job was launched
by this helper's implementation worker.
