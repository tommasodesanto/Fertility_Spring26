# Original-target revised-timing continuation

The [cluster driver](../../../../code/cluster/soft_timing_calibration/continuation/README.md)
owns this isolated ten-chain Torch continuation. [`start_plan.json`](start_plan.json)
pins the ten selected starts, source receipts, exact target/weight contract,
original bounds, and economic-change disclosure. `deployment/` contains the
immutable upload archive and receipt. `collection/` will hold independently
validated terminal results when the search ends.
