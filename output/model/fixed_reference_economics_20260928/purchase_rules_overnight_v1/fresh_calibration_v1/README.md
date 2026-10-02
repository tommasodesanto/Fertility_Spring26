# Fresh 80% hard and quarter calibration

The new independent 24-chain Torch array is **19040483**. All 24 slots were verified RUNNING after submission. An initial read found slot 0 in objective call 2 after one completed full-equilibrium evaluation and slot 12 initializing, with no launcher failure receipt. The earlier verified hard and quarter incumbents, 97.0112198128 and 51.5560360491, remain the comparison points. No new chain is a calibrated result until its fresh selected-point native postcheck passes.

Read [DESIGN.md](DESIGN.md) for the six verified centers, 24 deterministic starts, unchanged economics/targets/bounds, source pins, budgets and stop rules. [submission_receipt.json](submission_receipt.json) records the job and preflight gates. The staged root is `/scratch/td2248/projects/purchase_fresh_calibration_v1`; results are under `results/slot_N/`. One read-only status snapshot can be collected with:

```sh
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor.py
```

The monitor writes `monitor/` locally and never mutates the cluster jobs. The search source is frozen at submission; do not restage it into the running root.
