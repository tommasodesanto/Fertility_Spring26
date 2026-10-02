# Fresh 80% hard and quarter calibration

The independent 24-chain Torch array **19040483** is complete. The 15:29 New York [monitor snapshot](monitor/STATUS.md) records 24 terminal chains and 24 passed fresh selected-point postchecks. Hard has 909 valid scored cases, 27 budget-uncomputed cases and zero numerically inadmissible cases; quarter-saving has 995, 21 and zero respectively. The economic rules, targets, weights, ten free coordinates and bounds did not change.

| 80% purchase rule | Verified slot | Loss | Previous verified loss | Full target fit | Estimates and bounds |
|---|---:|---:|---:|---|---|
| Hard | 8 | 88.5884027814 | 97.0112198128 | [14 rows](monitor/hard_verified_slot8_target_fit.csv) | [31 rows](monitor/hard_verified_slot8_parameters.csv) |
| Quarter-saving | 19 | 48.3199378291 | 51.5560360491 | [14 rows](monitor/quarter_verified_slot19_target_fit.csv) | [31 rows](monitor/quarter_verified_slot19_parameters.csv) |

Material misses remain. Hard mean rooms are 6.188 versus 5.729, first-birth room increase 1.098 versus 1.465, ownership ages 30–55 0.589 versus 0.676, and early fertility 0.520 versus 0.810. Quarter-saving values are 6.100, 1.224, 0.641 and 0.521 against the same respective targets. Hard \(h_P=2.6\) contacts its upper bound; quarter-saving \(h_P=2.5720975\) remains below 2.6. Both fertility curvature parameters are flagged near a bound in both arms, but neither is at an endpoint. These are numerically verified experimental selected points, without optimizer-convergence certificates or paper-baseline adoption. The [collection receipt](collection/completed.json) authenticates all 24 native selected-point postchecks and SHA-256-verifies both winners' complete root/repeat packets, including 17 standard plots each: [hard slot 8](collection/hard/slot_8/selected_root) and [quarter slot 19](collection/quarter/slot_19/selected_root). The monitor is paused after completion. Previous policy outcomes have not been recomputed at these fits.

Read [DESIGN.md](DESIGN.md) for the six verified centers, 24 deterministic starts, unchanged economics/targets/bounds, source pins, budgets and stop rules. [submission_receipt.json](submission_receipt.json) records the job and preflight gates. The staged root is `/scratch/td2248/projects/purchase_fresh_calibration_v1`; results are under `results/slot_N/`. One read-only status snapshot can be collected with:

```sh
code/model/.venv/bin/python output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/monitor.py
```

The monitor writes `monitor/` locally and never mutates the cluster jobs. The search source is frozen at submission; do not restage it into the running root.
