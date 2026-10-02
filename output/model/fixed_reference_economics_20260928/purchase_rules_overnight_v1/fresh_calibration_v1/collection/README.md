# Fresh calibration native collection

All 24 Torch slots ended successfully and passed selected-point native postchecks. The winning selected points are hard slot 8 (loss 88.58840278139655) and quarter slot 19 (loss 48.31993782905433). These are experimental fits, not optimizer convergence certificates or adopted paper baselines.

`collect.py` reads the frozen Torch results and copies files into this directory only. It authenticates the campaign design, submission, source and runner pins, each slot’s terminal and stage contracts, numerical postcheck, target/weight contract, bounds, equilibrium closure, 14 target rows, 31 parameter rows, and 17 standard plots. It compares every copied file with its remote SHA-256.

Reproduce, when Torch access is available:

```sh
python3 output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/fresh_calibration_v1/collect.py
```

[`completed.json`](completed.json) records the two winners and all 122 copied-file hashes. Its SHA-256 is `5b9b15f11ede097be1b54eff3eba71a5350c309474d598c0aef142e99d33fdc6`. [`all24_postchecks.json`](all24_postchecks.json) retains the compact terminal, search, contract, input and selected-postcheck receipts for all 24 slots (SHA-256 `227aec8eef545d99e6b2e2ac1c8a375eec5177d49e5109701b3c456d0d81de87`).

| Arm | Selected root and standard plots | Native repeat | Final repeat | Files | Bytes |
|---|---|---|---|---:|---:|
| hard | [`selected_root`](hard/slot_8/selected_root/) | [`selected_repeat`](hard/slot_8/selected_repeat/) | [`selected_repeat_final`](hard/slot_8/selected_repeat_final/) | 61 | 73860041 |
| quarter | [`selected_root`](quarter/slot_19/selected_root/) | [`selected_repeat`](quarter/slot_19/selected_repeat/) | [`selected_repeat_final`](quarter/slot_19/selected_repeat_final/) | 61 | 74948263 |

The full target and parameter tables are in each `selected_root/target_fit.csv` and `selected_root/parameters.csv`; the standard 17 PNGs are in each `selected_root/standard_diagnostics/`. The selected native arrays are retained under each `selected_repeat/stage/solution_arrays.npz`.
