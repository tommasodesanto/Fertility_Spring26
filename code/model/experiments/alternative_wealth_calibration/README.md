# Alternative timing with model-matched wealth target

This isolated local experiment reuses the reviewed soft-financing calibration
driver and the post-interest transaction overlay. It changes one empirical target
value: pooled 2005/07 PSID aggregate net worth excluding business/farm equity,
other real estate, and vehicles, divided by gross head-plus-spouse labor earnings.
The ratio is 4.45838713455674. The source is
`output/model/wealth_numerator_match_20261002/wealth_numerator_match_results.csv`
(`ratio_c_noveh`). The catch-all “other assets” component remains included.

The original `wealth_earnings` target is 6.92658379107299. Its numerical
objective weight, 7.595098472533724, is retained to isolate the target-value
change; the new ratio has no estimated standard error. The bequest/wealth target
is unchanged, and its denominator compatibility remains outstanding. No model
observer, economic primitive, financing rule, floor, entry distribution,
earnings process, numerical gate, or other target changes. The active target
count remains ten scored moments for ten free coordinates; informative rank is
not certified.

`calibrate.py` preserves each native `target_fit.csv` as the original-contract
diagnostic. It writes `target_fit_new_contract.csv` and
`wealth_rescore_receipt.json` for the experimental score. The receipt records
both target fingerprints and verifies the exact one-row loss identity. Final
native verification compares the original reports and the rescored target
tables across exact repeats. `supervise.py` queues ten distinct starts in two
waves of at most five with a memory cap, 2-hour-45-minute per-chain budget,
250-call cap, 1800-second native reserve,
one-thread environment, and no automatic retry.

The run packet and explicit launch receipt are under
`output/model/fixed_reference_economics_20260928/alternative_wealth_local_20261003_v1/`.

## Widened-beta Torch variant

`cluster_calibrate.py` is a separate, pinned ten-chain variant for Torch. It
starts at the numerically verified new-wealth local chain 02 point and changes
only the **search bound** for annual discount factor $\beta$ from $[0.940,0.990]$
to $[0.930,0.990]$ relative to the local experiment. All ten coordinates
remain free during optimization. The other nine starting coordinates equal the
verified point; the ten annual-$\beta$ seeds are recorded in the run packet's
`start_plan.json`. It requires that exact start file and its SHA-256 on every
invocation. The objective-call cap is 500 per chain, subject to a six-hour wall
limit that reserves the final 1800 seconds for native verification. These are
search-design changes, not changes to the target, observer, economics, or
numerical acceptance gates. The active 10 scored moments for 10 free
coordinates still have unverified informative rank.

The cluster variant's input and smoke receipts are under
`output/model/fixed_reference_economics_20260928/alternative_wealth_cluster_20261003_v1/`.
Direct Linux/Apptainer execution uses `cluster_calibrate.py`; local macOS
validation uses `local_entry.py` to load the frozen source overlay. A two-call
native smoke and a zero-solve child preflight at $\beta=0.930$ passed before
Torch staging. Neither is a production calibration result.
