# Fixed-β sensitivity: preparation stopped before a GE solve

The proposed diagnostic would hold the verified alternative-timing,
model-matched-wealth chain 2 point fixed except for annual β at 0.95 and 0.94.
The candidate binding uses the established four-year conversion
`P.beta = beta_annual ** P.period_years`; ψ and all other parameters remain fixed.
The native evaluator would re-solve the equilibrium price. The request for 0.95
is in `beta_095_v2/request.json` and the isolated driver is `run_fixed_beta.py`.

No new numerical sensitivity result was obtained. The first 0.95 attempt
stopped before a model solve because the wrapper checked β on the unbound seed
object. The guard was corrected to check `inputs.bind(P, point, bounds,
"floor")`. The second 0.95 attempt passed that check but stopped at the
unchanged native reserve gate before its initial price evaluation. The native
budget sets `stage_deadline_seconds = 300` in
`output/model/publication_refactor_20260929/single_market_verification_v1/source/driver.py`;
`run_phase_b` requires more than `300 + 400 = 700` seconds remaining for
selected reporting and exact repeat. The authorized seven-minute (420-second)
OS and native cap cannot satisfy that gate. No 0.94 run was launched, no gate
or time budget was changed, and no automatic retry is scheduled. See
`recovery_receipt.json` for the two attempt statuses and source paths.

An earlier β sweep exists at
`output/model/wealth_diag_slot19_20261002/beta_sweep/progress.csv`, including
0.95 and 0.94. It uses the older quarter-saving purchase rule, old wealth
target, and another parameter point. Those results are historical illustration,
not matched evidence for this calibration.
