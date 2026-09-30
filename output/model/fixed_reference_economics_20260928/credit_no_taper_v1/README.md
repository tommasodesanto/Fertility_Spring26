# Isolated renter no-taper estate bound — targeted checks passed

## Current author-selected contract: scalar renter credit and repayment on sale

The author now selects mortgage repayment on sale into renting and an explicit
constant `unsecured_credit_limit`, with initial planned value zero. This replaces
the rollover proposal described below as the specification to prepare. The
isolated sources, overrides and verification are in
[`fixed_credit_contract_v1/README.md`](fixed_credit_contract_v1/README.md).
The frozen reference remains **2007 stationary reference — block0506, September 28 verified export**.
No positive credit magnitude, entry change, recalibration or new GE has been selected.

## Author-requested 30-minute comparison, September 29 evening

The new comparison uses the retained estimated parameters of **2007 stationary reference — block0506, September 28 verified export**; no refit or reference switch. `credit_rule_quick_v1/` compares common-price household outcomes; `credit_rule_ge_quick_v1/` attempts closed stationary GE for each rule with fixed preferences including child-benefit psi. The two rules are rollover of inherited renter debt until mortality repayment (`ours`) and zero renter saving debt at every age plus debt-clearing owner-to-renter sales (`author`). Both retain buyer LTV and incumbent-owner grandfathering. Positive unsecured borrowing is not introduced.

Torch zero-lifecycle compiled smoke **18837749** passed in 11 seconds. Parallel arrays **18837765** (common price) and **18837766** (GE) were submitted once with after-ok dependency. Immutable source packets are under `/scratch/td2248/projects/fixed_reference_credit_rule_quick_20260929/source_credit_rule_quick_v1` and `/scratch/td2248/projects/fixed_reference_credit_rule_ge_quick_20260929/source_packet`. Hard delivery deadline is September 30 **01:00:25 UTC** (September 29 21:00:25 New York); no extension or automatic retries. A numerical root without an exact repeat remains preliminary; failure or incomplete candidates must not be called GE.

GE retains the absolute housing-supply curve and its scale parameter, with population clearing housing and price solving actual birth renewal. This differs from fixing physical housing stock. The inherited estate-counterparty and finite-grid limitations remain. Checkpoints stay on Torch.

Recovery record: GE v1 stopped before any lifecycle solve on a dictionary/set type error. PE v1 completed two lifecycle solves but stopped while writing inherited-state diagnostics into the read-only reference directory, before saving aggregates. These failed outputs and immutable sources are preserved. GE v2 was prepared but never submitted. The bounded logging recovery is **PE v2 job 18838216**, with two additional lifecycle solves at most, and **GE v3 job 18838220**, retaining the original GE solve cap and absolute deadline. Both redirect only the diagnostic destination and save `raw_lifecycle.json` before reconstruction. No feasibility tolerance, projection, entry position, economic primitive, or reference source changed. The original deadline was not reset.

Final result: [compact comparison, full tables and diagnostics](credit_rule_quick_v2/collected/README.md). Rollover passed common-price gates; strict zero debt rejected approximately 0.00802% of entrants. Both GE cases remain uncomputed. All jobs are terminal. The remaining sections describe the earlier standalone overlay test; their uncomputed list applies to that original test, not to the later comparison.

Reference: **2007 stationary reference — block0506, September 28 verified export**.

This packet does not change `code/model/`, the frozen Torch source, its checkpoint,
or any other packet. `patch_driver.py` first verifies the supplied frozen hashes for
`parameters.py`, `solver.py`, and `kernels.py`, then writes exactly one changed file
to a new `overlay_v1/` directory. It fails if that directory already exists.

Economic change when `P.renter_no_taper_estate_bound=True`: only when
`P.lambda_d == 0`, a renter with current unsecured balance \(b\) faces
\(b'\ge\min(b,0)\) if the current choice has no death probability and \(b'\ge0\) if
`j == J-1` or (`use_age_survival` and `survival_probs[j] < 1`). Thus existing
negative renter debt may roll over before death risk, but there is zero new
unsecured borrowing and no age-42--62 taper. Positive credit is explicitly
rejected; it remains an author choice. Lifetime repayment/nonnegative renter
estates are retained. Owner, buyer, incumbent, estate, purchase-cash, natural-
credit, earnings, fiscal, target, grid, and timing logic are unchanged.

The native scalar floor is retained: setting the next-decision weight to one and
its cap to zero gives `min(b,0)`; setting both to zero gives zero. The terminal
appended entry remains zero. The default flag is absent/off, so the frozen
age-tapered schedule is unchanged; it is also listed as a dynamic override key
so parameter rebuilds retain it.

## Torch verification and limits

Lead-reviewed Torch job **18835438** completed successfully in eight seconds.
All **11** targeted checks passed (zero model solves); compact receipts and log
are in `verification_v1/`. The effective isolated `parameters.py` SHA256 is
`9c6def300f76b2d5ac55c392e8a595fca881b1d78b15c1cba96016a30a3b83b9`.
The reference manifest and all three frozen source hashes were checked again
after the test and remained unchanged. Lead inspected the actual generated
diff: only the default-off flag, override registration, and renter schedule
construction differ. Solver and kernel source bytes are unchanged.

The checks used pure Python mode (`NUMBA_DISABLE_JIT=1`). They validate the
changed parameter builder and native floor functions, not a compiled lifecycle
solution, production evaluator integration or a full baseline repeat. The
source driver/test/launcher are pinned in `source.sha256` and read-only on Torch.

Original one-shot launch recipe (already completed; **do not resubmit**):

```sh
sbatch output/model/fixed_reference_economics_20260928/credit_no_taper_v1/run.sh
```

The launcher requests 1 CPU, 16 GiB, and 5 minutes; it has no retry behavior. It
performs zero lifecycle, distribution, or GE solves and reads no checkpoint. It
writes a compact overlay manifest with exact frozen/effective/test hashes and a
PASS/FAIL test receipt. The tests call the changed `build_debt_caps` and the native
`solver.renter_borrowing_floor`, covering absent-flag identity, zero/negative debt
at every age, first and later mortality, terminal death, positive-credit rejection,
owner boundaries, and static scalar plumbing to both native saving-stage paths.

`overlay_v1/` contains only `parameters.py`; it is not production-runtime acceptance.
The targeted test uses a fresh Python process, registers that overlay as
`intergen_eqscale_seq_optimized.parameters`, and only then imports the native solver.
Any future Torch runtime check must use that same fresh-process registration before
any runtime package import, or a separately authenticated isolated package copy. It
must not copy the repository or any checkpoint, and production integration remains
unvalidated by this packet.

Uncomputed: production evaluator integration, native baseline repeat,
behavioral response, lifecycle solution, stationary distribution,
equilibrium prices, target fit, recalibration, and all 17 standard diagnostic plots.
No production run or reference switch has occurred. The frozen reference is
preserved. Positive unsecured credit remains a separate open author choice.
