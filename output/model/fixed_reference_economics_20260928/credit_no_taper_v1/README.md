# Isolated renter no-taper estate bound — targeted checks passed

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
