# Matched existing-owner credit check — preparation only

No solves, source freeze, numerical approval or model adoption performed here.
Driver: `code/model/tools/run_e5f_due_stayer_matched_check.py`, owned exclusively
by the matched-check preparer. Other workers own the native core and accounting.

## Explicit comparison

Use original **de_0093** parameters/checkpoint, baseline price and wealth grid.
Do not substitute the newer 41.112 calibration candidate. Hold child-benefit
level fixed; no fertility normalization or parameter search. The only economic
change is `native_due_stayer_credit=True` for the DUE arm. Purchases keep their
baseline borrowing rule; native solvency credit must remain off.

The baseline arm must reproduce the old full 14-target table and all saved
numeric arrays; both arms output all 31 parameter rows after checking unchanged
parameter objects. The diagnostic retains strict operator, probability, policy,
household-budget, pension and estate-funding checks. It produces the same 17
standard plots. A failed gate saves a failure, never a fabricated fit.

Both arms hold the same price. DUE housing demand therefore need not clear.
This is a **partial-equilibrium diagnostic**, not a new equilibrium. The separate
observer path reports the true housing residual and reuses the frozen scoring
function. It does not patch or relax the production reporter, which correctly
requires market clearing. No population rescaling is used to conceal the residual.

## Required preparation and remaining blockers

1. Finish lead review of core and accounting changes and freeze a fresh immutable
   source snapshot. The preparer hashes the complete Python model tree and saves
   a source copy; execution checks every pin. Do not prepare while peers edit.
2. Authenticate the new origin-specific estate audit and dated budget helpers.
   The old frozen estate audit cannot be reused: the runtime deliberately rejects
   DUE with that contract. The plan explicitly pins the replacement, which must
   declare the DUE capability.
3. **Purchase audit remains unresolved.** The runtime's existing purchase audit
   is frozen. This driver requires an explicitly pinned audit with the same
   five-argument `audit_purchase_accounting` interface and a declared
   `SUPPORTS_NATIVE_DUE_STAYER_CREDIT=True` capability. A lead review must establish
   that buyer and stayer saving/consumption rules are handled correctly. Missing
   or origin-blind audit fails before any solve. No capability is invented here.
4. Finalize the separate price-fall experiment. A 10% fall is an explicit experimental magnitude
   that the lead may choose or replace within the existing matched-test authority;
   no additional user permission is required. It must use exactly the inherited pre-fertility distribution with no
   frontier projection, and strict infeasibility must be reported. The draft
   contains a strict pre-gate, but **price-fall execution is deliberately blocked**
   until the expected future price path and dated normalization-row observer are
   specified. A stationary-policy solution combined with an inherited population
   must not silently mix stationary and dated fertility measurements.
5. Set `execution_authorized=true` and an explicit absolute end only after review.
   Keep at most 300 seconds per arm and 1,200 seconds globally. A parent supervisor
   must own child processes and enforce these deadlines even inside compiled code;
   the driver's signal alarm is a secondary guard, not sufficient process-owner
   enforcement. No launcher is supplied and none was run.

Preparation command, after source freeze, using the main project's Python:

```sh
python code/model/tools/run_e5f_due_stayer_matched_check.py \
  --prepare --output NEW_PREPARATION_DIRECTORY
```

This writes `plan_draft.json` with execution disabled, complete source/reference
pins, explicit audit placeholders and source snapshots. A separate reviewed plan
must retain the original prepared evidence and pin its own driver bytes. Arm API:

```sh
python code/model/tools/run_e5f_due_stayer_matched_check.py \
  --plan REVIEWED_PLAN --arm baseline --output NEW_BASELINE_DIRECTORY
```

Run baseline first; lead checks exact 14/31 tables and all policy arrays before
running `--arm due` in a fresh interpreter. Price-fall arms are not runnable yet.
Compile check passed. No runtime import, model solve or numerical result claimed.

## Preparation correction

Imported current-source receipt keys may be absolute or relative. The runner now
normalizes absolute keys relative to its authenticated source root before checking
the plan inventory; paths outside that root fail. Compile check passed, no solves.
The strict no-projection fix is underway in main and is not yet merged into this
isolated worktree. Its final source identity must be reconciled before execution.
