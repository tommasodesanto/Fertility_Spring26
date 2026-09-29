# Optimized two-birth diagnostic

The author approved one fully solved two-birth experiment on September 28.
Reference: **2007 stationary reference — block0506, September 28 verified
export**. The preceding fixed-policy replay is in `../two_births_v1/`; its
0.698 age-25 count and 2.362 lifetime count are not results for this model.

## Experimental economic contract

After one successful birth during a four-year period, a household with room
below the existing maximum of three children may choose one additional
attempt. There is no additional opportunity after waiting or failure, and
at most two new births may occur in the period. Housing, consumption and
saving respond to the resulting number of children.

The new conditional wait/try decision uses the existing later-birth choice
scale. Its expected maximized value enters the preceding birth decision.
Thus the experiment adds a choice/shock opportunity as well as relaxing
spacing; it is not the previous fixed-policy mechanical experiment. The
second conception draw is conditionally independent and uses the existing
age-specific probability. These within-period assumptions are experimental,
not adopted. They introduce no new estimated parameter.

The ten reference calibration coordinates stay fixed. Child benefit is
renormalized to completed fertility 2.1 and demographic renewal is enforced.
Housing prices and the inherited pension balance are solved normally. Income,
entry wealth/income, housing/fiscal primitives, targets, weights, bounds and
scientific tolerances are unchanged. The empirical/model within-cell age
projection remains linear in pre/post stocks and does not date the two births
separately. First-birth housing and recent-parent observers must count the
correct families, including families whose first success leads to two children.

## Isolation and verification

No shared source file is edited. Patch builders generate full source copies
inside the Torch job's output, with original/effective hashes and diffs. Only
that worker process installs changed function bodies and the extended owned
policy cache. Original ancestry authentication remains distinct from the
additional effective-source manifest. Checkpoints stay on Torch outside Git.

Before the experimental evaluation, synthetic Bellman/flow/observer/cache
tests and a full reference Bellman/distribution replay with the option disabled
must pass. The replay uses the saved reference price and child benefit and
compares every core value/policy/population array. The normalized experimental
evaluation then retains all market, probability, household budget, transaction,
estate, value, pension and renewal gates. It reports all 14 target rows, all
31 parameter rows and the 17 standard plots. No automatic promotion follows.

## Bounded execution plan

One Torch smoke job: synthetic tests plus one fixed-price reference solve,
20-minute internal cap. One separately authorized diagnostic evaluation after
that smoke: at most 23 stationary solves, 70-minute internal total cap,
one model worker. The reference normalization used six stationary solves and
971 seconds; the changed model's time and normalization difficulty are unknown.
The hard cap replaces any assumption of a speed gain. There is no calibration
search, continuation of an expired job, transition or automatic retry.

The supervisor writes a heartbeat every 15 seconds and native normalization
writes each stationary solve. Unknown failures, scientific gate failures or
deadlines stop the experiment and preserve the failure. A passed point still
requires review of all targets and plots before any decision about further
recalibration. Source changes after the smoke require a new version and smoke.

All implementation changes received lead review and an independent integration
review. The final case also audits extra-probability menu sums, occupied dead
menus, and equality plus independent ownership of solution/policy/parameter
caches. The standard fertility probability plots retain their original meaning:
they show the outer attempt, not the added conditional opportunity.

Smoke **18753550** stopped before any model solve: 14 of 15 synthetic tests
passed; a later canonical import bypassed the already-patched private recent-
parent observer module. `failed_smoke_v1/` preserves the driver and exact log.
The installer now loads canonical modules first and patches all loaded aliases.
Review also identified the live fertility observer's frozen source path; its
bytes match the authenticated current definition (SHA `ca120b5bddc7cf6a8ad48237c821320432c90aeb737b6905ad156d6883474208`).
The installer now patches that exact live function and records every source
path/module name. No target, definition or tolerance was changed for this repair.

Queued smoke **18753837** was cancelled before starting to include both alias
fixes in one test. Revised smoke **18753907** runs the same 15 synthetic tests
and one full flag-off reference replay. At that stage no optimized evaluation
had launched.

**Revised smoke passed:** job 18753907 completed in 105.1s, all 15 tests and
one full reference-price solve. All 12 value/policy/population arrays match
exactly (maximum difference zero). Receipt and effective-source diffs are in
`smoke_v3/`; the receipt SHA is
`c998e63b3655e5c83320a77c689ee871bdb31e6473e3653fa16ef9cf55f60498`.

After review of that receipt and every source pin, the single experimental
evaluation was submitted as **Torch 18754308**, currently pending. `launch.json`
records the changes, budget and exact source pins. No additional calibration,
transition or automatic promotion follows. Latest/best summaries on Torch
explicitly report no completed experimental point until one passes.
