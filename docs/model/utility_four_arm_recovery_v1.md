# Utility comparison recovery proposal, version 1

**Preparation only; no launch is approved or implemented.** The main research
task requested this separate recovery design after all four searches in
`run_001` stopped during the initial population. It does not alter array
18567879, its source snapshot, contract, records or eight-hour clock. Those
controllers finish their existing repetitions and exports. The experimental
utility choices and adopted common pension update remain exactly as disclosed
in [the original design](utility_four_arm_preparation.md).

The proposed recovery completes the same initial proposal bank and, if its
checks pass, the same three finite differential-evolution generations. It
reuses authenticated completed evaluations and never retries a failed or
timed-out point. The necessary changes are to failure handling, time budgets
and coordinated search bookkeeping, not the objective or model mathematics.

## What the failures establish

Both floor evaluations at `initial_0005` reached the native feasibility gate:
mass above \(10^{-12}\) occupied states without a feasible value at age 34.
The recorded renters had negative budget slack. The proposed policy rejects
these attempted evaluations under the **same** gate. Small mass is not an
excuse to accept them. Nor does the census prove that the parameter vector
is infeasible at every equilibrium price or fertility-benefit coefficient:
the failure can occur in an intermediate solve. A primitive or model bug is
not ruled out. This is a reviewed computational rejection with no scored loss.

Continue only for the reviewed native exception, stage and complete finite
renter-deficit census, authenticated against source, contract, arm and point.
Entry failures, owners, unfamiliar stages/states, nonnegative slack, malformed
or potentially truncated censuses, inconsistent mass totals and unknown
exceptions remain fatal. Smoke or selected-repeat failures also remain fatal.
The [independent review](../../output/model/utility_comparison/recovery_v1/feasibility_review.json)
sets out the exact evidence and negative test cases.

The original share `initial_0012` traces explicitly identify the wall-cap
alarm as the cause, later wrapped in `SystemError`. They completed ten and
eleven native stationary solves before the interrupted solve. This does not
justify treating arbitrary `SystemError` messages as timeouts. A prospective
controller must own the deadline and record termination/reaping outside the
numerical kernel; the worker's inner alarm must not erase that evidence.
Other legacy failures need their own positive evidence before reuse.

## Proposed execution policy

1. A successful candidate must still pass every existing source, target,
   normalization, equilibrium, fiscal, entry, feasibility, purchase-accounting,
   value/probability and output-integrity check. Exit code zero or a partial
   checkpoint is insufficient. Nothing in this proposal weakens those checks.
2. A reviewed feasibility rejection consumes its planned slot with no loss.
   A controller-proven candidate timeout also consumes its slot but remains
   **censored: its scientific loss is unknown**. These are separate statuses.
   Neither can enter the best-point ranking or selected repetitions.
3. Unknown failures stop new dispatch. Already-running siblings finish within
   their original deadlines. A shared four-arm barrier blocks the next common
   phase after any fatal failure or missing/unrun result. This avoids letting
   unaffected arms advance through additional generations alone.
4. Keep the original stop above 50% inadmissible proposals. Additionally stop
   above 50% combined rejected and timed-out slots, reporting the categories
   separately. This proposed computational guard is at least as restrictive
   as the old one and leaves at least twenty scored slots in a forty-slot
   barrier. It is an explicit search-policy choice for lead review.
5. Use typed outcomes in DE bookkeeping. A timeout is not the old `None`
   representation for a scientifically rejected evaluation. Finite planned
   coordinates can still provide proposal geometry. A successful new trial
   may replace an unscored slot, but that is not evidence of an objective
   improvement over its unknown parent. Unrun/missing results cannot cross
   the barrier. Donor streams, bounds, mutation 0.7 and crossover 1 stay fixed.
6. Deduplicate by arm, scientific source/target identity and the full parameter
   vector. Reuse a prior success through its original receipt; a duplicate
   failed or timed-out point keeps its linked rejected/censored outcome and
   is counted as a skipped duplicate with no new attempt or score. Do not
   solve it again or invent a replacement slot. Such an explicit reference
   is distinct from an unexplained missing or unrun result.

These rules produce matched **planned** coverage and common stopping phases.
They cannot guarantee equal numbers of successful evaluations. Final reports
must retain failures, time censoring and all incomplete/unrun slots; a best
loss is not evidence of an equally explored or globally optimized comparison.

## Reuse and finite cost

The [bounded inventory](../../output/model/utility_comparison/recovery_v1/observed_reuse_inventory.json)
records existing controller outcomes and receipt hashes. It is preliminary:
the share controllers were still finishing, and this inventory does not
independently rehash large checkpoints. A launchable reuse manifest must be
made only after the existing run is terminal and checked against the original
approved contract
`c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a`.

Reuse the four exact smoke pairs, the inherited anchors and every other valid
completed point. Authenticate original receipts, checkpoint bytes, complete
target and parameter tables, arm/point identity and all scientific pins on
Torch. Record the old contract as ancestry in a distinct recovery contract;
never rewrite an old receipt with a new fingerprint. New controller code must
not change the economic source, target system, observers or numerical gates.
If that identity cannot be established, reuse is not approved. No extra anchor
solve is needed merely because orchestration changed.

Each floor arm has 29 never-dispatched initial proposals left. The share arms
already dispatched all 39 new proposals; their duplicate smoke seed supplies
slot zero. Thus the ceiling is **58 new initial objectives, 480 genuinely new
DE objectives and at most eight selected repetitions: 546 new objectives**.
An unchanged selected point can reuse its existing two exact repetitions.
Failures and duplicates can reduce work; they do not create replacements.

Keep the existing normalization algorithm, target 2.1, tolerance 0.0005 and
maximum 23 stationary solves. Propose a 9,000-second process cap using an
explicit planning allowance of 360 seconds per native solve and 720 seconds
for setup, observers and scoring: \(23\times360+720=9{,}000\). The largest
logged completed solve in the two diagnosed share timeouts was 313.080
seconds; 360 and 720 are assumptions, not measured upper bounds. Check the
solve-count limit and absolute deadline before every native call. No partial
normalization is scored, and no wider fertility tolerance is allowed.

With ten workers per arm, three remaining initial waves and twelve DE waves
take at most 37.500 hours at that cap; one concurrent repeat wave plus a
15-minute export allowance gives **40.250 hours** of wave accounting. Propose
a **fresh 42-hour ceiling** to leave 1.750 hours for verification, scheduling
and bounded shutdown. This is a new budget for review, not an extension of
tonight's run or a promise that all numerical cases converge. At most 12,558
native stationary calls fit the finite case count. The objective caps account
for 1,365 CPU-hours before supervision overhead; reserving all forty CPUs for
the full 42 hours would
allocate 1,680 CPU-hours and 20,160 GiB-hours. Readiness has a separate
30-minute limit. A read-only check of partition `cs` reports `MaxTime=UNLIMITED`;
account/QOS admission and availability still require prelaunch checks. The exact arithmetic is in the
[budget proposal](../../output/model/utility_comparison/recovery_v1/budget_proposal.json).

A shorter fresh window must explicitly reduce the number of generations or
increase reviewed resources before launch. The current evidence cannot
support promising this full worst-case schedule in another eight hours.

## Concrete prototype and review boundary

The separate
[`e5f_utility_recovery_policy_v1.py`](../../code/model/tools/e5f_utility_recovery_policy_v1.py)
prototype and its focused fault-injection tests exercise classification,
process ownership, continued dispatch of distinct allowed candidates and
fatal stopping without importing or solving the model. It is not wired into
the frozen runner or controller and grants no launch permission.

Before deployment, lead review must approve the narrow rejection predicate,
censored-outcome semantics, combined availability guard and fresh budget.
Production integration still needs authenticated artifact reuse; typed DE
selection, scientific-key duplicate prevention and four-arm barriers; preserved full success validation; external
timeout ownership; and synthetic tests of those integration points. Tests
must demonstrate that missing/altered fingerprints, partial or late outputs,
smoke/repeat failures and unrelated kernel/accounting errors cannot become
accepted points. A new immutable snapshot and reviewed contract must bind
these changes. No numerical run or submission is authorized by this packet.
