# Utility comparison recovery: eight-hour option (plan revision 2)

**Preparation only; no launch is approved or implemented.** The main research
task requested this separate recovery design after all four searches in
`run_001` stopped during the initial population. It does not alter array
18567879, its source snapshot, contract, records or eight-hour clock. Those
controllers finish their existing repetitions and exports. The experimental
utility choices and adopted common pension update remain exactly as disclosed
in [the original design](utility_four_arm_preparation.md).

The recommended option for review is a **fresh eight-hour limit**, including
verification, repetitions, reports and visual review. It starts from three
retained feasible points per arm and allows **two rounds of twenty new points
per arm**, with paired local and global moves. It is a smaller search with a
clear coverage tradeoff, not completion of the original 546-objective recovery
scenario. That scenario's 42-hour calculation is retained only as an
**unapproved upper-bound illustration**, not the default proposal or an
authorized extension. No failed or timed-out point is retried.

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
4. Retain the original operational stop above 50% inadmissible proposals,
   applied prospectively to each twenty-proposal round. This is a conservative
   controller choice, not a scientific acceptance condition. The additional
   combined rejection-plus-timeout threshold proposed for the 42-hour scenario
   is not proposed here. Timeouts remain separately visible and unscored.
5. Use typed outcomes and choose adaptive centers only from successful scored
   cases. A timeout is not a rejected scientific evaluation or evidence of a
   worse loss. Unexplained missing/unrun results cannot silently cross the
   common round barrier. This small design uses direct local/global moves,
   not three generations of DE or an implicit claim that those ran.
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

## Retained points and two matched rounds

The [bounded inventory](../../output/model/utility_comparison/recovery_v1/observed_reuse_inventory.json)
records existing controller outcomes and receipt hashes. It is preliminary:
the share controllers were still finishing, and this inventory does not
independently rehash large checkpoints. A launchable reuse manifest must be
made only after the existing run is terminal and checked against the original
approved contract
`c3fa4d5b7925a54a747c72511538b8182029e1493e4d02e5ef703a30fcecc28a`.

Reuse the four exact smoke pairs and the common feasible initial slots
`0000`, `0001` and `0002`: these are successful in all four arms in the bounded
inventory. Keep the additional completed share-arm points as cumulative
evidence, but exclude those unmatched points from the primary matched-start
selection and adaptive centers. This avoids giving the share arms a larger
starting search in the primary comparison. Authenticate original receipts, checkpoint bytes, complete
target and parameter tables, arm/point identity and all scientific pins on
Torch. Record the old contract as ancestry in a distinct recovery contract;
never rewrite an old receipt with a new fingerprint. New controller code must
not change the economic source, target system, observers or numerical gates.
If that identity cannot be established, reuse is not approved. No extra anchor
solve is needed merely because orchestration changed.

Each round contains twelve local moves and eight global draws per arm. In
round one, take two opposite pairs of local perturbations around each of the
three common retained anchor sets. In round two, take six opposite pairs
around each arm's best eligible point from that common bank and round one.
Use the existing transformed coordinates, reflected original bounds and local
standard deviation 0.04. Draw eight additional full-range points in each
round. Pin the complete random streams before launch, with proposed seed
20260926. Common coordinates use common draws; floor and share coordinates
use their respective matched family streams.

Round-two centers may differ across arms. That is limited recalibration;
paired perturbations do not mean identical physical parameter vectors. The
common retained anchors separately show conditional utility comparisons.
Each arm has at most forty new points and the same two rounds. At most two
selected repetitions per arm bring the total to **168 new objectives**,
including 160 search points. Reuse existing exact repetition evidence if the
selected original is unchanged. Failures and duplicates consume their slots;
there are no replacement draws, extra rounds or automatic retries.

## Eight hours, including verification and reporting

Keep the existing normalization algorithm, target 2.1, tolerance 0.0005 and
maximum 23 stationary solves. Propose a **4,200-second candidate cap**. In the
saved timing inventory, successful objectives used four to thirteen native
solves and at most 3,018.968 seconds. The diagnosed 3,100-second timeouts had
completed ten or eleven solves; the additional 1,100 seconds allows several
more solves at their observed rates. This is a computational allowance, not
a guarantee: successful-only timing omits the unresolved tail. A case that
still exceeds the cap is censored, never scored from a partial normalization.
Check the inherited solve-count limit and absolute deadline before each call.

With ten workers per arm, each twenty-point round needs two waves. The four
search waves account for 280 minutes at the cap. Use one absolute clock,
beginning when the first recovery controller is ready:

| Window after start | Work and cutoff |
| --- | --- |
| 0–30 minutes | Four-arm readiness and authenticated reuse verification; abort if incomplete. |
| 30 minutes–6 hours | Two rounds; 280 minutes of capped waves plus 50 minutes for coordination and shutdown. |
| 6 hours–7 hours 10 minutes | Both selected repeats concurrently within each arm, if needed. |
| 7 hours 10 minutes–8 hours | Full reports, collection and qualified visual review. |

Readiness is included in the eight hours. All stage and case deadlines are
absolute, with no extension. Early completion creates slack, not permission
for an extra round. If verification, valid-case coverage or reports remain
incomplete, report that explicitly at the cutoff.

The maximum is 3,864 new stationary calls from the finite case count, not a
forecast. Candidate caps account for **196 CPU-hours** before supervision;
reserving forty CPUs for all eight hours would allocate **320 CPU-hours** and
3,840 GiB-hours at the original 120 GiB per arm. This buys a matched small
search, not 546 objectives or full reoptimization. The machine-readable
[eight-hour option](../../output/model/utility_comparison/recovery_v1/eight_hour_option.json)
records the timing evidence, counts and arithmetic. The earlier
[42-hour calculation](../../output/model/utility_comparison/recovery_v1/budget_proposal.json)
is historical, unapproved and superseded as the default proposal.

## Scientific requirements and controller choices

Scientific acceptance still requires the same economic specification,
targets, weights, observers, fertility target/tolerance, fiscal and numerical
gates, and source/target/point/checkpoint integrity. All full target-fit and
parameter-bound tables remain required. No deadline permits accepting a
failed gate or an incomplete normalization.

The eight-hour ceiling, 4,200-second cap, number and mix of moves, shared
barriers, candidate-local failure policy and 50% rejection stop are controller
choices. The inherited 23-call guard is an algorithmic safeguard, not an
economic target. Reusing verified smoke/repeat evidence is an execution
choice. Changes to these choices need a new reviewed contract, but need not
change scientific acceptance. The eight-hour option explicitly changes
search coverage and orchestration while preserving the science.

## Concrete prototype and review boundary

The separate
[`e5f_utility_recovery_policy_v1.py`](../../code/model/tools/e5f_utility_recovery_policy_v1.py)
prototype and its focused fault-injection tests exercise classification,
process ownership, continued dispatch of distinct allowed candidates and
fatal stopping without importing or solving the model. It is not wired into
the frozen runner or controller and grants no launch permission.

Before deployment, lead review must approve the narrow rejection predicate,
censored-outcome semantics and this smaller eight-hour design. Production
integration still needs authenticated artifact reuse, the finite paired-move
bank, typed round selection, scientific-key duplicate prevention and four-arm
barriers; preserved full success validation; external timeout ownership; and
synthetic tests of those integration points. Tests
must demonstrate that missing/altered fingerprints, partial or late outputs,
smoke/repeat failures and unrelated kernel/accounting errors cannot become
accepted points. A new immutable snapshot and reviewed contract must bind
these changes. No numerical run or submission is authorized by this packet.
This revision changes only the plan and cost documentation. It adds no agents,
model work or heavy tests; the existing monitor and result collection continue.
