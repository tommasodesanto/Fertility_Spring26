# Bounded completion of the unattempted global points

This is a derivative of the pinned one-birth Estate-A 64-point Sobol plan,
SHA-256 `44ffc96620deab7577b8ecaaa744e41782b9294e643d04d2d5bafa7a7ce7e01b`.
It changes only the stage's case-selection and numerical-rejection handling.
The original scrambled points (seed 20261004), linear/log transforms, source
inventory, 14-row target and weight contract, ten free bounds, 32-lifecycle
cap, and native selected/repeat gates remain unchanged. No new point is drawn.

The terminal collection of array 19183972 shows 13 checkpointed points and
15 fatal attempted points. This stage selects precisely the 36 **never
attempted** indices, excluding both groups. Nine one-core tasks each own four
specific indices:

| Task | Original Sobol indices |
| ---: | --- |
| 0 | 1, 2, 3, 5 |
| 1 | 6, 7, 9, 10 |
| 2 | 11, 13, 14, 15 |
| 3 | 18, 19, 22, 23 |
| 4 | 29, 30, 31, 34 |
| 5 | 35, 38, 39, 42 |
| 6 | 43, 46, 47, 49 |
| 7 | 50, 51, 53, 54 |
| 8 | 55, 61, 62, 63 |

The original attempted-fatal indices were 0, 4, 8, 12, 17, 21, 28, 33, 37,
41, 45, 48, 52, 59, 60. The stage will assert the complete index partition
against the preserved collection before it submits.

The only numerical orchestration repair is classifying recognized
`uncomputed_price_unbracketed` and typed inherited-distribution infeasibility
as rejected cases with no finite valid loss, then continuing the independent
points. Unknown exceptions and source/contract drift remain fatal. Each case
gets at most 20 minutes; each task has a 90-minute wall and checkpoints after
every case, plus latest and best summaries. The entire stage is one array
`0-8%9`, with no retries, extension, or replacement points. A two-case
zero-solve exact-loop mock and exception-classification test must pass before
submission. The native final acceptance and exact-repeat checks are unchanged.

Production array **19196205** was submitted exactly once at 11:27 PM New York,
with plan SHA-256
`5c738c7ee4c1c331eea2511fb8288a24af8d53acf02aa13ae0c80b3d8d9ce5a9`
and stage-manifest SHA-256
`79f82cc2f35a40af02359d92bfd594a634a43e2f21c9fffbe261ef065ff7aa79`.
The two-case exact-loop launcher smoke passed with Sobol indices 1 and 2,
zero lifecycle solves and exit 0; recognized-error and source-pin checks also
passed. The [saved submission receipt](../../../../output/model/experiments/birth_count_choice/estate_a_global_completion_20261005_v1/deployment/production_submission.json)
is retained with the complete stage plan, manifest, and source pins.

The production array depends on **both** already-running overnight arrays
19194495 and 19194496 reaching terminal state (`afterany`). This prevents
contention and does not change either overnight search's own deadline, call
cap, or recovery rights. The independent global-completion cutoff is October 5,
2026, 2:00 PM New York (18:00 UTC). If dependency release leaves too little
time or the stage ends incomplete, save the completed and remaining index
lists and stop; there is no automatic budget expansion. Until results pass
fresh native selected-point and exact-repeat checks, any best loss is
provisional. A subsequent refinement requires separate review of feasible,
distinct seeds and its own explicit finite budget.
