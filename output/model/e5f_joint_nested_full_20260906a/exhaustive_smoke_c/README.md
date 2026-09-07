# Exhaustive-saving correctness verification

The current full-loop job is 17087058 in the frozen Torch snapshot
`/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907c`.
It is a four-history/eight-policy-date verification, not a calibration search.
Root CALIBRATION_STATUS.md owns current status. No dependent long job exists.

- `contract.json` pins source, targets, both seeds, loop budgets and gates.
- `saving_source_manifest.json` records the 55 checked scientific/helper files.
- `integration_verification.json` links the independent review and prior fixed-price proof.
- The completed fixed-price comparison is `../exhaustive_smoke/saving_integration/`.
- The review is `../saving_integration_review.md`.
- The previous discussion PDF remains unchanged and does not report this new smoke.

The experimental source is committed as 3e48fe3 on codex/joint-nested-full.
Twenty-six local checks pass; compiled checks and a second-process cache reload
run before the cluster model. The default-off old model reproduces ten arrays
exactly. The new integrated optimizer matches the audited independent full
lifecycle optimizer on 15 of 16 arrays exactly; renter housing differs only
by 1.776e-15. All 17 standard PNGs match the already visually audited packet.

The provisional complete-history planning estimate is 1,600 seconds, two cases
at once. The whole smoke has a 90-minute limit, including four short policy
branches. Its measured runtime must determine any subsequent parallel budget.
Do not interpret the inherited 360-case maximum as a promise that the full
search fits in the overnight wall-clock cap.

Inspect queue and progress:

```sh
ssh torch 'squeue -j 17087058'
ssh torch 'cat /scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907c/output/model/joint_nested_overnight/smoke/search_state.json'
```

Collect without altering the frozen source:

```sh
rsync -a torch:/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260907c/output/model/joint_nested_overnight/ output/model/e5f_joint_nested_full_20260906a/exhaustive_smoke_c/
```

Verify every original case receipt and artifact hash with the experimental
adapter, compare both anchors using compare_reference, require all 17 PNGs to
match exactly, inspect all policy-date packets, then require the complete
policy_loop_verification.json and equilibrium_receipt.json. A partial smoke
never authorizes the larger search. The existing finite monitor owns collection,
verification and the concise readiness update; production remains unchanged.
