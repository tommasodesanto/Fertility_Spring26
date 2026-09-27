# Stationary calibration code integration

The author has held calibration launch pending a report of implemented/tested
coverage and further instructions. No parameter search or new household solve
was started in this implementation session. Provisional economic choices are
distinct from unfinished code validation.

## Implemented

- Native optional child-benefit curvature and compensated first-child housing
  shares, applied before the ordinary family-type compression. The legacy
  comparison wrapper cannot also be active. Defaults preserve old arrays.
- Optional one-market warm price search with bounds, at most eight proposals,
  the original refinement fallback and shared finalization. Fiscal certification
  precedes committing the warm state. The state resets for each parameter point.
- Explicit normalization first step, default 0.25; original solve cache retained.
  No saving-kernel changes or weakened numerical gates.
- Separate age-25 capped children-ever-born observer and AHS actual-room mean;
  existing capped-nine family-room definitions are unchanged.
- Explicit fourteen-row scoring registry and estate-funded positive entry ledger;
  debts remain separate and remaining positive net estates enter a sink. The
  net valuation and physical/financial counterparties remain provisional.
- Controller-owned deadlines, canonical recorded normalization inputs, complete
  target and parameter table repeats, frozen selected export, source/target pins,
  bounded proposal rounds, progress/checkpoint files, and a production-approval
  gate. Search is Gaussian/coordinate proposals around the best point, not DE.

## Verification and source scope

All execution occurred on Torch. The working snapshot is a copy of the frozen
September 25 source with only the reviewed integration edits. Current main has
other changes and was not substituted wholesale. The original frozen bundle is
untouched. The first staged copy is retained; the second copy preserves it by
replacing edited files rather than modifying shared hard-linked inodes.

The first component suite passed 71 tests after correcting one test's exact
comparison with a rounded decimal. The first runtime suite exposed a missing
staged estate-account dependency; staging the existing authenticated helper
resolved it. Four runtime/estate tests and thirteen controller tests then passed,
and job 18623058 passed the assembled zero-solve preflight. Its review-only
contract fingerprint was
`7b23ee7cc5acf7a79cb687376f52fc17e5c4b29cbfcbf167e3da578a9deb8ca9`.
Subsequent review fixes require the second snapshot's receipts; that first
contract is not launch authority.

Final preparation check: Torch job **18623087** passed **90 tests** (72 component
checks plus 18 runtime/controller checks) and the assembled zero-solve preflight.
The second snapshot is
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/calibration_code_integration_20260927_v2`.
Its review-only contract hash is
`399abb6e9eab0d447dca627f94e3de4e6a8006d8920d02241fadac21f6c6ebae`.
`component_tests.log`, `runtime_controller_preflight.log` and `preflight.json`
retain the receipts. `contract.json`, `objective.json` and
`target_provenance.json` are local copies for inspection; their executable paths
refer to Torch. No real objective, complete native repeat, speed benchmark or
search has been run with this integration. No job remains running.

The Claude review is retained in `claude_review.txt`. Lead adjudication:

1. Return the same JSON-safe receipt that is saved, avoiding a later raw-object
   serialization failure in the controller.
2. Remove competing child timeout checks; preserve the reviewed controller's
   ownership and fatal-integrity precedence.
3. Reject a nonpositive **final** normalized benefit with a dedicated type.
   Do not reject negative trial intercepts: they may bracket a positive root.
   An inconsistent solve count or failed pension account remains fatal.
4. Allow only the annual-discount-factor round trip a tight tolerance; repeated
   computed parameter tables still require exact agreement. This was already
   corrected during the parallel review.
5. Exclude administrative production blockers from the numerical identity, and
   require explicit production authorization plus an empty blocker list to run
   a search. Economic/target/source changes still invalidate the smoke identity.

The review's statement that the initial benefit came from a linear-benefit
diagnostic is incorrect: the retained parameter table records curvature 0.14
and exponent 0.86. Its claim that every loss difference is bounded by numerical
tolerances is also not established; full warm/cold comparisons remain required.

## Still required before a calibration launch

- Author's further instructions, including the proposed weight of 100 for early
  fertility and curvature bounds [0,0.8]. No existing target is dropped.
- Full warm/cold objective comparison, including all target/parameter rows and
  the stable seventeen diagnostics; no new-utility speedup is yet measured.
- Full independent repeated objective and exact production subprocess loop,
  then verification of checkpoints and selected export on real outputs.
- Size the bounded search from measured total objective time. The prepared
  ceiling is ten single-thread workers, at most thirty rounds of ten points,
  6.5 hours of search, one hour for repeats and half an hour for export.

The prepared specification has fourteen target rows: thirteen scored moments
and completed fertility 2.1 as a separate normalization; nine free structural
parameters plus the normalized child-benefit level. Counting does not establish
local identification. Tomorrow's review covers tenure scale, interest/credit,
rental support, income/entrant wealth, estate and wealth measurement, weights
and identification, plus the separately deferred owner-grid and conception
schedule tests.
