**Review update, September 6 at 23:15 UTC.** The full overnight search
remains stopped for the author's two-hour review. Rebuilt smoke `17076426`
completed all four certified histories. Both anchors exactly reproduce the
preceding scientific results: twelve fit rows, parameters, 253 numeric history
entries and seventeen PNGs. The exact original-population gate replay passed.
A separate policy-reader schema error (`calendar_year` versus the historical
writer's `period` and `years_from_start`) was found independently by the lead
and a bounded reviewer and corrected. The reviewer verified the original
population and 2019-end-queue handoff and the stated closure/copying semantics.

Policy-stage job `17079223` then passed the baseline at 2023 and 2027 but
stopped at the supply expansion's 2023 budget audit. Violating mass was
`1.9390478108e-7`, exceeding the unchanged `2e-10` gate. Source inspection and
state ledgers identify the previously documented owner consumption-reporting
floor: the optimizer uses positive budget-feasible consumption but its output
raises some values to 0.04. One state reports 0.04 while its budget supports
0.0153663. This is not a new shock-law or housing-market nonexistence result.
The baseline's maximum market residual is `2.413e-5`; only that two-date policy
path is certified so far. LTV and tax paths have not yet run.

A narrow isolated repair reconstructs owner consumption from the unchanged
saving choice and budget only on feasible solved branches, in joint mode.
It does not change the optimizer, housing, saving, transaction maps or gates.
Twenty local tests pass. Diagnostic job `17080030` makes two single-price
replays, before and after this repair, and requires identical values, choices,
distributions, prices, quantities and all seventeen standard plots, while the
budget audit passes. No full smoke or search is implied by that diagnostic.
New exclusive snapshot `Fertility_Spring26_joint_nested_full_20260906e` has
scientific bundle `85450db0d7611f7206fba933a74c0f962c18990917a926f9bb2888057494ff39`;
contract SHA `59fb3a15911cf8692b30402d3d867b145c0c851479cac59f0fef720194b4e8c4`.
The prior policy-only snapshot d and its failure receipts remain preserved.
A subsequent bounded full-loop verification, if this diagnostic passes, must
finish within the review window. Full calibration remains held until discussion.
