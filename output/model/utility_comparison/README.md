# Four-arm utility comparison

The author delegated tonight's experimental choice; the reviewed fixed-rent
utility normalization is chosen for the comparison. The design and assumptions are in
[`docs/model/utility_four_arm_preparation.md`](../../../docs/model/utility_four_arm_preparation.md).

`preparation/` retains the first adapter-only check. The integrated
`preparation_v2/validation_receipt.json` records Torch job `18567413`: 54 tests,
four native arm preflights, two structural seeds per arm, zero equilibrium
solves. `preparation_v2/budget.json` records the eight-hour, 652-objective
finite design. Its `contract.json` is the checked preparation snapshot and
deliberately refuses objectives; a separate approved launch contract is
required. Subsequent source edits require a new snapshot, preserving old pins.

Remote preparation root:
`/scratch/td2248/projects/Fertility_Spring26_native_financing_20260919a/utility_four_arm_preparation_20260925_v2/`.
The current checked packet is `checks_v2/`. Large checkpoints remain remote.

`pdf_layout_validation/` contains a visibly labeled historical layout fixture,
its rendered pages and dependency receipt. This packet is not a new fit or
calibration result. The first fixture attempt (`18567452`) stopped because
ReportLab was missing; its partial outputs are retained. Only the reporting
environment is corrected before a fresh layout check.

All targets, weights and measurement definitions inherit the immutable
September25 nightpair. The common adopted fiscal change and experimental
preference changes are separately disclosed in the design. Pending estate,
entry-wealth and mortality proposals are not inserted into these arms.
