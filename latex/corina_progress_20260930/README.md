# Corina progress presentation

Separate adviser update covering **September 17–October 1, 2026**, prepared for Corina
Boar. This dated presentation was explicitly requested and does not replace the
continuing [JMP Slides](../JMP_slides/JMP_slides.tex).

- [Presentation PDF](corina_progress.pdf)
- [Editable standalone Beamer source](corina_progress.tex)
- [Source and evidence bundle](corina_progress_source_bundle.zip)
- [Slide-by-slide provenance](source_map.md)
- [Verification receipt](verification.json)

## Verified scientific facts

The displayed calibration is the selected-repeat-verified parenthood-only
Stone–Geary point, overnight chain 7, case 0173_nm, collected October 1 at
09:41 New York. Its authoritative packet is
[ROOT/REPEAT evidence](../../output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/).
The latest author decision selects this preference specification; parameter
estimates remain provisional and optimizer convergence is not certified.

The packet contains all 14 fit rows: ten scored moments, three validation rows
and a birth-renewal check. Ten coordinates, including the child-benefit level,
are jointly free; housing supply scale is fixed at 6.293507689200028. The
selected search uses the `early4x` weighting profile: search objective
54.52095879194978, versus 31.284007255664957 under the original base weights.
These are two objectives at the same point, not an improvement comparison with
the old reference. Full fit contributions, parameter bounds and near-bound
flags are retained in the selected packet's ROOT tables.

Material utility uses a housing floor whenever children live at home, constant
consumption/housing shares, the power equivalence scale and nonlinear child
benefit. The compensated-share factor is inactive. The selected housing floor
is 2.3, exactly its upper search bound. Children ever born at age 25 are
0.5312175503600473 versus 0.8095276384290021 in the data. The selected packet
lacks the motherhood-share and children-per-mother observers; no earlier
point's decomposition is attributed to this winner.

The selected closure uses house price to clear birth renewal and population to
clear absolute housing supply. Entrant wealth uses nonnegative five-bin ratios
with positive nodes rescaled to preserve the original mean; the common annual
real interest rate is 2%, with no unsecured borrowing. These entry, credit,
closure and preference changes differ from the frozen September 28 block0506
reference. That reference and previous experiments remain historical evidence.

## Review and changes

The ten frames explain reviewed work and author decisions, the selected model,
its complete calibration, preference alternatives and remaining questions:

1. Model review and revisions
2. Household decisions
3. Children and housing needs
4. Demographic and market equilibrium
5. Housing around the first birth
6. Calibration strategy and inputs
7. Stone–Geary calibration
8. Calibrated parameters
9. Preference alternatives
10. Credit, initial wealth and remaining fit

The preference alternatives include compensated housing shares and the requested
equivalence-scale comparison. A requested comparison is not a completed estimate;
no new result is claimed without its own verified packet. The fit uses
`Moment | Target | Model`; the parameter frame reports all ten searched estimates
in three columns. Complete restrictions, bounds and bound flags stay in the
supporting CSVs. Visible source footers, internal failure history and a dedicated
computation section are omitted. Private provenance is in
[source_map.md](source_map.md).

The PSID figure and empirical provenance are retained. Earlier source/QA receipts
and their slide numbering describe preserved versions, not this rewritten deck.
Final verification checks the new equations, candidate identity, objectives,
full tables and every rendered page. No model, test suite, benchmark, cluster
job or job stop was run for this documentation update. The main JMP deck and
manuscripts are outside this task's edit scope.

## Preservation and build

The received source, PDF, source map, QA receipt and original source bundle are
preserved in `archive/received/`. `archive/transfer_manifest.json` records the
incoming hashes. The original task directory is also preserved:
`/Users/tommasodesanto/Documents/Codex/2026-09-30/task-4/corina_progress_20260930/`.
Its transient build products were moved out of this new active folder to the
review's temporary build directory.

Compile the source twice with `pdflatex -interaction=nonstopmode -halt-on-error`,
directing auxiliary files to a temporary build directory. The final document
must contain exactly ten pages and ten frames with no overlays. The final
verification receipt records compilation, rendered-page inspection, source
hashes, preservation checks and Library delivery availability.

The calibration slide places all three validation moments in a separate bottom
table; the ten scored moments and completed-fertility renewal check remain above.
All fourteen target/model pairs are unchanged.
