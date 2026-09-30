# Corina progress presentation

Separate adviser update covering **September 17–30, 2026**, prepared for Corina
Boar. This dated presentation was explicitly requested and does not replace the
continuing [JMP Slides](../JMP_slides/JMP_slides.tex).

- [Presentation PDF](corina_progress.pdf)
- [Editable standalone Beamer source](corina_progress.tex)
- [Source and evidence bundle](corina_progress_source_bundle.zip)
- [Slide-by-slide provenance](source_map.md)
- [Verification receipt](verification.json)

## Verified scientific facts

The reference is **2007 stationary reference — block0506, September 28 verified
export**. The complete [reference fit](evidence/reference_target_fit.csv) contains
14 rows: ten scored moments, three validation rows and the separate completed-
fertility normalization. The [31-row parameter table](evidence/reference_parameters.csv)
includes estimates, restrictions, bounds and bound flags. The objective is
19.581310760138322. Full corresponding tables for the later one-birth and
two-birth experiments remain in `evidence/`; neither experiment is adopted.

The preference equations distinguish children at home from lifetime births.
The code stores the one-child benefit directly in `psi_child`, giving
\(v(m)=\psi m^{1-\gamma}\), with zero benefit at \(m=0\). The first-child share
change lowers the consumption share, and its material-utility compensation uses
a fixed reference rent. Dependency departure and the 16/20-year adult-entry
queue are separate objects.

The young-age miss chiefly concerns children among mothers in the experimental
comparison. The weak curvature direction is local to a search-round center;
it is not a statistical nonidentification result. The completed September 29
price experiment shows that relaxed credit raises fertility levels without
attenuating the local response to house prices and mapped rents.

The matched grid comparison keeps diagnostic debt allowance \(D=0.53\) in both
arms and takes 468.606 versus 187.254 seconds. Small aggregate differences do
not certify local policy or transition accuracy. The three subsequent pilots
use \(D=0.25,0,0\), different entrant-wealth rules and a common 2% annual real
interest rate. Their dated launch/status evidence is separate from completed
scientific results. Historical shock fits yielded no estimates; diagnostic
18820811 was cancelled by the author, as recorded in the incoming handoff.

## Review and changes

The deck now opens with household preferences and demographic accounting, then
proceeds through the empirical evidence, reference fit, model diagnostics and
next calibration. The reading order is: Household preferences; Demographic
accounting; Housing around the first birth; Calibration inputs and empirical
counterparts; 2007 stationary reference fit; Young mothers and the fertility
intensive margin; Mortgage access and unsecured liquidity; House prices, credit,
and births; Stationary computation and grid resolution; Entry wealth, credit,
and the next calibration. The former general overview frame was removed to keep
the deck at ten frames. The preferences and accounting material formerly shared
one frame and now open the deck as two separate frames.

An independent read-only numerical review authenticated all 15 copied source
artifacts, checked all 14 reference-fit rows and all chart coefficients and
intervals, and verified the other displayed numbers. The lead checked the
economic interpretation and preference implementation. The sources cover early
September 17–19 mechanism work, September 23–28 input/specification revisions,
the September 24 empirical rebuild, September 28–29 fertility experiments, and
September 29–30 price, computation and entry work.

The review found a dated status update and three presentation issues: crowded
negative event-time labels, tight table-column spacing, and the debt allowance
appearing before its units were defined. The review's references to slides 4, 6,
9 and 10 use the pre-reordering numbering; the corresponding current frames are
3, 5, 9 and 10. The 14-row fit remains in the author's three-column format.

No model, calibration, test suite, benchmark, new cluster job, or job stop was
run for this review. The main deck and both manuscripts retain their existing
contents; this update does not certify synchronization of those documents.

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
