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
\(\psi m^{1-\gamma}\), written directly inside utility on the slide, with zero benefit at \(m=0\). The first-child share
change lowers the consumption share, and its material-utility compensation uses
a fixed reference rent. Dependency departure and the 16/20-year adult-entry
queue are separate objects.

The reference price clears aggregate housing demand against elastic supply,
with household population normalized to one. The child-benefit level is
separately normalized to completed fertility 2.1. Equal retiree pensions
balance PAYGO. This is distinct from the later experimental closure that uses
price for birth renewal and population for absolute housing supply.

The young-age miss chiefly concerns children among mothers in the experimental
comparison. The weak curvature direction is local to a search-round center;
it is not a statistical nonidentification result. The three entrant-wealth
pilots use \(D=0.25,0,0\), different entrant-wealth rules and a common 2% annual
real interest rate. All three completed and remain unadopted; the deck reports
their designs and completion only. The authoritative completion record is
[terminal results](../../output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/RESULTS.md).

## Review and changes

The latest author request calls for a matter-of-fact adviser narrative: reviewed
work and changed decisions, a detailed explanation of the model and calibration,
and brief selected problems. The ten-frame order is:

1. Model review and revisions
2. Household decisions
3. Children and housing demand
4. Demographic and market equilibrium
5. Housing around the first birth
6. Calibration strategy and inputs
7. 2007 stationary calibration
8. Calibrated parameters
9. Fertility at young ages
10. Credit and further calibration

The fit retains all fourteen rows in three columns. The parameter frame reports
all ten searched estimates and the normalized child-benefit level; complete
bounds, restrictions and near-bound flags remain in the supporting CSV. Visible
source footers, the early sandbox table, price-elasticity table, computation
section and internal failure history are omitted. Private provenance remains
in [source_map.md](source_map.md), and prior evidence and receipts are retained.

An independent read-only numerical review authenticated all 15 copied source
artifacts, checked all 14 reference-fit rows and all chart coefficients and
intervals, and verified the other displayed numbers. The lead checked the
economic interpretation and preference implementation. The sources cover early
September 17–19 mechanism work, September 23–28 input/specification revisions,
the September 24 empirical rebuild, September 28–29 fertility experiments, and
September 29–30 price, computation and entry work.

Earlier reviews and their slide numbers describe the preserved versions; they
should not be read as the numbering of this rewritten deck. The final review
must check the new model equations, full fit and parameter frames, and every
rendered page against the rewritten source.

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
