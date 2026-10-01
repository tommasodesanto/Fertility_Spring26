# Housing-unit equation audit

Requested October 1, 2026: read-only unit-invariance audit and one-to-two-page advisor explanation.

- [Two-page advisor note](housing_unit_normalization_note.pdf).
- [Editable LaTeX](../../../../latex/housing_unit_normalization_note.tex).
- [Detailed source map and proof](../../../../docs/model/housing_unit_normalization_audit.md).
- [Standalone arithmetic record](arithmetic.json): 1,260 scalar comparisons; maximum relative discrepancy 1.911e-14. Six inspected engine files match experiment source pins.
- `arithmetic_recipe.py`: reproduction of scalar arithmetic using saved JSON/CSV and standard-library math; no model imports or solves.

Reference identity: September 28 adopted block0506 shares reference. Separately evaluated input contract: October 1 09:41 NY experimental floor chain7 case0173_nm, D=0, not adopted. D=.25 and D=.53 cases are not used.

Established: economic-equation equivalence under a complete factor-of-ten conversion. Not certified: a rescaled numerical equilibrium, policy arrays, historical conversion choices or empirical validity of demand/supply elasticities and financial levels. No model code or calibration data changed.

PDF export: `pdflatex -interaction=nonstopmode -halt-on-error -output-directory=output/model/fixed_reference_economics_20260928/room_unit_audit_v1 latex/housing_unit_normalization_note.tex` (twice). Built-in LaTeX compiler diagnostics are also recorded during creation. Only this documentation build is run; no model tests or jobs.

## External review

- [Full Claude Fable 5.1 Max review](claude_full_review.md), completed after explicit user authorization.
- [Execution receipt](claude_review_receipt.json).
- [Lead verification and remaining limits](claude_review_response.md).

The note was tightened following review; the arithmetic recipe remains an illustrative equation audit, with the renter closed-form coverage and historical-conversion checks outstanding.
