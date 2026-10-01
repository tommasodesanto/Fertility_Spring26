# Housing-unit equation audit

Requested October 1, 2026: read-only unit-invariance audit and one-to-two-page advisor explanation.

- [Advisor note](housing_unit_normalization_note.pdf), updated with a worked dollar conversion, the distinct reference rents, and the supply coefficient consistent with benchmark population one. The note separates the exact accounting guarantee from empirical fit and recalibration.
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

## October 1 empirical supply-level check

[Price-level readout](price_level_readout.md) records fresh national2007 AHS price/room aggregates and a PSID gross-earnings money conversion, separately for adopted block0506 and experimental floor chain7 0173. Adopted implied rent is2.59% below the contract-rent stock ratio; owner price0.31% above the self-reported-value stock ratio. The external stock/rent anchor gives H0=6.06484 versus current6.29351. Experimental rent/value gaps are−11.31%/−8.67%. These are provisional point diagnostics, with income-date, quality and uncertainty limits stated; no model solve or target changes. [Comparison record](price_level_comparison.json), [AHS aggregates](ahs2007_price_quantities.json), [PSID age-cell earnings](psid2007_earnings_by_age.csv) and the retained AHS recipe provide the arithmetic evidence.
