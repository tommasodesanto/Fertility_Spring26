# Lead review of specification disclosure

Reference: **2007 stationary reference — block0506, September 28 verified export**.

Two cheap bounded passes are preserved in `audit_findings.md` and
`coverage_findings.md`. They inspected active flags and code across household
preferences, housing/tenure, credit, earnings/entry, birth/dependency timing,
survival/estates, fiscal rules, renewal, supply and measurements. This is a
specification inventory, not a numerical certification or a proof of no bugs.

## Lead corrections to raw worker findings

- The native renter taper applies to renter debt. The first table's suggestion
  that it governs owners' “unsecured excess” imports a legacy rule into the
  native branch and is incorrect. Native buyers and incumbent owners have
  separate bounds; the source checks in the parent README remain authoritative.
- The second table's household-conversion row incorrectly calls the serialized
  `entrant_conversion_factor=0.5` operative. The active split-birth queue uses
  adjusted births **divided by 2.1 once**; its factor 0.5 splits the result
  equally between 16- and 20-year entry queues. It is not another division by
  two. Lead verified `adult_entry.py:32–37, 57–73` and `solver.py:6126–6129`.
- September 14's displayed 17.9% tax belongs to a historical model. The pension
  rule changed subsequently; comparison to today's 8.028% is not evidence of
  a hidden feature in the original presentation. Preserve historical slides.
- The all-positive-estate versus child-directed bequest-target mismatch is
  already explicit in the reference receipt and advisor Doc. It remains a
  measurement problem, not another newly discovered undisclosed restriction.
- Feasibility sentinels are numerical conventions. Their existence alone is
  not evidence of economically material contamination.

## What this establishes

The renter repayment taper was omitted from the actual September 14 deck and
advisor Doc, and its resolution was not tracked to an explicit author decision.
The author now selects removal; isolated implementation is in `../credit_no_taper_v1/`.
Neither cheap pass established another hidden active restriction comparable
to that taper. Several disclosed approximations remain unresolved, especially
recipient-consistent bequest measurement, the entry wealth/income construction,
and the separate full transition. These must not be represented as settled.

Remaining audit limits: the empirical construction of the fixed 160×15 entry
array was not rederived; the separate dated transition was not verified; no
occupied-state numerical extraction, solve, grid sensitivity or causal test
ran. Additional unknowns are not automatically bugs. Positive unsecured credit
remains an open author choice, separate from taper removal.

The audit compares active economic rules and disclosure. It does not mean that
every source line or all 1,241 pinned files were independently reviewed.
