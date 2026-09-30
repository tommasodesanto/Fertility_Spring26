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

## Short list for the author — September 29 evening

1. **Urgent credit decision:** remove the arbitrary renter age-42--62 repayment
   taper. The isolated patch retains zero new unsecured credit and nonnegative
   estates. Native compiled household validation is the next gate; positive new
   unsecured credit remains a separate, unchosen specification.
2. **Bequest measurement:** reconcile all positive model estates with the
   child-directed empirical target before treating the associated preference
   estimate as recipient-consistent. This is an already disclosed mismatch.
3. **Entry wealth/income:** validate the construction and sensitivity of the
   maintained joint distribution; no construction error is established here.
4. **Dependent-child departure and adult entry:** previously discussed at length
   with the author, not a new discovery. Retain the independent departure rule
   and separate 16/20-year entry queue as explicit approximations for reconsideration.
5. **Behavioral diagnostics:** high-wealth ownership declines, the age-30 housing
   downturn, retirement wealth profiles and conditional versus occupied policy
   interpretation remain questions, not certified findings.
6. **Transition:** the separate 2023 transition must be computed and validated
   before using it for policy responses. Stationary or prescribed-price results
   do not supply that missing path.

This list separates unresolved measurement, maintained assumptions and numerical
validation from confirmed coding errors; it is not a claim that all items are bugs.

## Renter rollover and primary-paper comparison — September 29

The tested rollover proposal is an experiment, not the DUE renter contract or an adopted baseline. Relative to **2007 stationary reference — block0506, September 28 verified export**, it removes the 42–62 taper and permits renters to choose \(b'\ge\min(b,0)\) before mortality risk. Thus full principal may remain outstanding, with interest paid through the budget; the floor becomes zero at positive mortality (first at age 66 here) and terminal decisions. The strict experiment instead imposes \(b'\ge0\) and nonnegative raw sale balances. No parameters were re-estimated, and neither experiment supplied a completed GE.

| Primary paper/version | Verified credit contract | Relation to renter rollover |
|---|---|---|
| [DUE, February 16, 2025](https://www.andrii-parkhomenko.com/files/Dynamic_Urban_Economics.pdf), §2.1.7, pp.11–12, equations 2.2–2.4 | Origination collateral constraint; incumbent owners need not delever after price declines. Renters cannot borrow. Substituting zero owned housing into the post-adjustment constraint requires a nonnegative balance after sale. | Protects incumbent owner debt; does not grant a former owner unsecured mortgage rollover. |
| [Boar–Gorea–Midrigan, author PDF linked to ReStud 2022](https://drive.google.com/file/d/1NBdnPJ-YW4MbHJ9R8Y1eQczP1fh_feAF/view), September 2021 PDF, §3.1, pp.9–14, equations 1–7 | Separate liquid account with \(a'\ge\bar a<0\) and mortgage balance. Thirty-year fixed-rate mortgage with minimum payments; LTV/PTI at origination. Selling into renting repays mortgage principal and interest and leaves no mortgage in the renter's next state. | Explicit bounded unsecured credit is permitted; carrying the former mortgage as unconstrained renter debt is not its rule. Numerical liquid-credit limit was not verified; do not import the 0.036 value from older two-author versions. |
| [Sommer–Sullivan–Verbrugge, published JME 2013](https://kamilasommer.net/RentPriceRatio.pdf), p.858, equations 3–4 | Interest-only/HELOC mortgage representation; no required principal reduction while retaining the property. Sale or becoming a renter requires full mortgage payoff. | Full principal rollover while owning is a maintained mortgage assumption; it does not persist as renter mortgage debt. |
| [Boar, Dynastic Precautionary Savings, December 2019 NBER version](https://www.nber.org/system/files/working_papers/w26635/w26635.pdf), pp.31,35,37 and footnote 58 | Bond-only model, no housing choice. Zero borrowing in baseline; age-22 entrants start with zero assets. Sensitivity allows working-age debt up to 18.5% of average annual income and zero retirement borrowing. | Provides zero-debt and explicitly bounded credit benchmarks, with a compatible entry contract; it cannot justify a housing-sale debt rule. |

Lead checked the DUE equations against the retained primary text and independently read the Boar PDFs after cheap extraction. The strict experiment's negative-debt renter entrants must be reconciled with the selected credit contract before a valid strict GE. An explicit unsecured-credit cap would be a separate economic choice, not a mechanical consequence of removing the age taper. No new code, calibration, or numerical experiment followed this literature check.
