# Borrowing resolution evidence receipt

**Scope:** read-only source and decision audit, 2026-09-26. No model imports, solves, tests, empirical runs, frozen-source edits, canonical updates, or manuscript edits.

**Conclusion:** default code implements renter unsecured rollover and an owner collateral floor plus rollover of debt in excess of collateral. `lambda_d=0` removes new unsecured credit but not inherited debt. Current income enters purchase eligibility once in the isolated adapter retained by the Sept. 25 frozen pair. A falling price does not force an owner to pay down the existing balance; selling transfers net proceeds and a residual shortfall to renter wealth. Default owner stayers can reborrow against available collateral, so “no refinancing” is false for this code path. The age taper's economic rationale and terminal creditor settlement are not established. Keep all frozen runs unchanged.

## Source evidence

- Active formula and age grid: `code/model/intergen_eqscale_seq_optimized/parameters.py:718-755`, `solver.py:141-170`; defaults at `parameters.py:171-174`.
- Cash/downpayment, secured debt floor, net sale equity, owner cost: `solver.py:2579-2598`. Forward sale-shortfall map and purchase wealth: `solver.py:4085-4148`. Tenure gate source: `kernels.py:719-760`.
- Current-income adapter and guards: `code/model/tools/e5f_earnings_wealth_contract.py:243-280`. It shifts the purchase-stage down payment and collateral gate by (y/R), leaves transaction wealth unchanged, and relies on the flow budget to count income once.
- Frozen Sept. 25 pair: `.../target_review_v1/payroll_tax_review/overnight_pair_20260925_preparation/README.md` says current-income eligibility was retained; its frozen source/result pointers and hashes are listed there. Sept. 26 four-arm frozen contract pointer and SHA are in `CALIBRATION_STATUS.md:8-21`. The source snapshot `tmp/utility_overnight_20260923_v1/source` includes the adapter, while its generated runtime diff is separate. No frozen bundle was opened for mutation.
- Adopted direction versus recommendation: `docs/model/ACTIVE_DECISION_LEDGER.md:78-80, 157-160` records current-income financing as author-agreed and full amortizing loan states as an explicit extension. `docs/model/POST_PRESENTATION_ISSUES.md:346-348` records the existing taper but not its rationale. `docs/model/structural_model_review_chatgpt_20260915.md:188-208` is a dated recommendation, not author adoption: distinguish origination from incumbent debt, replace or explain taper, and identify sale-shortfall/terminal creditors.
- Estate dependency received from the separate estate owner: `code/model/tools/e5f_estate_receiver_adapter.py:147ff` clips negative net donor wealth at zero in the recipient case, while its test does not certify creditor/debt writeoff accounting. The estate task owns the recipient rule and negative-estate measurement; this task does not infer costless debt cancellation.
- Prior independent source audit: `docs/model/e5f_full_code_correctness_efficiency_review_20260905.md:116` reports the cash test and collateral boundary agree algebraically, including mover sale proceeds, and checks both sides of the purchase threshold.

## Hand arithmetic independently checked

**Adapter limitation:** e5f_earnings_wealth_contract.py:263-265 replaces owner-kernel rollover inputs with zero, while kernels.py:1238-1244 uses them to preserve debt beyond the collateral floor. This prevents new-purchase temporary debt from rolling forward but also removes incumbent-underwater owner rollover after a price decline. It is unverified outside the frozen stationary use; future price-transition/inherited-underwater use needs a transaction-specific floor and separate stayer/origination tests.

At (P=1,h=4,\phi=.8), required cash is .8 and secured floor is (-3.2). If the owner holds (b=-3.2) and price falls to .75, the floor becomes (-2.4), excess unsecured debt is (-.8); at (s=1,D=0), the floor remains (-3.2), so there is no margin call. A sale nets ((1-.06)(.75)(4)=2.82), leaving renter wealth (-.38). At age weight .5, the same residual-debt floor is (-.19). With (b=.5) and (y/R=.3), current income exactly bridges the .8 down payment; the purchase adapter admits it without changing (x=b-4=-3.5), and (Rx+y=R(-3.2)), so income is counted once.

## Primary literature checked

- Kaplan, Mitman, and Violante (2020), published-paper PDF: [The Housing Boom and Bust: Model Meets Evidence](https://violante.economics.princeton.edu/sites/g/files/toruqf5621/files/documents/kaplan-et-al-2020-the-housing-boom-and-bust-model-meets-evidence.pdf). Used only for the distinction between origination LTV restrictions and longer-term debt/default/refinance institutions; it does not identify this model's age taper.
- Hurst and Stafford (2004), author-hosted paper: [Home Is Where the Equity Is](https://erikhurst.com/wp-content/uploads/2020/02/hurst_stafford_jmbc_final.pdf). Supports treating equity extraction/refinancing as an economically meaningful choice.
- Campbell and Cocco (2015), *Journal of Finance*: [A Model of Mortgage Default](https://onlinelibrary.wiley.com/doi/10.1111/jofi.12252). Example of a fuller mortgage/default framework, not a required minimum for this project.

**Validation:** no tests or solver were run by design. The accompanying report specifies the smallest post-decision accounting tests and a Sunday completion criterion.
