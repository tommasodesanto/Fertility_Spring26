# Corina adviser progress deck, September 17–30, 2026

This separate standalone presentation for Corina Boar has ten frames, with no
overlays. It follows the author's latest request for reviewed work and changed
decisions, detailed model and calibration exposition, and brief selected
problems. There are no visible `Source:` footers or dedicated computation
section. This file records private provenance. Paths below are relative to
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/`.

## Reference and interpretation

The scientific reference is **2007 stationary reference — block0506, September
28 verified export**, authenticated in
`output/model/fertility_identification_20260928/fixed_reference_manifest.json`
and its README. Its unchanged fourteen-row fit and thirty-one-row parameter
table are copied to `evidence/`. The weighted objective is 19.581310760138322.
The later E01 one-birth and E02 two-birth selected experiments are distinct,
unadopted candidates, with objectives 7.826226594410982 and 7.8420175378092205.
Never compare these objectives with September 17–19 sandbox objectives.

In the retained reference, house price clears aggregate housing demand against
elastic housing supply and household population is normalized to one. The
calibration separately adjusts the child-benefit level to completed fertility
2.1. Equal retiree pensions balance PAYGO. The later entrant-wealth pilots instead
use price to clear birth renewal and population to clear absolute housing supply;
that closure must not be assigned to the September 28 reference.

## Slide-by-slide sources

1. **Model review and revisions.** `CALIBRATION_STATUS.md`, the September 23–28
   accepted decisions, `docs/model/accepted_input_reconciliation_20260926.md`,
   and `docs/model/calibration_identification_review_20260929.md`. Early work
   remains documented by `docs/model/situation_report_20260918.md`. Distinguish
   adopted changes from provisional accounting and experimental calibrations.

2. **Household decisions.** The reference manifest's
   `actual_serialized_parameters` and
   `code/model/intergen_eqscale_seq_optimized/solver.py:278–295,3349,3486–3551`.
   Households differ in age, wealth, persistent earnings, tenure, children ever
   born and children at home. They choose birth attempts, housing/tenure,
   consumption and saving. Earnings are stochastic endowments, not a labor-supply
   choice. The reference permits one birth opportunity per four-year cell;
   success moves (n,m) to (n+1,m+1). The two-birth option belongs only to E02.
   Numeric period lengths belong in calibration, rather than model exposition.

3. **Children and housing demand.**
   `code/model/intergen_eqscale_seq_optimized/child_preferences.py:11–14,27–51`
   and `solver.py:2579–2607`; reference `parameters.csv` and manifest.
   Material utility is CRRA of A(m)c^alpha(m)s^(1-alpha(m))/e(m), with
   e(m)=((2+0.7m)/2)^0.7 and alpha(m)=alpha0−Delta_alpha*1{m>0}.
   Direct child benefit is psi*m^(1−gamma), zero at m=0, written directly inside the utility equation on the slide. The stored `psi_child`
   is the one-child benefit, not the equivalent CRRA coefficient.
   A(m)=K(alpha0,r*)/K(alpha(m),r*) with
   K(a,r)=a^a*((1−a)/r)^(1−a). Compensation holds optimized material utility
   constant when the expenditure share changes at a fixed reference rent,
   conditional on expenditure and equivalence scale. It does not compensate
   the cost of children or preserve utility at every endogenous rent.
   Bequest utility and the first-birth fixed utility cost are separate objects.

4. **Demographic and market equilibrium.** Reference
   `resume_v1/selected_export/primary/receipt.json` fields `fiscal_rule`,
   `normalization`, `adult_entry_gate` and `estate_funding`;
   `solver.py:1478–1492,2074–2086`; September 24/26 decisions in
   `CALIBRATION_STATUS.md`. Housing supply is H0*(r/r_bar)^xi, with owner user
   cost supplying r. Dependency departure and adult entry are separate clocks;
   entry splits birth vintages at 16/20 years and applies the 1/2.1 conversion
   once. Positive net estates fund positive entrant assets and the remainder
   enters a nonutility sink. This provisional stationary ledger does not trace
   own-family inheritance or certify counterparties/physical settlement.

5. **Housing around the first birth.**
   `code/data/psid_followup_mar2026/output/sa_rooms_first_birth_v2/A2h/{coefficients,fit_receipt}.csv`
   and `sa_first_birth_outcomes_v3/A2h_{own,moved}/fit_receipt.csv`.
   Adults were reference person/spouse in the −3/−2 baseline window; the
   estimator uses first biological births, confirmed-childless controls, IW
   weights, person/year effects, age/education covariates and person-clustered
   uncertainty. Bin definitions are in `sa_rooms_first_birth_v2.do:170–180`.
   The omitted baseline is hollow at zero; whiskers are 1.96 times raw SE,
   pointwise intervals. At years 3–4, rooms rise 1.465293280235685 (SE .050270378),
   ownership rises .2688572025306991 and moving falls .2114666960744194;
   fractions convert to percentage points by multiplying by 100. Each outcome
   uses its own complete-case sample. The 1.025-room row-by-row household design
   is a different group. These associations are not an exogenous fertility shock.

6. **Calibration strategy and inputs.** Author-adopted national housing inputs,
   bequest target and pension rule in `CALIBRATION_STATUS.md`, reconciled in
   `docs/model/accepted_input_reconciliation_20260926.md`; reference contract
   `output/model/overnight_calibration_20260928/contract_v1/contract.json`.
   AHS rooms=5.729434240102641; CPS pension/earnings=.2294460118659327;
   payroll tax=.08028070961950022; annual depreciation=.01416143718381309;
   annual property tax=.010598360773872594; SCF estate/wealth=.007291023472616158.
   Pension/earnings compares household means including zeros, not individual
   replacement rates. Bequest recipient scope, income proxies and the model
   first-birth housing observer retain documented empirical counterpart gaps.

7. **2007 stationary calibration.** Exact fourteen-row source:
   `output/model/fertility_identification_20260928/resume_v1/selected_export/primary/target_fit.csv`.
   All rows appear in `Moment | Target | Model` format. Percentage and pp rows
   multiply source fractions by 100. One-child share conditions on mothers ages
   40–44. Three zero-weight rows are validation. Ten moments are scored and
   completed fertility is separately normalized. Full precision, gaps, weights,
   loss contributions and roles remain in `evidence/reference_target_fit.csv`.

8. **Calibrated parameters.** Exact source is the same reference's
   `parameters.csv`, copied to `evidence/reference_parameters.csv`. The slide
   reports all ten searched coordinates plus normalized `psi_child` in three
   columns. Bounds, external restrictions and near-bound flags remain in the
   full CSV. Both fertility taste scales are near their lower bounds under the
   inherited 1%-of-range screen, but neither equals its lower bound.

9. **Fertility at young ages.**
   `output/model/fertility_identification_20260928/two_stream_overnight_v1/comparison_v1/age25_decomposition_inputs.csv`
   and `morning_readout_v1/RESULTS.md`. Use `Original selected` for E01 and
   `Two-birth selected` for E02, not the separate block0506 `Reference` row.
   The 27.1% gap closure is
   (.6060544766097269−.5304463498285714)/(.8095276384290021−.5304463498285714).
   Both selected candidates passed repeats and remain unadopted. E02 adds an
   independent taste opportunity and common event-time proxy, so it does not
   isolate spacing alone. Any curvature discussion is grounded in
   `docs/model/calibration_identification_review_20260929.md`: the Jacobian is
   at a search-round center, not a fresh selected-point identification test.

10. **Credit and further calibration.** Current completion authority:
    `CALIBRATION_STATUS.md` and
    `output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/RESULTS.md`.
    The completion receipt is copied to `evidence/pilot_completion_verification.json`.
    All three arms completed and remain unadopted. The deck displays their
    designs and completion, without a new pilot-results table: independent
    empirical five-ratio/current-earnings entry with D=.25; zero entrant wealth
    with D=0; and nonnegative five-ratio nodes with positive nodes rescaled to
    preserve the original mean, D=0. The transformation applies to the five-bin
    approximation, not raw survey observations; factor=.3632385158888715.
    All use a common 2% annual real rate, no premium, fixed supply scale and
    child benefit, and search nine remaining coordinates against ten scored
    targets. Their price-renewal/population-housing closure is experimental.
    Earlier launch/status receipts remain historical, unchanged evidence;
    their running statements are superseded by the terminal results.

## Retained supporting evidence, no longer displayed

The early sandbox mechanism table is sourced by
`docs/model/situation_report_20260918.md`, section 1. The completed prescribed-price
elasticities remain in
`output/model/fixed_reference_economics_20260928/elasticity_v1/recovery_v1/{README.md,collected_v1/elasticities.csv}`.
Matched saving and grid comparisons remain in
`output/model/publication_refactor_20260929/small_credit_replication_v1/README.md`
and `grid_resolution_v1/credit053_v2/runner/README.md`. The fixed-parameter
five-ratio entry comparison remains in
`output/model/fixed_reference_economics_20260928/entry_ratio_comparison_v1/README.md`.
Their copied supporting artifacts and all receipt history are retained unchanged;
removing their frames does not change or adopt those experiments. The archived
received `source_map.md` preserves the original detailed provenance.

## Build and preservation

The presentation is in `latex/corina_progress_20260930/`; the received version
is retained in `archive/received/`, with the original task location preserved.
Compile twice into a temporary build directory and verify exactly ten pages and
frames with no overlays. Render and inspect every rewritten page for equations,
legibility, clipping, axes and labels. Final build/preservation claims belong
in the lead's updated `verification.json`; this documentation edit does not
itself certify the rewritten PDF. No numerical work is required for this rewrite.
