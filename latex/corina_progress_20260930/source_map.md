# Corina adviser progress deck: selected Stone–Geary specification

This separate ten-frame adviser presentation preserves the empirical review
covering September 17–30 and incorporates the October 1 selected-repeat packet
and latest author decision. It has no overlays, visible `Source:` footers,
internal failure history or dedicated computation section. This file records
private provenance. Paths are relative to
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/`.

## Selected point and objective contract

Let **selected packet** below mean
`output/model/fixed_reference_economics_20260928/utility_floor_psi_v1/deployment/monitor_snapshot/verified_global_20261001T0941NY_chain7_0173/`.
`ROOT/target_fit.csv` and `ROOT/parameters.csv` are the displayed scientific
sources, copied locally as `evidence/floor_target_fit.csv` and
`evidence/floor_parameters.csv`. Local `floor_input_contract.json`,
`floor_closure.json`, `floor_repeat_verification.json` and
`floor_evidence_identity.json` preserve the supporting contract and identity. `verification.json` certifies selected ROOT/REPEAT identity for all
14 target rows, 31 parameter rows and 17 standard-plot hash matches. These
checks establish repeatability, not optimizer convergence or economic
plausibility. The author now selects
the parenthood-only Stone–Geary specification. The saved packet predates that
choice and says `experimental_not_adopted`; author selection of a specification
does not certify parameter convergence or production adoption of the point.

The selected `input_contract.json` uses `weight_profile=early4x`, with the
`early_fertility` weight multiplied by four. Its search objective is
54.52095879194978. ROOT's original-base-weight contributions sum to
31.284007255664957. Distinguish both objectives; neither is a controlled loss
comparison with September 28 block0506. Ten coordinates are free, including
psi_child; H0 is fixed at 6.293507689200028. The first CSV row retains the
historical identifier `initial_normalization`, but economically it checks birth
renewal under this closure: psi is no longer separately normalized.

## Economic changes relative to block0506

- **Author-selected preferences:** parenthood-only Stone–Geary housing floor,
  no compensated-share factor A(m), constant consumption/housing shares;
  nonlinear child benefit and the power equivalence scale remain.
- **Joint estimation:** psi_child is free alongside nine other coordinates,
  rather than adjusted separately to completed fertility; H0 is now held fixed.
  Its [.01,.5] psi interval is a diagnostic search restriction, not an
  externally identified economic bound.
- **Entry/credit working contract:** negative five-bin wealth/income ratios
  become zero and positive nodes are scaled by .3632385158888715 to preserve
  the original mean. This transforms bin means, not raw survey observations.
  Renter saving floor is zero, with no unsecured borrowing; the common annual
  real rate remains 2%. Corrected homeowner sale/no-taper and mortality-repayment
  rules are retained from the isolated credit workflow.
- **Equilibrium working contract:** price clears birth renewal; population
  clears absolute housing supply. This differs from block0506's unit-population
  housing-price clearing with separately normalized child benefit.
- **Weighting/search:** the selected case uses early-fertility ×4. Target values
  remain unchanged. Estimates, identification strength and optimizer convergence
  remain provisional; repeated solutions do not certify search convergence.

Evidence: selected `input_contract.json`, `ROOT/closure.json`, `verification.json`,
ROOT tables, and the corresponding sections of `CALIBRATION_STATUS.md`.

## Slide-by-slide sources

1. **Model review and revisions.** Latest author decision, canonical
   `CALIBRATION_STATUS.md`, and selected input contract. Previous reviewed
   choices remain documented by `docs/model/accepted_input_reconciliation_20260926.md`
   and `docs/model/calibration_identification_review_20260929.md`. Describe
   specification selection separately from numerical adoption of estimates.

2. **Household decisions.** Retained lifecycle architecture in
   `code/model/intergen_eqscale_seq_optimized/solver.py:278–295,3349,3486–3551`
   and selected fixed-reference bindings. Households differ in age, wealth,
   persistent earnings, tenure, children ever born n and children at home m.
   They choose birth attempts, tenure/housing, consumption and saving; labor
   earnings are stochastic endowments. The one-birth-per-period architecture
   is retained, not the separate two-birth experiment. Keep numerical timing
   constants in calibration rather than model exposition.

3. **Children and housing needs.** Selected input contract and native utility
   bindings. Material utility is CRRA of
   c^alpha * [eta_t*(h−h_P*1{m>0})]^(1−alpha)/e(m), with
   e(m)=((2+.7m)/2)^.7, eta_t=1 for renters and the owner service premium for
   owners. Add psi*m^(1−gamma), with zero child benefit at m=0.
   The floor applies to parenthood, not each additional child; alpha is constant
   and A(m) is inactive. Bequest utility and first-birth cost remain separate.
   Candidate parameters must be taken from ROOT/parameters.csv, not the seed
   utility receipt: h_P=2.3 at the selected upper bound, not 2.080890957425898.
   `utility_verification.json` records the input-contract seed. The runner's
   candidate binding applies the selected coordinates before model evaluation;
   its utility_checks document the floor, constant shares and inactive factor.
   Exact implementation: `output/model/fixed_reference_economics_20260928/utility_floor_round2_v1/runner.py:214–220,285–304`,
   with `code/model/refactor_lab/engine/child_preferences.py` for nonlinear
   child benefits. Seed utility checks must not replace candidate-table identity.

4. **Demographic and market equilibrium.** Selected `ROOT/closure.json` and
   input contract; reference receipt's inherited demographic/fiscal primitives.
   House price clears birth renewal and population scales housing demand to
   absolute supply. H0 is fixed. Equal retiree pensions balance PAYGO at the
   inherited tax rule. Dependency departure and adult entry are separate
   clocks, with the birth-to-household conversion applied once. Positive net
   estates fund positive entrant assets, leaving a nonutility residual sink;
   creditor/recipient and physical-settlement limitations remain provisional.

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

6. **Calibration strategy and inputs.** Selected input contract, ROOT tables,
   author-adopted input decisions in `CALIBRATION_STATUS.md`, and
   `docs/model/accepted_input_reconciliation_20260926.md`. Ten jointly free
   coordinates include psi; H0 is fixed. Target values are unchanged; selected
   search weighting is early4x. Empirical inputs remain AHS rooms
   5.729434240102641, CPS pension/earnings .2294460118659327, payroll tax
   .08028070961950022, annual depreciation .01416143718381309, annual property
   tax .010598360773872594 and SCF estate/wealth .007291023472616158.
   Pension/earnings compares household means including zeros, not individual
   replacement rates. Recipient scope, income proxies and the model event-study
   counterpart gaps remain documented limitations.

7. **Stone–Geary calibration.** Selected ROOT/target_fit.csv: all 14 rows,
   ten scored, three validation and one unscored birth-renewal check. Percentage
   and pp entries multiply fractions by 100. One-child share conditions on
   mothers ages 40–44. The CSV retains original base weights/contributions;
   distinguish its 31.284007255664957 sum from search loss 54.52095879194978.
   Selected early fertility is .5312175503600473 versus .8095276384290021.
   Do not present that miss as a winner-specific young-mother decomposition:
   ROOT/REPEAT observer reports needed for motherhood and conditional children
   are absent, as recorded in `CALIBRATION_STATUS.md` and
   `utility_floor_psi_v1/mechanism_responses_v1/matched_credit_diagnostic/WINNER_AGE25_LOCAL_CHECK.md`.

8. **Calibrated parameters.** Selected ROOT/parameters.csv supplies all ten
   jointly free estimates: beta_annual, chi, first_birth_fixed_cost, kappa_fert,
   kappa_fert_continuation, theta0, h_P, child_benefit_curvature,
   tenure_choice_kappa and psi_child. H0 is fixed. The selected h_P=2.3 equals
   its upper bound. Full restrictions and near-bound flags remain in the CSV;
   those flags also mark both fertility dispersions and theta0. Psi is estimated,
   not normalized. ROOT/REPEAT parameter SHA256 is
   `7949b86204796b6dd6c8db78988003c0bdd2a5077ce3d9f2a6eeb4fb90cb1dd5`.

9. **Preference alternatives.** Previous compensated-share design and its
   fixed-reference-rent interpretation remain in
   `code/model/intergen_eqscale_seq_optimized/child_preferences.py:11–14,39–51`
   and historical reference tables. The fixed-parameter `floor`, `no_A` and
   `constant_alpha` comparisons are recorded in
   `output/model/fixed_reference_economics_20260928/parenthood_floor_quick_v1/`
   and the September 30 utility-comparison sections of `CALIBRATION_STATUS.md`.
   The softened-threshold test is historical, fixed-parameter evidence in
   `output/model/soft_housing_probe/README.md` and the September 26 status
   section; it does not supply this selected calibration. The author requests
   an equivalence-scale comparison; this is a question/preparation until a matching verified results
   packet exists. Do not claim that its estimates or fit are already completed.
   Previous E01/E02 or earlier floor-point decompositions cannot be transferred
   to the selected floor point.

10. **Credit, initial wealth and remaining fit.** Selected input/closure and
    the three completed unadopted entrant-wealth pilots in
    `output/model/fixed_reference_economics_20260928/entry_calibration_pilot_v1/RESULTS.md`.
    Distinguish their alternative laws/limits from the selected nonnegative-mean,
    no-unsecured-borrowing contract. The latest questions concern initial wealth,
    preference alternatives and remaining young-fertility/housing misses.
    Do not substitute undocumented mechanism or recalibration results.

## Retained historical supporting evidence

September 28 block0506 and its 19.581310760138322 objective remain frozen,
historical scientific evidence in
`output/model/fertility_identification_20260928/resume_v1/selected_export/primary/`
and copied `evidence/reference_*` CSVs. Those files do not supply the displayed
selected floor fit. E01/E02 results remain in the two-stream overnight packet;
prior copied fit/parameter and age-25 tables remain historical. Early sandbox,
price elasticities, saving/grid timing and five-ratio sensitivity provenance
remain in `archive/received/source_map.md` and their original evidence files.
All prior receipt history is preserved unchanged. No visible source footers
are needed to preserve these private records.

## Build and preservation

The deck is in `latex/corina_progress_20260930/`; received versions are preserved
under `archive/received/`. Compile twice into a temporary build directory and
verify ten pages/frames without overlays. Inspect every rewritten page for
source identity, equations, objective labels, fit/parameter tables and clipping.
Final build/preservation claims belong to the lead's verification receipt;
this documentation edit itself does not certify the PDF. No model execution
is required or performed for this rewrite.
