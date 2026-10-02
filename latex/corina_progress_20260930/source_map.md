# Corina adviser progress deck: October 2 purchase-rule comparison

This separate fourteen-frame presentation incorporates the fresh-postchecked
October 2 hard and quarter-saving calibrations. It preserves the empirical
review and selected model exposition. There are no overlays, visible `Source:`
footers, internal failure history or dedicated computation section. Paths below
are relative to `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/`.

## Selected sources and paired contract

The authoritative readout is
`output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/collection/readout/`.
Displayed targets and estimates come from `{hard,quarter}/selected_root/{target_fit.csv,parameters.csv}`;
closure comes from each arm's `closure.json`. Local copies are
`evidence/overnight_{hard,quarter}_{target_fit.csv,parameters.csv,closure.json}`.
`overnight_selection_summary.json` records the hard Torch restart chain11 and
quarter local chain54 selection; `overnight_evidence_identity.json` authenticates
all six copies and records selected-root/repeat byte equality. Fresh native
postchecks establish numerical repeatability, not optimizer convergence,
strong identification or baseline adoption.

Losses 97.01121981277964 (hard) and 51.55603604909936 (quarter) use the same
original weights, with target fingerprint
`db60605ef444b747b2ed7f9482b6b8c90ef5e8b1d75c49a3ac1ff6f4275c9ba1`
and weight fingerprint
`2391cd2d4a39a6669a405be34ff14116ad314354353173cb619f7fa7c66043b0`.
The paired comparison is within that contract. Do not compare either loss to
the earlier 31.284 fixed-H0 point or 19.581 block0506 as a calibration improvement:
the closure, eligibility, bounds and weighting changes are bundled.

## Economic changes relative to the October 1 fixed-H0 point

- Population is normalized to N0=1; price clears birth renewal and H0 is derived
  from housing demand. H0 is 7.546002524007095 / 7.288573389887633, rather than
  fixed at 6.293507689200028 with endogenous population scale.
- Original weights replace the earlier selected early-fertility ×4 profile.
  Target values and definitions remain pinned. Ten parameters, including psi,
  remain jointly free for ten scored moments; three validation and one renewal
  row are additionally reported.
- Purchase eligibility is hard closing or quarter-saving credit timing, with
  a buyer's end-of-period net-estate solvency floor in the saving choice.
  No unsecured renter borrowing, 80% mortgage financing, 2% annual real rate,
  nonnegative mean-preserving entry and 120×9 grid are retained.
- The parenthood housing-floor upper search bound is 2.6 instead of 2.3.
  Selected hP=2.562319632969801 / 2.49169824815624 is below that upper bound.
  The parenthood-only floor, constant shares, nonlinear benefit and power
  equivalence scale remain; compensation A(m) is inactive.

Evidence: purchase-rule `README.md`, `plan.json`, selected parameter/closure
records and current `CALIBRATION_STATUS.md`. These are experimental fitted
candidates; specification work does not certify adoption of their estimates.

## Slide-by-slide provenance

1. **Model review and revisions.** Latest author decisions,
   `CALIBRATION_STATUS.md`, and purchase-rule README. Describe reviewed work
   and selected preference architecture without calling either fitted purchase
   rule a paper baseline.

2. **Household decisions.** Retained lifecycle architecture in
   `code/model/intergen_eqscale_seq_optimized/solver.py:278–295,3349,3486–3551`
   and the isolated purchase-engine bindings. Household differences include
   age, wealth, persistent earnings, tenure, children ever born n and children
   at home m. Choices are birth attempts, housing/tenure, consumption and saving;
   earnings are stochastic endowments. The one-birth-per-period architecture
   remains, unrelated to the quarter-saving credit fraction.

3. **Children and housing needs.** Parenthood-only Stone–Geary utility is
   CRRA of c^alpha*[eta_t*(h−hP*1{m>0})]^(1−alpha)/e(m), plus
   psi*m^(1−gamma), with zero benefit at m=0. Here e(m)=((2+.7m)/2)^.7,
   eta_t=1 for renters and chi for owners. There is no additional-child floor
   loading, child-dependent share or compensation A(m). Bequest utility and
   first-birth fixed cost remain separate. Native bindings are documented by
   `output/model/fixed_reference_economics_20260928/utility_floor_round2_v1/runner.py:214–220,285–304` and inherited child-benefit
   code; selected estimates must come from the overnight ROOT parameter tables,
   not earlier seed-only utility-verification receipts.

4. **Demographic and market equilibrium.** Overnight selected `closure.json`
   and purchase-rule README: price clears birth renewal, household population
   is normalized to one, H0 is derived from housing demand. Equal retiree
   pensions balance PAYGO at the inherited tax rule. Dependency departure and
   adult entry remain distinct clocks, with birth-to-household conversion
   applied once. Estate-funded positive entrant assets and residual sink retain
   the documented provisional creditor/recipient/physical-settlement limitations.

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

6. **Calibration strategy and inputs.** Purchase-rule plan/README and selected
   tables pin ten jointly searched coordinates and original ten scored weights.
   Author-adopted empirical inputs remain documented in
   `docs/model/accepted_input_reconciliation_20260926.md` and status: AHS rooms
   5.729434240102641; CPS pension/earnings .2294460118659327; payroll tax
   .08028070961950022; annual depreciation .01416143718381309; annual property
   tax .010598360773872594; SCF estate/wealth .007291023472616158. H0 is derived,
   not fixed or a searched coordinate. The model/empirical counterpart gaps
   for first-birth rooms, recipient scope and income definitions remain.

7. **Identification.** All ten free coordinates are beta_annual, chi,
   first_birth_fixed_cost, kappa_fert, kappa_fert_continuation, theta0, hP,
   child_benefit_curvature, tenure_choice_kappa and psi_child. The ten scored
   moments are childlessness, one-child share among mothers, first-birth age,
   wealth/earnings, bequest/wealth, mean rooms, ownership ages30–55, first-birth
   rooms, recent-parent ownership and early fertility. Treat these as joint
   identifying information: wealth/tenure/fertility moments depend on multiple
   parameters, not a certified one-to-one mapping. Three validation rows and
   birth renewal are separate. H0 is derived and theta1 externally restricted.
   Historical main JMP identification lines432–458 use the old H0/theta1-free
   contract; they do not describe this search. Curvature and tenure dispersion
   are free here. Main JMP slides and both manuscripts remain unchanged; this
   presentation does not claim they are synchronized. No fresh Jacobian rank
   or statistical identification certificate is inferred from fit.

8. **Stone–Geary calibration.** Both selected ROOT target_fit.csv files provide
   all fourteen rows in `Moment | Target | Hard | Quarter-saving` format. Ten
   scored and one birth-renewal check remain above; three validation rows are
   separated below. Fraction-to-percent/pp conversions multiply by 100.
   One-child share is conditional on mothers ages40–44. The historical CSV
   identifier `initial_normalization` checks birth renewal; psi is jointly free,
   not separately adjusted to fertility2.1. Full gaps, weights, contributions
   and roles are retained privately. No earlier young-mother decomposition is
   assigned to these new fits.

9–11. **Figure slots.** Household policies and lifecycle paths; housing supply
   and borrowing constraints; transition dynamics. These are titled placeholders
   requested by the author, not inserted figures or quantitative transition
   results. The selected arms' standard plots remain in the original readout.
   All seventeen per arm were visually reviewed in
   `collection/readout/visual_review/REVIEW.md`; that review records off-support
   boundary behavior, finite/gate screens and plotting-scale limitations.

12. **Parameter estimates.** Both selected ROOT/parameters.csv files provide
    all ten searched estimates. H0 is derived at N0=1; theta1 is externally
    fixed. hP bounds are [.1,2.6], and neither estimate contacts the upper bound.
    Both fertility taste scales are flagged near-bound under the inherited
    1%-of-range screen on [.02,50], but neither equals .02. Complete bounds,
    restrictions and statuses remain in supporting CSVs.

13. **Credit timing.** Purchase-rule README and isolated quarter
    `engines/quarter/refactor_lab/engine/kernels.py:912–920` document
    hard A>=.2Q versus quarter-saving A+.25[b'−(A−Q)]>=.2Q.
    A is financial wealth plus net sale proceeds, Q chosen home's price, and
    b' end-of-period financial balance. The bracket is net saving relative to
    the post-purchase balance A−Q. This is a fraction of within-period net saving,
    not quarterly fertility or quarterly period length. Budget/interest timing,
    mortgage ceiling and buyer-end net-estate solvency remain separate rules.

14. **Calibration comparison.** The two fresh-postchecked fits can be compared
    under their common pinned target/weight contract. Main remaining misses
    include rooms, first-birth housing and early fertility; full paired fit
    tables provide the evidence. There is no optimizer-convergence, policy-response
    or paper-baseline-adoption claim. Policy preparations or running jobs are
    not numerical outcomes for this slide.

## Historical evidence and preservation

All prior evidence and receipts remain unchanged: September28 block0506,
E01/E02, fixed-H0 October1 floor, earlier pilot/price/grid comparisons and their
copied evidence are historical, not sources for the current displayed fit.
The old detailed maps remain in archived received versions. This separate deck
leaves the main JMP slides/manuscripts unchanged and keeps private provenance
without visible source footers.

Compile twice into a temporary build directory; verify fourteen pages/frames
without overlays and inspect every changed page for equations, source identity,
full paired tables, credit definitions, figure-slot labels and clipping. Final
build claims belong to the lead's QA receipt. This documentation update performs
no model, test, benchmark, cluster, policy or Git work.

## October 2 soft comparison and notation revision

The deck now includes the best saved soft point: chain 16, case 0046_nm, loss 23.078309294160004. Its full fourteen-row fit and thirty-one-row parameter tables are copied unchanged to evidence/overnight_soft_*.csv; the saved best receipt is overnight_soft_best.json. Targets, weights and roles match hard exactly. Soft is exploration-unverified and was stopped before final postchecks; hard and quarter-saving were postchecked. Hard/quarter also add a buyer ending-net-estate restriction, so these are separately recalibrated specifications, not a pure timing experiment.

Notation follows the older deck where possible: u(c,s;m), tenure subscripts, xi for first_birth_fixed_cost, and kappa_C for kappa_fert_continuation. This changes notation only.
