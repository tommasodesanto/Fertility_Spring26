# Corina adviser progress deck, September 17–30, 2026

Exactly ten Beamer frames and ten PDF pages, with no overlays. This is a separate
progress deck for Corina Boar. The main deck was read only as a style reference.
All paths below are relative to `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/`.
The editable source is standalone and needs no external figures.

## Reference and interpretation

The retained scientific reference is **2007 stationary reference — block0506,
September 28 verified export**. Its full fourteen-row fit and thirty-one-row
parameter table are copied unchanged into `evidence/`. Its weighted objective is
19.581310760138322. The one-birth and two-birth candidates on slide 7 are later,
unadopted experiments, with objectives 7.826226594410982 and 7.8420175378092205.
Their full fit and parameter tables are also supplied. The presentation's fit
table deliberately uses the author's `Moment / Target / Model` format.
Supporting CSVs retain gaps, weights, contributions, restrictions and bound flags.
Never compare these objectives with September 17–19 sandbox objectives.

The fixed reference normalized completed fertility to 2.1 by changing the child
benefit during calibration. Fixed-reference price and credit experiments hold
that saved benefit fixed. The September 30 price/population equilibrium diagnostics
instead clear actual birth renewal through prices and physical housing supply
through population. Fertility near 2.1 in these runs is not independent validation.

## Slide-by-slide sources

1. **Fertility and Housing Tenure Choice.** Two-week scope and synthesis.
   `docs/model/situation_report_20260918.md`, the September 23–30 sections of
   `CALIBRATION_STATUS.md`, and `docs/model/calibration_identification_review_20260929.md`.

2. **Mortgage access and unsecured liquidity.** First table in section 1 of
   `docs/model/situation_report_20260918.md`. These are fixed-benefit early-sandbox
   diagnostics. The sandbox baseline is not the exact September 14 replay.
   The five-year unsecured allowance is a large mechanism probe, not the later
   0.25 annual-income pilot. The rental wedge, mortgage amortization, parent-age
   exit and child earnings-penalty proposals in the early review were not adopted
   wholesale. The slide reports their economic lesson without implying adoption.

3. **Household preferences and demographic accounting.**
   `output/model/fertility_identification_20260928/resume_v1/selected_export/primary/parameters.csv`,
   `code/model/intergen_eqscale_seq_optimized/child_preferences.py`, and the
   September 24 dependency/entry and September 26 estate-funding decisions in
   `CALIBRATION_STATUS.md`. The model stores `psi_child` as the one-child benefit
   b in `B(m)=b*m^(1-gamma)`. Slide notation writes this as
   `v(m)=psi*m^(1-gamma)`, so slide psi is the stored one-child benefit, preserving
   b for financial wealth and B for bequests. The equivalent CRRA coefficient is `(1-gamma)*b`.
   Compensation is at a fixed reference rent, not the endogenous market rent.
   Fifteen persistent income states describe the retained reference, not a claim
   that the earlier AR(1)+iid battery identified a universal earnings process.
   Independent dependency departure and the adult-entry queue are distinct.
   The birth-to-household conversion 1/2.1 applies once. Estate funding keeps
   entrant assets fixed and rejects a shortfall; it does not trace own-family
   inheritance or grant welfare to the residual sink.

4. **Housing around the first birth.** Raw coefficients and receipt:
   `code/data/psid_followup_mar2026/output/sa_rooms_first_birth_v2/A2h/coefficients.csv`
   and `fit_receipt.csv`. The group is all current adults who were reference
   person/spouse in the −3/−2 baseline window, with first biological births,
   confirmed-childless controls, IW weights, person/year effects, age/education
   covariates and person-clustered uncertainty. The omitted baseline is shown
   hollow at zero. Error bars are 1.96 times the raw SE and are pointwise, not
   simultaneous. The title and x-axis use first birth, not household formation.
   Ownership and moving changes come from
   `code/data/psid_followup_mar2026/output/sa_first_birth_outcomes_v3/A2h_own/fit_receipt.csv`
   and `A2h_moved/fit_receipt.csv`. Each uses its own complete-case sample.
   The row-by-row household design gives 1.025 rooms, while the baseline-status
   design gives 1.465293280235685. Those are different groups. The active target
   rounds the latter to 1.465. Rooms/moving used rebuilt dated items; ownership
   was contemporaneous. These associations are not an exogenous fertility shock.

5. **Calibration inputs and empirical counterparts.** Author-adopted national
   housing inputs (September 23), bequest target (September 24) and pension rule
   (September 25) in `CALIBRATION_STATUS.md`; reconciled in
   `docs/model/accepted_input_reconciliation_20260926.md`. The AHS quantity is
   5.729434240102641, the CPS pension/earnings ratio 0.2294460118659327, payroll
   tax 0.08028070961950022, annual depreciation 0.01416143718381309, annual
   property tax 0.010598360773872594, and SCF estate/wealth target
   0.007291023472616158. Pension/earnings compares household means including
   zeros, not an individual earnings replacement rate. The bequest target uses
   mortality-weighted positive NETWORTH−TRUSTS, not an exact published replication.
   The stationary first-birth rooms observer does not implement the empirical
   event-study estimator. The active bequest observer and the target differ in
   recipient/estate scope. Some wealth and income counterparts use proxies.

6. **2007 stationary reference fit.** Exact source:
   `output/model/fertility_identification_20260928/resume_v1/selected_export/primary/target_fit.csv`.
   All fourteen rows appear. Percentage and percentage-point rows multiply
   the source fractions by 100. One-child share is conditional on mothers ages
   40–44. Three rows have zero weight and are marked validation. Completed
   fertility is a separate normalization. The reference has ten searched
   coordinates and ten scored moments; the benefit is separately normalized.

7. **Young mothers and the fertility intensive margin.** Raw age-25 values in
   `output/model/fertility_identification_20260928/two_stream_overnight_v1/comparison_v1/age25_decomposition_inputs.csv`.
   Complete selected results in `two_stream_overnight_v1/morning_readout_v1/RESULTS.md`.
   The 27.1% closes `(0.6060544766097269−0.5304463498285714)/(0.8095276384290021−0.5304463498285714)`.
   Both selected candidates passed repeats, but neither supersedes block0506.
   The weak-curvature conclusion is the round-center Jacobian analysis in
   `docs/model/calibration_identification_review_20260929.md`, not a fresh
   selected-point check or statistical nonidentification proof. The two-birth
   option adds an independent taste opportunity and a common event-time proxy.
   Matching motherhood instead would remove the intensive margin on which the
   model fails. The age clock was checked; relabeling age 25 as 26 does not fix it.

8. **House prices, credit, and births.** Completed September 29 job 18815133:
   `output/model/fixed_reference_economics_20260928/elasticity_v1/recovery_v1/README.md`
   and `collected_v1/elasticities.csv`. Table entries are centered log slopes
   between 0.99 and 1.01 of the reference house price, with mapped rent changing.
   Immediate responses use the same inherited households; completed fertility
   recomputes the lifecycle cohort. The relaxed regime removes artificial
   renter, purchaser and incumbent restrictions while retaining lifetime
   repayment and net-estate solvency. No benefit renormalization, GE, transition,
   fixed-stock assumption or causal credit-channel decomposition is implied.
   The 87.57% first-birth contribution is the +1% reference impact comparison.
   Only one step size is certified. Estate-counterparty limitations remain.

9. **Stationary computation and grid resolution.**
   `output/model/publication_refactor_20260929/small_credit_replication_v1/README.md`
   reports scalar/indexed full workflow 549 versus 417 seconds at identical D=.14,
   160×15, corrected credit and birth-renewal/population-housing closure, six
   lifecycle calls each, with exact outputs. This is not the historical 15-minute
   credit run. The matched grid run is completed18883994:
   `output/model/publication_refactor_20260929/grid_resolution_v1/credit053_v2/runner/README.md`.
   Full workflows 468.606/187.254 seconds, seven calls each, D=.53 in both arms.
   Ownership changes +0.2005pp and price +0.0755%. The nine-state joint-entry
   projection preserves wealth marginals but reduces wealth-income covariance
   5.93%. Weighted ownership discrepancy >5pp covers 3.1521% of household mass,
   from `collected/supplemental_discrepancies.json`. Weights are post-fertility,
   pre-tenure. Fine policies are income-interpolated; large local discrepancies
   combine grid and interpolation effects. No production grid promotion or
   transition/policy-accuracy certificate is asserted. Singleton-market cleanup
   passed equivalence at D=.14 without a separate measured speedup.

10. **Entry wealth, credit, and the next calibration.** Live first section of
   `CALIBRATION_STATUS.md` and owner task `01a0f334-33d4-7d11-855f-157792c5f87f`,
   checked 2026-09-30 20:54:57 UTC. The three 120×9 one-hour exploratory arms use
   D=.25/0/0 at common annual real rate 2%. No rate premium is adopted. The
   nonnegative arm transforms five ratio nodes and rescales positive ones by
   0.3632385158888715, not raw survey observations. Supply scale and benefit are
   fixed, nine other historical coordinates searched against ten scored rows.
   These are preparation settings, not results or production adoption. No pilot
   submission was recorded at this timestamp; update that dated statement only
   from an actual launch/result receipt.
   Completed conventional five-ratio sensitivity at160×15,D=.53 is
   `output/model/fixed_reference_economics_20260928/entry_ratio_comparison_v1/README.md`,
   job18888956. Its native projection clips0.0183%draw mass at the lower support,
   raising mean from0.186519679 to0.186633617. No forward relocation is allowed.
   Main moment differences are small at fixed parameters, not after re-estimation.
   Prior corrected D=0 was blocked by genuinely negative entrant cells, not
   permission to forgive debt or delete mass.
   Historical shock jobs18801439/18801451 failed before any estimated shocks:
   `output/model/fixed_reference_transition_20260928/four_shock_v1/launch_v3/`.
   The requested path consists of successive surprises with inherited states
   and both entry queues, not an announced four-shock path. Follow-up diagnostic
   18820811 was cancelled per the author's instruction, as supplied by the parent;
   no estimates are claimed. This cancellation is authoritative session context,
   while the local status's older running description is stale.

## Build and preservation

Compile twice with `pdflatex -interaction=nonstopmode -halt-on-error corina_progress.tex`.
The deliverable PDF has ten pages. All pages must be rendered and checked for
legibility, clipping, equations, plot axes and labels before delivery.
No models, tests, benchmarks, scientific code edits, cluster jobs or job stops
were performed by this slide task. Source and build products stay in the
separate task workspace because the project directory is read-only for this task.
