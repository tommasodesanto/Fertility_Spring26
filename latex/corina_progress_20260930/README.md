# Corina progress presentation

Separate adviser update covering **September 17–October 2, 2026**, prepared for Corina
Boar. This dated presentation was explicitly requested and does not replace the
continuing [JMP Slides](../JMP_slides/JMP_slides.tex).

- [Presentation PDF](corina_progress.pdf)
- [Editable standalone Beamer source](corina_progress.tex)
- [Source and evidence bundle](corina_progress_source_bundle.zip)
- [Slide-by-slide provenance](source_map.md)
- [Verification receipt](verification.json)

## Verified scientific facts

The displayed results compare the best saved soft point with the October 2
fresh-postchecked hard-closing and quarter-saving calibrations. Authoritative full packets are in
[overnight readout](../../output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/collection/readout/),
with unchanged local copies under `evidence/overnight_*`. Each has all fourteen
fit rows, thirty-one parameter rows and seventeen standard plots. The copied
identity receipt records byte-identical selected-root/repeat tables and closure.
Estimates remain provisional: postchecks do not certify optimizer convergence
or paper-baseline adoption.

All three use the same original target values and weights: soft loss
23.078309294160004, hard loss
97.01121981277964 and quarter-saving loss 51.55603604909936. The loss uses a common objective, but hard and quarter-saving also impose
a buyer net-estate floor; this is not a pure timing experiment. They are not improvement comparisons
with the earlier 31.284 or 19.581 objectives, whose economic/weighting contracts
differ. Ten parameters, including the child-benefit level, are jointly free.
Ten moments are scored, three are validation, and one checks birth renewal.
Housing supply scale is now derived at population normalized to one: H0 is
7.546002524007095 (hard) and 7.288573389887633 (quarter). Price clears birth renewal.

The parenthood-only Stone–Geary specification retains constant expenditure
shares, the power equivalence scale and nonlinear child benefit, with no
compensated-share factor. The housing-floor search bound is now [0.1,2.6];
estimates 2.562319632969801 and 2.49169824815624 are not at the upper bound.
Nonnegative mean-preserving entrant wealth, the common 2% annual real rate and
no unsecured renter borrowing are retained. Both purchase rules enforce buyer
end-of-period net-estate solvency. These changed closure, purchase-rule and
bound settings must be distinguished from the earlier fixed-H0 floor point.

The hard purchase test is A >= .2Q. The quarter-saving test is
A + .25[b'−(A−Q)] >= .2Q, where A is financial wealth plus net sale proceeds,
Q the chosen home's price and b' the end-of-period financial balance. The
quarter refers to the eligible fraction of net saving, not quarterly fertility.
No new policy-response result is claimed.

## Review and changes

The fourteen-frame narrative explains the model, joint identification, complete
paired calibration, purchase-credit timing and remaining fit:

1. Model review and revisions
2. Household decisions
3. Children and housing needs
4. Demographic and market equilibrium
5. Housing around the first birth
6. Calibration strategy and inputs
7. Identification
8. Stone–Geary calibration
9. Household policies and lifecycle paths (figure slot)
10. Housing supply and borrowing constraints (figure slot)
11. Transition dynamics (figure slot)
12. Parameter estimates
13. Credit timing
14. Calibration comparison

Identification discusses all ten parameters and all ten scored moments as joint
sources of information, not a one-parameter/one-moment mapping. H0 is derived
and theta1 externally restricted; child-benefit curvature and tenure taste
scale are free. The old main JMP identification text describes a different
parameter contract and remains unchanged. The fit uses
`Moment | Target | Hard | Quarter-saving`; all fourteen rows remain, with three
validation rows separated below. Full gaps, weights, contributions, restrictions
and bound flags remain in supporting CSVs. Three titled figure slots remain
pending; no transition estimate is implied by a placeholder.

The prior preference-alternatives section is replaced by credit timing. The
PSID figure and private empirical provenance are preserved. No visible source
footers, internal failure history or dedicated computation section are added.
All seventeen standard plots per selected arm have a saved visual review at
[plot review](../../output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/collection/readout/visual_review/REVIEW.md).
It documents boundary/support limitations; visual screening is not a convergence
or policy-acceptance certificate. Earlier receipts describe preserved versions.
Final deck verification checked the candidate identities, equations,
all paired table cells and all fourteen rendered pages. No model, test suite, benchmark, cluster
job or job stop was run for this documentation update. The main JMP deck and
manuscripts remain outside this task's edit scope.

## Preservation and build

The received source, PDF, source map, QA receipt and original source bundle are
preserved in `archive/received/`. `archive/transfer_manifest.json` records the
incoming hashes. The original task directory is also preserved:
`/Users/tommasodesanto/Documents/Codex/2026-09-30/task-4/corina_progress_20260930/`.
Its transient build products were moved out of this new active folder to the
review's temporary build directory.

Compile the source twice with `pdflatex -interaction=nonstopmode -halt-on-error`,
directing auxiliary files to a temporary build directory. The final document
must contain fourteen pages and fourteen frames with no overlays. The final
verification receipt records compilation, rendered-page inspection, source
hashes, preservation checks and Library delivery availability.


Three titled illustration slots follow the calibration table (pages 9–11):
household policies and lifecycle paths; housing supply and borrowing constraints;
and transition dynamics. Figure selection and insertion remain pending.
No fixed slide-count limit remains.

## October 2 soft comparison and notation revision

The deck now includes the best saved soft point: chain 16, case 0046_nm, loss 23.078309294160004. Its full fourteen-row fit and thirty-one-row parameter tables are copied unchanged to evidence/overnight_soft_*.csv; the saved best receipt is overnight_soft_best.json. Targets, weights and roles match hard exactly. Soft is exploration-unverified and was stopped before final postchecks; hard and quarter-saving were postchecked. Hard/quarter also add a buyer ending-net-estate restriction, so these are separately recalibrated specifications, not a pure timing experiment.

Notation follows the older deck where possible: u(c,s;m), tenure subscripts, xi for first_birth_fixed_cost, and kappa_C for kappa_fert_continuation. This changes notation only.

October 2 prose cleanup: shortened review and model exposition, removed repeated workflow summaries, and retained the full identification and fit tables. Equations and numerical results are unchanged.

Calibration setup and identification are combined in one parameter/target slide (13 frames total); recomputed moments are marked. Regression slide notes the earlier mistake without changing estimates.

## GE financing comparison and impact response

Replaced the fixed-price policy slide with permanent steady-state and temporary-impact GE slides. The endpoint comparison is copied unchanged to evidence/permanent_steady_state_comparison.json. H0 and preferences remain fixed within each arm. Permanent steady states pass, but their dated transitions fail. Temporary 48/64-date paths pass; the displayed impact uses 48 dates: hard first-birth flow -0.5416896477%, hazard -0.08487845 pp; quarter flow -0.4305795241%, hazard -0.06709008 pp. Evidence: output/model/fixed_reference_economics_20260928/purchase_rules_overnight_v1/mechanism_deployment/README.md and accepted_hard_temporary_h48_0722.json / accepted_quarter_temporary_h48_0752.json. These use the earlier overnight fits, not the later fresh-search winners. No new model runs or recalibration.
