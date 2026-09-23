# First-birth rooms estimand audit

## Additional saved-coefficient aggregation, September 23

`common_cohort_contrast.py` examines a concrete way to remove one demonstrated
normalization dependence. It first reproduces the current 0.720246 target from
the saved cohort coefficients and original aggregation weights. It then keeps
the 27 birth cohorts observed at both -1 and +3, takes the contrast within each
cohort, and uses the same cohort weights at both endpoints. Fixed pre-birth
weights give 0.811397 rooms; fixed post-birth weights give 0.797759. The common
cohorts retain 80.89% of fitted pre-birth weight and 88.31% of fitted post-birth
weight. These are candidate aggregations of an existing regression, not new
regressions or adopted targets.

An additive cohort-specific normalization cancels in each within-cohort
contrast. The script verifies this algebra and the source hashes. A first
positive-variance guard stopped on cohort 1971; its observed -1 coefficient is
the sole omitted baseline rather than missing endpoint support. The corrected
guard verifies that case explicitly and retains the cohort. The independent
`common_cohort_contrast_review.json` checks this treatment and the calculation.

The change from 0.720246 also changes cohort composition and weights; it cannot
be attributed solely to normalization. This construction does not establish
causal identification, individual balance, parallel trends, or equivalence to
the model's counterfactual or the author's -2 reference preference. The full
coefficient covariance is unavailable, so no standard error or calibration
weight is supplied. The current target and running objective are unchanged.

This is a read-only audit of the saved PSID event-study outputs, the active
model observer, and a sample-only PSID birth-history count. It ran no regression,
model solve, or remote checkpoint fetch, and it changed no live target or
observer. Reproduce the saved coefficient calculations with:

```bash
python3 output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/audit_rooms_estimands.py
Rscript output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/reconstruct_second_birth_shares.R
```

The script writes `room_estimand_results.csv` and `source_hashes.json` here.

## Empirical construction and numerical receipts

The frozen target is the Sun–Abraham contrast
\(\hat\beta_{+3}-\hat\beta_{-1}\), with event time \(-2\) omitted. The
underlying weighted event study uses person and survey-year fixed effects, age
and education controls, person-clustered standard errors, the PSID longitudinal
weight `IW`, and confirmed childless women as controls. Its estimation sample
has 49,457 woman-household-years and 4,112 women; it uses a single woman per
single-FID dwelling-year. The target receipt records \(\hat\beta_{+3}=1.2017450879\),
\(\hat\beta_{-1}=0.4814988255\), covariance
\(\widehat{\mathrm{Cov}}(\hat\beta_{+3},\hat\beta_{-1})=0.0047806127\), target
0.7202462623815278, SE 0.0852600513385958, and weight 137.5652749002964.
Reconstructing the contrast SE from the rounded event CSV's two marginal SEs
and the saved covariance gives 0.0852600477794, a 3.6e-9 rounding difference
from the full-precision receipt. The fit log shows the covariance extraction
from Stata's full `e(V_iw)` matrix; the compact CSV retains coefficient
marginal variances only.

The saved outputs also support \(\hat\beta_{+3}=1.2017451\) with marginal SE
0.1025444, \(\hat\beta_{+4}=0.8311618\) with marginal SE 0.0953112, and the
point estimate mean of \(\hat\beta_{+2},\hat\beta_{+3},\hat\beta_{+4}\),
0.9196888. They do **not** preserve the pairwise covariance terms among +2,
+3, and +4, so the mean's covariance-based SE is unavailable. `room_estimand_results.csv`
marks it unavailable rather than using an independence approximation.

The saved PSID correction folder contains aggregate event-curve, cohort-support,
and timing-validation receipts, but no person-level estimation sample with
ordered biological birth dates. I therefore reconstructed the builder's
pre-regression sample directly from the local 5.9 GB PSID shelf, using its exact
age, sex, relation, positive-weight, death-year, one-FID-dwelling, room-wave
alignment/missing-code, education, household-year de-duplication, birth-history,
and first-birth-after-entry restrictions. The resulting 49,872 rows, 4,527
women, and 2,486 treated women exactly match the target receipt's pre-regression
input counts. No regression was fitted.

For the exact event-time subset with valid rooms and covariates, the documented
share with a second biological child by \(+3\) is 677/1,344 (50.37%); by \(+4\)
it is 807/1,329 (60.72%). These counts include a second child recorded in the
same birth year, which may include a multiple birth but cannot distinguish one
from another annual-year tie. Requiring a strictly later distinct recorded
biological birth year gives 665/1,344 (49.48%) by \(+3\) and 791/1,329
(59.52%) by \(+4\). The singleton-pruned reconstruction leaves these exact
event-window denominators and rates unchanged. The pruned estimation-sample
reconstruction has 2,458 treated women; 1,108/2,458 (45.08%) have a documented
later distinct biological birth year by +3, and 1,331/2,458 (54.15%) by +4.
For comparison, among all 2,486 pre-regression treated women, the respective
documented shares are 46.06% and 54.95% for a second child, or 44.77% and
53.74% for a later distinct birth year. Counts by denominator and history
coverage are in `second_birth_event_shares.csv`; sample and identity receipts
are in `second_birth_sample_receipt.json`.

The target's saved event-study dataset has no person-level `e(sample)` marker,
and its text log reports the aggregate counts but no singleton-drop tally.
However, an independent iterative two-way-FE singleton reconstruction from the
exact pre-regression rows finds 415 one-row person IDs and zero singleton years.
Dropping them once leaves 49,457 rows, 4,112 people, 2,458 treated people and
1,654 confirmed never-treated controls, exactly matching the saved fit's
observation, person, and control counts. The event-time +3/+4 group sizes also
match before and after this pruning. This is strong aggregate evidence for the
sample reconstruction; because the individual marker is unavailable, it does
not prove that each row's membership is identical to the saved `e(sample)`.

The local shelf provides `RELCHI1ID`–`RELCHI20ID`, each labeled “Ind's child #,
unique ID.” In the treated women's full panel history, the audit found no
timed biological child record with a missing child ID, no person-slot linked to
multiple child IDs, no child ID moving across slots, no conflicting birth years
for an ID, and no difference between slot-based and child-ID-based first,
second-child, or later-distinct birth years. These checks verify the slot
identity assumption for the records used in this calculation; the output
retains child-ID comparisons so the result is not inferred from first-birth
minimum reproduction alone.

`RELCHIREP` is labeled “reported number of children, with or without records.”
A missing second biological child-year can therefore mean unrecorded timing;
the documented shares should not be read as proof that all such births are
absent. The table also records the number of event-time women whose reported
child count is at least two and the number with an explicitly untimed
biological slot. The count of reported children is not treated as a
second-birth event or used to impute timing.

## Model observer contract

`code/model/tools/e5f_initial_housing_observer.py::_stationary_birth_diagnostic`
calls `begin_dated_first_birth_housing_branch(..., origin_period=0)` and
`finish_dated_first_birth_housing_branch(..., destination_period=1)` in
`code/model/tools/run_e5f_transition_calibration.py`.

- **Origin and treated cohort:** the stationary pre-fertility distribution is
  selected at childless states with settled readiness. The sequential observer
  weights a successful first birth by age-specific fecundity times attempt
  probability for model age cells satisfying `A_f_start <= j+1 <= A_f_end`.
  The joint-choice path uses its corresponding first-birth factorization.
- **Control:** the same origin states are cloned at equal mass and remain
  childless by construction. The treated branch records the first child.
- **Horizon and aging:** both branches advance one model period, which is four
  years under the active model contract. The origin-date survival, location,
  tenure, saving, income-transition, and child-aging kernels carry them forward.
  The aggregate Census age bridge is explicitly excluded. The stationary
  diagnostic reuses the same stationary `evaluation` at origin and destination.
- **Destination choices:** the branches are feasibility-gated using destination
  policies. At the destination the treated branch can have a continuation
  birth; the control is explicitly held childless. Destination location/tenure
  choices and housing policies determine realized rooms. The observer reports
  treated and control mean rooms, their difference, origin/destination mass,
  and aggregate treated continuation-birth mass.
- **Normalization and weighting:** each destination mean is total realized
  rooms divided by its branch mass. Equal origin masses and equal destination
  masses are asserted. The result is weighted by the model's stationary
  successful-first-birth flow and its induced state composition, not by the
  Sun–Abraham cohort weights.

The target and model object are thus demonstrably different in horizon and
starting point: the target compares event time +3 with event time -1 (a
four-calendar-year contrast) relative to the omitted -2 period, while the model
starts before the fertility/tenure decision and measures at the next four-year
date. Controls also differ (confirmed never-treated people versus exact
same-state no-birth counterfactual), as do weights. At the model destination the
treated branch may have a second birth. These are estimand-alignment differences
by construction. This audit finds no coding defect in those explicit mechanics;
it also does not prove that they explain a particular share of the numerical
gap.

## What saved B exposes and what remains unresolved

The literature memo reports B's scalar response as 1.341457. The local
`earnings_entry_battery_v1/final_readout/B/selected/native_summary.json` exists,
but its compact summary format has no matched-branch decomposition. The nearby
`selected_B_checkpoint_audit.json`
concerns the separate earnings-entry battery and identifies its native
checkpoint only by a remote `/scratch/...` path; it is not a local B
housing-branch checkpoint and its saved summaries do not contain branch age,
tenure, or second-birth arrays. No remote checkpoint was downloaded.

The observer's saved result schema, when a branch result is collected, exposes
only aggregate origin mass, destination mass, treated/control rooms, and total
treated continuation-birth mass. It does not itself expose the requested
origin-age-by-tenure response table or counterfactual response with all
second-birth households removed. The 1.341457 scalar cannot be decomposed from
that scalar alone. No 0.2-room decision threshold is supported by the receipts.

Accordingly, the audit cannot attribute the full \(1.341457-0.720246=0.621211\)
gap to either observer alignment or economics. The observer has a horizon,
baseline, treatment-control, second-birth, and weighting mismatch with the
empirical contrast; its magnitude is not yet quantified. The parenthood housing
floor being at its bound is a parameter-bound fact, not evidence that the floor
causes the residual. Economic fit, estimand mismatch, and possible bugs remain
separate propositions; this audit establishes the mismatch, identifies no
observer coding bug, and does not establish the economic mechanism.

## Minimal next test

1. If an exact local saved B evaluation becomes available, run a read-only
   postprocessor against it that tabulates branch origin mass and response by
   origin age and tenure, and rebuilds the destination comparison after removing
   treated states with a continuation birth. Preserve the same mass normalizers
   and report the discarded mass. This needs the saved state/policy arrays but
   no equilibrium solve.
2. Separately, use the local PSID shelf to reconstruct only treated-woman
   eligibility and count second biological births by first-birth year +3 and
   +4 under the target's exact sample rules. Record each denominator and wave
   availability. This is a sample-only data audit, not a regression.
3. If covariance-based uncertainty for the +2..+4 average is needed, rerun the
   already specified regression only after authorization or recover the full
   `e(V_iw)` matrix from an existing saved Stata artifact. The current target
   receipt and coefficient table are insufficient.

Hashes for the reviewed source/receipt files and audit scripts, including the
active observer and branch implementation, are recorded in `source_hashes.json`.
