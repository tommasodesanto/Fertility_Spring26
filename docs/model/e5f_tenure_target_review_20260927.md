# Tenure-choice scale: existing target review, September 27

Read-only review of saved empirical builders, tables and existing model evidence;
no empirical reconstruction or model solve. The only new artifact is this note.
The author now requests internal calibration of the tenure-choice scale. This
supersedes the older Google-ledger description of the scale as provisionally fixed
at 0.005. Keep the current recent-parent ownership gap in the main candidate;
this note does not change the running objective or adopt a replacement target.

## Recommendation for the 22:00 discussion

Use the existing recent-parent ownership gap as a **provisional joint restriction**
on tenure dispersion, the ownership preference and the first-child housing
loading. Do not describe it as uniquely identifying the tenure scale. Retain the
ownership level and housing-response rows alongside it. The best existing
conditional-wealth-gradient comparison is the saved PSID initial-renter transition
gradient below; use it as validation before proposing to add or substitute it.
Its model counterpart has not been authenticated on the current specification.

A larger tenure shock can flatten sorting, but equilibrium wealth, prices and
fertility also change. Neither the ownership level nor a single parent gap is a
one-to-one measure of the shock scale. Joint sensitivity and parameter tradeoffs
must be checked in the new search rather than inferred from parameter labels.

## Maintained recent-parent ownership row

The pinned fourteen-row target provenance is
`output/model/calibration_code_integration_20260927/target_provenance.json`.
The target is 0.127608 (12.761 percentage points), defined as

\[
\Pr(\text{owner}\mid NCHILD>0,ELDCH<4)
-\Pr(\text{owner}\mid NCHILD=0).
\]

Sample: national ACS 2005–06 household heads ages 30–55, housing structures
UNITSSTR 3:10, positive HHWT; no metro restriction. This is a weighted difference
of group means, not a regression conditional on income or wealth. There are no
fixed effects or causal interpretation. Builders named by the contract are
`output/model/e5f_matched_pf_20260909a/design_research/housing/inspect_early_housing.py`
and `summarize_early_housing.py`; saved source is
`specification_followup/housing_profiles_v1/full/target_recomputed.json`.
National sampling uncertainty has **not** been estimated. Weight 27055.823 is
inherited from the former 42-metro bootstrap, not a national inverse variance.

The current model observer compares actual births from previously empty-dependent
homes with currently empty homes, including former parents, in a synchronized
post-fertility snapshot; head-age overlap is uniform within the four-year age
cell. It is explicitly an approximate ACS observer, not a reconstruction of
annual survey interviews. The contract records this mismatch rather than hiding
it. These qualifications matter if this row gains responsibility for another
free parameter.

Historical evidence in
`docs/model/e5f_ssj_experiments_final_report_20260916.md:211–242` shows the old
family ownership gap falling from 0.161 to 0.121, 0.062 and 0.029 as the tenure
scale rises from 0.005 to 0.01, 0.02 and 0.05. This establishes sensitivity in
that old model. It does not establish current identification or rule out current
recalibration: utility, targets and reference parameters have since changed.
That experiment also shows why raising taste dispersion solely to cure numerical
threshold jumps is an economic change.

## Best existing conditional-gradient validation

Builder:
`code/data/psid_followup_mar2026/build_tenure_liquid_wealth_gradient.R`.
Saved moments and quintile cells:
`code/data/psid_followup_mar2026/output/tenure_liquid_wealth_gradient/`.

The top-minus-bottom initial liquid-financial-wealth/income quintile difference
in subsequent ownership is **0.124114**, person-bootstrap SE **0.037782**,
percentile 95% interval **[0.045835, 0.199467]**, with 2,933 people.
The bottom and top quintile ownership shares are 0.315858 and 0.439973.

Exact construction matters. PSID records span 1984–2019; reference persons with
valid tenure, positive person weights and no post-death observation are retained.
For each person, the script averages observations at ages 25–30, requires mean
ownership below 0.5 and family income above 1000, and constructs initial wealth
from savings plus financial funds less other debt, divided by family income.
It trims the weighted first and ninety-ninth percentiles before defining quintiles.
The later outcome equals one when mean ownership at ages 31–35 is at least 0.5.
Thus “owner by 35” is a shorthand for **majority ownership during observed ages
31–35**, not ever buying before the thirty-fifth birthday. Initial renters can
have some owned observations because selection uses an average below 0.5.
Weights are the individual's mean weight over initial observations.

There are no age, income, year or household fixed effects: conditioning is on
initial renter status and wealth/income quintile. The 500-draw bootstrap resamples
people and recomputes quintile cutoffs, conditional on the original trimming.
The reported quintile shares are not monotone across all five bins; ties also
produce unequal weighted bin sizes. Do not call this a smooth estimated slope.

This is useful because it measures sorting into ownership by **pre-purchase**
resources, avoiding mechanical contemporaneous wealth changes at purchase. It
still reflects credit, preferences, income persistence and saving behavior jointly;
it cannot identify tenure dispersion alone. It spans a different historical window
from the 2007 baseline, and its financial-fund measure includes assets such as IRAs.
A current model comparison must reproduce the observation windows, initial-renter
selection, asset concept, income denominator, weights and quintile construction.
No such current-source matched observation is established by this review.

## Existing residual-dispersion candidate, not a ready substitute

The separate builder `code/data/psid_followup_mar2026/build_tenure_residual_variance.R`
and saved `output/tenure_residual_variance/` measure four-year-ahead ownership
prediction error (the Brier score, mean squared error of predicted probabilities).
Its five-fold person-cross-fitted value is 0.117113, conditional bootstrap SE
0.002102, using 32,378 person-years from 8,775 people (initial years 1984–2015).
Initial reference persons are ages 25–55; the outcome follows the same person
four years later even if reference-person status changes. Covariates are current
tenure, income, age, liquid wealth/income, children, marital status and year.

This is conceptually closer to residual choice dispersion, but is not pure taste
noise: omitted persistent heterogeneity and prediction-model error also contribute.
Its SE holds fitted probabilities fixed rather than re-estimating the prediction
model in every bootstrap. The saved “model feasible” alternative is 0.117613
(full sample) or 0.118086 (cross-fitted), has no reported SE, and uses the obsolete
one-shot 0 / 1–2 / 3+ child grouping. It is not ready to replace a current target.
Historical July notes explicitly required a matched simulated prediction exercise
before hard-target use and kept the wealth gradient as validation. Do not silently
revive the older target designation as an accepted September specification.

## Decisions still required

1. Keep the recent-parent gap for the present internally calibrated tenure-scale
   candidate, with the current measurement and inherited-weight caveats visible.
2. At the discussion, evaluate whether tenure scale, ownership preference and
   first-child housing loading remain separately informative using the complete
   moment system; do not declare identification from one sensitive row.
3. If an additional empirical discipline is needed, first compare the saved PSID
   wealth gradient with an exactly matched model statistic. Reconstructing its
   sample or adopting a residual-dispersion row requires an explicit new target
   contract; no change to the live fourteen-row system is made here.

Google context was read from the existing local readback
`/tmp/codex-ledger-final-readback-20260926d/document-text.md`; the lead confirmed
that the latest direct read has the same September 27 05:13 modification time.
The latest author instruction on internal calibration takes precedence. No Google
write occurred, and no new empirical numbers were generated.
