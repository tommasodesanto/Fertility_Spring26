# Second-birth housing proxy design

This note specifies a bounded extension of the ACS housing analysis. The
estimand is a matched repeated-cross-section housing contrast around a
**coresident second-birth proxy**. It is not an actual second-birth event,
because ACS has no longitudinal person identifier or complete birth history.
The proxy must never be described as completed fertility or as a causal
second-birth effect.

## Source and author alignment

The source is the local extract27 contract, subject to the existing exact-key
and source-fingerprint gate. A child row is linked to a candidate mother when
`MOMLOC` equals the mother's `PERNUM` **within** `(YEAR,SAMPLE,SERIAL)`; a
`PERNUM` value is never joined across households. `NCHILD`, `ELDCH`, and
`YNGCH` remain source diagnostics. The linked roster is used to recover the
second-oldest age; no `MOMRULE` is invented.

The author code is the template, not an exact second-birth estimator:

- `matching.R:35--75` uses exact demographic bins and `Match` with replacement,
  all ties, and matching weights.
- `matching.R:95--109` creates numeric matching variables and restricts the
  first-birth age-at-event to 25--45.
- `matching.R:132--169` takes an event-0 newborn pool (`eldch == 0`) and makes
  shifted age/year observations for negative event times.
- `clean_acs.R:33--40` defines `age1b = age - eldch` and `ch1yob = year -
  eldch`, so this extension uses the same integer-age clock.
- `matching.R:315--335` assigns event time from child age and caps the event
  window.

The extension has two deliberate deviations from the first-birth matching
template. Its controls are mothers with one coresident child, and the donor's
first-child age is added as a lifecycle-alignment variable. It also uses the
second-oldest linked child to define the proxy clock. None of these variables
is a housing or outcome matching variable.

## Proxy and event clock

For each candidate mother, sort valid linked own-child ages in descending
order, $a_1 \ge a_2 \ge a_3$. The proposed clock is

\[
t = a_2,\qquad e = YEAR-a_2,\qquad A_e = AGE-a_2,
\qquad g = a_1-a_2.
\]

Here $e$ is the inferred proxy-event year, $A_e$ is mother age at that proxy
birth, and $g$ is the first-to-second-child age gap. The arithmetic follows
the author's `year - eldch` convention; it carries the same possible one-year
survey-age timing error.

The strict primary diagnostic requires:

1. `NCHILD == 2`, exactly two linked own children, $a_1 \ge 1$, and
   $a_1>a_2\ge0$;
2. if `FERTYR` is observed at $t=0$, it equals 2; missing `FERTYR` is
   retained under the existing contract but must be a separately reported
   support category;
3. $25 \le A_e \le 45$, which extends the author's age-at-first-birth band
   by applying it to age at the proxy birth; and
4. $g\ge1$, excluding age ties that could be twins or births too close to
   distinguish in integer ages.

The strict sample is therefore conditional on no observed third coresident
child. This improves interpretation of post-event housing outcomes but can
select on later fertility and household retention. A separate descriptive
sensitivity may allow `NCHILD >= 2` and use $a_2$ from the full roster. It must
report `NCHILD == 2` and `NCHILD >= 3` separately: later births contaminate the
latter's post-event housing path. Two or more age-zero linked children, or a
tied second-oldest age, are flagged as twin/close-birth ambiguity rather than
silently treated as one second birth.

The event-0 rule is the existing contract's `NCHILD==2`, `YNGCH==0`,
`ELDCH>=1`, with observed `FERTYR==1` excluded and missing/unknown status
retained as flagged support. For $t>0$, `FERTYR` is not required; the
second-oldest linked child supplies the clock. An interview with
no older child (for example ages `[0,0]`) cannot identify a second birth. Rows
where `NCHILD` and the linked-child count disagree are link-quality failures
for the strict sample and are reported, not repaired from `ELDCH` or `YNGCH`.

## Donor construction and support

The donor construction starts once from the event-0 newborn anchor pool. It
does not create a new anchor or matching run for every post-event row, which
would introduce future information and multiplicity. For an anchor with event
year $e$, mother age $A_e$, and gap $g$, each negative target
$k\in\{-5,\ldots,-1\}$ is

\[
YEAR=e+k,\qquad AGE=A_e+k,\qquad
\text{donor's sole-child age}=g+k.
\]

The donor is a one-coresident-child mother in that target cell, matched on the
frozen author variables (year, age, gender, education, marital status, race,
and state) plus the target sole-child age. Matching is with replacement and
all ties, preserving the author matching-weight logic. Donors are not matched
on rooms, bedrooms, ownership, rent, value, labor outcomes, or other
post-treatment characteristics. Carrying education and marital status from
the event-0 anchor follows the author design but is a strong stability
assumption; those variables can themselves change around fertility. The
results should show support and balance rather than treat this as proof of
identification.

The support geometry is binding. A donor target is impossible when $g+k<0$,
because a one-child donor cannot have a negative child age. A full `-5`
through `+10` curve therefore requires $g\ge5$ for every event-0 anchor. This
is a proposed **gap-conditioned diagnostic restriction**, not an author-
original restriction, and it limits generalizability. The author reference
period `-2` additionally requires $g\ge2$; a gap-one anchor cannot contribute
to the reference cell and must be excluded from that comparison, with its
missing support reported. A general $g>0$ curve is possible only with
event-specific support, a changing composition across negative times, and no
claim of one common pretrend or common estimand.

Post-event rows are mapped independently from their second-oldest age, using
$t=a_2$ and $e=YEAR-a_2$. They do not redefine the event-0 anchor or donor
pool. They remain repeated cross-sectional observations, not follows of the
event-0 women.

## Outcomes, weights, and inference

`ROOMS`, `BEDROOMS`, and `OWNERSHP` are interview-level outcomes. The rooms
top-code and year comparability rules in the housing contract remain in force;
bedrooms and ownership are separate outcomes. The vendor cleaning code sets
`wgt = perwt` (`clean_acs.R:217--224`). The matching code then saves
`wgt_original = wgt` and, for matched controls only, sets `wgt =
wgt_match*wgt` (`matching.R:307--313`); parent rows retain their original
`wgt`, while `match_wgt` separately records the matching quantity. The adapter
must preserve `wgt_original`, `wgt_match`, and `wgt` and reuse the same
row-specific role. It must not claim a universal `PERWT`-times-matching
formula, which could double-count weights. A household-level descriptive
level uses one household row and `HHWT`; the two units must not be mixed.

Donor reuse is expected and must be visible in match counts and weights.
Inference should follow the author's matching/event-study variance path, with
household clustering and survey strata/cluster metadata documented where
available. The current source payload does not establish access to survey
replicate weights, so replicate-weight design uncertainty cannot be claimed.
Missing rooms and missing donor outcomes remain missing. Report support by
event time, gap restriction, `NCHILD` group, link-quality status, `FERTYR`
status, housing validity, and donor reuse.

The existing PSID files can provide an external comparison of pre-event rooms
and response direction. They do not identify this ACS proxy and must not be
rerun or used to relabel it as an actual second birth.

## Executable diagnostic builder

The bounded roster/proxy builder is
`code/empirical/acs/kleven_pseudo/build_second_birth_proxy.R`. It consumes a
normalized `data.table`, requires explicit `fertyr_codes` and
`match_covariates`, and returns raw input rows, household-scoped links,
`mother_rows`, strict `anchors`, strict `post_rows`, one-child donors, negative
donor targets, an exclusion-count audit, and the frozen configuration. It does
not match, estimate, assign generated weights, or load real ACS data.

The hand-built deterministic test is
`code/empirical/acs/kleven_pseudo/test_second_birth_proxy_builder.R`; run it
from that directory with:

```text
Rscript test_second_birth_proxy_builder.R
```

The source-level tests are supplemented by the NE-only empirical support
diagnostic in Torch job `18078758`, whose compact receipts are under
`code/empirical/acs/kleven_pseudo/output/second_birth_proxy_diagnostic_20260920/`.
It reports roster, link, gap, event-time, and donor support without matching,
housing estimation, or causal interpretation.

The current housing estimator consumes these strict `post_rows`: its primary
post population remains mothers with `NCHILD==2` and exactly two valid linked
children at each retained row. Consequently, an observed third coresident
child later in the clock is selected out of the primary post path; this is a
composition restriction, not evidence that no third birth occurred. Any wider
`NCHILD>=2` sensitivity must be labeled separately and must not be pooled into
the strict result.
