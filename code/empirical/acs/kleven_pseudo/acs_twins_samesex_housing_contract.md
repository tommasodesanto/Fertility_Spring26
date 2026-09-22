# ACS twins-like and same-sex housing extension contract

Status: 2026-09-21 execution authorization received for a bounded national
Twin1/SameSex2 fixture-test, smoke, and single production run (see
`run_twins_samesex_iv.R`, `twins_samesex_iv_lib.R`,
`test_twins_samesex_iv.R`). The completed national first-birth rooms/ownership
results remain unchanged and are not rerun or relabeled by this work.

## Feasibility evidence already available

The verified extract27 source audit contains 2,784,052 Northeast person rows,
2,784,052 unique `(YEAR,SAMPLE,SERIAL,PERNUM)` keys, 722,337 valid `MOMLOC`
links, and exact key/value concordance for 2005--2019 ACS 1-year support. The
strict existing second-birth roster diagnostic found 1,439,607 female rows with
valid age, 132,622 strict eligible mother rows, 6,872 event-0 anchors, 1,771
anchors supporting the full -5:-1 donor window, 6,325 supporting reference
event -2, 70,063 post rows, 211,056 one-child donors, and 21,919 donor targets.
It did not construct a sex-composition instrument or run matching/housing
estimation. These are Northeast feasibility counts, not national counts.

The existing corrected PSID re-audit is a design benchmark only: the twins
same-birth-year proxy has 34 positive mothers in a 1,974-mother rooms sample
(weighted first stage 0.3668, F=61.8), while same-sex first two has 916 positive
groups in 1,835 mothers (weighted first stage 0.0721, F=7.2). The ACS counts
for twin-like positives and same-sex groups are therefore still `not measured`;
the first implementation must emit them before any estimate.

## Source and identity gates

The current extract field inventory is recorded in
`source_audit_extract27.R:47-50`: `YEAR SAMPLE SERIAL PERNUM HHWT PERWT
STATEFIP PUMA GQ OWNERSHP ROOMS BEDROOMS MOMLOC POPLOC NCHILD NCHLT5 ELDCH
YNGCH RELATE SEX AGE MARST FERTYR RACE EDUC`. The exact bridge key is
`(YEAR,SAMPLE,SERIAL,PERNUM)`; `PERNUM` must never be joined across households.
`build_second_birth_proxy.R:135-166` provides the existing full-household
`MOMLOC == mother PERNUM` join, valid-child-age checks, self-link rejection,
and link-quality audit.

`MOMLOC` is a household relationship pointer and does not establish biological
motherhood by itself. The source gate must retain `RELATE`, require the author
child/mother relationship rule where available, and report social, step, and
adoptive or otherwise non-biological links. `NCHILD`, `ELDCH`, and `YNGCH` are
diagnostics; they cannot repair a failed roster link.

The current extract does not include `BIRTHQTR`; the source contract and
`psid_first_birth_acs_mapping_20260920.md` require its presence/absence to be
recorded separately and prohibit substituting it for `hispan` or changing the
frozen author matching cells. Inspect the author raw schema for `BIRTHQTR`
before requesting any amended source. Even if present, a respondent birth
quarter does not establish a linked child's birth quarter. With the current
roster, a same-age or same-age-zero pair is a **twin-like proxy**, never a
verified biological twin indicator. Do not call it twins without a child
birth-date/quarter field and a valid source concordance.

## Two distinct estimands

### Twin-like first-birth margin (`Twin1`)

The proxy is constructed separately at EACH observed oldest-linked-child age
0:5 in the repeated cross-section (age 0, 1, 2, 3, 4, 5 pooled, not only age
0), among mothers with at least one valid linked child and oldest-child age in
that range, using same observed age (not a birth quarter, since extract27 has
no `BIRTHQTR`) between the first two oldest linked children and no older
linked child. The instrument label is `twin_like_proxy`; it is not a
biological twin treatment, does not infer which child was born first, and is
not the same object as the completed first-birth pseudo-panel's event-0 row —
this design does not follow a single mother longitudinally across ages, it
observes different mothers (or the same mother in different survey years) at
each cross-sectional oldest-child age. Report both the pooled 0:5 estimate and
separate event-age-3 and event-age-5 cells. Same-age contamination (three or
more same-age linked children at the oldest age) is measured and reported
separately, not folded into the positive count. An optional `Twin2` design
around the second-birth proxy is deferred until `Twin1` source and support
gates pass; it would be a distinct third-child margin.

The primary Twin1 readout is a reduced form for housing outcomes at each
observed post-instrument event age 0:5 (pooled and at event ages 3 and 5). An
IV/2SLS treatment may be considered only if a separately measured
additional-child treatment (`>=2` linked children at interview) and a
non-degenerate first stage exist, and is reported as an assumption-dependent
diagnostic with Anderson--Rubin uncertainty, not as identified.

### Same-sex first-two margin (`SameSex2`)

The proposed instrument is equality of the inferred oldest two linked-child
sexes, after a valid age-based order is established. ACS integer child ages do
not identify biological birth order when ages tie; those rows are excluded from
the primary group and reported as an order-uncertain sensitivity. The source
row must have at least two linked children. The primary eligible pool is
`NCHILD >= 2`, with no restriction to exactly two children: conditioning on
`NCHILD == 2` after the instrument is realized is post-instrument selection.
An exactly-two-child sample is a labeled sensitivity only, never the primary
same-sex estimate.

The SameSex2 event clock is the inferred second-child age `a2`, with event year
`YEAR-a2`, and the readout is housing at fixed post-realization times 0:5.
The three-plus-child treatment is `NCHILD >= 3` only when the composition link
and event-time construction support it; do not condition on future child count
or housing. This is a second-birth-to-third-child estimand and must not be
merged with Twin1. The existing coresident proxy documentation
(`second_birth_proxy_design.md`) remains a supplemental matched-pseudo-panel
design, not evidence of a true longitudinal second birth.

## Outcomes, weights, and controls

The outcomes are `ROOMS` (valid 1:27 and 30, cap at 9, codes 0/28 missing),
`OWNERSHP` (1 owner, 2 renter), and `BEDROOMS` as a diagnostic (cap at 5).
Use original maternal-household rows first, one row per source household and
mother identity, with the source person/household weight specified in the
receipt and household clustering. A matched pseudo-panel extension may use the
author matched weight `wgt` only after the composition link and source identity
are proven; it must not treat the 14-million-row first-birth pseudo-panel as
14 million independent instrument units.

The primary control set is limited to variables determined before the relevant
instrument realization: flexible maternal age-at-event terms, race, survey
year, and event-age indicators. Current residence state, education, and
marital status may themselves respond to fertility (interstate moves,
enrollment and marriage timing can follow a birth); any state-adjusted or
current-demographic specification is a labeled conventional-adjustment
sensitivity, not a pretreatment control. Never match or condition on rooms,
bedrooms, ownership, post-event child count, or any other post-treatment
housing/fertility object.

## Identification and inference plan

Estimate the reduced form first for each outcome, with first-stage and treatment
definitions reported alongside it. For SameSex2, the first stage is the effect
of same-sex composition on `NCHILD >= 3` by +5; for Twin1 it is the effect of
the twin-like proxy on the separately defined additional-child treatment.
Report coefficient, robust or source-household-clustered SE, support N, positive
instrument count, and first-stage F. If an IV is attempted, report the
Kleibergen--Paap/weak-IV robust confidence set or Anderson--Rubin interval, not
only conventional 2SLS intervals.

The housing exclusion restriction is doubtful for SameSex2 because child sex
composition can directly change room sharing and housing demand. A positive RF
is therefore not evidence of a valid IV. A same-sex sample restricted to exactly
two children is a post-instrument selection and is not an exclusion test.
Primary interpretation should remain RF/descriptive unless lead review accepts
the first-stage and exclusion discussion.

## Bounded implementation and stop rules

Before a national run, use the existing small Northeast packet or a synthetic
fixture to test: full-key uniqueness; MOMLOC/RELATE link categories; child SEX
availability; age-tie/order flags; `NCHILD >= 2` and `NCHILD >= 3` support;
ROOMS/OWNERSHP/BEDROOMS coding; weight positivity; and one row per maternal
household anchor. The receipt must emit `usable_N`, `Twin1_proxy_positive_N`,
`SameSex2_eligible_N`, `SameSex2_positive_N`, treatment means, first-stage
coefficient/F, and counts by age-tie, missing-sex, and MOMLOC relationship.

Stop before submission if the full key is not unique, the source/sample
fingerprint differs, linked-child sex is unavailable, biological/order status
cannot be distinguished and no proxy label is accepted, or positive support is
not measurable. Report RF only when the first stage is weak; do not submit an
IV run when the weak-IV interval is uninformative or when the housing exclusion
concern is unresolved. No raw-panel laptop read, rematch, or national job is
allowed until lead reviews the fixture receipt and the exact resource estimate.
