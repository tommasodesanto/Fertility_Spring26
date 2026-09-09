# Households with children: what the current experiment establishes

September 9, 2026. Bounded review requested by Tommaso, starting 08:38 EDT. This is an accounting and measurement review of the completed rebated-property-tax experiment. No calibration, target, numerical gate or policy specification was changed.

**Subsequent author correction:** the intended final policy exercise is perfect foresight with the developed person/head demographics. The recent temporary-equilibrium tax packet did not carry that work forward. A verified paired perfect-foresight result existed under the September1 sensitivity; a later parameter configuration had an unresolved baseline convergence failure. The September3 conversation explicitly made the person-demography/H128 specification primary. Resolving this handoff now takes priority over the measurement sequence proposed below. See [perfect_foresight_result_audit.md](perfect_foresight_result_audit.md). The statement that2023 lies on the historical calibration path does not establish perfect-foresight household expectations.

**The model can count its households with dependent children consistently. The unresolved problem is how that category maps into observed families and into future resident persons and household formation.** Those are economically substantive assumptions, not just alternative labels for the same number.

## 1. Four different objects

A birth is a flow during a period. Households with children are a stock at a date. Total households count separate household decision units. Resident population counts people, including children and adults who are not household heads.

For example, ten first births among 100 existing childless households leave the economy with 100 households, ten households with children and ten additional children. A second birth normally adds a child without adding another household with children. A later birth can also return a household whose previous children have left to the household-with-children category.

The current state records lifetime parity, \(n\), separately from the number of children currently in the dependent state, \(m\). Consequently, a household with \(m=0\) may be a never-parent or an empty nester. These groups should not both be called childless when discussing lifetime fertility.

If \(g_t(s)\) is household mass across economic and family states, the three internal stocks are

\[
H_t=\sum_s g_t(s),\qquad
F_t=\sum_s\mathbf 1\{m(s)>0\}g_t(s),\qquad
C_t=\sum_s m(s)g_t(s).
\]

Here \(H_t\) counts households, \(F_t\) counts households with dependent children, and \(C_t\) counts the model's explicit dependent-child units. Neither \(H_t+C_t\) nor \(2H_t+C_t\) has been established as resident population. The couple normalization in utility is not a measured number of resident adults per ACS household head.

## 2. What last night's experiment says

Both policies rebate their own property-tax revenue equally per household. The baseline tax is 1% annually and the reform tax is 2%. The table uses the model's dependent-child category and expresses counts per 100 households in the common 2023 economy; these are not US population totals.

| Year | Baseline households with children | Reform households with children | Reform effect on their number |
|---|---:|---:|---:|
| 2023 | 36.4149 | 36.4488 | +0.0929% |
| 2043 | 25.6914 | 25.8189 | +0.4960% |
| 2063 | 19.5240 | 19.6991 | +0.8972% |

In 2063 the corresponding shares of contemporaneous households are 26.3531% and 26.5366%, a difference of **0.1835 percentage points**. The count effect is larger than the share effect because total households are also 0.1995% higher under the reform. Among young households, using the established model age nodes 26, 30 and 34, the number with dependent children is 1.4451% higher under reform in 2063.

The positive policy difference coexists with a large baseline decline: dependent-child households fall **46.38%** from 2023 to 2063, while total households fall **25.91%**. The paths are closed household-renewal diagnostics, not forecasts of American families. Matching completed fertility does not by itself validate this household composition or transition.

On impact, per 100,000 initial households, the reform produces about **42.73 additional explicit births and 33.84 additional households with dependent children**, with no additional total households. The latter consists of 28.66 additional first births and a net 5.18 returns from positive parity with no current dependent child. These are net policy differences, not identified individual switchers. Explicit births here exclude the separate 3+ top-code adjustment used in headline birth totals.

The source partition identities and all 44 case/date/age-group records reconcile. The audit checked 94 input files and 750 conditions, including the later entry-handoff check; the lead independently rechecked their hashes, stock partitions and headline effects. See [counts_readout.md](counts_readout.md), [counts.csv](counts.csv) and [counts_audit.json](counts_audit.json).

## 3. Why the child count is a proxy

Each dependent child leaves independently with probability \(2/9\) per four-year step. This gives 18 years at home **in expectation**, not departure on the eighteenth birthday. The state has no actual child ages, sibling birth dates or newborn protection in the next calendar advance.

Conditional on the household surviving, this means 22.2% leave the dependent category by the first four-year advance, while 28.5% remain after twenty years. These properties were built into the explicitly authorized August 5 approximation. They are not evidence that a deterministic age-18 rule was implemented incorrectly. They do mean that the category cannot be interpreted literally as children under 18.

This matters for economics: the rule determines how long children affect housing demand and consumption needs, how siblings overlap at home, and when households become empty nesters. An average duration of 18 years does not determine the correct age profile of households with children. We have not quantified how much changing this approximation would alter the room response or the policy effect; it would be premature to attribute the first-birth rooms miss to this rule.

There is also a top-code limitation. The model explicitly represents up to three children, while the 3+ fertility category is measured as 3.602359 children. The additional 0.602359 enters adjusted birth reporting and future renewal when a household reaches the top category. It does not add explicit dependent states, housing needs or dated births. A consistent resident-child ledger must address that bridge explicitly.

## 4. Two empirical group comparisons need reconciliation

The empirical target estimates are traceable; the issue is their model counterparts.

| Target | Empirical group | Current model group |
|---|---|---|
| Ownership gap, heads ages 30–55 | Oldest resident own child under four, versus no resident own children | Any dependent child, versus no previous birth |
| Rooms gap, heads ages 30–55 | Three or more resident own children versus one or two, requiring at least one under 18 | Three dependent children versus one or two |

The ownership comparison differs on both sides: recent parents are not all parents, and households without resident children include empty nesters. The rooms comparison uses the ACS number of resident own children of any age. A head with children aged 22, 19 and 10 enters its empirical 3+ group, although only one child is under 18.

The ownership discrepancy was already documented in the July code audit and remained an outstanding reconciliation in the September 7 status. The scoped history check found awareness, but no explicit decision accepting that exact approximation. Authorization of independent child maturation is a separate decision. This review is not evidence that the author never discussed the target elsewhere.

There is no direct target for the overall number or share of households with children in the current twelve-moment calibration. CPS completed fertility concerns children ever born to women; it does not independently validate today's stock of household heads with resident children. See [empirical_alignment.md](empirical_alignment.md) for the builders, samples, definitions and source anchors. Nothing here drops or reweights a target.

A small new tabulation from the pinned ACS cache makes the definition issue concrete:

| 2023 household-head age range | Any resident own child | At least one resident own child under 18 | Oldest resident own child under four |
|---|---:|---:|---:|
| 18–85 | 39.56% | 26.92% | 3.63% |
| 25–34 | 31.59% | 31.55% | 11.15% |
| 30–55 | 57.34% | 48.39% | 5.46% |

These are HHWT-weighted shares in the 42-metropolitan-area MMS household-head cache, restricted to owner/renter households with positive rooms. They are **not national estimates**. All child-presence classifications are observed in these samples. The model's 2023 all-age dependent share is 36.41%, and its young-node share is 47.66%; these are useful warning comparisons, not newly certified calibration gaps. Child definitions, geography/population composition and the young annual-age window need alignment before interpreting their differences as model error.

The cache hash matches its August 17 receipt; the extraction also reproduces both existing empirical ownership and rooms target values. [acs_validation.csv](acs_validation.csv) contains all 24 aggregate rows, including pooled 2012–2023 counterparts; [acs_validation_receipt.json](acs_validation_receipt.json) records denominators and limitations. One cache read used one thread and took about 12 seconds. The 9.92 GB raw extract was not reread. No uncertainty estimate or new target weight was constructed.

## 5. Children leaving home and new households currently follow different rules

Current births enter a separate queue that creates adjusted births divided by 2.1 new households twenty years later. For example, 21 adjusted births schedule ten entrants. The queue does not use the dependent-child departures described above. A child can therefore disappear from the dependent state long before its associated birth cohort contributes to household entry, or remain in that state afterward. The present model has no intervening ledger of non-head adults linking the two.

The household identity itself is sound: surviving households plus entrants determine next period's household count. The unresolved object is the translation from actual persons and their family membership into household heads. The initial population uses empirical household age masses; its continuation switches to the birth queue. In the current baseline, the youngest household cell rises from 0.014337 in 2023 to 0.063851 in 2027, a factor of 4.45. The latter exactly equals the scheduled entrants inserted in that advance. This is a verified handoff across different rules, not evidence that current births immediately create households or that the arithmetic is broken.

Total households initially rise 1.23%, peak in 2027, and then decline at every date. Households with dependent children fall 6.39% during 2023–2027 and decline throughout the reported path. That distinction further limits an interpretation of the large long-run decline as a smooth demographic forecast.

We are not starting the population work from zero. Tommaso explicitly requested a credible persons forecast on August 26. Existing code under `code/model/demographic_transition/` tracks annual persons by age and sex, survival and migration, and maps them into heads using fixed 2023 ACS age-sex headship rates. A separate experimental driver already invokes that system; last night's driver does not.

That existing system is reusable, but its fixed headship rates are an accounting assumption, not estimated policy-sensitive household formation. Its economic bridge rescales households within age cells while preserving their conditional family composition. It therefore does not automatically repair the dependent-child state or the empirical family-group comparisons. Older population numbers used different experiments and must not be attached to the new tax results.

## 6. Recommended next decision

**Discussion clarification, September 9:** Tommaso confirms that the expected 18-year duration was chosen as an approximation to children living at home. Retain that approximation while resolving measurement. The earlier recommendation to make under-18 children the primary definition was a proposal, not an adopted author decision. “Children at home” remains the intended economic object; show any-age resident children and under-18 resident children as explicit empirical comparisons, without silently equating either to recent parenthood. State household-head ages and geography, and use household weights.

First align and validate those measurement definitions. If the child-age profile is materially wrong, test a child-age representation in an isolated version; an age of the youngest child alone is insufficient to know when each older sibling leaves. Preserve the current specification as the comparison and assess runtime before a full recalibration. Merely changing the label or adding an aggregate population series will not repair this mechanism.

For resident-population reporting, reuse the existing annual person/headship implementation. An accounting-only calculation using saved births can answer a conditional population question; feeding revised household counts back into housing markets requires a new equilibrium path. Keep that distinction explicit.

The immediate priority is thus a **consistent family definition and validation**, followed by a decision on the child-age and household-formation approximations. Another large calibration search with these mappings unresolved would leave the central question unanswered.

## 7. Is 2023 currently a steady state?

**No.** The active calibration already places 2023 on the simulated 2007–2023 path. Each candidate first constructs an old stationary economy with completed fertility normalized to 2.1. Its age distribution is then reweighted to the observed 2007 household-head distribution. Preferences change over 2007, 2011, 2015, 2019 and 2023; wealth and family states propagate, and the dated household totals and age marginals are imposed from Census/ACS. The ownership target is measured on the resulting 2023 cross-section. The first-birth housing target uses the separate 2019–2023 matched branch.

This does not impose a stationary 2023 distribution. It does retain static expectations: at each date, households choose using current prices and primitives as if those will persist. Actual market prices change along the computed path. A fully anticipated future-price transition is a separate specification. The post-2023 switch from imposed demographic age masses to the birth queue is another distinct assumption.

The historical fit was difficult, but the current method does not force 2023 to be stationary to improve that fit. Nor are all empirical targets 2023 observations: the ownership target pools 2012–2023 data. Its time window should be reconciled with model measurement alongside the parent/control definitions. Do not silently substitute a 2023-only estimate.

The bounded follow-up reproduces a pooled recent-parent ownership gap of 16.7662pp, versus 15.2117pp for any resident own child and 15.0155pp for at least one under-18 own child, all against the same no-resident-child control. The latter two are diagnostic alternatives, not adopted targets, and neither fixes the model's never-parent comparator by changing the scalar alone. See [ownership_target_followup.md](ownership_target_followup.md) for dates, definitions, the corresponding 2023-only rows and positive source verification.

Source check: frozen `run_e5f_transition_calibration.py` lines 2031, 2241, 2344, 2387 and 2600; `run_dynamic_population_transition.py` lines 469 and 1783. The selected summary also explicitly describes remaining targets as the 2023 transition cross-section.

## Evidence and limits

The verified claims form separate pieces:

| Object | Established | Limit |
|---|---|---|
| Current simultaneous fertility nests | Implemented and passed recorded numerical/choice checks | Not verified through the PF transition adapter |
| Current calibration | Reproduced historical temporary-expectations fit against the retained objective | Not a jointly calibrated PF history; parent-group measurement remains unresolved |
| Earlier perfect-foresight work | September1 paired policy sensitivity passed; code and results preserved | Different configuration; later September3 baseline did not converge |
| Latest rebated-tax paths | Both temporary-equilibrium paths pass their declared numerical gates | Neither the current PF exercise nor a resident-population or welfare result |

The operational divergence predates the final simple nest: September5,17:13 New York, the task chose the older policy route and called the person-demographic implementation a separate branch; September6 retained temporary expectations explicitly. The latest extension was submitted September8 at23:24:47 New York under sourcee483254a. The documented rationale was a reproduced pipeline, comparability and bounded runtime. No new PF failure demonstrated a need to abandon PF under the new nest. See [the exact conversation chronology](perfect_foresight_handoff_audit.md). A launch statement that assumptions stayed fixed referred to continuity within the diagnostic pipeline, not to the intended PF/person-demographic specification.

The numerical interpretation uses the frozen `tmp/e5f_fertility_nest_compute_20260907a` source, not unrelated working-tree changes. [code_audit.md](code_audit.md) and [code_audit.json](code_audit.json) provide exact paths and hashes. The current transition contract hash is `ce80b6ad241bec9556d2f3d3cdccb9f89e0ccfdd5b2e80ddcd167ed5156b5b68`.

The existing [morning PDF](../../../pdf/rebated_property_tax_morning_20260909.pdf) contains the complete twelve target fits and eleven parameter estimates/bounds. This review qualifies their group interpretation; it does not alter the saved objective or verified arithmetic. No raw-data re-estimation, model solve, causal claim about the effect of changing maturation, or new national demographic forecast is claimed here.
