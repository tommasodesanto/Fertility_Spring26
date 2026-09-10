# Initial economy and historical shock: decision note

September 10, 2026. Bounded primary-literature check and read-only reconstruction of prior implementation. Recommendation for author discussion; no target, parameter, model or running-job contract changed.

**Subsequent same-day empirical check:** the [concrete parameter map](parameter_target_audit/README.md)
supersedes the limited availability assessment below. The active CPS builder is
2024-only, but raw June 2004/2006/2008 CPS is already local and preliminary early
moments have now been computed. Older natality raw files also exist beyond the
1987-start cache. Initial PSID wealth is directly recoverable from saved annual
components. A fully matched 2007 metro housing target remains unavailable as a
ready input. The new proposal anchors the initial fertility level with a period
rate and fits an announced preference decline for conditional policy scenarios;
it does not claim to identify the decline's cause or activate a new objective.

## Recommendation and tomorrow's claim

Retain the initial-economy-plus-anticipated-transition design. An accurately fitted initial economy followed by an estimated low-dimensional shock is a defensible eventual specification. However, a complete 2007 empirical recalibration is not available as a ready switch. For tomorrow, retain the existing approximate initial economy and describe the current exercise as a conditional historical transition at inherited parameters. Do not describe it as an accurately calibrated 2007 economy, an estimated annual fertility history, a completed PF recalibration, or a new policy result.

If completed quantitative results are required tomorrow, that is a separate feasibility constraint: the current price root and horizon checks remain unfinished. Literature precedent establishes the legitimacy of a design, not certification of these results. Do not switch to a stationary 2023 economy solely to conceal that limitation.

The substantive interpretation should be a calibrated fertility-demand shift. Its fitted size summarizes changes not separately modeled; fitting fertility does not establish their cause. Policy counterfactuals hold that fitted shift fixed. They must not re-estimate it to undo the policy's fertility response.

## What the relevant papers actually do

**Sommer, Sullivan and Verbrugge (2013), JME, sections 3.5 and 6.** They fit structural parameters to baseline household moments, then study an unexpected permanent change in interest rates and down-payment requirements along a perfect-foresight transition. Their baseline moments use several empirical vintages, rather than every statistic from one exact year. This supports our solution architecture and explicit use of stable pooled restrictions; it does not justify calling later housing levels observed 2007 facts. [Author-hosted published paper](https://www.kamilasommer.net/RentPriceRatio.pdf).

**Kaplan, Mitman and Violante (2020), JPE, section III.B and footnote 20.** Their housing model uses observed income/credit conditions and calibrated beliefs about future housing demand. Belief-process parameters are disciplined by expectations evidence and boom/bust statistics; the authors explicitly allow a residual interpretation beyond income and credit. This is precedent for a fitted unobserved driver with a transparent interpretation. Their stochastic expectations and episode construction are different from our fully announced deterministic path; they do not validate our specific fertility wedge. [Author-hosted published paper](https://violante.economics.princeton.edu/sites/g/files/toruqf5621/files/documents/kaplan-et-al-2020-the-housing-boom-and-bust-model-meets-evidence.pdf).

**Sommer (2016), JME, introduction.** The fertility model is calibrated to a later cohort and externally measured income risk, then compared with an economy using earlier risk estimates. This is a comparison between steady states, not a fitted historical transition. It shows that calibrating a modern reference economy is a legitimate different exercise; one must not present it as the historical path into that economy. [Author-hosted published paper](https://www.kamilasommer.net/Fertility.pdf).

The distinction across these examples is between an observed driver supplied to the model and an unobserved driver inferred from the outcomes. Our current fertility shifter is the latter. These examples are methodological precedents, not evidence that taste changes caused the US fertility decline.

## How a genuine two-stage calibration would work

Let theta denote the parameters common across dates, and psi the fertility-preference intercept. First estimate theta and the initial psi level from initial-era household/fertility moments plus explicitly maintained structural micro restrictions. Solve the stationary distribution jointly with initial market prices. Give every free parameter informative moments or an explicit external restriction; replacing the present target system by the scalar fertility normalization would be underidentified.

Next hold those parameters and the inherited household distribution fixed and estimate one shock amplitude delta. The existing fixed-shape specification is

\[
\psi_t=\psi_0+\delta g_t,\qquad
 g_t=\min\{1,\max\{0,(t-2007)/16\}\}.
\]

The path is known in 2007; its level changes gradually until 2023 and is constant afterward. Announcement timing and implementation timing are different objects. One amplitude does not identify an unrestricted annual shock history.

For each candidate delta, solve the entire anticipated equilibrium price path. Compare the resulting fertility observations with their exact data counterparts and choose delta by weighted minimum distance. This is a scalar outer search around a full equilibrium solve, not a sequence of myopic annual inversions. Birth timing, childlessness and completed parity provide additional fit restrictions; one scalar need not match all of them exactly. The assumption of a flat post-2023 preference level affects current decisions and needs to remain visible.

Choose either a consistently measured cohort fit or build correctly exposed age-specific birth rates for a period-history fit. Do not use the current completed-parity target as if it were annual TFR. Holdout housing and population outcomes are useful only after preserving identification and distinguishing imposed demographics from model predictions. In particular, the existing hard housing/child-cost restrictions cannot simply be demoted to validation without replacement identifying moments.

## What our previous work did

The August 16 diary records a dated calibration that normalized old completed fertility (then 2.12), imposed the observed 2007 householder-age distribution and measured the twelve-row objective in 2023. The current normalization is 2.1. Common parameters and the shock coordinate were estimated jointly against later moments, rather than fitting all initial-period moments first. Earlier dated runs solved individual market dates without the current PF continuation values.

The current PF experiment inherits the structural coordinates and shock amplitude. It changes the anticipated equilibrium calculation and has not yet re-estimated them. Historical household totals and age composition are conditioning inputs through 2023, so their historical fit is not an independent success of the estimated fertility shock.

Sources: [August 16 diary](../../../SESSION_DIARY.md); [initial-state implementation](../../../tmp/e5f_matched_pf/code/model/tools/e5f_matched_pf_initial_state.py); [full current fit and restrictions](HORIZON100_PROGRESS.md); [target mapping](overnight_target_mapping_review.md). The dated diary entry begins at line1553 in the reviewed version.

## Why a complete 2007 refit is new work

- The authoritative CPS builder uses June2024 women40-44. An initial-period completed-cohort target requires an earlier supplement and an explicit cohort choice.
- The timing cache begins in1987. It misses early births for cohorts already around age40 in2007; even its1970-74 reference cohort is explicitly marked incomplete at young ages.
- The active42-MMS ACS cache is pinned to2012-2023. A2007 version needs earlier geography treatment and fresh target uncertainty. An existing national housing-path series has different geography and a2007 room-topcode caveat; it cannot simply replace the current metro target values.
- Reviewed pooled PSID/event-study restrictions may be retained as explicit structural approximations if adopted; they are not automatically2007 observations. The Sun-Abraham regression does not need reopening.

Relevant builders: `code/data/cps_fertility/README.md`, `code/data/nchs_natality_timing/timing_target_metadata.json`, `code/data/moment_standard_errors/build_active_acs_room_target_receipt.R`, and `code/data/Spatial_aggregate_withmicrodata/build_national_householder_housing_path.py`.

A useful immediate decision is therefore to adopt the initial-economy-plus-one-shift architecture while being explicit that tomorrow's numerical evidence, if still at inherited parameters, is conditional and provisional. A later fully initial-era calibration needs its own target contract and identification map. No silent date relabeling, target dropping or new shock process is required for this methodological decision.
