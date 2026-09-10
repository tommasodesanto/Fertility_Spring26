# Can the pre-2007 fertility evidence support an approximate stationary benchmark?

Read-only design research, September 10, 2026. No calibration, target/weight change, model change, cluster action or git action. This note supports an architecture decision; it does not make one.

**Answer: an explicitly approximate initial benchmark is defensible, but the data do not establish literal stationarity. Period TFR2.0605 and early completed-fertility1.85661 can coexist in actual data because they describe different women and histories. They cannot simply be assigned to the same stationary model statistic.** The proposed early childlessness/one-child shape constraints can coexist algebraically with a2.0605 fertility mean, but fitting them can force a larger high-parity tail than the observed old CPS distribution.

## What each observation describes

- **Period TFR2.0605:** equal-year average of published national2003--2006 age-specific fertility rates, summed over reproductive age and converted to births per hypothetical woman. It is not the realized lifetime births of the women who were40--44 then. The underlying rates are2.0425,2.0455,2.0535,2.1005. This pool avoids the2007 peak2.1220 but is not perfectly flat. Definition and values: [NCHS Births: Final Data for2007, printed pp6 and25, Table4](https://www.cdc.gov/nchs/data/nvsr/nvsr58/nvsr58_24.pdf).
- **CPS pooled CEB1.856608:** population-weighted mean of min(live births ever had,5) among women40--44 in June2004 and June2006. Those women were born approximately1959--1966, depending on exact survey date/birthday; their birth histories mostly precede the proposed initial period. Uncapped mean1.878384. The age40--44 endpoint is near-completed fertility, not exact age50 completion. Childlessness.198279 and exactly-one conditional on motherhood.213655 use those same women and histories.
- **Initial period first-birth timing25.976264 and share30+.249278:** counts of observed first births at all maternal ages12--49 pooled2003--2006, transformed to the established four-year midpoint labels. These are distributions among current first births. They are not the timing history of the old CPS mothers. Their count weighting also reflects the current female age composition; a model must use its matching dated births distribution, not an arbitrary cohort hazard-weighted mean.

The early CPS data use IPUMS FREVER and FRSUPPWT. [IPUMS documents the fertility supplement weight](https://cps.ipums.org/cps-action/variables/446065), its final-person-weight relationship and four implied decimals in raw extracts. The adjacent local loader supplies the exact schema. Existing extraction receipts hash the selected June partitions. Full official-raw crosswalk/source receipt and uncertainty remain outstanding.

## Stability across nearby observations

| CPS window, women40--44 | CEB capped5 | Childless | Exactly one among mothers |
|---|---:|---:|---:|
|2004|1.871292|.192839|.215110|
|2006|1.841533|.203864|.212141|
|2008, after announcement|1.871399|.177733|.223720|

These movements are moderate for completed fertility and the one-child share, but no sampling uncertainty has been constructed. Do not call differences statistically significant or declare stability from two pre-announcement surveys.2008 is a sensitivity, not an unambiguous pre-shock observation.

The first-birth timing cache supplies a stronger warning against a generic postponement explanation: midpoint mean age falls from26.0474 in2003 to25.8616 in2006 and25.8398 in2007; the age30+ share falls from.257806 to.238562 and.235928. Across1999--2008 the corresponding means first rise, then reverse, then recover slightly. Thus the selected initial window is reasonably close in average timing, but not literally time-invariant. This is first births only, with coarse midpoint bins and population count weighting; it does not identify all-order tempo effects.

One additional bounded read of the already-hashed June2006 partition found:

| Age band | N | CEB capped5 | Childless |
|---|---:|---:|---:|
|35--39|4908|1.850229|.189445|
|40--41|2110|1.832847|.201039|
|42|966|1.866480|.198961|
|43--44|2196|1.839093|.208813|

These are different birth cohorts, not successive observations of the same women. Sampling and cohort composition prevent interpreting them as a single stationary lifecycle. The available first-birth age/year fields are blank for all4215 mothers40--44 in that2006 partition, so they cannot provide a shortcut to old cohort timing. For mothers40--44 the observed top-parity mean E[min(N,5)|N>=3]=3.474114; the model's inherited3+ representative count3.602 is a separate normalization that also needs observation-unit care.

## What literal stationarity would imply

Let f(a,t) denote births per woman-year at age a and date t. The period index is

\[
TFR(t)=\int f(a,t)\,da,
\]

while completed fertility for birth cohort b is

\[
CF(b)=\int f(a,b+a)\,da.
\]

If the age-specific fertility schedule is constant over time and the population, exposure, birth-count and age-range definitions match, both equal the integral of that common schedule. This is the identity relevant to a fully stationary benchmark. It does not directly equate TFR to a capped stock observed at ages40--44, and it does not make a current period and old completed cohort empirically interchangeable.

The numeric discrepancy is0.203892 births between2.0605 and capped1.856608. Removing the CPS cap explains0.021776, leaving0.182116. At the frozen2003--2006 age-specific rates, the entire age40--49 fertility integral is only0.047875 births. This is a deliberately generous amount of remaining reproduction for someone already40; it is **not** a bound on actual future births under changing rates. It shows that capping and incomplete completion alone are not an adequate mechanical reconciliation under the proposed frozen schedule. Cohort histories, population/source differences and actual nonstationarity remain substantive.

The current proposal does not hard-target the old CEB mean. Its subset is not algebraically contradictory: with P0=.198279 and P1/(1-P0)=.213655, P1=.171292 and P(N>=2)=.630429. If completed mean were mechanically2.0605, E[N|N>=2] would need2.996701. The old capped CPS value is2.673284. Both are feasible parity means, but the difference predicts where a stationary compromise may appear: more high-parity fertility. This calculation is a diagnostic illustration; exact model age/topcode operators must replace it in a production assessment.

Also, the existing later CPS capped mean1.918425 is **higher**, not lower, than the early1.856608. The observed period decline and the comparison of these completed cohorts point in different directions. They should not be narrated as the same decline.

## What demographic methodology supports

1. **Bongaarts and Feeney(1998), On the Quantum and Tempo of Fertility.** The original authors explain that changing birth timing can move period rates independently of fertility quantum; postponement depresses period rates and advancement raises them. Their adjustment rests on explicit assumptions about the source of fertility change. I verified the Population Council working-paper abstract and its link to the published PDR24(2):271--291 version; the full download was not accessible in this check. This supports distinguishing timing and level, not inserting a mechanical adjustment into our calibration. [Author institution source](https://knowledgecommons.popcouncil.org/departments_sbsr-pgy/248/).
2. **Kohler and Ortega(2002), Demographic Research6(6):91--144.** Section1, printed pp92--93, explains why even tempo-adjusted period fertility does not identify completed fertility of actual cohorts without assumptions about future timing and quantum. Section3.1, pp102--103, defines age/parity hazards using women of the relevant parity as exposure. Their framework reinforces using matched period and cohort observers and warns that an unconditional birth-order rate is not a parity-transition hazard. [Full primary paper](https://www.demographic-research.org/volumes/vol6/6/6-6.pdf).
3. **Schoen(2004), Timing Effects and the Interpretation of Period Fertility.** The February PAA paper, pp3--6, defines period TFR as a synthetic-cohort object and compares it with actual cohort experience; the section The Bongaarts--Feeney Timing Adjustment, pp6--9, discusses why timing adjustments do not automatically recover completed fertility. The accessible source is the author's conference version, not the final Demography41(4):801--819 typeset article. [Primary conference paper](https://paa2004.populationassociation.org/papers/40293).

We should not assert that postponement explains this particular old CPS/period gap merely because tempo effects exist. The direction of the early first-birth timing change, differing cohorts, topcoding and exposures must be investigated empirically. The papers justify transparency, not automatic target reconciliation.

## Candidate initial designs and their costs

**A. Period-oriented synthetic initial benchmark.** Match early period TFR and the early dated first-birth distribution, with a declared approximately stationary environment. Use early cohort parity/childlessness as additional restrictions or explicit approximation checks, with actual age-window observers and visible fit gaps. This addresses the requested later decline in current birth rates. It does not claim to reconstruct every initial cohort's historical fertility. If both old parity shape and period level are hard, check the high-parity compromise and local identification rather than assuming the combination fits.

**B. Cohort-stock-oriented benchmark.** Anchor the initial level to the old completed-parity stock and shape, then use period TFR as a separate validation object. This better represents the older observed cohorts, but starts from a lower fertility level and does not itself match the2.0605 initial period index. Initial cohort timing needs the older NCHS raw builder extension. It is an alternative target contract, not a free relabeling.

**C. Historically heterogeneous initial state.** Retain a current period fertility environment while carrying observed old cohort parity histories in the initial state. Then period2.0605 and old stock1.85661 need not describe the same stationary distribution. This more faithfully separates the objects, but requires a specified empirical initial-state construction and changes the interpretation away from a fully stationary initial distribution. It is not accomplished merely by reweighting the old model's age margins.

The author/lead must select the scientific approximation; this note does not adopt A, B or C or remove any target.

## Smallest remaining work after that decision

For A, no additional large fertility data download is needed: formalize the existing CPS sample/capping/weights, compute joint uncertainty for old childlessness/parity, choose window-based or sampling weights for new period timing, and verify the model's female exposures, birth-order/topcode accounting and dated age windows. Produce a single initial-fit table that shows period and cohort observations separately, including the old completed mean even if diagnostic. For B, add the pre1987 NCHS age/order/sample checks and old cohort operator. For C, an explicit empirical conditional parity/wealth-state construction is additional work and should not be promised overnight.

## Evidence

Existing source pins and candidate numbers: `../fertility/fertility_availability.json`, `../fertility/cps_initial_window_candidates.csv`, `../fertility/nchs_period_timing_candidates.csv`, `../nchs_period_tfr_sources.csv` (paths relative to this note's parent need one level up as shown in repository navigation). New research receipts: `fertility/cps2006_age_shape_check.json` and `fertility/initial_fertility_coherence_arithmetic.json`. No plots or illustrations were generated.
