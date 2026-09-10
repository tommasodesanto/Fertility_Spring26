# Adversarial check: initial stationary calibration, then one preference shift

September 10, 2026. Bounded review using only the already checked fertility receipts and mathematical distinctions. No new data scan, model run, target change or architecture adoption.

**Verdict:** the evidence does not rule out a diagnostic initial stationary fit. It does rule out promising an accurate simultaneous fit of all early stock and flow observations before testing their model counterparts. A joint PF fit with the same stationary initial-distribution restriction does not remove that problem.

## Strongest objection

The initial CPS completed mean, childlessness and one-child share describe older cohorts' accumulated histories; initial period first-birth timing describes current births. An environment that is approximately constant around 2005--2007 need not have been constant over the preceding twenty-five years that produced those cohorts. Approximating both by the same stationary conditional parity distribution is therefore a substantive historical restriction, not just choosing a convenient benchmark year.

There are two different versions of the proposed objective:

1. **Early completed mean/P0/P1 and period timing, without imposing the early period-TFR level:** no simple accounting contradiction has been established. A stationary solution may fit these objects tolerably. But its initial implied period TFR must be reported against the observed 2.0605; it cannot silently acquire that observed initial rate. Fitting a later period rate then estimates a decline from the model's initial fertility level, which may differ materially from the observed initial-to-late decline.
2. **The same objective plus early period TFR 2.0605:** this adds real tension. Early capped CPS CEB is 1.856608 and uncapped CEB 1.878384. Matching age ranges, counts, population and exposures is essential before invoking equality, but the gap is not plausibly removed just by capping and a short remaining reproductive tail under the frozen early rates. The entire 40--49 integral under those rates is only 0.047875. An optimizer may compromise, or fit through unusual late/higher-parity fertility; successful numerical convergence is not proof the stationary approximation is credible.

The earlier illustrative calculation using only P0 and P1 showed that a 2.0605 completed mean would require about 2.997 births among women with at least two, versus the old capped observed 2.673. That is not proof of impossibility. Once the old completed mean itself is included, however, we cannot present the higher-parity adjustment as simultaneously matching that mean. The exact age40--44/topcode observer must show the remaining gap. The model's 3+ representative count 3.602 also differs from the 2006 observed capped top-group mean 3.47411.

**Thus a diagnostic is justified; a presumption of an accurate initial fit is not.** Its purpose should be to identify whether common stationary policies can approximate the early observations without distorting the observed age/parity pattern. A weighted compromise across incompatible restrictions must remain visible. New weights should not be enlarged merely to make that compromise look statistically acceptable.

## Does joint PF estimation fix it?

No, if it retains the same stationary pre-announcement initial state and the same early measurement contract. With a time-invariant fertility schedule, matched period and lifetime cohort fertility are integrals of the same schedule. Re-estimating common parameters using later observations changes the schedule and the compromise; it does not remove this restriction or manufacture different early cohort histories.

Joint estimation can improve statistical discipline of common parameters and reveal which parameters the later transition informs. That is a valid reason for it, distinct from fixing initialization. If it obtains a better total objective by worsening the initial fit, that tradeoff should be shown explicitly. To remove the historical restriction would require a different initial-state/prehistory construction, not merely a joint optimizer.

A further timing point matters: under full announcement in 2007, households can change their 2007 decisions in response to later conditions even though \(\psi_{2007}\) equals the old intercept. Therefore observations from 2007 are not automatically pre-announcement stationary choices. Calling the benchmark “2005--2007” is acceptable only as an explicit initial-era empirical approximation, with 2007-contaminated flow observations handled separately or the announcement/calendar convention stated. The recovered early CPS pool is 2004+2006, and the checked period timing/rate pool is 2003--2006; neither becomes a 2005--2007 estimate by relabeling it.

## Minimum conditions for recommending the bounded initial diagnostic

- Define the actual initial windows and the period, cohort, age and parity operators. State whether initial period TFR is targeted or only a check.
- Include a separate full initial fit panel for CEB, P0, P1|mother, period timing and period TFR, with no silent substitution of one fertility definition for another. Inspect late-age births and the high-parity tail if these reconcile aggregate targets.
- Treat the new initial objective and later validation set as an explicitly new scientific contract; retain the old production results and identify any former hard moment moved to validation. The lead must assess replacement identification across all common parameters.
- Freeze common parameters for the one-shock stage only after the initial fit has earned that interpretation. The later shock estimate is conditional on that baseline; it cannot repair a poor initial distribution.

Existing evidence: `early_fertility.md`, `fertility/initial_fertility_coherence_arithmetic.json`, `fertility/cps2006_age_shape_check.json`, and the parent parameter-target audit's original fertility receipts. This check makes no parameter or target choice.
