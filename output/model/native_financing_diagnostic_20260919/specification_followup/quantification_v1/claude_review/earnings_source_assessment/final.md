# Earnings-source assessment (Fable, 2026-09-21)

## Recommendation

**Baseline source/concept.** Household gross labor earnings of the PSID reference person plus spouse, the object the model already prices (working-age income is gross earnings times the age profile times $z$ times $(1-\tau_{pay})$, `code/model/intergen_eqscale_seq_optimized/solver.py:6830-6846`, and the wealth/earnings target grosses back up by $1/(1-\tau_{pay})$, `solver.py:327-335`). Process: age profile plus persistent AR(1) plus iid transitory, no permanent type, estimated on the project's own PSID panel, but only after three repairs: (i) a sourced measurement-error correction to the transitory variance, (ii) HSV progressivity compression $(1-\tau)$ applied to both innovation standard deviations to represent the tax-and-transfer insurance the flat payroll tax cannot deliver, (iii) the age profile taken from the same sample (the current points have no source, `docs/model/structural_model_review_fable_20260915.md:850`; defaults at `parameters.py:864`).

**Why not the literature import.** Every published candidate measures a different object from the model's. Floden–Linde is an individual hourly wage; Sommer–Sullivan is a round-number household productivity convention; BGM is equivalized disposable income including transfers and pensions. Importing any of them requires undocumented conversions that are larger than the estimation uncertainty in our own sample. The author's doubt about the pipeline is best answered by a bounded validation of that pipeline, not by switching to a wage process.

**Robustness alternative.** The already-estimated fixed-effect version (fixed effect plus AR(1) plus transitory, objective 14.3 versus 39.5 without the fixed effect), labeled "permanent heterogeneity robustness," using the existing E6b machinery.

**Cost of no types, stated plainly.** The no-type AR(1) is rejected by the long-lag autocovariances: it underpredicts lags 16–28 by 0.04–0.07 (roughly 2–3 bootstrap s.e.) and absorbs about 0.39 of fixed log variance into a near-unit-root state ($\rho_a = 0.970$, stationary persistent variance 0.69). Economically, heterogeneity households know at entry is being modeled as slowly unfolding risk, which overstates precautionary saving and uncertainty at young ages. This is a convention choice for state-space economy, not a fit-based choice. Report the misfit openly.

## Evidence table

| Source | Income concept, sample | Stochastic components | Compatibility with model object |
|---|---|---|---|
| Floden–Linde 2001, Table IV p.421 ([PDF](https://martinfloden.net/files/Floden%20M%20-%20Linde%20J%20-%20RED%202001.pdf)) | Individual hourly wage of PSID heads in the labor force 1988–92, SEO excluded, N=1789 (p.416–417); pre-tax; permanent part removed by 1988 regression on age, sex, education, occupation ($\sigma_\psi^2=0.1175$, p.420) | $x=\psi+z+\xi$, $z$ AR(1): $\rho=0.9136$ (0.0090), $\sigma_\varepsilon^2=0.0426$ (0.0048), $\sigma_\xi^2=0.0421$ (0.0039) treated as measurement error | Wage, not household earnings; no transitory risk; five-year window; fixed effect explicitly present in the data model |
| HSV 2017 QJE | Household pre- to post-government income, PSID 2000–06 + TAXSIM | $\tau_{US}=0.181$ (s.e. 0.002); CBO 0.200 | Pinned July 23 from the published paper (`memory/AGENT_MEMORY.md:1676-1679`); paywalled, not re-verified here |
| Sommer–Sullivan 2018 AER, Table 1 p.253, text p.254 ([PDF](https://kamilasommer.net/Taxes.pdf)) | Household labor productivity, pre-tax (model has payroll tax 0.076 and a bracket income tax) | $\rho_w=0.90$, $\sigma_w=0.20$, set "in the range" of Card, HSZ, HSV 2010; seven states; "based on SSV 2013" | A convention, not an estimate; no transitory; verified |
| Sommer 2016 JME (FEDS 2014-32 pp.11, 16–17, [PDF](https://www.federalreserve.gov/pubs/feds/2014/201432/201432pap.pdf)) | Household wage, annual, age profile from 2004 CPS family earnings over couple hours; labor supply endogenous | AR(1) $\rho=0.95$, $\sigma_\varepsilon=0.21$ (upper end, Meghir–Pistaferri), iid $\sigma_\nu=0.17$; $\varepsilon_1=0$; no fixed effect | Closest fertility precedent for the functional form; values are chosen, not estimated |
| Boar–Gorea–Midrigan (NBER w23345 pp.8, 17–18, [PDF](https://www.nber.org/system/files/working_papers/w23345/w23345.pdf)) | Disposable income: wages + SS + pensions + UI + transfers minus federal/state income tax, CPI-deflated, OECD-equivalized; PSID 1999–2007; quarterly model | $y=\lambda_t z e$, $z$ AR(1), $e$ iid, no type; moments: var, autocov at 2 and 4 years, s.d. of 2-year growth; transitory s.d. scaled down 55% for measurement error | Right functional form and estimation template; wrong income concept (equivalization interacts with child costs; includes pensions the model provides). Table 2 numbers quoted in the Sept 20 review were not re-verified |
| Own PSID packet 2026-07-27 (`code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/`) | RP+spouse gross labor earnings, ages 25–60, 1984–2019, 89,936 person-years, age/year FE residualized | With FE: FE var 0.393, $\rho_a=0.886$, persistent var 0.332, transitory 0.310, objective 14.3. Without FE: $\rho_a=0.970$, persistent 0.693, transitory 0.344, objective 39.5. No measurement-error adjustment | Exactly the model's object and the wealth/earnings denominator; sample starts at 25 while the model enters at 18 |

## Historical timeline (author decisions versus implementations)

- **July 16 (author):** D1 = literature AR(1), Sommer–Sullivan $\rho=0.90$, $\sigma=0.20$; own PSID estimates relegated to appendix evidence; "not a license to shrink income risk" (`docs/model/m5_recalibration_contract_20260716.md:11-30, 109-112`).
- **July 22 (author directive):** finalize the process "properly, with the proper variance"; open items were concept, source, and 5 versus 7 states (`docs/model/eqscale_calibration_reconciliation_20260722.md:647-670`).
- **July 23 (author-confirmed):** Floden–Linde persistence and innovation variance, HSV $(1-\tau)$ scaling, 5-state Rouwenhorst (`eqscale_calibration_reconciliation_20260722.md:597-606`; `docs/prompts/HANDOFF_fable_eseries_20260723.md:151-164`). Wired in `externals.py:6-32` with annual innovation s.d. $0.2064\times0.819=0.169$.
- **July 27 (implementation, no author approval located):** the PSID decomposition was built "for E6b" to impose a fixed-effect variance externally (packet README lines 21–22; `memory/daily/2026-07-27.md:41-55`). E6b multiplies the FL×HSV chain by three permanent levels (`e6b_profile.py:42-51`) and is wired into the frozen baseline via `e5f_income_entry_profile.py:30`. The July 28 memory records a recommendation "for author review"; the Sept 14 ledger asks to reconstruct the adoption history, so I treat E6b adoption as undocumented.
- **Sept 14 (author):** prefers removing the permanent component (`docs/model/ACTIVE_DECISION_LEDGER.md:301-335`).
- **Sept 19 (implementation):** the diagnostic candidate "reuses directly" the nested no-FE fit of the E6b packet (`earnings_candidate/README.md:9-28`). No author decision to switch the source from literature to own PSID exists; the switch happened because that packet was the only in-house no-type estimate at hand.
- **Sept 21:** ledger E2 calls own PSID the "established source" (`ACTIVE_DECISION_LEDGER.md:126-139`). That overstates the record: the two explicit author source decisions were literature imports.

## Verdict on the July choice

The July 23 pair was a coherent repair of an incoherent pairing (M5 persistence 0.9602 with a 0.20 innovation), and the HSV rationale is sound: a flat payroll tax shifts levels and cannot compress log risk. But it fails the author's Sept 21 measurement principle on three counts. First, it is a wage of an individual head, not household earnings; with inelastic labor supply and no spouse in the model, wage risk is not the risk households face. Second, it carries no transitory component at all, because externals.py drops Floden–Linde's $\sigma_\xi^2$ as measurement error; the frozen baseline therefore has zero short-run earnings risk, the margin most relevant to down-payment saving. Third, its stationary variance (0.17 after HSV) is far below any household-earnings estimate, so the permanent component was later needed to recover wealth dispersion. The July choice should be retired, not because a run failed, but because its object does not match.

## Importing a published process: adjustments and double counting

- **Labor supply:** the model has none. Any wage-based process (FL, Sommer 2016, SS 2018) must be converted to earnings; no primary source for that conversion exists in the repo. Household earnings estimates avoid it.
- **Taxes and transfers:** model income is gross earnings after a flat payroll tax (0.179, `tmp/paper_baseline_sep14/PAPER_BASELINE.md:15`) with PAYGO pensions. Applying HSV $(1-\tau)$ to innovation s.d.s does not double count the payroll tax, because a proportional tax leaves log risk unchanged. Double counting would occur only if BGM's disposable-income variances were used, since they already net out income taxes and transfers, or if HSV's level factor $\lambda$ were also applied. Do not stack BGM variances with HSV.
- **Household composition:** BGM equivalizes by household size, so income mechanically falls at each birth. The model prices children through spending and space, so equivalized income would double count child costs. Use unequivalized earnings.
- **Age profile:** neither FL nor HSV supplies one. Estimate it on the same PSID sample, including ages 18–24 currently outside the 25–60 window.
- **Measurement error:** the candidate's transitory variance 0.344 is unadjusted; BGM scale the transitory s.d. by 0.45 citing Krueger–Perri. Use a sourced correction, and do not tune it to the housing fit.

## Four-year mapping

- **Endpoint mapping** ($\rho_4=\rho_a^4$, $\sigma_{\eta,4}^2=\sigma_{\eta,a}^2\sum_{k=0}^{3}\rho_a^{2k}$, `local_panel.py:1063-1070`) is exact for the latent state observed every fourth year. It needs no citation; it is algebra. It is not a statement about period income.
- **Period income is a four-year sum** in the budget (`solver.py:329-332` divides by period years). The variance of a sum along the path is smaller than the endpoint's; the September 20 receipt shows the 15-state chain's level variance 0.962 versus 1.205 (endpoint) and 1.155 (exact four-year average) (`CALIBRATION_STATUS.md:269-276`).
- **Supported precedent:** De Nardi 2004 forms five-year PSID cells and estimates at that frequency (verified by the Sept 21 ledger entry, `ACTIVE_DECISION_LEDGER.md:188-199`). Doepke–Kindermann 2019 fix wages and support nothing about conversion. I found no primary source for the candidate's iid-averaging formula $\log[1+(e^{V}-1)/4]$ or for the moment-matched proxy; both are approximations.
- **Recommendation:** estimate the period process directly on four-year average earnings where the annual PSID era permits (1984–1996 gives complete four-year cells), and use simulation of the annual process aggregated to four-year sums as the implementation, accepting it only if it reproduces the direct period covariances at lags 0–2.

## Minimal next steps before calibration

**Author decisions (three):**
1. Confirm the income concept: household gross labor earnings with HSV compression of risk, flat payroll tax on levels.
2. Confirm no permanent type in the baseline and the fixed-effect process as the named robustness process, accepting the disclosed long-lag misfit.
3. Confirm four-year period-average earnings as the measurement target for the period process.

**Numerical validation (bounded, no calibration):**
1. Pipeline check that addresses the author's doubt: the empirical lag-1 autocovariance (0.658, n=40,297, annual era only) is below lag-2 (0.673), which a stationary AR(1)+iid cannot produce; re-estimate separately for 1984–1996 and 1999–2019 and confirm the persistence is stable. Cross-check: the no-FE persistence near 0.97 appears in both in-house pipelines (July 16 Block 1: 0.9749; July 27: 0.9703), which is mild evidence the pipeline is not aberrant.
2. Sourced measurement-error adjustment, applied before aggregation.
3. Same-sample age profile including 18–24.
4. Direct four-year period moments from the annual era versus simulated aggregation; then a grid-resolution check at fixed entrant composition (E5/E6).

## Status of claims

**Verified:** all local paths above; FL Table IV; SS 2018 Table 1; Sommer 2016 equations and values; BGM pp.8, 17–18; both PSID parameter and autocovariance tables. **Inference:** why own estimates entered the candidate; the E6b adoption gap. **Unresolved:** HSV $\tau=0.181$ and BGM Table 2 values (not re-verified); a primary source for any annual-to-period closed form beyond the AR(1) identity; the correct measurement-error variance for PSID annual household earnings.
