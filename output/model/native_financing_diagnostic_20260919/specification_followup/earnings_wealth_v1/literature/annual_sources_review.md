# Annual earnings-risk sources: bounded literature review

Date: 2026-09-21
Scope: annual US earnings processes relevant to a life-cycle/OLG model; no code, model run, or parameter adoption follows from this note.

## Executive assessment

The strongest plain-vanilla precedent for this project is Sommer (2016): a one-year household model with a persistent AR(1) earnings component and a separate iid transitory component, without an independent permanent earnings type. Its process is transparent and operational, but its ρ, σ values are imposed from the literature rather than estimated for the model's household sample.

Hubbard, Skinner, and Zeldes (1995) provide a closely related saving/insurance implementation. They use education groups as initial lifetime-income heterogeneity and estimate an AR(1) residual plus a transitory component that combines iid variation and measurement error. This is useful for model architecture, but it does not identify real transitory risk separately from measurement error.

Meghir and Pistaferri (2004) is the key warning against treating the annual process as homoskedastic AR(1)+iid. They use PSID labor income, separate permanent and transitory components, allow an MA(q) transitory process, explicitly model measurement error, and find education-, time-, age-, and person-specific variance dynamics. Their preferred permanent component is a unit-root martingale, so it is not a direct justification for a stationary AR(1) benchmark.

It does not establish that Sommer's parameter triplet (0.95, 0.21, 0.17) should be adopted unchanged, nor that a permanent fixed type should be added without a separate identification decision. A fresh annual estimate is needed if the project wants its household earnings concept, sample, and risk decomposition to be empirically authoritative.

## Sommer (2016), Journal of Monetary Economics

Source: Kamila Sommer, “Fertility Choice in a Life Cycle Model with Idiosyncratic Uninsurable Earnings Risk,” *Journal of Monetary Economics* 83 (2016), pp. 27–38, DOI 10.1016/j.jmoneco.2016.08.002. [Author PDF](https://www.kamilasommer.net/Fertility.pdf).

The model period is one year; households start at age 18, work through age 65, and live to age 80 (p. 29). The wage is household wage income in a unitary husband–wife model:

\[
\log w_t=\log w_0+h(t)+\epsilon_t+\nu_t,
\qquad
\epsilon_t=\rho\epsilon_{t-1}+\psi_t.
\]

Here \(h(t)\) is the deterministic age profile, \(\nu_t\sim N(0,\sigma_\nu^2)\) is a transitory shock received each period, and \(\psi_t\sim iid(0,\sigma_\epsilon^2)\) is the persistent innovation; the initial persistent state is set to \(\epsilon_1=0\) (p. 29, eq. in Section 3.1). There is no separate permanent fixed earnings type in the stated process.

Table 2 (p. 31) sets \(\rho=0.95\), persistent innovation standard deviation \(\sigma_\epsilon=0.21\), and transitory standard deviation \(\sigma_\nu=0.17\). These are standard deviations, not variances. The paper says these are selected from published ranges: \(\rho\) and \(\sigma_\nu\) are set at the middle of available estimates, while \(\sigma_\epsilon=0.21\) is chosen near the upper end to represent the 1980s–1990s cohort (pp. 31–32).

Observed wage calibration is household-based but the age profile is constructed from 2004 CPS husband–wife families: family labor income is the sum of both spouses' yearly earnings, divided by the couple's total hours (p. 32, footnote 14). The retirement transfer is 40% of household earnings in the last working period (p. 32, footnote 15). The paper's fertility targets come from NLSY79; the risk parameters themselves are external literature values rather than a new PSID household-process estimate. The model starts childless with zero assets and limited credit (p. 29).

Useful exact statement (p. 29): “The persistent shock, \(\epsilon_t\), also received each period, follows a first-order autoregressive process.”

Interpretation: Sommer is a clean annual benchmark for the project's intended no-fixed-type architecture. It does not validate the numerical values for a gross-labor-earnings process in this project, and it does not separately identify measurement error.

## Hubbard, Skinner, and Zeldes (1995), Journal of Political Economy

Source: R. Glenn Hubbard, Jonathan Skinner, and Stephen P. Zeldes, “Precautionary Saving and Social Insurance,” *Journal of Political Economy* 103(2) (1995), pp. 360–399, DOI 10.1086/261987. [Author-hosted PDF](https://business.columbia.edu/sites/default/files-efs/imce-uploads/szeldes/pdfs/precautionary.pdf); [journal page](https://www.journals.uchicago.edu/doi/10.1086/261987).

The empirical earnings object is log earnings after controlling for a cubic age profile and year effects; the paper's earnings measure includes unemployment insurance and subtracts taxes. It separates households into three education categories. The published model uses common preferences across lifetime-income groups, while education proxies lifetime-income heterogeneity.

The process reported on p. 379 is:

\[
y_{it}=Z_{it}\beta+u_{it}+v_{it},
\qquad
u_{it}=\rho_eu_{i,t-1}+e_{it},
\]

where \(e_{it}\) is a white-noise innovation and \(v_{it}\) combines iid transitory earnings variation with measurement error. The model therefore has an AR(1) persistent residual, an iid/transitory residual, and education-group initial heterogeneity. It does not justify interpreting the whole persistent component as a permanent type.

The paper's practical limitation is identification: the \(v\) object is not pure real transitory risk because measurement error is folded into it. Its net-after-tax/including-UI earnings concept also differs from the project's gross labor-earnings object. It is a useful benchmark for precautionary saving and social insurance, but the parameters should not be transplanted without reconciling the income definition.

Useful exact statement (p. 379): “The variable \(v_{it}\) is a combination of independently and identically distributed transitory variation in earnings and measurement error.”

## Meghir and Pistaferri (2004), Econometrica

Source: Costas Meghir and Luigi Pistaferri, “Income Variance Dynamics and Heterogeneity,” *Econometrica* 72(1) (2004), pp. 1–32, DOI 10.1111/j.1468-0262.2004.00476.x. [Stanford PDF](https://web.stanford.edu/~pista/meghir.pdf).

Data are PSID family and individual files, waves covering 1967–1992 (paper pp. 2–3). The earnings variable is the labor portion of money income from all sources—wages, bonuses, overtime, commissions, professional practice, and labor portions of farm/business income (p. 3). It is an individual/head earnings concept, not the project's household gross labor-earnings aggregate.

The paper's conditional-mean controls are education-specific: year dummies, a quadratic in age, race, region, and SMSA residence (p. 6). Its baseline decomposition is more general than AR(1)+iid:

\[
u_{it}=r_{it}+e_{it}+p_{it},
\qquad
p_{it}=p_{i,t-1}+\zeta_{it},
\]

where \(r_{it}\) is classical iid measurement error, \(e_{it}\) is transitory, and \(p_{it}\) is a permanent martingale component (p. 5). The transitory component is allowed to follow an MA(q) process rather than iid (p. 5). Thus their permanent component is a unit root, not a stationary AR(1).

They explicitly caution that the variance of real transitory risk cannot be separated from measurement-error variance and the MA coefficients using the basic autocovariances alone (p. 11). With an MA(1), they report bounds; using an external measurement-error share of 25%, they obtain transitory variances of 0.0548, 0.0267, and 0.0049 for high-school dropouts, high-school graduates, and college graduates, respectively (pp. 12–13). These are variances, not standard deviations. The permanent innovation variance is time-varying and education-specific; the paper reports no single plain-vanilla \(\rho,\sigma_\eta,\sigma_\varepsilon\) triplet suitable for direct transplantation.

Their variance dynamics are ARCH(1) with education-specific year effects, life-cycle effects, and person-specific fixed effects in the variance (p. 13). This is evidence against assuming iid risk is constant across age, time, education, or people.

Useful exact statement (p. 2): “We find strong evidence of state dependence in the variance of both permanent and transitory components.”

## Storesletten, Telmer, and Yaron (2004)

Sources: “Cyclical Dynamics in Idiosyncratic Labor Market Risk,” *Journal of Political Economy* 112(3) (2004), pp. 695–717, DOI 10.1086/383105, [journal page](https://doi.org/10.1086/383105); and “Consumption and Risk Sharing over the Life Cycle,” *Journal of Monetary Economics* 51(3) (2004), pp. 609–633, [journal page](https://doi.org/10.1016/j.jmoneco.2003.06.005).

The verified published summary for the JPE paper reports annual idiosyncratic-risk persistence of approximately \(\rho=0.95\), with conditional shock standard deviation increasing from 0.12 at a macro peak to 0.21 at a trough. The paper estimates these objects from household-level PSID labor earnings and conditions on age/education and macroeconomic histories. The JME paper's published summary emphasizes a highly persistent life-cycle component with autocorrelation between 0.98 and one and separates risk realized before labor-market entry from risk realized during working life.

These papers are important annual empirical benchmarks, but they are not clean sources for the project's exact income object: the JPE measure is household labor earnings residual risk, while the JME paper's decomposition includes pre-labor-market heterogeneity and life-cycle risk. The verified public summaries do not establish a single iid-transitory standard deviation, a measurement-error correction, or an actual code-level initialization rule for the project's process. They should therefore discipline plausible annual persistence and cyclical risk ranges, not determine a full production process by themselves.

## Implications for the project's annual specification

1. A defensible plain-vanilla benchmark is Sommer's architecture: deterministic age profile plus stationary persistent AR(1) innovation and separate iid transitory shock, with no independent permanent type. Record all three risk parameters as standard deviations/coefficients with their units.
2. Do not label the persistent state a “permanent type” unless an author decision adds fixed heterogeneity. HSZ's education groups are observed initial heterogeneity; Meghir–Pistaferri's person-specific variance effects are not the same object as a fixed earnings level.
3. Do not treat HSZ's transitory residual as pure real risk: it includes measurement error. Meghir–Pistaferri show that even a richer panel cannot identify real transitory variance without an external measurement-error restriction or multiple indicators.
4. The project's gross labor-earnings concept differs from Sommer's household wage calibration and HSZ's after-tax, UI-adjusted earnings. A fresh data estimate is required before calling any literature parameter the project's empirical estimate.
5. The literature review supports a plan of annual data construction, explicit sample and household-unit definitions, education/age controls, separate real transitory and measurement-error treatment where possible, then a transparent annual calibration. It does not support silently importing \((0.95,0.21,0.17)\) as a settled target or adding a fixed permanent component merely to fit dispersion.
