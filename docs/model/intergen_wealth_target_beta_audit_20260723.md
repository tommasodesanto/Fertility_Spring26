# Matched wealth target and conditional-beta audit

Date: 2026-07-23

Scope: the timing-repaired **current one-shot model**, not the E-series. The
conditional profile replaces only the aggregate wealth row. The other thirteen
moments, weights, free nuisance parameters, and external restriction
\(\theta_n=0\) are unchanged. It is a diagnostic profile, not a promoted
calibration.

## 1. Target definition

The hard target is now

\[
 \frac{\sum_{a=18}^{85} w_i\,\text{NETWORTHR}_i}
      {\sum_{a=18}^{65} w_i\,\text{EARNINDRRC}_i}
 = 6.873077.
\]

It uses one PSID reference-person family record per wave, waves 2005--2019,
RP/spouse combined gross labor earnings (`EARNINDRRC`), and the longitudinal
weight `IW`. The sample contains 49,550 family-years and 10,432 persons. A
999-draw bootstrap clustered by reference-person ID gives standard error
0.398836 and percentile interval [6.189879, 7.692717].

Gross labor earnings are preferable to total family income (`INCFAMR`) here:
they map directly to the model's labor-income process, whereas `INCFAMR`
includes pensions, transfers, asset income, and other components. Holding the
sample fixed, the `INCFAMR` denominator would give 5.569710; that is recorded
only as a definitional sensitivity.

The model statistic now matches the same object:

\[
  \frac{\sum_{\text{all living households}}(b_t+pH_t)}
       {\sum_{\text{working households}} y_t^{\mathrm{gross}}}.
\]

Cross-sectional wealth uses the coherent beginning-of-period balance sheet.
Mortgage debt is carried in \(b_t\), so \(b_t+pH_t\) is net worth. Working
states 18 through 62 represent empirical ages 18--65; the terminal state 82
represents ages 82--85.

The old borrowed 6.9 after-tax target is retired. The fact that the newly
measured target also rounds to 6.9 is coincidence, not validation of the old
denominator. At the previous timing-repaired winner, the historical
after-payroll ratio is 6.148, while the matched gross/gross ratio is only
5.048.

## 2. Lifecycle robustness profile

The same waves and definitions yield:

| Reference-person age | PSID ratio | Bootstrap s.e. |
|---|---:|---:|
| 26--35 | 1.248257 | 0.075186 |
| 36--45 | 2.670562 | 0.218991 |
| 46--55 | 4.453772 | 0.391760 |
| 56--65 | 9.643537 | 0.908924 |

These four ratios are **robustness checks only**, in both the one-shot and
sequential code paths. They are not calibration targets. Because the model
uses four-year age states, boundary states are prorated by their overlap with
each empirical bin; for example, the 34--37 state contributes one half to
26--35 and one half to 36--45.

## 3. Literature comparison

[Borella, De Nardi, Yang, and Torres Chain
(2026)](https://users.nber.org/~denardim/research/NBERwp33874.pdf) estimate an
annual discount factor of 0.9958 (s.e. 0.00036), not 0.9900. Their estimate is
supported by 334 moments, including average and median wealth profiles by age
and marital status from ages 28 to 84. Their Appendix H varies each structural
parameter and shows that beta shifts wealth broadly over the lifecycle,
whereas bequest motives have more age- and household-type-specific effects.
This establishes that a high annual beta can be defensible, but it is not a
direct validation of a value near 0.999 in the present, much thinner target
system.

## 4. Conditional annual-beta profile

Ten bounded Torch searches reoptimized the other thirteen parameters at five
fixed annual betas. Every chain completed 360 evaluations, passed the strict
equilibrium gate, and produced two exact strict repeats. The lower loss from
the two chains in each cell is:

| Annual beta | Strict loss | Wealth / gross earnings | Living-old p90/p50 | TFR |
|---:|---:|---:|---:|---:|
| 0.9800 | 0.494155 | 3.179018 | 2.037876 | 2.031583 |
| 0.9900 | 0.380447 | 4.117083 | 1.899402 | 2.042949 |
| 0.9950 | 0.316187 | 4.640850 | 1.932108 | 2.078991 |
| 0.9990 | **0.295647** | **5.163442** | 1.853158 | 2.081957 |
| 0.9995 | 0.296005 | 5.060266 | 2.026796 | 2.079868 |

The objective has a high-beta plateau rather than a sharply identified point:
0.999 beats 0.9995 by only 0.000357. The nearest cell to Borella et al.'s
estimate, 0.995, costs 0.020540 relative to the profile minimum; 0.99 costs
0.084800.

Even the best cell reaches only 5.163 against the 6.873 wealth target, a 25
percent shortfall. Its lifecycle ratios are 0.910, 2.421, 4.545, and 6.483.
Thus the model approximately matches ages 46--55 but under-accumulates
especially at ages 56--65, where the data ratio is 9.644. These age moments
remain untargeted.

At the best cell, the living-old wealth-dispersion row contributes 0.213960
to the loss and the aggregate wealth row contributes 0.061873. Together they
account for 93 percent of the objective. The remaining rows fit much more
closely. The profile therefore exposes a wealth-structure problem, not merely
an unconventional discount factor.

Complete per-cell target fits, parameter tables, and strict-chain checks are
in
`output/model/intergen_new_moment_beta_recent_gross_20260723/report/`.

## 5. Decision

Do not cap beta at 0.98 or 0.99: those restrictions materially worsen the
matched objective and do not solve the wealth structure. Also do not present
0.999 as a precise estimate. The defensible statement is that the current
one-shot model selects a very patient region and still cannot jointly generate
the level and late-life dispersion of wealth.

The next diagnosis should use the lifecycle profile to isolate the missing
saving mechanism. Candidates include richer earnings risk and persistent
heterogeneity, retirement and pension incentives, and a better disciplined
intergenerational-saving channel. This diagnosis belongs in the reconciled
model rather than in another reweighting of the current fourteen rows.
