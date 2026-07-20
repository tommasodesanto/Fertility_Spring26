# Re-audit and study guide for the July 17 calibration note

Date: 2026-07-20

Circulated artifact: `latex/quantification_rewrite_review.pdf`

Circulation date: 2026-07-17
Calibration audited: M5 (`intergen_income_disciplined_recalibration_20260716`)

## Bottom line

The circulated note is numerically faithful to the canonical M5 solve. All 14
reported parameter estimates and all 15 model moments reproduce the strict,
exactly repeated Torch winner, whose weighted loss is $9.044422$ and whose
equilibrium residual is $1.48 \times 10^{-5}$. The note also disclosed the two
central substantive limitations that were known at circulation: the low income
innovation variance and the excessive old-age ownership rate.

The re-audit nevertheless changes how the note should be defended. Four issues
need to be stated more sharply in a meeting:

1. **The calibration is only weakly locally identified.** The SMM-weighted
   Jacobian has effective rank $9/14$ at relative tolerance $10^{-2}$ and
   $12/14$ at $10^{-3}$, with condition number $2925.9$. The weakest
   direction jointly loads on $(\bar h_0,\theta_1,\theta_0,\kappa_T)$. The note's
   moment-to-parameter language should therefore be read as economic intuition,
   not separate identification.
2. **Two headline targets are source-unconfirmed.** The repository's July 6
   reproduction audit found no builder or saved source artifact for completed
   fertility (1.918) or completed childlessness (0.188). The note attributes
   them to the CPS, but the survey year, age window, children-ever-born variable,
   and weights are not pinned in the repository. They are not known to be wrong;
   they are not presently reproducible.
3. **The scalar loss is not a conventional efficient-SMM statistic.** The
   objective mixes inverse-variance weights for some late-life moments with
   judgmental/legacy weights for most other moments. It is useful for ranking
   candidates under the same target system, but it has no clean chi-squared
   interpretation and should not be compared across specifications with changed
   moments or weights.
4. **The architecture is more stylized than the note makes explicit.** Family
   size is chosen in one shot, all dependent children share one geometric
   maturity clock, and the housing-finance block uses a collateral/down-payment
   shortcut rather than a mortgage contract. There is no mortgage coupon,
   amortization, refinancing, payment-to-income limit, default, or continuous
   housing-equity state. Bequest utility values $b+PH$ at death, with mortgage
   debt already inside $b$, but without a selling-cost deduction.

The defensible description is therefore: **M5 is a reproducible provisional
calibration that matches several housing-response and wealth-level objects, but
it is not yet a fully identified or fully sourced structural estimate. Its main
quantitative use is conditional analysis of family housing needs, tenure, and
the down-payment margin.**

## What Corina actually received

The five-page PDF was compiled at 17:32 on July 17. It differs in two places
from the `.tex` source that was edited later that evening:

- The circulated PDF says that the median total-wealth-to-income ratio for
  households with two or more children exceeds the one-child median by
  (0.101). The later source adds the essential qualification that the
  difference is statistically insignificant. The person-bootstrap standard
  error is (0.563), so the empirical difference is only (0.18) standard
  errors from zero.
- The later source adds a measurement qualification to the estate-tail
  discussion: PSID nonhousing wealth contains business wealth, retirement
  accounts, and other assets absent from the model.

The PDF, rather than the current `.tex`, is the evidentiary record of what was
circulated.

## The model in plain language

### Lifecycle

- A period is four years; model ages are $18,22,\ldots,82$, with the final
  decision period interpreted as living through age 85.
- Retirement begins at age 66.
- Survival is certain before retirement. Thereafter, four-year survival
  probabilities are constructed from the combined male and female survivor
  counts in the 2023 SSA period life table.
- Dependent children mature with probability $2/9$ each four-year period.
  This is a geometric duration with mean $4/(2/9)=18$ years, not a
  deterministic 18-year child-age profile.

### Fertility

- The M5 model is a **one-shot completed-family-size model**. A childless
  household chooses its completed family-size category rather than moving
  through separately modelled first- and second-birth hazards.
- The reported object called `tfr` in code is
  $2 \times$ mean completed household fertility after the fertility window.
  The note appropriately labels it completed fertility; it is not a period TFR.
- Childlessness is the post-fertility parity-zero stock, not the absence of a
  child in the current household.
- Children affect utility directly through $\psi_{child}$, consumption
  requirements through $\bar c_n$, and housing requirements through a first-
  child jump $\bar h_{jump}$ plus $\bar h_n$.

### Housing and tenure

- Renters choose continuous housing services up to six rooms.
- Owners choose a discrete rung from $\{2,4,6,8,10\}$ rooms.
- Owners receive a reduced-form service premium $\chi_O$.
- Housing services clear against an upward-sloping supply curve with elasticity
  $1.75$ and internally estimated scale $\bar H$.
- A purchase requires cash of at least $(1-\phi)PH=0.2PH$. In the code,
  $\phi=0.8$ is the financed share, not the down-payment share.
- The model charges housing user cost each period and imposes a collateral
  threshold, but it does not track a standard amortizing mortgage contract.

### Saving and bequests

- Liquid wealth $b$ may be negative and includes secured debt.
- At death, the bequest object is $b+PH$. Because mortgage debt is in $b$,
  this is net of debt but gross of selling costs.
- The bequest utility is a normalized, child-blind De Nardi-style warm glow in
  M5. $\theta_0$ controls strength, $\theta_1$ shifts the wealth dependence,
  and $\theta_n=0$ is imposed externally.
- M5 matches the median estate level and the share of older households with
  nonhousing wealth above annual income, but it does not match the upper tail or
  the late-life tenure path.

## Claim-by-claim audit

| Claim in circulated note | Verdict | What to remember |
|---|---|---|
| Four-year periods; entry 18; retirement 66; maximum age 85 | Verified | Implemented as states 18 through 82, with the last four-year period extending through 85. |
| Post-retirement survival from 2023 SSA life table | Verified | Code uses combined male/female survivor counts and four-year survival ratios. |
| Expected dependent-child duration 18 years; maturity probability (2/9) | Verified with qualification | It is a geometric shared clock, not age-resolved children. |
| (sigma=2) | Verified | Maintained restriction. |
| (i=2%), (delta=1.1%), sale cost (6%), payroll tax (17.9%), (phi=0.8) | Verified | Greaney et al. report these benchmark values/sources; their property tax is (0.71%), whereas this note separately maintains (1%). |
| No payment-to-income limit | Verified | The optional screen exists but is off in M5. |
| Annual income process (rho_z,sigma_z)=(0.960,0.0645), five-state Rouwenhorst | Verified | Persistence is defensible; innovation risk is deliberately too low and is not a PSID estimate. |
| Income variance is provisional because higher risk creates empty budget sets | Verified and strengthened after circulation | The July 18 frontier shows failure already at (sigma_z=0.09) with the M5 subsistence bundle. |
| Owner sizes ({2,4,6,8,10}); renter cap six rooms | Verified | Implemented in the M5 override. The two-room owner rung is effectively dormant near the childless housing floor. |
| Supply elasticity (1.75) from Saiz | Verified | Saiz reports a population-weighted metropolitan elasticity of 1.75. Applying one national elasticity to this one-market model remains a maintained simplification. |
| Entrant wealth is the PSID childless-renter 18--24 distribution | Verified | Weighted quintile-bin means are mapped into model wealth units at entry. |
| The 2+-minus-1-child estate difference is (0.101) | Numerically verified but incomplete | Bootstrap SE (0.563); not statistically different from zero. The later source correctly added this. |
| 14 parameters are estimated with 15 moments | Count verified | Overidentified by one only mechanically; effective local rank is materially lower. |
| Specific moments are “particularly informative” for named parameters | Directionally plausible, too strong as identification language | Parameters are jointly disciplined; do not claim one-to-one identification. |
| All 15 displayed model moments | Verified | They match the canonical strict M5 record. |
| Completed fertility and childlessness come from the CPS | Source-unconfirmed | Exact values exist as constants, but the CPS extract/builder and sample definition are missing. |
| The model “closely matches” renter mean rooms and the owner-renter gap | Reasonable in raw units, weaker statistically | Gaps are about (2.14) and (1.34) available bootstrap SEs, respectively. |
| Largest miss is old-age ownership | Verified | It contributes (5.759), or (63.7%), of the total weighted loss and is (15.1) available bootstrap SEs above the target. |
| Missing health/LTC risk and age-dependent burdens can explain excessive retention | Plausible diagnosis, not established by M5 | These mechanisms are absent and are natural candidates; M5 does not identify which one is responsible. |
| Estate p90/p50 is (1.75) vs (3.45) | Verified | The target is untargeted; M5 lacks the wealth tail. |
| Quantitative exercises should be conditional on limitations | Correct | This is the right interpretation of the circulated results. |

## Full M5 parameter table

No reported estimate is flagged as near its search bound by the collector. That
does not imply precise identification.

| Parameter | Estimate | Search bound / restriction | Status |
|---|---:|---:|---|
| $\beta_{annual}$ | 0.991232 | [0.8000, 0.9995] | Estimated |
| $\alpha_c$ | 0.591226 | [0.02, 0.98] | Estimated |
| $\bar c_0$ | 1.259594 | [0, 2] | Estimated; equals 0.3149 of mean annual income after period conversion |
| $\bar c_n$ | 0.392945 | [0, 3] | Estimated |
| $\bar h_0$ | 0.391542 | [0.05, 5.8] | Estimated |
| $\bar h_{jump}$ | 1.602118 | [0, 8] | Estimated |
| $\bar h_n$ | 0.164365 | [0, 5] | Estimated |
| $\psi_{child}$ | 0.197914 | [-3, 3] | Estimated |
| $\kappa_F$ | 2.043386 | [0.02, 50] | Estimated |
| $\chi_O$ | 1.113266 | [0.1, 5] | Estimated |
| $\bar H$ / `H0` | 8.645578 | [0.2, 80] | Estimated supply scale |
| $\theta_0$ | 0.311765 | [0, 8] | Estimated |
| $\theta_1$ | 0.397282 | [0.02, 16] | Estimated |
| $\kappa_T$ | 0.010017 | [0, 0.12] | Estimated; jointly weak with housing/bequest block |
| $\theta_n$ | 0 | External restriction | Fixed; empirical 2+-minus-1 gap (0.101), SE (0.563) |

## Full M5 target-fit table

The loss is $\sum_m w_m(\hat m_m-m_m^{data})^2$. The weights are a mixture of
inverse-variance and judgmental weights, so loss contributions are useful for
diagnosis but not for a formal chi-squared test.

| Moment | Data | Model | Gap | Weight | Loss contribution |
|---|---:|---:|---:|---:|---:|
| Completed-fertility-equivalent | 1.918000 | 2.034398 | +0.116398 | 20.000 | 0.270972 |
| Completed childlessness | 0.188000 | 0.188735 | +0.000735 | 20.000 | 0.000011 |
| Ownership, ages 30--55 | 0.575472 | 0.658378 | +0.082906 | 100.000 | 0.687337 |
| New-parent minus nonparent ownership, ages 30--55 | 0.167662 | 0.219444 | +0.051782 | 45.000 | 0.120662 |
| First-child housing response, rooms | 0.664435 | 0.680817 | +0.016383 | 14.000 | 0.003757 |
| Childless-renter mean rooms, ages 30--55 | 3.805288 | 3.963278 | +0.157990 | 6.000 | 0.149765 |
| Childless-owner share with at least six rooms, ages 30--55 | 0.596131 | 0.760146 | +0.164015 | 25.000 | 0.672524 |
| Young childless-renter liquid wealth / annual gross income | 0.179226 | 0.327928 | +0.148702 | 12.000 | 0.265348 |
| Childless owner-renter mean rooms gap, ages 30--55 | 2.418762 | 2.322768 | -0.095993 | 12.000 | 0.110577 |
| Ownership, ages 25--34 | 0.341166 | 0.274442 | -0.066725 | 80.000 | 0.356173 |
| Rooms gap, 3+ versus 1--2 children, ages 30--55 | 0.367700 | 0.607477 | +0.239778 | 8.000 | 0.459947 |
| Median total estate / annual income, ages 76--84 | 6.501316 | 6.430926 | -0.070389 | 18.586 | 0.092087 |
| Share with nonhousing wealth at least annual income, ages 65--75 | 0.608333 | 0.610017 | +0.001684 | 9435.187 | 0.026744 |
| Ownership, ages 65--75 | 0.764261 | 0.953987 | +0.189727 | 160.000 | 5.759384 |
| Mean occupied rooms, ages 18--85 | 5.779970 | 5.672629 | -0.107341 | 6.000 | 0.069133 |
| **Total** |  |  |  |  | **9.044422** |

### How to read the fit

- The model genuinely does well on completed childlessness, the first-child
  housing response, the two targeted late-life wealth levels, and aggregate
  rooms, conditional on accepting their target definitions and weights.
- It misses the tenure block: ownership is too high at prime and old ages, too
  low for ages 25--34, and too strongly associated with family formation.
- It exaggerates large-home selection among childless owners and the room gap
  between large and small families.
- The old-age ownership row alone supplies (63.7%) of total loss. The next
  largest contributions are prime-age ownership (7.6%), large-home share
  among childless owners (7.4%), and the 3+-versus-1--2 rooms gap
  (5.1%).
- Available bootstrap SEs make the same problem look more severe: old-age
  ownership is (15.1) SEs high; the large-home share is (9.4) SEs high;
  the family ownership gap is (8.0) SEs high; prime-age ownership is (4.1)
  SEs high; and young ownership is (3.7) SEs low. These are descriptive
  ratios, not a joint test, because the active weight matrix is not the full
  covariance matrix.

## Identification: the precise answer

If asked “are the 14 parameters identified by the 15 moments?”, the best answer
is:

> The system is overidentified by one in the mechanical moment-count sense, and
> the exact winner is numerically reproducible. But the local sensitivity matrix
> is ill-conditioned: its effective SMM-weighted rank is 9 at a 1% relative
> threshold and 12 at a 0.1% threshold. So the data jointly discipline the broad
> blocks, but several parameter combinations are weakly separated. I would not
> make claims of precise or one-to-one identification from this calibration.

The two source-unconfirmed fertility targets matter for that answer. If they
were removed without replacement, only 13 empirically verified hard moments
would remain for 14 free parameters, making the system underidentified even by
count. They must be sourced and retained, replaced with equally informative
fertility moments, or paired with an additional external restriction.

## Limitations omitted or understated in the circulated note

1. **One-shot fertility.** The model does not contain separately timed first-
   and second-birth decisions. The rooms gap and first-child housing response
   discipline housing demand around family formation; they do not identify a
   sequential second-birth hazard.
2. **Geometric shared child clock.** All children age out together under a
   memoryless maturity probability. Child spacing and age-specific child costs
   are absent.
3. **Mortgage shortcut.** The model has a down-payment/collateral constraint but
   no amortizing contract, payment schedule, refinancing, default, or lock-in.
4. **Gross terminal house valuation.** Terminal bequest utility includes (PH)
   and does not subtract the sale cost. This convention was accepted
   provisionally, not empirically validated.
5. **Late-life liquid-wealth path.** At the M5 solution, mean liquid wealth turns
   negative by the age-74 state and falls further at ages 78 and 82 while
   ownership remains above 90%. The model matches selected old-age wealth
   targets largely through housing and collateralized debt, not through a
   realistic portfolio path.
6. **Mixed weights and missing uncertainty.** The circulated table omits target
   SEs, weights, loss contributions, parameter bounds, and parameter uncertainty.
7. **National one-market interpretation.** A single housing market and one
   supply elasticity abstract from the spatial variation central to the broader
   project.

## What changed after circulation

These results should be presented as subsequent evidence, not silently folded
into what Corina saw.

### Income-feasibility frontier

The July 18 fixed-parameter frontier shows that M5 is closer to the feasibility
boundary than the footnote suggested. With the M5 Stone--Geary bundle and no
safety net, annual innovation risk of (0.09) already produces positive-mass
infeasible states at age 22. A measured transfer floor does not help at the M5
bundle because the guarantee is below the required bundle. Lowering annual
(bar c_0) from (0.315) to (0.10) is what moves most of the frontier; the
floor adds the final step from (sigma_z=0.18) to (0.20).

### Transfer-floor probe

At fixed M5 parameters, the combination of a low bundle, a measured debt-blind
means-tested floor, and (sigma_z=0.20) is feasible. It raises estate p90/p50
from (1.75) to (3.07), with receipt by (0.91%) of households and outlays
of (0.054%) of gross income. It simultaneously raises young liquid wealth to
(1.33) versus a (0.179) target and changes several fertility and tenure
moments sharply. This is evidence about the trade-off, not a recalibration or a
replacement for M5.

### Equivalence-scale plus sequential-fertility experiment (E1)

E1 removes the Stone--Geary intercepts, admits (sigma_z=0.20), and introduces
sequential first- and second-birth attempts. After recalibration on the unchanged
15 moments, its loss is (12.608), worse than M5's (9.044). It resolves the
old-age ownership level (0.710) versus (0.764) and produces estate
p90/p50 (3.67) versus (3.45), but misses fertility and the first-child
housing response: completed fertility is (2.474), childlessness (0.112),
and the first-child room response (1.091). The fitted direct child-utility
parameter is near zero. E1 is an informative experiment, not the new baseline.

## Study questions and short answers

### What is the paper's central quantitative mechanism?

Children raise minimum housing needs. Because owner housing is lumpy and a
purchase requires a down payment, high housing costs and limited liquid wealth
change tenure, housing consumption, saving, and the attractiveness of family
formation. The calibration tries to discipline that joint mechanism with
fertility, housing-response, tenure, and wealth moments.

### Which parameters are external and which are estimated?

Timing, survival, (sigma), finance/tax rates, the housing menus, the supply
elasticity, the income process, entrant wealth distribution, and
(theta_n=0) are external or maintained. Fourteen parameters governing
discounting, consumption/housing requirements, fertility utility and
dispersion, owner utility, supply scale, bequests, and tenure dispersion are
estimated jointly.

### Why is (bar c_0) so consequential?

It is a hard Stone--Geary minimum. At M5, (bar c_0=1.2596) per four-year
period, or (0.3149) of mean annual income before adding minimum housing. It
creates strong saving/fertility discipline but makes realistic low income
states infeasible without insurance or a redesigned preference block.

### Why does the model overpredict old-age ownership?

Owners face transaction costs, value the house in terminal bequests, and lack
health/LTC expense risk, age-dependent maintenance burdens, default, and a
realistic mortgage/downsizing margin. M5 establishes the symptom, not which
missing mechanism is causal.

### Why can the estate median fit while the estate tail fails?

The bequest and saving parameters can locate the center of the distribution,
but the low-risk income process and limited portfolio heterogeneity generate too
little dispersion. Matching a median does not identify the upper tail.

### What does (theta_n=0) mean?

Bequest utility does not directly scale with the number of children. Children
still affect estates indirectly through housing needs, tenure, and saving. The
PSID 2+-minus-1-child median estate gap is (0.101) with SE (0.563), which
supports treating a direct child-count shifter as weak, but it does not prove
the parameter is exactly zero.

### Can the M5 loss be compared with E1's loss?

Yes, cautiously, because E1 intentionally uses the unchanged 15 moments and
weights. M5 fits those targets better (9.044<12.608). The architectures
differ, so the comparison is a specification-fit comparison, not a likelihood
ratio test or a welfare ranking.

### What can be said safely about identification?

The broad blocks are jointly disciplined, but the parameter vector is weakly
separated locally. Do not say that one moment “identifies” one parameter or that
the numerical estimates are precise.

### What should not be claimed tomorrow?

- Do not call (1.918) a period TFR.
- Do not claim the fertility targets are reproducibly sourced until the CPS
  builder is pinned.
- Do not describe the model as containing a sequential second-birth hazard.
- Do not describe the housing-finance block as an amortizing mortgage model.
- Do not interpret the loss as a chi-squared statistic.
- Do not present E1 or the transfer-floor probe as the new benchmark.
- Do not claim separate identification of (bar h_0,theta_0,theta_1,kappa_T).

## Suggested meeting opening

> I re-audited the five-page calibration note we circulated. The numerical
> tables reproduce the exact M5 solution, and the two limitations we disclosed
> -- low income risk and excessive old-age ownership -- are real. The audit also
> showed that I should be more careful about three things: two fertility targets
> still need a pinned CPS source, the 14 parameters are only weakly separated by
> the local Jacobian, and the baseline is a one-shot fertility model with a
> collateral shortcut rather than a full mortgage contract. Since circulation,
> we have run specification experiments that clarify these tensions, but M5
> remains the benchmark. I would like to use today's discussion to decide which
> limitation is most important for the paper's next version.

## Productive questions for Corina

1. Is the one-shot completed-family-size architecture acceptable for the paper's
   current mechanism, provided timing claims are removed, or is sequential
   fertility essential now?
2. Should the next baseline prioritize realistic income risk and insurance, or
   first repair the old-age ownership/mortgage-exit margin?
3. Is the current empirical target set too broad relative to what the stylized
   one-market model can credibly identify?
4. Would she prefer externally fixing the weak bequest/tenure combination, or
   adding new moments/mechanisms before treating those parameters as estimated?
5. Which quantitative claims in the paper actually require the late-life
   portfolio block to be credible, and which survive if that block is explicitly
   quarantined?

## Evidence used in this re-audit

- Circulated PDF and its July 17 Git version.
- Canonical M5 target-fit, parameter, acceptance, lifecycle, and identification
  artifacts under `output/model/intergen_income_disciplined_recalibration_20260716/`.
- Target-object ledger and July 6 reproduction/standard-error audit under
  `code/data/moment_standard_errors/output/`.
- Active model implementation and M5 runner under
  `code/model/intergen_housing_fertility/` and
  `code/model/tools/run_intergen_bequest_exit_chain.py`.
- July 18 income-feasibility and transfer-floor artifacts.
- July 19 E1 collection report.
- Primary external sources: [2023 SSA period life table](https://www.ssa.gov/oact/STATS/table4c6.html);
  [Greaney, Parkhomenko, and Van Nieuwerburgh (2025)](https://www.nber.org/system/files/working_papers/w33512/w33512.pdf);
  [Saiz (2010)](https://academic.oup.com/qje/article-pdf/125/3/1253/5373851/125-3-1253.pdf);
  [Sommer and Sullivan (2018)](https://www.aeaweb.org/articles?id=10.1257%2Faer.20141751).
