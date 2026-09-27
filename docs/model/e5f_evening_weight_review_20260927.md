# Evening calibration: weighting and identification review

September 27, 2026. Bounded methodological review; no model evaluations or
changes to targets, code, or running contracts. This note proposes diagnostic
weights for the author-authorized six-hour Torch search. It does not claim
that any diagonal weights estimate an optimal covariance matrix.

## Target contract and parameter count

The instructed system has ten scored moments, ten searched coordinates
(including the newly free tenure taste scale), and one internally normalized
child-benefit level. Completed fertility 2.1 supplies that normalization, so
there are eleven fitted quantities and eleven restrictions by count. This is
just identified by count, not an empirical demonstration of identification.
The three other displayed moments—first births at age 30 or older, the rooms
gap between larger and smaller families, and older wealth dispersion—remain
validation outcomes. Keep all fourteen rows in every readout, with explicit
scored/normalization/validation labels. Validation contributions must be blank
or zero under the new objective, rather than carrying inherited scores.

The numerical targets and inherited weights below are copied from
`output/model/supervised_calibration_20260927/final_readout/target_fit.csv`.
Production must retain the full-precision authoritative target/weight contract;
this rounded display is not a replacement source.

| Scored moment | Target | Inherited weight | Block | Relative-error scale |
|---|---:|---:|---|---|
| Childlessness, ages 40–44 | 0.198 | 35,532.304 | Fertility | Exact target |
| Exactly one child, among mothers | 0.214 | 26,952.821 | Fertility | Exact target |
| Mean age at first birth | 25.976 | 139.828 | Fertility | Exact target |
| Children ever born at 25, capped at three | 0.810 | 100.000 | Fertility | Exact target |
| Aggregate wealth / annual earnings | 6.927 | 7.595 | Wealth | Exact target |
| Annual bequest flow / wealth | 0.007 | 5.165e6 | Wealth | **0.010** |
| Mean occupied rooms | 5.729 | 128.021 | Housing | Exact target |
| Ownership, ages 30–55 | 0.676 | 2,339.362 | Housing | Exact target |
| Housing response to first birth | 1.465 | 137.565 | Housing | Exact target |
| Recent-parent ownership gap | 0.128 | 27,055.823 | Housing | Exact target |

## Three separate weighting lanes

Let $g_i(\theta)=m_i(\theta)-t_i$ be the model-minus-target gap, evaluated after
the unchanged internal fertility normalization. Let $w_i^0$ denote the retained
inherited weight of a scored moment. All sums below contain only the ten
scored moments. No lane changes the model, bounds, empirical definitions or
normalization. Each lane chooses its own winner; compare all winners under
all three weight systems and always include the common primary score.

1. **Primary: inherited active weights.**
   $$L_P(\theta)=\sum_{i=1}^{10}w_i^0g_i(\theta)^2.$$
   Keep the weights of retained moments exactly unchanged. Removing the three
   validation contributions creates a new objective; its loss is not directly
   comparable with the old thirteen-scored-moment loss of 42.282. Re-score the
   starting point under this new contract before reporting improvement.

2. **Diagnostic: identity on safely scaled relative gaps.**
   Define $s_i=\max\{|t_i|,d_i\}$ and
   $$L_R(\theta)=\sum_{i=1}^{10}[g_i(\theta)/s_i]^2.$$
   Thus raw-gap weights are $1/s_i^2$. Use fixed floors $d_i=0.1$ for
   proportions/ownership gaps and children ever born, $d_i=1$ for age, rooms,
   room responses and wealth/earnings, and $d_i=0.01$ for annual bequests/wealth.
   At the current targets only the bequest floor binds: its weight is 10,000,
   rather than dividing by a tiny or possibly zero target. These floors are
   frozen before search; never divide by the current model value or current
   residual. “Identity” here means identity after the declared row scaling,
   not identity on heterogeneous raw units. Relative age errors are naturally
   inexpensive because age is around 26, an intentional diagnostic tradeoff
   that must not be described as statistically neutral.

3. **Diagnostic: equal block averages of inherited standardized errors.**
   Let fertility contain the four fertility rows, housing the four housing
   rows, and wealth the two wealth rows in the table. For block $b$, let
   $n_b$ be its row count. Set
   $$L_B(\theta)=\frac13\sum_{b\in\{F,H,W\}}\frac1{n_b}
       \sum_{i\in b}w_i^0g_i(\theta)^2.$$
   The implemented raw-gap weights are $w_i^0/12$ for fertility and housing
   and $w_i^0/6$ for wealth. This equalizes block *averaging coefficients*,
   not realized contributions or statistical information. It preserves the
   inherited relative importance of rows within each block. Do not normalize
   by current or starting residual magnitudes: near-zero fits would generate
   arbitrary extreme weights. The overall factor 1/3 does not change a lane's
   minimizer but must be fixed in the fingerprint and reported scores.

The primary lane remains authoritative for selection unless the author adopts
a different system. The two alternatives diagnose objective sensitivity, not
provide extra votes for a candidate. The six-hour controller must share the
same hard deadline (approximately 22:12 EDT if that is the recorded start plus
six hours), with time reserved for repetitions and export. Pin exact launch
and end epochs; do not extend a lane merely because another is unfinished.

## What the existing Jacobian does—and does not—identify

The daytime Jacobian is at the old nine-coordinate, fixed-tenure-scale anchor,
with the original thirteen scored moments. It has numerical rank nine and
condition number 262.387 in its documented dimensionless coordinates. It does
**not** establish rank ten for this new system, particularly after DUE changes
and three scored rows are removed. No new rank computation was performed in
this review.

The weakest old joint direction is dominated by child-benefit curvature and
first-birth fixed cost (unit-vector coefficients -0.856 and -0.473). The removed
family-rooms moment had a substantial weighted response, approximately -0.647,
along this direction. Demoting it removes some information precisely where the
old system was relatively weak. The retained childlessness, one-child share,
mean first-birth age and early-fertility moments must now identify the initial
and continuation fertility scales, curvature and fixed cost after the benefit
level is normalized. This may work, but requires a fresh local sensitivity
check rather than a count argument. The benefit normalization also moves with
all coordinates; do not hold it fixed when forming calibration derivatives.

The bequest parameter had the smallest old individual weighted column norm,
1.750. Dropping older wealth dispersion leaves aggregate wealth and bequest
flow to carry more of the discount-factor/bequest distinction. This does not
prove underidentification, but a single near-flat bequest profile can make an
apparently square system fragile. Retain parameter profiles and bound flags.

The recent-parent ownership gap is potentially informative about the tenure
scale because it measures a change in ownership probabilities across groups,
not just the overall ownership level. It is not an exclusive identifying
moment: the ownership utility shifter chi changes latent tenure advantages,
and the first-child housing loading changes parents' desired housing and
financing exposure. The four housing moments jointly supply variation across
levels, services and parent differences, but a new tenure-scale column is
necessary to distinguish these mechanisms. Diagnose near-collinearity among
chi, the housing loading and tenure scale before saying the gap identifies
one of them. Inspect ownership probabilities near indifference, not just
aggregate ownership.

## Tenure-scale bound

The proposed interval $[0.001,0.1]$ is a reasonable **exploratory** positive
interval for log sampling. It contains the inherited 0.005 and proposed 0.05,
allowing a fivefold reduction below the inherited value and a twentyfold
increase above it. It spans two orders of magnitude; “moderate” should refer
to an exploration budget, not to a small economic change. The scale has utility
units, so its admissibility does not follow from the number 0.1 alone.

Sample uniformly or propose locally in $\log\kappa_{tenure}$, and retain both
0.005 and 0.05 among deliberate initial anchors. Keep sigma and the utility
normalization fixed across lanes. The lower endpoint approaches sharp tenure
choice and may aggravate grid kinks; the upper endpoint can flatten ownership
responses and become substitutable with chi. Boundary selection is a request
for a profile/bound review, not automatic permission to widen it. These are
exploratory bounds, not externally estimated restrictions.

Finally, the old beta/H0/chi central differences showed finite-step nonlinear
responses. Even a new full-rank matrix should not trigger an unrestricted
Newton update without narrower-step checks in those columns. Use search
improvement and validated full moments as evidence; do not label a finite-budget
termination as optimizer convergence.
