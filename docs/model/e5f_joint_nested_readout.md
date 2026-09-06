---
title: "Simultaneous fertility and housing choice"
subtitle: "First experimental results | read page 1 first"
author: "Prepared for Tommaso De Santo"
date: "6 September 2026"
---

# The useful result

**Your proposed simultaneous nested-logit structure is implementable. The first experiment passes its checks. The immediate quantitative issue is the shock scales, rather than a large difference between the two tested choice rules.**

I evaluated 13 scale/correlation combinations under joint choice and a sequential control: **26 completed cases**. Both use the same four plans: rent/wait, rent/attempt, own/wait, own/attempt. Prices, the inherited 2023 population and future values are held at the retained benchmark.

- **Changing the rule has a tiny effect here.** Across the panel, the largest joint-versus-sequential difference is about **one birth unit per million households** over a four-year period. The largest ownership difference is **0.0000054 percentage points**. This is evidence about this particular common menu and frozen continuation, not proof that shock timing never matters.
- **Changing the scales has a large effect.** Keeping tenure dispersion at **0.005** makes current births essentially zero. At outer scale **2.5**, with independent alternatives, birth units are **9.785 per 100 households**, versus **8.478** in the retained reference; ownership is **42.29%**, versus **63.56%**. These are all-age current quantities, not the calibration's prime-age ownership target.
- **The next useful step is a fully solved experimental lifecycle model.** Replace future choice values as well as today's rule, clear the housing market, and then evaluate all twelve original moments. This first test does not establish an improved calibration or resolve the first-birth rooms target.

![Solid joint-choice and dashed sequential-control lines nearly coincide. Different colors change within-tenure dispersion relative to the outer scale.](output/model/e5f_joint_nested_experiment_20260906a/panel_comparison.png){width=100%}

Production remains unchanged. Pages 2-3 state the experiment and all its results; pages 4-5 preserve the complete reference calibration record.

\clearpage

# What was actually changed

A household observes all current taste components and chooses a tenure and a birth attempt jointly. The taste components have a nested extreme-value distribution: alternatives sharing tenure have correlated tastes. **The nesting describes correlation, not a chronological order of decisions.**

The conditional value $q(\tau,x)$ is the best deterministic housing-product value within tenure $\tau$ at household state $x$, using the original budgets, transaction costs, borrowing limits and baseline future values. There is exactly one location. The birth attempt succeeds with the original age-specific probability $\pi$; $x^+$ adds one birth and child at home. The first-success cost $F$ is zero for subsequent births. The two plan values within each tenure are

\[
Q_{\tau,0}=q(\tau,x),\qquad
Q_{\tau,1}=(1-\pi)q(\tau,x)+\pi[q(\tau,x^+)-F].
\]

**The first experiment commits tenure across success and failure.** Housing size, consumption and saving may respond to conception within that tenure. Original renter housing is continuous and the five owner sizes remain available, but owner-product tastes are removed from this experimental block. These are explicit additional restrictions; comparison with production is not a pure timing decomposition.

The outer tenure scale $\kappa$ and dissimilarity $\lambda$ determine the within-tenure scale $\lambda\kappa$, with $0<\lambda\leq1$. Mean-zero-shock expected values satisfy

\[
S_\tau=\lambda\kappa\log\sum_a\exp\{Q_{\tau,a}/(\lambda\kappa)\},
\qquad V=\kappa\log\sum_\tau\exp\{S_\tau/\kappa\}.
\]

The joint probability is the derivative of $V$ with respect to the corresponding plan value. At $\lambda=1$ it is flat four-alternative logit. As $\lambda$ approaches zero, fertility becomes deterministic conditional on tenure, linking this architecture conceptually to the simplified theory's renter and owner maximization problems. The lifecycle and fertility state spaces still differ. The joint-law foundation is standard [nested GEV; Train (2009), chapter 4](https://eml.berkeley.edu/books/choice2nd/Ch04_p76-96.pdf).

The sequential control instead first integrates over tenure at scale $\kappa$, then over attempts at scale $\lambda\kappa$, using exactly the same plan values. It changes the shock law along with revelation; for $\lambda<1$ it is not a reversed nested-GEV law.

**Verification:** seven independent mathematical/accounting tests; two exact reproductions of twelve reference arrays; 26 verified case checkpoints; four exact smoke-to-panel matches; flat-logit equality also checked on full saved arrays. Every case has valid probabilities, conserved population and birth stocks/flows, and zero flagged occupied wealth-value decreases. The original seventeen reference graphs and all 26 supplemental choice plots were inspected.

The known benchmark reporting-floor budget exception remains: $2.914\times10^{-9}$ of household mass. No case increases budget-violating mass at any destination state in the additional-exposure audit. The failed initial absolute-budget smoke and its resolution are preserved. This is not a claim of perfect benchmark feasibility or a full equilibrium solution.

\clearpage

# Complete experimental panel

$\kappa$ is the outer tenure scale; $\lambda\kappa$ is the within-tenure fertility scale. All new scales are diagnostic, not estimated. The table shows **joint-choice** outcomes. Birth units include the original adjustment for the top-coded 3+ family category, over one four-year model period. The machine-readable cases also retain explicit births without this adjustment.

```{=latex}
\begin{center}\small
\begin{tabular}{rrrrrr}
\toprule
$\kappa$ & $\lambda$ & Birth units / 100 HH & Own (\%) & Rooms / HH & Market gap (\%)\\
\midrule
0.005 & 0.05 & 0.000000 & 63.123 & 6.30289 & -0.229 \\
0.005 & 0.5 & 0.000000 & 63.123 & 6.30289 & -0.229 \\
0.005 & 1 & 0.000000 & 63.123 & 6.30289 & -0.229 \\
0.05 & 0.05 & 0.000000 & 49.162 & 5.87401 & -7.018 \\
0.05 & 0.5 & 0.001532 & 49.162 & 5.87401 & -7.018 \\
0.05 & 1 & 0.038342 & 49.163 & 5.87405 & -7.017 \\
0.5 & 0.05 & 0.001507 & 42.792 & 5.65830 & -10.433 \\
0.5 & 0.5 & 1.063243 & 42.792 & 5.65972 & -10.410 \\
0.5 & 1 & 2.358372 & 42.790 & 5.66231 & -10.369 \\
2.5 & 0.05 & 0.391825 & 42.293 & 5.63709 & -10.768 \\
2.5 & 0.5 & 5.934513 & 42.290 & 5.64927 & -10.576 \\
2.5 & 1 & 9.784676 & 42.287 & 5.66115 & -10.387 \\
0.5 & 1e-06 & 0.000000 & 42.792 & 5.65830 & -10.433 \\
\midrule
\multicolumn{2}{l}{Retained reference} & 8.477610 & 63.562 & 6.31729 & -0.001\\
\bottomrule
\end{tabular}\end{center}
```

The signed market gap is demand minus fixed supply, divided by supply. The larger-scale cases reduce demand by roughly 10% at the old price. They therefore cannot be treated as equilibria. Values reported as zero below the displayed precision are economically negligible; at the near-deterministic limit, explicit current births are exactly zero.

The following differences are **joint minus sequential**, on the same menu. The flat-logit cases ($\lambda=1$) agree to floating-point precision. Small-scale near-zero birth flows and the deterministic-limit case also have negligible or zero differences.

```{=latex}
\begin{center}\small
\begin{tabular}{rrrrr}
\toprule
$\kappa$ & $\lambda$ & Birth units / 100 HH & Ownership (pp) & Rooms / HH\\
\midrule
0.05 & 0.5 & 0.00000148 & -0.00000031 & -0.00000001 \\
0.5 & 0.05 & 0.00000321 & -0.00000000 & 0.00000000 \\
0.5 & 0.5 & 0.00007388 & 0.00000535 & 0.00000232 \\
2.5 & 0.05 & 0.00010634 & -0.00000002 & 0.00000140 \\
2.5 & 0.5 & 0.00004259 & 0.00000101 & 0.00000977 \\
\bottomrule
\end{tabular}\end{center}
```

**Interpretation for further work.** The retained first- and subsequent-birth scales are 2.168 and 1.736, both larger than the externally fixed tenure scale 0.005. That ordering cannot be carried into these tenure nests. This is a specification and identification issue. Making tenure dispersion estimable would add a free parameter unless another restriction replaces it; the full moment-to-parameter mapping and Jacobian must be reassessed before calibration. A count of moments alone does not establish informative identification.

Future values in this experiment still anticipate the old shock process. The low-current-birth outcomes can therefore reflect waiting for future options under that process. They do not establish that low-dispersion simultaneous choice is incapable of fitting fertility after the whole model is solved again.

\clearpage

# Reference calibration: every target

The retained September 4 `task_010` has loss **30.4829667**, with eleven estimated parameters and twelve target rows. This is the reference behind the frozen 2023 state, not a fit of simultaneous choice. The newer overnight candidate remains unpromoted and is not substituted here. Gaps equal model minus target; each loss contribution is the displayed weight times the squared gap, using unrounded values.

| Moment | Target | Model | Gap | Weight | Loss |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.922907 | 0.004907 | 1425.739 | 0.0343 |
| Childlessness | 0.188000 | 0.189407 | 0.001407 | 17180.744 | 0.0340 |
| Mean age, first birth | 26.044627 | 26.256032 | 0.211404 | 44.444 | 1.9863 |
| First births at age 30+ | 0.260327 | 0.237579 | -0.022748 | 10000.000 | 5.1748 |
| Rooms response, first birth | 0.720246 | 0.436418 | -0.283828 | 137.565 | 11.0821 |
| Rooms, 3+ minus 1-2 children | 0.367700 | 0.403441 | 0.035741 | 2958.515 | 3.7793 |
| Family ownership gap | 0.167662 | 0.161699 | -0.005963 | 14229.591 | 0.5060 |
| Ownership, ages 30-55 | 0.575472 | 0.544488 | -0.030984 | 1207.846 | 1.1596 |
| Mean occupied rooms | 5.779970 | 6.317291 | 0.537321 | 11.973 | 3.4568 |
| Wealth / annual earnings | 6.873100 | 6.932652 | 0.059552 | 6.288 | 0.0223 |
| Annual bequests / wealth | 0.008800 | 0.008433 | -0.000367 | 5165289.256 | 0.6971 |
| Old wealth/income p90 / p50 | 3.448111 | 3.236511 | -0.211600 | 56.960 | 2.5503 |

The last row uses ages 76-84. The first-birth room response remains the empirical event-study target, with its original contract and weights. No target is dropped, demoted, reweighted or replaced by a current-flow statistic in this experiment.

The all-age ownership rate on page 3 is 63.56% in the reference. The ownership calibration row above uses ages 30-55 and is 54.45%. Their difference is a sample definition, not a failed reproduction. Likewise, four-year birth flows are distinct from completed fertility and from the matched historical room response.

Exact unrounded records are retained in `reference_target_fits.csv` and `reference_parameters.csv` in the experiment's result folder. The original target fingerprint is

`3726c17e62c8233ce62d5f4c95f44fd2cc2ea6cfa3d2492795461b4569300497`.

\clearpage

# Reference calibration: every parameter

These estimates define the original continuation values and conditional housing/saving problems. The experimental current-date taste block uses the diagnostic scales on page 3. It does not re-estimate any preference, housing, fiscal or demographic parameter.

| Parameter | Estimate | Bounds or restriction | Status / bound flag |
|---|---:|---|---|
| Annual discount factor | 0.995117 | [0.94, 0.9995] | Estimated; no flag |
| First-birth dispersion | 2.168173 | [0.02, 50] | Estimated; no flag |
| Subsequent-birth dispersion | 1.736471 | [0.02, 50] | Estimated; no flag |
| Owner-service premium | 1.043472 | [0.1, 5] | Estimated; no flag |
| Housing-supply intercept | 14.562959 | [0.2, 80] | Estimated; no flag |
| Bequest strength | 0.528428 | [0, 8] | Estimated; no flag |
| Bequest wealth shifter | 0.107249 | [0.02, 16] | Estimated; raw-scale flag |
| Per-child room floor | 0.282210 | [0.1, 1.8] | Estimated; no flag |
| First-success fixed cost | 4.559138 | [0, 8] | Estimated; no flag |
| First-child room jump | 0.364931 | [0, 0.5] | Estimated; no flag |
| Child-value change, 2007-2023 | -0.328714 | [-1.5, 0.2] | Estimated; no flag |
| Child-value intercept, 2007 | 0.288017 | Externally normalized | Not free |
| Child-value intercept, 2023 | -0.040697 | Derived from intercept + change | Not free |
| Tenure dispersion | 0.005000 | Externally fixed | Not free |
| Housing-supply elasticity | 0.630000 | Externally fixed | Not free |

The bequest wealth shifter's stored flag uses raw parameter units. Its normalized logarithmic search coordinate is 0.251237, so it is not at the lower edge in the coordinate actually searched. All other estimated parameters have no stored near-bound flag.

**Evidence and reproduction.** The design, exact source reconstruction, failed development smoke, independent review and its resolution, contracts, all case states and the collection checks are indexed in

`output/model/e5f_joint_nested_experiment_20260906a/README.md`.

Final Torch smoke **17068184** and panel **17068310** completed with exit code zero. The panel required no additional Bellman solve. Including the two development smokes, total allocated-job elapsed time was **6 minutes 35 seconds**, within the 55-minute first-stage cap. No cluster job remains running for this experiment.

The final scientific bundle is identical to the retained reference. The experimental driver and collector are separate tools. Reproducing the full lifecycle and historical calibration under simultaneous choice remains the next substantive experiment, not a completed result of this packet.
