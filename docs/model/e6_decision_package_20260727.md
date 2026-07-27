# E6 calibration decision package

**Status: working package, not an adoption recommendation.** E6a is certified.
E6b and E6a+E6b production refits are running. E6c is governed by the
preregistered trigger below and has not been implemented. This document becomes
final only after every launched collector is certified, the E6c gate is
resolved, and the recommended configuration reproduces its full graph packet.

## 1. Decision question and fixed contract

The standing comparison is whether a small, default-off extension improves the
certified E5b fit enough to justify its economic and computational cost. The
author, not this exercise, decides whether any extension becomes the maintained
model.

Every reportable comparison uses the same twelve moments, target values,
weights, and ten baseline free parameters. The following external restrictions
do not move:

- the financed housing share is `phi = 0.80`;
- tenure choice is deterministic;
- the child-dependent bequest shifter is zero;
- the annual discount factor is at least `0.94`;
- no median statistic is a hard target;
- no housing gate is introduced.

Only strict, twice-repeated collector winners are reportable. Losses are
comparable because the target system and weights are identical. A diagnostic
graph packet is usable only after an independent strict solve reproduces all
twelve certified model moments within `1e-6`.

## 2. Ranked evidence available so far

| Rank among certified rows | Configuration | Signed loss | Change from E5b | Strict solve time | What it buys | Main price |
|---:|---|---:|---:|---:|---|---|
| 1 | E6a terminal fecundity tail | 312.094 | -73.049 | 9.25 s | Childlessness rises from 0.0802 to 0.1026; mean first-birth age moves toward target | 30+ share barely changes; old wealth dispersion remains near 2.0 |
| 2 | E5b baseline | 385.143 | -- | 12.45 s | Fits eight of twelve signed rows closely | Misses childlessness, timing shape, aggregate wealth, and old wealth dispersion |

E6b and E6a+E6b are omitted from the ranking until their collectors finish.
Loose search cases are not evidence.

## 3. Certified E5b baseline

E5b has loss `385.1428801590` with residual `1.11e-6`. Its four material
misses are:

| Moment | Target | E5b | Gap | Loss contribution |
|---|---:|---:|---:|---:|
| Childlessness | 0.1880 | 0.0802 | -0.1078 | 199.70 |
| First births at 30+ | 0.2701 | 0.3323 | 0.0623 | 60.59 |
| Aggregate wealth / annual gross labor earnings | 6.8731 | 5.4548 | -1.4183 | 12.65 |
| Old total wealth p90/p50 | 3.4481 | 2.0684 | -1.3797 | 108.43 |

The autopsy establishes three distinct failures. First, the model lacks a
permanent source of earnings dispersion. Second, the hand-set biological tail
allows too many late first births. Third, the entry margin has no mechanism
that produces the observed 25--30 timing shape. Separately, housing prices have
negligible influence on births; the author reserved any housing-gate decision.

Certified baseline artifacts:
`output/model/eqscale_seq_e5b_recalibration_20260725/report/`.

## 4. E6a: externally fitted terminal fecundity tail

E6a preserves the fitted four-year conception schedule through age 40 and
decays success to `0.03` immediately before the hard close at 45. The decay is
external demographic evidence, not a free calibration parameter. The gate is
default-off and is used identically in household optimization and population
flow.

Eight of eight production chains are eligible. The selected chain 3 winner has
loss `312.0938643710`, residual `4.73e-6`, and exactly equal strict repeats.
All ten estimates satisfy their restrictions. Two estimates lie within two
percent of the width of their broad search ranges: the continuation noise scale
is near its lower bound, and the bequest curvature shifter is near its lower
bound.

Relative to E5b:

| Moment | E5b | E6a | Model change | Absolute-gap change | Loss-contribution change |
|---|---:|---:|---:|---:|---:|
| Completed fertility | 1.9036 | 1.8926 | -0.0110 | +0.0110 | +0.6266 |
| Childlessness | 0.0802 | 0.1026 | +0.0225 | -0.0225 | -74.5212 |
| Mean first-birth age | 25.6374 | 25.5185 | -0.1190 | -0.1190 | -1.0179 |
| First births at 30+ | 0.3323 | 0.3309 | -0.0014 | -0.0014 | -2.7366 |
| First-child rooms response | 0.6572 | 0.6987 | +0.0415 | +0.0271 | +1.0182 |
| Larger-family rooms contrast | 0.3722 | 0.3576 | -0.0147 | +0.0056 | +0.2427 |
| Family ownership gap | 0.1659 | 0.1647 | -0.0012 | +0.0012 | +0.0790 |
| Ownership | 0.5638 | 0.5775 | +0.0137 | -0.0097 | -0.1608 |
| Mean occupied rooms | 6.0755 | 6.1728 | +0.0973 | +0.0973 | +0.8019 |
| Aggregate wealth / earnings | 5.4548 | 5.6111 | +0.1564 | -0.1564 | -2.6357 |
| Bequest flow / wealth | 0.009078 | 0.008983 | -0.000095 | -0.000095 | -0.2267 |
| Old wealth p90/p50 | 2.0684 | 2.0339 | -0.0344 | +0.0344 | +5.4814 |

Economics: E6a truncates late opportunities and thereby raises the mass that
finishes childless. It does not create the missing 25--30 timing shape, and it
cannot create an earnings or wealth tail.

The independent diagnostic solve reproduces all twelve certified moments and
the scalar loss exactly. The stable 17-plot packet passes visual inspection.
Artifacts:

- `output/model/eqscale_seq_e6a_recalibration_20260727/report/`
- `output/model/eqscale_seq_e6a_diagnostic_packet_20260727/`

## 5. E6b: measured permanent earnings levels

The PSID measurement uses 89,936 person-years for 12,508 reference persons
ages 25--60, 1984--2019. A long-lag minimum-distance decomposition of log
household labor earnings estimates:

| Component | Variance |
|---|---:|
| Fixed effect | 0.3931 |
| Persistent AR(1) | 0.3319 |
| Transitory | 0.3098 |

Annual persistence is `0.8863`. The fixed-effect variance has bootstrap
standard error `0.0333` and 95 percent interval `[0.3300, 0.4528]`; alternative
lag and weighting choices give `0.365--0.419`. This replaces the invalid
shortcut that treated most between-person variance as fixed.

E6b maps the measured fixed variance into three mean-one Gauss--Hermite levels,
with multipliers `[0.2775, 0.8220, 2.4347]` and weights
`[1/6, 2/3, 1/6]`. Their product with the existing five-state process produces
15 earnings states. Permanent levels do not change after entry. The structure
is standard lifecycle heterogeneity in the tradition of Huggett (1996),
Storesletten, Telmer, and Yaron (2004), and Kaplan and Violante (2010); it is
not a novel household type mechanism.

The fixed-parameter diagnostic shows exactly what the state buys before
refitting:

- old wealth p90/p50 rises from `2.07` to `3.91`;
- aggregate wealth / earnings rises from `5.45` to `5.62`;
- childlessness barely moves, from `0.080` to `0.084`;
- the childlessness gradient across permanent levels has the wrong sign;
- the family ownership gap falls from `0.166` to `0.095`.

Thus E6b is a wealth-dispersion mechanism, not the missing childlessness
mechanism at fixed parameters. The strict production refit remains pending.
No free parameter was added.

Empirical packet:
`code/data/psid_followup_mar2026/output/psid_income_fixed_effect_md_20260727/`.

## 6. E6c preregistered decision rule

The signed 30+ first-birth share has target `0.270062` and declared standard
error `0.008`. E6c activates only if the certified E6a+E6b estimate remains
more than two standard errors (`0.016`) from target **and** the supplemental
age-bin diagnostic still displays the missing 25--30 shape.

If triggered, E6c is the handoff's binary readiness state. An exogenous,
age-dependent event moves an entrant irreversibly from unsettled to settled;
first-child entry is available only when settled. The hazard has exactly two
free coordinates, location and spread, so the system has twelve free
parameters and twelve signed moments. A strict local finite-difference check
must show rank two for the two timing rows with respect to the two hazard
coordinates before the refit can be called identified.

No E6c code is written before the trigger. If either trigger condition fails,
the architecture is not added.

## 7. Search and reporting safeguards

- Every refit uses two baseline seeds, three `start_mix=0.10` seeds, and three
  `start_mix=0.25` seeds.
- Each exact loop was smoke-tested before its production array.
- Every chain writes completed cases, best-so-far state, and a final summary.
- Collectors require strict convergence and two exactly equal tight repeats.
- Every certified report includes all twelve target rows and all free
  parameters with bounds.
- The E6a and E6b launches inherited a mistakenly declared annual-discount
  search lower bound of `0.80`. All certified/current best estimates satisfy
  the standing `0.94` external restriction. The code was repaired before the
  combined run; this discrepancy remains disclosed rather than rewritten.

## 8. Author decisions reserved

The final package will leave these decisions to the author:

1. whether any E6 variant becomes the maintained specification;
2. whether to introduce a housing gate capable of making births materially
   price-sensitive;
3. whether to revise any target or weight after reviewing the completed
   mechanism evidence.

This exercise makes no target-system change, housing-gate choice, or
policy-experiment claim.

## 9. Attempts that did not solve the stated problem

- Refitting only the conception levels at ages 30, 35, and 40 did not remove
  the terminal excess.
- The E6a tail lowers the signed loss but does not create the 25--30 timing
  shape or the wealth tail.
- Permanent earnings levels repair wealth dispersion at fixed parameters but
  do not generate the empirical childlessness gradient and damage aggregate
  ownership margins.
- Reweighting diagnostics do not change the diagnosis: childlessness and old
  wealth dispersion dominate the loss because mechanisms are missing, not
  because the signed weights are unusual.

## 10. Items required before this package is final

- certify the E6b collector;
- certify the E6a+E6b collector;
- apply the preregistered E6c trigger;
- if E6c triggers, implement, identify, smoke, refit, and certify it and any
  required combination;
- regenerate the ranked comparison tables with every certified variant;
- select a recommendation and include its complete twelve-row target table,
  complete free-parameter/bounds table, and independently verified standard
  graph packet;
- run the completion audit against the handoff contract.
