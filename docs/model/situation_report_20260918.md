---
title: "Situation report, September 18, 2026"
subtitle: "Where the model stands after the overnight tests, and what to decide"
date: "September 18, 2026, 10:00"
---

# The short version

The baseline is solid: the paper's code is on `main`, tagged, and the
September 14 initial state replays exactly. The sandbox solves one change at a
time in a few minutes on that code. Overnight it produced the first full set
of one-change comparisons, and they change what the paper can claim.

1. **The mechanism on the slide does not operate at this calibration.**
   Removing the down payment raises ownership by ten points and leaves
   fertility and its timing unchanged; giving households unsecured credit
   raises fertility and brings births forward but leaves ownership unchanged.
   Making space scarce (rooms at the data's level) does not change this. At
   these parameters, births respond to liquidity and to the price of space,
   not to the housing-collateral requirement, which governs ownership.
2. **Children at home are net costs at the retained parameters.** The value
   of a household falls at every state when a child stays one more period.
   Fertility is held up by the taste shocks and by the child-preference level
   being set to hit 2.1. This is the P1 problem in its sharpest form, and it
   means any change that lengthens dependency or adds a child cost must be
   compared with that level re-normalized, not fixed.
3. **The maturation fix works.** A parent-age exit law with a newborn
   exemption makes children at home by parent age match the ACS at both ends,
   where the current law is 40 percent short at ages 30–42 and keeps adult
   "dependents" past 58.
4. **Two switches move the fit a lot, in the right places.** A size-dependent
   rental wedge in place of the cap and the owner premium cuts the loss from
   2607 to 402 (ownership overshoots and would be calibrated down). A
   20 percent earnings cost of children at home cuts it to 1746 and moves
   first births later and the rooms response toward its target.
5. **Earnings:** one consistent PSID decomposition brings the old-age wealth
   tail to its target (3.38 against 3.52); the Boar–Gorea–Midrigan process
   without permanent types loses the tail (2.9) and gives too few childless
   households. Permanent types are what the tail buys.

# 1. The mechanism, with the evidence

All numbers below are steady states at the retained parameters. "Fixed ψ"
holds the child-preference level at its retained value; "ψ root" re-sets it so
completed fertility is 2.1. The sandbox's baseline is close to, not equal to,
the slides' levels (its ownership is 0.46 against the slides' 0.54; the
reason is a step in the original launch script the sandbox reproduces only
approximately). Differences between columns are the object; levels are not.

**Frictionless benchmark (fixed ψ).**

| Moment | Baseline | No down payment (φ = 1) | Unsecured credit line (five years' earnings) |
|---|---:|---:|---:|
| Completed fertility | 1.872 | 1.866 | 1.948 |
| Childless 40–44 | 0.239 | 0.241 | 0.225 |
| Mean first-birth age | 26.96 | 26.97 | 25.79 |
| Ownership 30–55 | 0.459 | 0.557 | 0.450 |
| Wealth / earnings | 5.18 | 5.07 | 4.43 |
| First-birth rooms response | 1.07 | 1.15 | 1.17 |

**Who is constrained.** The down payment binds for 4.5 percent of
family-forming households (ages 26–38, zero or one child at home). With the
supply scale lowered until mean rooms equal the data (5.56, price up
4 percent), that share falls to 2.8 percent, and fertility is flat with and
without the down payment (1.840 against 1.835). Scarcer space lowers fertility
a little through cost (1.872 to 1.840). The author's prior that the rooms
misses hide the mechanism is rejected on this test.

**Why.** A renter can have six rooms and the parents' space floor is 2.3
rooms, so a family with one or two children fits in a rental; the owner-only
sizes start at eight rooms and matter for three-plus families. What blocks an
early birth is that a renter cannot borrow at all and must save the child's
cost in advance. That is a precautionary channel, and it is real, but it is
not the collateral story.

**What survives.** Easier mortgages raise ownership and not births; easier
unsecured credit raises births and not ownership; housing costs lower
fertility through the price of space. The property-tax channel to young
families runs through ownership and the allocation of large homes, so its
fertility effect will be small unless the space margin binds. The slide's
mechanism sentence has to change; the model, the calibration machinery and
the experiment stand.

# 2. Decisions, with their evidence

Every row of the table below has completed fertility 2.1 (ψ re-normalized).
Baseline loss 2607, ψ 0.190. Losses use the sandbox's twelve-row objective and
are comparable across rows, not with the slides' 179.

| Change | Loss | ψ | Childless | First-birth age | Ownership 30–55 | Old p90/p50 | Rooms resp. | Recent-parent gap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Baseline | 2607 | 0.190 | 0.201 | 26.1 | 0.455 | 4.32 | 1.12 | 0.47 |
| Rental wedge, no cap, χ = 1 | 402 | 0.205 | 0.208 | 26.3 | 0.730 | 3.83 | 0.97 | 0.28 |
| Earnings penalty 20 % while children home | 1746 | 0.315 | 0.228 | 27.1 | 0.416 | 4.46 | 0.85 | 0.39 |
| Mortgage block (origination-only, 11 % amortization) | 2050 | 0.189 | 0.201 | 26.1 | 0.349 | 5.05 | 1.08 | 0.42 |
| Boar–Gorea–Midrigan earnings, no types | 2100 | 0.156 | 0.144 | 27.0 | 0.499 | 2.92 | 1.24 | 0.41 |
| E6b earnings, one decomposition | 2325 | 0.180 | 0.184 | 26.6 | 0.477 | 3.38 | 1.14 | 0.44 |
| Tenure shock κ_H = 0 | 2379 | 0.190 | 0.202 | 26.1 | 0.451 | 4.30 | 1.12 | 0.45 |
| Parent-age maturation | 2464 | 0.167 | 0.210 | 26.5 | 0.449 | 4.36 | 1.39 | 0.45 |
| Modest unsecured credit line | 2539 | 0.187 | 0.202 | 25.8 | 0.446 | 4.32 | 1.15 | 0.46 |
| Concave child benefit (log) | 2111 | 0.199 | 0.210 | 25.9 | 0.449 | 4.15 | 1.10 | 0.43 |
| Estate transfer to ages 45–65 | 2582 | 0.181 | 0.193 | 26.4 | 0.446 | 4.07 | 1.06 | 0.46 | | | | | | | |
| All switches together | 962 | 0.480 | 0.223 | 27.6 | 0.426 | 4.44 | 0.97 | 0.27 |
| Rental cap at 8 rooms (no other change) | 489 | 0.189 | 0.202 | 26.1 | 0.300 | 4.38 | 1.36 | 0.23 | | | | | | | |

Targets: childless 0.198, first-birth age 26.0, ownership 0.648, p90/p50
3.52, rooms response 0.72, recent-parent gap 0.16.

**F3, how children leave home.** Recommend the parent-age law with the
newborn exemption (exit 0.05 per period for young parents, rising from parent
age 34 to certain exit at 62). Evidence: with ψ re-normalized, the ACS profile
matches at 22–34 and is zero after 62; the two rooms moments rise (a refit
would lower the space floor). The switch is built, tested and verified line
by line; off by default.

**P1, the child term.** Two facts: children at home are net costs at every
state, and the intensive margin is carried by the taste scales. Candidates:
keep linear and say so; a concave benefit; or the earnings penalty
(S3), which at 20 percent halves the loss, moves first births later and the
rooms response toward target, and needs ψ to double. The concave benefit
(log) cuts the loss to 2111 with the first-birth age exactly on target (25.9
against 26.0) and fewer large families, at ψ 0.199. Recommend both: the
concave benefit for the intensive margin and S3 with the penalty taken from
the child-penalty literature for the income gradient, keeping ξ only if
childlessness still needs it. Author decision.

**H3/H6, the cap and the premium.** The size-dependent wedge is the single
largest improvement and it also disciplines the ownership level from the
rent–price ratio rather than from a taste parameter. Its intercept is too
strong at 0.02 (ownership 0.73); it would be calibrated. Recommend adopting
the wedge and retiring χ, with the AHS renter share by size as its moment.

**H1, the mortgage block.** Amortization with origination-only borrowing
lowers ownership by ten points and raises old-age wealth dispersion; it does
nothing to fertility. It makes the model's finance side standard and the
Coven comparison like for like. Recommend adopting it and letting the wedge
intercept absorb the ownership level. A modest unsecured credit line brings
first births forward by a third of a year and is the cheapest change that
improves the timing rows; recommend it in the baseline at the
Kaplan–Violante size.

**E1, earnings.** The permanent types are what the old-age wealth tail buys;
without them (Boar–Gorea–Midrigan) the tail is 2.9 against 3.5 and
childlessness falls to 0.14. The one-decomposition E6b process puts the tail
at 3.38. Recommend retaining the types with the E6b decomposition, one tax
treatment throughout.

**F2, the tenure shock.** Inert in the steady state; loss slightly better at
zero. Recommend zero, with the kink handled numerically.

**P5, estates.** At fixed ψ the transfer to ages 45–65 barely moves anything
(loss 2599) and costs an hour per solve; the ψ-root result is on Torch.
Recommend the net-of-selling-cost valuation regardless; the receiver is a
welfare-accounting choice more than a fit lever.

# 3. What ran overnight, what failed, what is still running

- Built and verified, all default off, full suite 227 tests: parent-age
  maturation (with an exact second housing-stage solve after my rejection of
  a shortcut); child earnings penalty; mortgage block; rental wedge; estate
  receiver. Muse sessions died silently three times on long tasks; the work
  was recovered from disk each time and every switch was re-run through the
  suite. Two sandbox bugs found and fixed on the way: a φ override silently
  reset by the package, and a spec loader without list support.
- Solved locally: the frictionless pair, the scarce-space pair, the
  who-is-constrained tables, the ACS profiles, eight fixed-ψ single switches,
  and eight ψ-root single switches (see the table).
- On Torch: ten ψ-root jobs (ids 17938995–17939004) and the smoke were
  queued from 01:30 in a congested queue (about 2,800 jobs pending) and at
  11:00 all show CANCELLED with zero elapsed time: they never ran and were
  cancelled by a user or an administrator, not by this session. The batch
  folder on scratch was missing one folder the launcher expects; it was
  repaired at noon, the smoke passed in 39 seconds, and the two rows that
  depend on it (estate transfer and all switches with ψ re-normalized; each
  estate iteration re-solves the equilibrium, so they are too slow for the
  laptop) were resubmitted as jobs 17949060 and 17949061 and completed in
  2.2 and 8.8 hours; their rows are in the table. Estate transfer: loss 2582
  (2607), first births later by a third of a year, old-age tail 4.07 (4.32),
  nothing else moves. All switches together: loss 962, ψ 0.480 (two and a
  half times the baseline root), first-birth age 27.6 against 26.0, first
  births 30+ 0.32 against 0.25, wealth/earnings 4.4, rooms 5.25 (below the
  5.56 target for the first time), three-plus rooms gap 0.58 (target 0.35).
  The switches do not add up: the penalty and the maturation law each need a
  higher ψ and together push timing past the target, while the wedge alone
  fits better (402) than all of them together.
- Added September 19: rental cap raised from six to eight rooms, nothing
  else changed. Fixed ψ: completed fertility 1.876 against 1.872, childless
  0.238 against 0.239, first-birth age unchanged, ownership 0.31 against
  0.46, recent-parent gap 0.21 against 0.46. ψ root: loss 489, ψ 0.189
  (baseline 0.190). On the paper's model the cap is an ownership lever, not a
  fertility lever; the July 1 audit claim (TFR +0.071 from the cap) was
  measured on the earlier package and does not carry over.
- Not done: nothing was fitted, no shock estimation, no transition, no
  production default changed, per the standing rule. The recovery batch
  (job 17858740) stays stopped.

# 4. Baseline and provenance

`main` equals the tag `paper-baseline-2026-09-14` in the model package and
carries the exact solver and the pension helpers; Codex replayed the initial
recipe on Torch and reproduced all 13 moments, 17 parameters and the loss
179.298 (`output/model/paper_baseline_sep14/replay_20260917/README.md`). A
byte copy of the presented code is under
`calibration_archive/legacy_presentation_20260914/`. The sandbox uses the
recipe's own assembly functions; two parameter bugs were fixed, and its
baseline still differs from the slides' state in two wealth-distribution rows
because the recipe starts from a checkpoint that only exists on the cluster.
Rule: the cluster replay is the exact tool for levels; the sandbox is the
quick tool for differences. A stray git ref `refs/remotes/origin/main 2` blocks
automatic repacking; left for the Codex session.

# 5. Next steps, in order

1. Done: both Torch runs collected (rows above).
2. Decide F3, P1, H3/H6, H1, E1, F2 from the table; each has a recommendation
   above.
3. With the decisions made, one refit under the same target contract, on
   Torch, then the frictionless and constrained diagnostics again on the
   refit, so the mechanism statement is made once on the model the paper
   will keep.
4. Rewrite the mechanism paragraph of the slides and the paper around what
   the model does: liquidity and the price of space move births; the down
   payment moves ownership and the allocation of large homes.
5. Only then the history refit and the policy transition.

Everything above is also in the follow-up log of
`docs/model/structural_model_review_consolidated_20260915.md`, with paths.
