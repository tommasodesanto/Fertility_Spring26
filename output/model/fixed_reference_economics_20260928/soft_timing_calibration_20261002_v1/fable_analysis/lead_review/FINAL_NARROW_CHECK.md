# Final narrow check of revised memo

October 3, 2026; zero solves. Source: `fable_analysis/ECONOMIC_MEMO.md`
after the first revision, saved native arrays for selected original chain 15
and alternative chain 13. No memo or analysis code edited here.

1. Memo lines 80–82 say the only general-equilibrium financing results are
   historical permanent steady states. Its own table at line 75 lists
   accepted temporary dated 100% financing paths, with first-period birth
   flows −0.54% and −0.43%. Suggested text: “The accepted permanent
   steady-state financing comparisons are the historical hard/quarter runs;
   accepted temporary dated paths at those older points are also available.
   Neither is a soft-rule policy result.”
2. Memo lines 85–87 say the held-coordinate timing swap leaves “every
   fertility row within 0.003.” Its own table line 77 reports mean age at
   first birth 26.03→26.07, a roughly 0.038-year change. Suggested text:
   “Age-25 children ever born changes 0.5276→0.5249 (−0.0027), while mean
   age at first birth rises by about 0.038 years; both are small relative to
   the 11.2 percentage-point ownership change.” Do not combine moments with
   different units under one numeric bound.
3. Memo lines 59–60 and 190 say the model “matches” the share of age-25
   mothers. It is 0.441959 versus CPS 0.457254, a 1.53 percentage-point
   gap (`analysis/out/native_grid_analysis.json`, original `age25`).
   “Is close to” describes the evidence more accurately.
4. Memo lines 157–159 say first-birth attempt probability in the four lowest
   income states is zero “at every wealth node.” This is false on the saved
   native grid. In the original arm at age 22, low-income states 1–4 have
   full-grid maxima about 0.936 each; maxima among nodes with reconstructed
   pre-fertility renter mass >1e−14 are 0.837, 0.846, 0.855 and 0.864,
   respectively. These high-wealth nodes have very small mass. Mass-weighted
   means are 0.0000002, 0.0000034, 0.0000917 and 0.00649, matching the
   JSON's `young_childless_renters[22].attempt_prob_by_income_state`.
   Suggested text: “The four lowest income states have near-zero
   mass-weighted attempt rates, although wealthy occupied nodes within them
   can have high attempt probabilities. The feasibility split mixes income
   and wealth and is not causal.” The alternative arm shows the same pattern.

The first three items are wording corrections; item four also corrects a
false native-policy claim. These do not disturb the verified age-25 stock
decomposition.
