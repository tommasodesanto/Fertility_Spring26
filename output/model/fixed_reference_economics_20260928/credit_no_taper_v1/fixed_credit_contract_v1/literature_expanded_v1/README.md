# Unsecured credit: focused literature comparison

September 30, 2026. Read-only source review by two Sol 6.1 workers, with lead
review of the economic mapping. No model changes, recalibration, or solves.

- [Kaplan–Violante and Kaplan–Moll–Violante evidence](violante_evidence.md):
  KV2014 uses an income-dependent allowance of 74% of quarterly current income
  (18.5% annualized), from the median SCF credit-limit ratio; zero in retirement.
  KMV2018 uses a common allowance of one quarter of mean annual labor income,
  an external benchmark. Both distinguish borrowing and saving rates.
- [Additional comparators](comparators_evidence.md): Maxted–Laibson–Moll2025
  uses a constant one-third annual permanent-income allowance; SCF homeowner
  credit-limit ratios provide its empirical anchor. The manuscript-only
  Herkenhoff comparison has unresolved measurement/version issues and is not
  used as a numerical benchmark.
- [Published Boar–Gorea–Midrigan verification](../literature_published_v1/README.md):
  7.324% mean annual income in the official replication; a separate empirical
  rationale for this assigned value was not verified.

## Implication for our fixed reference

Reference: **2007 stationary reference — block0506, September 28 verified export**.
The manifest's serialized parameters imply a mean income-state multiplier of
1 and equal-weight working-age mean annual gross income of 1:
mean(income[0,:J_R]) / period_years / (1-tau_pay) times dot(z_weights,z_grid).
Thus the experimental D=.14 equals 14% of that annual gross-earnings
normalization. This is not the mean disposable income of all occupied households.
These income denominators differ across papers; the percentages are not directly
interchangeable. A debt stock is not multiplied by the four-year decision period.

The new SCF-based evidence favors prioritizing an empirically disciplined
positive allowance over choosing credit solely to repair entrant feasibility.
For a constant scalar D, use credit-capacity evidence for the relevant population
and state the approximation. Borrowing spreads require an explicit separate
decision: the current symmetric 2% annual rate differs from these benchmarks.
Zero credit remains a useful comparison; no parameter or entry transfer has
been adopted by this review.
