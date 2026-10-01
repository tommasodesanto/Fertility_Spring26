# Winner 31: saved-policy lifecycle readout

**Identity and scope.** This is the original winner 31 fixed-price experiment at
`q0=0.719168368828958`, using its original broad-reference PRE distribution
(SHA-256 `c6589a65ad74e624579c7abe2e4069ab4e07712f299496deefda6ec7a6d5e48a`),
saved baseline and broad-credit policies, and the broad-credit stationary PRE.
It is not the later mortgage-control baseline. All quantities are at the frozen
price, in four-year ages 18–42, for households with no children ever born and
none at home immediately before fertility. `receipt.json` authenticates the
source, parameter table, grid, exact baseline PRE replay, policy source, mass,
and first-birth totals. There were **zero household, lifecycle, or GE solves**.
The remote read-only reduction took 95 seconds and used one numerical thread.

`lifecycle_validated.csv` is the analysis table. `lifecycle.csv` is the raw
reducer output, `age_summary.csv` aggregates first-birth flows, and
`postprocess_checks.json` records independent CSV reconciliation. Scenario A
uses the reference policy on reference PRE, B uses broad credit on reference
PRE, and C uses broad credit on its own stationary PRE. Thus B−A measures the
policy response at fixed inherited state; C−B is stationary-distribution
accounting at fixed broad-credit policy. Neither is a causal saving experiment.

The first policy divergence is already at age 18: the first-birth rate rises
from 0.26413 (A) to 0.28091 (B) and the owner destination share from 0.1930
to 0.3373. A and B have exactly the same age-18 inherited distribution. The
first PRE-distribution divergence is age 22. Among inherited renters with no
children, the mean PRE financial balance is 0.3425 under B's reference PRE and
−0.0067 under C's broad-credit PRE; the inherited negative-balance share rises
from zero to 0.5416. Their first-birth rate under the **same credit policy**
falls from 0.25789 to 0.23397. Among inherited owners, mean PRE financial
assets shift from −1.7241 to −2.3673, and the first-birth rate falls from
0.58030 to 0.32662. The owner-origin group also grows from mass 0.00510 to
0.01277, so these conditional averages alone cannot isolate a wealth effect.

At ages 22–30, C−B first-birth flows total −0.004902159, or about 90.7% of
the full −0.005403706 first-birth distribution component. The contribution
within inherited-renter groups is −0.008453852 and within inherited-owner
groups +0.003551693, reflecting reallocation between groups as well as their
conditional rates. Age 18 has no distribution component. The full first-birth
flows reconcile to A 0.049912692675, B 0.052734988853, C 0.047331282986.
The finer age × income × inherited tenure/housing × location × number of
children × children-home decomposition in the parent README attributes the
full C−B first-birth difference to within-group financial-asset distribution
−0.007221411500, group mass +0.001817705633, and no unmatched groups. That
is an accounting localization, not a controlled wealth intervention.

Under the broad-credit policy, the mean chosen next financial balance and
nonhousing consumption are lower in C than B at age 22 for both inherited
renters (next balance −1.0953 versus −0.8656; consumption 1.4910 versus
1.5351) and inherited owners (−1.9690 versus −1.1243; consumption 1.7387
versus 2.3213). These are distinct populations evaluated under the same policy,
not a fixed-household saving response. The next-minus-current balance column
is an asset-change proxy, not an accounting measure of saving. Net negative
financial assets among *destination renters* reveal unsecured borrowing in
that tenure, whereas negative balances among owners combine mortgage and
other liabilities and cannot distinguish them. No state-specific lifetime
solvency floor was saved, so exact mass binding that floor is unavailable;
zero-balance and grid-minimum masses are descriptive only.

**Reference DUE stayer caveat.** The reference policy has a separate
incumbent-owner stayer `bp_pol_stay`/`c_pol_stay` array. The raw reducer used
the generic `bp_pol`/`c_pol` for all realized owner mass. Consequently the raw
reference-policy, inherited-owner rows do **not** correctly measure mean next
assets, consumption, owner next-negative mass, or owner grid-minimum mass.
`lifecycle_validated.csv` blanks exactly those six fields. First-birth flows,
PRE assets, inherited negative-balance shares, destination tenure, and renter
next-balance fields remain valid. All B and C balance/consumption fields remain
valid because broad credit disables DUE stayer treatment. Source semantics:
`code/model/tools/run_dynamic_population_transition.py` `evaluate_period`
and `code/model/intergen_eqscale_seq_optimized/solver.py`
`realize_stayer_cross_section`/forward transition.

This evidence shows that broadened credit immediately raises age-18 births
and ownership, while its stationary cohorts enter subsequent ages with more
negative financial balances and lower first-birth rates at fixed credit
policy. It does **not** establish that renter unsecured credit rather than
new-buyer gate removal or incumbent-owner treatment caused the reversal. A
decisive next experiment must separate those active rules while holding the
winner 31 economics, price, entry, and PRE reference fixed.
