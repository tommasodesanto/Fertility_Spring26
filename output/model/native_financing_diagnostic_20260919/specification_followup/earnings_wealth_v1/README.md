# Earnings, entry wealth and purchase timing

**September 21 evening: implementation and cheap checks completed; model runs
await the income-period decision. No new calibration result or submitted job.**

The author authorized implementation after returning home. This packet keeps
the September 14 reference intact and prepares three diagnostic arms: the
original model once, literature income once, and literature income plus
current-income purchase accounting twice. Structural parameters are fixed at
the original retained benchmark, not the later PSID-income refit. The existing
equilibrium, pension balance and fertility normalization are re-solved within
each full objective; this is not a parameter search or policy exercise.

## The decision found by the cheap check

The published Sommer (2016, Section 4.1/Table 2) annual process supplies
\(\rho_a=0.95\), persistent innovation standard deviation \(0.21\), and iid
standard deviation \(0.17\), with no permanent type. Applying that wage-risk
process to exogenous household earnings remains an explicit approximation.
The original deterministic age profile is held fixed for this diagnostic.

Exact four-year averages of this annual level process have variance
\(0.541047309\) and first two covariances \(0.447837257,0.351691726\).
Matching all three with an independent lognormal AR(1)-plus-iid period process
would require iid log variance \(-0.021990423\). The constructor rejects this;
no variance was clipped and no household solve was attempted.

The proposed alternative uses the conventional four-year persistent transition
\(\rho_4=\rho_a^4\) with stationary variance unchanged, and iid log variance
\(\log[1+(\exp(0.17^2)-1)/4]\). Its continuous level variance is 7.84% above
the exact four-year-average variance; lag-one through lag-four covariances
differ by less than 0.6%. It approximates the aggregate flow; it is not the exact
distribution of four years of total earnings. **Author choice is pending.**

[Complete aggregation and grid receipt](income_aggregation_check.json).
The prepared grid has 15 persistent by 3 iid nodes. Its first five level
covariances differ from the intended continuous endpoint process by less than
5%; its log covariances match that process to floating-point precision. The
5% threshold is an explicit diagnostic grid criterion, not evidence of
household-policy convergence. The smaller 5-node persistent grid fails this
criterion. Do not choose a coarse grid because its variance error happens to
offset the time-aggregation error. Further household-grid validation is unrun.

## Operative inherited objects and entry treatment

The actual reference has four-year periods, age 18 entry, 17 age nodes and
retirement at 66; \(R=1.02^4=1.08243216\) and \(\beta=0.99^4=0.96059601\).
The wealth grid remains 120 points from \(-12\) to \(30\). Taxes, normalized
age profile, housing requirement, owner services, bequests, pension closure,
all target definitions and weights are held at their inherited specifications.
Generic source defaults are not substituted for these operative values.

The old entry rule scales empirical wealth/income ratios by realized earnings.
Using it unchanged would let a new iid shock mechanically alter entry wealth.
The diagnostic instead preserves the reference wealth marginal and maps old
total-income ranks to new persistent-income ranks by probability-interval
overlap. Iid draws are independent of wealth conditional on persistent rank.
This coupling is an explicit diagnostic assumption, not a newly estimated
empirical joint distribution. The native-input check preserves the wealth
marginal to \(1.12\times10^{-16}\) and mean wealth to floating-point precision.
[Full entry distribution and operative values](entry_wealth_check.json).

## Purchase accounting implemented

Let \(b\) be financial wealth, \(S\) sale proceeds, \(pH\) the new house price,
\(y\) period disposable resources, and \(\phi\) the financed share. Preserve
the inherited transaction balance \(x=b+S-pH\) and conditional budget

\[
c+\text{owner costs}+b'=Rx+y.
\]

Current income can help meet the purchase condition
\(b+S+y/R\ge(1-\phi)pH\). It is not also added to \(x\): that would count
income twice. The discount \(y/R\) follows the existing timing convention in
which \(R\) applies to the transaction balance. This treats period resources
as available within the purchase decision; it does not model annual liquidity
or information within the four-year period.

Two associated numerical paths change together: the forward map uses actual
\(x\), without lifting it to the collateral limit; and the final owner saving
floor remains \(b'\ge-\phi pH\), without carrying the temporary purchase-stage
shortfall into the next period. The latter also changes how off-limit grid
nodes are treated; zero income nests the eligibility formula, not necessarily
every interpolated policy from the original solver.

An independent occupied-state audit rejects transaction clipping, purchase
eligibility violations and final mortgage-limit violations. The native dated
budget check remains active. Joint-choice, grant/waiver, PTI and unsecured-credit
variants are outside this adapter's supported scope. Reference and income-only
arms do not install the purchase changes.

## Validation, run contract and remaining work

The four focused test files cover process algebra/admissibility, entry marginal
preservation, native source-patch compilation, tenure eligibility, transaction
wealth, final mortgage limits, stop-on-failure, timeout partial records and
zero-solve preflight. The reference's native preflight checks all 641 source
files, objective/input hashes and the unchanged 17-figure contract.
[Verification receipt](verification.json).

[Prepared plan](plan.json): four full objective evaluations, at most 32 nested
stationary solves (24 expected), sequential cases with a 110-minute controller
cap. Each native case has 30 minutes and its wrapper 35 minutes. Expected
cluster runtime is roughly 40–70 minutes, extrapolated from previous native
timings and the larger income grid; it has not been measured on these arms.
The first failed case stops the loop. No silent retry, target change, search,
floor relaxation or policy run is authorized by this packet.

Before launch: resolve the income approximation; re-pin files after any edits;
stage an independent Torch bundle; pass all three native zero-solve preflights;
then submit this exact smoke loop. The local plan is deliberately blocked at
the income decision and is not a submitted batch.

After successful completion, collect all 13 target rows and 17 parameter rows,
the standard 17 figures, native exact repetitions and the additional accounting
receipts. Overlay actual beta upper bound 0.99. Native score fingerprints pin
the frozen base; the additional runtime contract and generated source diff
identify the changed economic implementation. No adoption or fit conclusion
is justified before those results are reviewed.

Active tools: `code/model/tools/build_literature_period_income.py`,
`e5f_earnings_wealth_contract.py`, `run_e5f_earnings_wealth_candidate.py`, and
`run_e5f_earnings_wealth_smoke.py`. Tests use the corresponding `test_*.py` files.
