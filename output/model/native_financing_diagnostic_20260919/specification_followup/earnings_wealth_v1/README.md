# Earnings, entry wealth and purchase timing

**September 21 evening: implementation and cheap checks completed; the author
requests a broader, top-five/top-field literature grounding and improved
four-year aggregation. No new calibration result or submitted job.**

## Concrete four-year proposal — September 21

Use Sommer's annual risk inputs (0.95, 0.21, 0.17), aggregate annual earnings
**levels** over four years, and approximate normalized block income with one
lognormal AR(1). An equal-weight fit of relative errors in stationary level
variance and four autocovariances, allowing a nonnegative independent iid
variance, gives period persistence **0.823078407** and innovation SD
**0.376839547**. The fitted independent iid variance is effectively zero
(8.4e-20, an active lower bound). Thus this is a recommendation for one period
income state, revising the earlier proposed separate four-year iid component.
Annual iid risk is still included in the aggregated target distribution.

The stationary log variance is 0.44027777; normalize exp(z) by its mean.
Continuous level-moment errors are +2.23%, -2.48%, -1.18%, +0.04%, +1.22%
for variance and lags 1–4 respectively. This is an approximation selected by
a declared income-moment criterion, not an exact temporal aggregation or
empirical parameter estimate. The annual covariance formula was independently
recomputed and three optimizer starts agree in objective within 1e-10.
[Complete calculation and restrictions](period_proxy_fit.json).

This calculation assumes stationary annual risk and excludes the deterministic
age profile. It does not copy Sommer's zero persistent state at entry. Using
her wage process for exogenous household earnings remains a diagnostic proxy.
Low-income tails, finite-grid accuracy, lifecycle/entry assumptions and household
responses are not yet validated. The existing adapter and run plan have **not**
been changed or launched; no model specification has been adopted.

## Literature review and revised recommendation

Prioritize journal quality and the relevant economic object over finding exactly
four-year periods. The following are verified precedents, not interchangeable
parameter estimates:

| Reference | Period and verified implementation | Use here |
|---|---|---|
| [De Nardi (2004), Review of Economic Studies](https://users.nber.org/~denardim/research/denardi.pdf), pp. 747, 752, 765–766 | Five-year OLG model; estimates the earnings process after aggregating PSID earnings into complete five-year cells. | Main precedent for matching income measurement to model frequency. Published Table 3 uses persistence 0.85 and variance 0.30; do not substitute values attributed to this source by later implementations. |
| [Bick (2016), JEEA; inspected working-paper version](https://mpra.ub.uni-muenchen.de/41757/1/MPRA_paper_41757.pdf), printed pp. 10, 44–45, Appendix C.1 | Three-year fertility/labor-supply model; constructs period income by summing monthly allocations of annual gross income, then estimates the process. Reports period persistence 0.882 and innovation SD 0.272. | Closely relevant family-model measurement method; these German spousal-income estimates are not US household parameters. |
| [Sommer (2016), JME](https://www.kamilasommer.net/Fertility.pdf), Section 4.1/Table 2 | Annual fertility model with persistent AR(1) plus iid wage risk: 0.95, 0.21, 0.17. | US fertility/risk architecture and candidate annual inputs; not a four-year conversion recipe. |
| [Doepke and Kindermann (2019), AER](https://faculty.wcas.northwestern.edu/mdo738/research/Doepke_Kindermann_AER_2019.pdf) | Three-year fertility decisions. | Supports multi-year family-model timing; does not supply our stochastic-income aggregation. |

A secondary exact-four-year example is [Kolasa (2024), JEDC; inspected Warsaw
working paper](https://www.wne.uw.edu.pl/download_file/4023/0), pp. 12, 20, 22.
It combines fertility timing with persistent and transitory earnings risk, but
Table 1 explicitly labels 0.9, sqrt(0.03), and 0.25 as **annual values**. The
reviewed text does not establish its complete four-year shock conversion.
It is not the main source under the author's journal preference, nor a ready-made
four-year parameter vector. No exact-four-year stochastic-earnings implementation
meeting the requested journal priority was verified in this bounded search.

**Lead recommendation:** retain a simple persistent-risk architecture without a
fixed permanent type, but construct the four-year earnings object before fitting
its finite process. Use four-year totals (or means with the explicit factor four)
consistently in data and model. First aggregate a sourced annual process; fit a
nonnegative, parsimonious period approximation and disclose errors in variance,
serial covariance and low-income tails. Direct estimation on matching four-year
data is the alternative when an adequate annual source cannot represent the
intended household-income concept. The simple process is an approximation, not
an exact decomposition of aggregated annual shocks. Do not require exact matching
of three covariances when that implies a negative variance, or add a permanent
type merely to avoid that algebraic rejection.

This supersedes the earlier request for an immediate yes/no adoption of the
endpoint approximation. The calculations below remain useful diagnostic evidence:
their 7.84% excess variance is not a demonstrated household-fit failure. The concrete continuous fit above is now available; tail, grid, lifecycle/entry
validation and household runs remain unrun.

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
distribution of four years of total earnings. **Retained diagnostic, not adopted;
aggregation is under revision as described above.**

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
