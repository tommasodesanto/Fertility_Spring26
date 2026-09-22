# Earnings, entry wealth and purchase timing

## Live review — September 22, 01:54 EDT

The primary-source reviews and direct four-year PSID point estimates are done.
The main candidate uses the earliest complete four-year block grid (1984-based
survey-year labels), with gross RP/spouse labor earnings and age/year fixed
effects. It has **persistent AR(1) plus iid transitory risk, no fixed permanent
type**. Three covariance moments identify three parameters; their exact fit is
not evidence that the specification is universally correct.

| Object | Four-year estimate | Interpretation |
|---|---:|---|
| Persistent coefficient | 0.776144670 | Per model period |
| Persistent innovation SD | 0.437435453 | Log earnings innovation, not stationary SD |
| Transitory SD | 0.164990611 | Period iid log earnings shock |
| Persistent stationary variance | 0.481262675 | Initial stationary distribution is externally imposed |

The other two estimable block alignments give persistence 0.7971 and 0.8135,
with positive transitory variances. A fourth alignment lacks the third moment
and is reported as underidentified. The main sample has 7,812 blocks and 1,449
eight-year pairs. The corrected 199-draw person bootstrap shows that the main
transitory estimate is uncertain, with its lower percentile at zero. The
100-replication synthetic exercise applies the same age/year regression,
missingness, weights and covariance observer; its small mean projection effects
are reported explicitly. [Complete data evidence](data_validation/README.md).

The 15 persistent by 3 iid grid reproduces all five log covariances to rounding;
level covariance errors are at most 4.61%. This is distribution validation,
not household-grid convergence. The analogous annual persistent coefficient is
0.93861 and innovation SD 0.23932 under an endpoint embedding. These conversions
are comparisons only: the period transitory SD cannot be mechanically called an
annual estimate. [Literature](literature/annual_sources_review.md),
[multi-year precedents](literature/multiyear_review.md),
[timing](literature/timing_review.md), [lead receipt](lead_review.json).

**Native testing found a specific issue.** Income-only and income-plus-purchase
both stop at the unchanged feasibility gate with about 3.16e-12 dead mass
(tolerance 1e-12). The difference is only 2.5e-18, isolating this failure from the
purchase adapter. The evidence points to inherited negative entrant debt under
new low-income support. That entry distribution was measured for ages 25–35 but
imposed at age 18. A distinct zero-assets-at-18 test, supported as an external
entry convention by Sommer and De Nardi, completed two stationary solves (246–255 seconds each) and was deliberately stopped before its full objective to set an adequate smoke budget. It is not yet adopted.
[Failure diagnosis](staging/feasibility_diagnosis.md).

**V2 smoke failure and V3 correction:** V2 completed six stationary equilibrium
solves, then failed the independent purchase-accounting gate before any scored
objective or complete repetition. Out-of-grid transaction mass was
9.173572702e-5; maximum wealth-map error was6.721741917. Threshold and final
mortgage-floor violations were zero. The wealth error equals net proceeds from
selling the largest home at the final price, consistent with upper-grid sales;
no solved policy arrays survived, so this state attribution is an inference.
The backward tenure kernels still clipped transactions despite the documented
support contract. V3 rejects interpolation outside either endpoint in both
compiled kernels and preserves native interpolation exactly within support.
Boundary tests check upper sales, lower purchases and exact endpoints, with both
argmax and logit choices. Native Numba compilation passes. This implements the
stated support contract; no accounting gate, income input, target or bound was
relaxed. Finite-grid robustness remains outstanding. The failed V2 evidence is
preserved in [reviewed correction](staging/purchase_support_correction.json).

**V3 failure and V4 numerical-domain correction:** V3 passed purchase accounting
(zero unsupported transaction mass, maximum map error3.55e-15), then failed the
occupied-value check. All25 drops are in the final wealth interval27.807→30,
owner product5 and income state44; affected lower-node mass9.18e-8 and maximum
value drop6.64e-4. It completed six stationary solves, zero scored objectives.
Rejecting transactions at the upper edge can remove a sale/downsizing option
when wealth rises. The regression test reproduces that loss and removes it by
extending support. V4 preserves all120 original wealth knots and appends40
geometric knots from30 to3000. No economic saving cap is added. The generated
probe changes only grid assignment and expected geometry to160 nodes; source,
target, normalization, accounting and value-quality gates remain intact.
The upper endpoint is a declared numerical choice, not a proven reachability
bound or convergence result.42 checks and a native zero-solve full-setup test
passed, including all641 source files and exact entry/grid agreement.
[Failure, correction and budgets](staging/wealth_grid_correction.json),
[full native setup](staging/v4_native_geometry_preflight.json).

Torch authentication has expired; no cluster job has been submitted. Under the
author's explicit open-laptop overnight authorization, the corrected frozen V4 bundle is
**running native smoke locally**: two exact anchor repetitions and one inward-beta
probe. Its successful receipt alone triggers the longer search. V1 stopped
before solving on a controller metadata mismatch; the correction and old receipt
are preserved in [startup repair](staging/startup_contract_repair.json). All
three V4 zero-solve native preflights passed, 42 unit checks passed, and native
progress is confirmed. [Execution/monitor registration](local_execution.json),
[immutable V4 plan](staging/frozen_v4_plan.json),
[source inventory](staging/frozen_v4_hash_manifest.json).

The longer search has four single-threaded workers, a fixed deterministic seed,
at most 64 joint proposals over the nine existing coordinates, and two exact
selected repetitions. The production controller has 18,800 seconds including 6,400
reserved for verification; the smoke supervisor has 10,000 seconds and the full
local supervisor 32,000 seconds (about 8.9 hours) including overhead. At most 69
full objectives /552 nested solves are authorized; time limits will likely bind
first. Estimated full objectives take 1,900–2,800 seconds from measured stationary
solve times scaled for the larger grid. Per-objective outer cap is3,200 seconds;
two-repetition cap is6,400 seconds. Every case records progress, and latest/best results are saved.
Unexpected failures stop the stage; no automatic retry or gate relaxation.
The follow-up checks every 20 minutes and stays quiet while healthy.

No full objective has completed yet. The final review requires complete 13-target
and 17-parameter tables, actual bounds, exact selected repeats, and 17 standard
plots, followed by the [parameter and literature audit](literature/parameter_validation_rubric.md).
The September 14 reference, target system and numerical gates remain unchanged.
The four-year information approximation, weak iid estimate, and stationary
initial income at18 remain explicit limitations even if the search succeeds.

## Authorized sequence recorded September 21 late night (status superseded above)

The author authorizes reviews, own-data estimation and/or transformation of
published estimates, implementation, tested longer calibration, and a final
parameter/plausibility review. The older endpoint-only run plan below is held;
the live September22 section above records execution. The older launch-status notes below are historical.

1. Verify annual and multi-year primary sources, actual gross-household earnings
   inputs, and purchase timing. Evidence is saved under [literature](literature/).
2. Construct observed complete four-year PSID gross-earnings blocks (no filled
   biennial years); report support, covariance fit and sensitivity to block
   alignment before choosing direct period estimates or transformed literature
   inputs. A binding zero transitory variance is a finding, not permission to
   quietly drop the intended component.
3. Freeze a single main specification, including its income concept, period
   information set, entry law, grid and purchase accounting. State approximation
   errors and any empirically unmeasured restriction. No preference parameter may
   compensate for a silent income-definition or target change.
4. Verify the income grid and budget math; run the exact calibration loop as a
   smoke with repeated anchor and a nearby proposal. Measure runtime and memory.
5. Run a finite, seeded joint calibration in the nine existing coordinates on
   Torch with the complete unchanged target system, actual bounds, progress,
   checkpoints, separate verification budget, and stop-on-contract-failure.
   Re-evaluate the selected point twice and retain all 17 standard diagnostics.
6. Review all targets and parameters, grid sensitivity, boundary behavior and
   literature comparability. Deliver a supported candidate or a precise failure;
   do not call a bounded search convergence, or equate a successful solve with a
   publishable calibration.

Planning budget: up to 35 minutes for direct data validation, 45 minutes for
controller preparation in parallel, followed by exact-loop smoke and a bounded
four-hour joint search if the specification and numerical gates pass. Actual
solve counts and allocations are set from the smoke, not promised in advance.

## Historical decisions and superseded preparation

**Launch on hold after author clarification:** acceleration did not authorize
selection of an arbitrary period approximation. Worker staging stopped; no
lead submission. The plan gate is pending again. The proposed overnight
specification below is retained as evidence, not an approved run.

**Live decision: retain persistent plus transitory income. The single-component
proposal below is withdrawn. A bounded three-arm overnight diagnostic is being
staged on Torch; no final earnings specification or calibration is adopted.**

## Live overnight specification

The classic decomposition is an age profile, a persistent AR(1), and independent
iid risk; see [Storesletten, Telmer and Yaron (2004), JME](https://doi.org/10.1016/j.jmoneco.2003.06.005).
Their richer specification also permits fixed heterogeneity. For the no-fixed-type
fertility application we use [Sommer (2016), JME, Section 4.1/Table 2](https://www.kamilasommer.net/Fertility.pdf)
as the annual parameter source: rho 0.95, persistent innovation SD 0.21, iid SD
0.17. The classic reference does not certify our four-year approximation.

The **diagnostic** four-year conversion retains independent persistent and iid
components. Persistent endpoint persistence is 0.95^4 = 0.81450625 and its
innovation SD is 0.21 sqrt(1 + 0.95^2 + 0.95^4 + 0.95^6) = 0.390176278.
The temporary log variance is log(1 + [exp(0.17^2)-1]/4), giving SD 0.085461555.
This matches the level variance of the average of four independent mean-one
lognormal temporary factors, approximating their average as lognormal. Combining
it with the persistent endpoint is an approximation to total period earnings,
not an exact annual-to-four-year aggregation. Its continuous variance exceeds the
exact four-year-average variance by 7.8387%; lag 1–4 covariances differ by <0.6%.

The model uses 15 persistent by 3 iid nodes, independent iid transitions, mean-one
income risk, and the retained age profile, taxes and period units. Stationary
entry risk and wage-to-exogenous-household-earnings mapping are diagnostic
assumptions. Entry wealth marginal is preserved under the documented rank
coupling. No permanent type or single-component collapse is introduced.

Tonight compares the reference once, new earnings once, and new earnings plus
current-income purchase accounting twice, with original structural parameters.
This is four full objectives, at most 32 nested stationary solves (24 expected),
110-minute controller cap, and stop at the first failure. It is the numerical
smoke/diagnostic itself; no dependent search or policy run is scheduled. Each
successful arm must retain full targets/parameters and the standard 17 plots.
The constructor, accounting and controller code/pins are unchanged from the
previous 27 passing focused tests. All three local zero-solve native preflights
now pass; Torch relocation/preflight and scheduler receipt are still required.

## Withdrawn single-component proposal — September 21

Use Sommer's annual risk inputs (0.95, 0.21, 0.17), aggregate annual earnings
**levels** over four years, and approximate normalized block income with one
lognormal AR(1). An equal-weight fit of relative errors in stationary level
variance and four autocovariances, allowing a nonnegative independent iid
variance, gives period persistence **0.823078407** and innovation SD
**0.376839547**. The fitted independent iid variance is effectively zero
(8.4e-20, an active lower bound). This was a proposed one-state
approximation; it was withdrawn after the author reaffirmed the intended
persistent-plus-transitory structure. It is not the overnight specification.
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
responses are not yet validated. The single-component fit has not been installed in the adapter or run plan;
no model specification has been adopted.

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

Before launch: retain the live diagnostic approximation above; re-pin files after any edits;
stage an independent Torch bundle; pass all three native zero-solve preflights;
then submit this exact smoke loop. The local plan authorizes this bounded diagnostic and is not itself a scheduler
submission. Record the scheduler receipt before claiming launch.

After successful completion, collect all 13 target rows and 17 parameter rows,
the standard 17 figures, native exact repetitions and the additional accounting
receipts. Overlay actual beta upper bound 0.99. Native score fingerprints pin
the frozen base; the additional runtime contract and generated source diff
identify the changed economic implementation. No adoption or fit conclusion
is justified before those results are reviewed.

Active tools: `code/model/tools/build_literature_period_income.py`,
`e5f_earnings_wealth_contract.py`, `run_e5f_earnings_wealth_candidate.py`, and
`run_e5f_earnings_wealth_smoke.py`. Tests use the corresponding `test_*.py` files.
