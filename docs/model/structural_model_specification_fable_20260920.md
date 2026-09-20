---
Provenance: Claude Max / fable (claude-fable-5-1); session fbfab0f8-fd81-4d5e-9453-23b6912efaf9; completed 2026-09-20T15:24:12.385254+00:00; reviewed-not-adopted.
---

# Provisional baseline recommendation

**Reading limits.** The memory files could not be read: the file tools reported the memory path as outside the restricted working directory, and no memory folder matched inside the repository, so I started from CALIBRATION_STATUS.md as CLAUDE.md prescribes. All listed sources were read, plus the local copy of Boar, Gorea and Midrigan NBER w23345, the earnings-candidate receipt, the retained fit and parameter tables, the sandbox BGM spec, and the debt-cap builder in the active development code. No jobs, edits, or subagents were used. Everything below is a recommendation; nothing is adopted.

**The recommendation.** Keep the paper's core: a four-year lifecycle household that chooses rent or own, rooms, saving and whether to try for a child, with a down payment at purchase, a parents' space floor, taste shocks on fertility, a warm-glow bequest, pay-as-you-go pensions and a static housing supply curve in the stationary state. Change four blocks that the evidence shows are internally inconsistent or unidentified: tenure access, the child cost device, the dependency process, and the mortgage contract. Adopt the persistent-plus-transitory earnings form conditionally, re-estimated on the project's own sample. Defer the estate receiver, stock-flow supply, the demographic closure, and every extra margin.

The economic logic is that the paper's claim concerns the price and financing of family-sized space. That claim is testable only if the baseline puts prospective parents in the right states, gives children a cost that does not fall with the number of children by construction, and prices the tenure margin instead of walling it. The current fit fails the first condition under both earnings processes, and both failures point at the tenure block.

| Row | Target | Original benchmark | New-income refit |
|---|---:|---:|---:|
| Ownership, heads 30–55 | 0.648 | 0.540 | 0.484 |
| Mean rooms, capped at 9 | 5.56 | 6.42 | 6.68 |
| Wealth / earnings | 6.15 | 4.89 | 6.72 |
| Old p90/p50, 76–84 | 3.52 | 4.50 | 4.28 |
| Mean first-birth age | 25.98 | 26.10 | 26.70 |
| Childless 40–44 | 0.198 | 0.196 | 0.225 |
| Loss, identical contract | | 179.3 | 353.7 |

Too few owners and too much space at the same time is what a six-room rental cap at one price per room produces: renters occupy family-sized space cheaply, so ownership must come from a taste premium and rooms overshoot. Both β and the first-child requirement sit at their upper bounds in both fits, and in the refit β is at its cap while wealth is above target, so β is being used for a timing row rather than identified by the wealth level. The size-dependent rental wedge nests the cap, has data of its own, and was the largest single sandbox improvement. That is why I put the tenure block first rather than the earnings process.

**Earnings, assessed against the primary source.** The referenced paper is Boar, Gorea and Midrigan, "Liquidity Constraints in the U.S. Housing Market," NBER w23345, published in the Review of Economic Studies 2022. On page 8 income is

$$y_{i,t}=\lambda_t\, z_{i,t}\, e_{i,t},\qquad \log z_{i,t+1}=\rho_z \log z_{i,t}+\sigma_z \varepsilon_{i,t+1},$$

with $e$ iid, both innovations standard normal, and no permanent type. The period is a quarter, page 17. Income is post-tax-and-transfer disposable income, OECD-equivalized, from the 1999–2007 PSID, with the transitory volatility scaled down by 55 percent for measurement error, page 18. Table 2, Panel B, reports annualized values.

| Object | BGM Table 2 | Native candidate |
|---|---|---|
| Income concept | post-tax-and-transfer, equivalized | gross head plus spouse earnings, flat payroll tax in model |
| Persistence, annual | 0.964 | 0.970 |
| Persistent stationary log variance | about 0.32 | 0.69 |
| Transitory log variance | 0.107 after measurement-error scaling | 0.344, unadjusted |
| Earnings-moment objective | not comparable | 39.5 versus 14.3 with a fixed effect |

So the candidate shares BGM's functional form but not its numbers or income concept, and the sandbox "BGM" spec uses BGM's numbers without the transitory shock. Neither is equivalent to BGM, whose process sits inside a model with home production, liquid borrowing and continuous rentals. Tommaso's preference is a defensible functional-form choice with strong precedent, including Kaplan, Mitman and Violante and Sommer and Sullivan. What it cannot be is a preference for whichever numbers produce a mechanism. Three conditions make it adoptable: re-estimate on the project PSID sample under one income definition; adjust the transitory variance for measurement error before aggregation, since the unadjusted variance triples BGM's and directly inflates precautionary saving, which is where the credit results live; and aggregate to four years exactly, by simulating the annual process and discretizing the period average, rather than the current coarse mapping. Report the long-lag autocovariance misfit openly. The old argument for types, the old-age wealth tail, no longer binds: the candidate over-produces the tail. The remaining argument is the autocovariance floor at long lags, which a single AR(1) cannot match. Keep the fixed-effect version as the designated robustness process.

# Decision table and the three decisions to take first

| Topic | Recommendation | Reason | Evidence; remaining uncertainty | State-space cost | Identifying data or restriction | Blocks freeze? Author choice? |
|---|---|---|---|---|---|---|
| E1 earnings | Change, conditional: age profile plus persistent AR(1) plus iid transitory, no type | Standard form; types over-produce dispersion in the candidate | BGM eq. p. 8 verified; candidate fit 39.5 vs 14.3; measurement error unadjusted | Same 15 states as now | External estimation on project PSID with error adjustment; age profile from same sample | Blocks. Author choice |
| P1/P3 child term | Change: keep linear $\psi m$, add external child earnings penalty $\tau_c(m)$ on working-age earnings; keep $\xi$ free | Marginal child cost currently falls in $m$; no income gradient | Sandbox penalty halves loss but 0.2 is too large; needs $\psi$ re-normalized | None | $\tau_c$ from child-penalty estimates times female earnings share; $\xi$ from childlessness | Blocks. Author choice |
| P1 concave benefit, SSK weighting | Defer | No external restriction; overlaps $\kappa_C$; SSK halves fertility | Sandbox S1, S2 | None | None available | No |
| F3/F4 dependents | Change: parent-age departure hazard with newborn exemption | Constant hazard gives 40% too few children at home at 30–42 and adult "dependents" after 58; orphans vanish by construction | ACS profile matched at 22–34; timing for early births imperfect | No new state; about +25% Bellman cost in fertile ages | Schedule set externally to ACS children at home by parent age | Blocks. Sign-off only |
| F3 child-stage state | Defer | Adds state; reopens closed decision | None | High | Same ACS profile | No |
| H1/H2 mortgage | Change: origination-only collateral test, amortization, no free cash-out; taper removed | Free re-borrowing makes the down payment a one-time fee and mutes cash-flow channels; KMV and Coven have amortization | Sandbox: ownership −0.11 at fixed parameters, fertility flat | None; debt already in signed assets | Amortization share from a 30-year contract | Blocks. Author choice |
| H1 unsecured credit | Retain $\lambda_d=0$ as an explicit restriction; credit stays a mechanism test | No identifying moment; cohort response unexplained and sign-unstable | Table below | None | KMV precedent: renters cannot borrow | No |
| H3/H6 tenure access | Change: wedge $r(h^R)=r+w_0+w_1\max\{0,h^R-\bar h_R\}$, cap removed, $\chi=1$ | Hard zero is false in ACS; cap plus premium is unidentified beyond their sum; largest fit lever | Sandbox loss 2607→402 at untuned intercept; ownership overshoots, so $w_0$ must be calibrated | Renter size grid extends to owner sizes; no new state | $w_0$ from ownership rate; $w_1$ from renter share by unit size, ACS | Blocks. Author choice |
| F2 tenure shock | Change: $\kappa_H=0$ | No moment; inert in levels; flips policy sign across vintages | Sandbox; transition kink handling untested | None | Deterministic choice, Greaney et al. | Blocks. Sign-off |
| P5 estate valuation | Change: net of selling cost | Code values the house gross; dying in it is cheaper than selling | Verified in review | None | Consistency | Blocks. Sign-off |
| P5 receiver, downsizing | Defer | Welfare accounting, not calibration | Small sandbox effect | Outer fixed point | Needed only for welfare claims | No |
| D1 accounting | Write now: implicit-couple unit, entrant rule, person identities | Prerequisite for population claims | Both reviews agree | None | Documentation | Freeze needs the note, not new structure |
| D2 closure | Defer to transition; choose inflow versus balanced contraction then | Blocks transition, not stationary calibration | No positive stationary population below 2.1 | None now | Terminal condition | No |
| D4 geography | Change to national; remeasure four ACS rows | Model has no location; fertility, wealth and earnings are national | Values after remeasurement unknown | None | Same ACS filters, national scope | Blocks. Author choice |
| H4 supply | Retain static curve for stationary state; stock-flow for transition | In steady state $H=I_0P^{\eta}/\delta$, so $H_0\equiv I_0/\delta$: calibration unchanged | Two elasticities coexist in code | None now | Depreciation, construction elasticity | No |
| Unintended births, formation, widowhood, health | Defer | No targets; none needed for the mechanism | Scaling of price elasticity by intended share | High | None | No |

The three decisions to discuss first, ranked:

1. **Tenure access and the tenure shock, H3/H6/F2.** This is the load-bearing block for both fit and mechanism, because the mortgage effect nearly disappears when the rental cap rises. Smallest evidence needed: none beyond the ACS renter share by rooms already in the tenure memo. The decision is a modeling stance: whether the paper's mechanism may rest on a size-dependent rental price rather than a hard availability wall.
2. **Child cost and earnings together, P1/P3/E1.** They meet at the income gradient of fertility, which the model cannot currently produce. Smallest evidence: the CPS children-ever-born gradient by household income from the existing supplement extract, and a two-column autocovariance table for the re-estimated process with and without a fixed effect, short and long lags.
3. **Geographic scope, D4.** It sets the target values, so it must precede the freeze. Smallest evidence: rerun the existing ACS builder nationally under the same filters and show the four rows side by side with the 42-metro values.

Maturation, net estate valuation and $\kappa_H=0$ need only sign-off. One prerequisite is not a decision but a fact: the overnight refit ran on the frozen September-14 source, which contains none of the switches. A baseline with any of these changes needs a new pinned source before calibration.

# Calibration contract, what survives, and the freeze

**Calibration contract.** Keep all thirteen rows and their weights. Add one row, the renter share of units with seven or more rooms, to identify $w_1$; $w_0$ replaces $\chi$ on the ownership row. Everything else changed is externally restricted, not estimated: the earnings process and age profile, $\tau_c$, the departure schedule, the amortization share, $\lambda_d=0$, $\kappa_H=0$, net estate valuation. The free vector stays at nine structural parameters plus the two wedge terms, with $\psi$ normalized separately, against thirteen scored moments, so identification is preserved. Rules: β is identified by the wealth level only, and a refit that leaves it at the cap is reported as misfit, not as an estimate; $h_P$ is identified by the first-birth rooms response, and a cap there means the rooms and recent-parent rows are in conflict; no row is dropped or reweighted without a named replacement. Validation profiles, not targeted: ownership by age, rooms by age and by number of children, children at home by parent age, wealth by age, first-birth hazard by age and wealth tercile, renter share by size, the marginal constrained share near 4 to 7 percent, and the fertility gradient by income.

Measurement reconciliations before the refit: the geographic scope of the four ACS rows; the first-birth-age convention, since the calibration observer and the cohort report differ by about two years on the same solution; the 2.1 normalization, which is a stationary construct and not the explicit cohort births near 1.86; the recent-parent group proxy, which still has no child ages under the new departure law; and one income definition across the earnings estimation, the wealth-to-earnings row and the model's after-payroll-tax income.

**Mechanism conclusions.**

| Family | Mortgage φ=1, cap 6: cohort births, ownership pp | Credit one annual earnings: flow, cohort | Credit four annual earnings: flow, cohort |
|---|---|---|---|
| Original | +0.23%, +10.6 | +3.22%, +1.11% | +6.89%, +2.92% |
| New-income pilot | +0.09%, +5.0 | +0.32%, −2.15% | +0.13%, −4.04% |
| New-income refit | +0.08%, +1.3 | +2.46%, −0.77% | +2.38%, −2.42% |

Source-code facts: mortgage relaxation moves deposit and collateral together; in the active code the unsecured cap is $\lambda$ times age-specific mean after-tax period income with a default taper from 42 to 62, in parameters.py, and whether the frozen source used that taper needs confirming; the estate is valued gross of selling cost; the value function falls when a child stays one more period at every state at the retained parameters. Mathematical implications: snapshot flow and lifetime cohort births are different objects by construction, so opposite signs are possible; the marginal child cost falls in $m$ under the deflating scale at $\sigma=2$ regardless of calibration; the stationary state is invariant to static versus stock-flow supply. Outcomes at particular calibrations: small mortgage effects in all three families, the rental-cap interaction, the 4.5 and 2.8 percent constrained shares, the switch table. Conjectures: the cohort sign reversal is debt overhang or a higher continuation value of childlessness when children are net costs; the plateau at one and four annual earnings is the natural borrowing limit or the asset-grid floor binding under the riskier process. The cheapest check for both is the age-by-age debt distribution in the saved arms, with no new solve.

Reassessing the old claims. "Children are net costs" is a value-function fact at one calibration, produced by the functional form and the $\psi$ normalization; its durable lesson is that every specification comparison must re-normalize $\psi$. "Unsecured credit raises fertility" does not survive as a general statement: it holds for snapshot flows, and for cohort births only under the original process. "Frictionless benchmark" should be retired; the arms relax financing while keeping the cap, the selling cost and one price per room. "Easier mortgages raise ownership, not births" survives as the best current evidence in all three families, but the ownership response shrinks from 10.6 to 1.3 points across families, and the joint refit with the wedge and the child cost remains untested, so the verdict on the deposit mechanism stays conditional. The rental-access interaction is the one result robust across families and it supports making the wedge the baseline object.

**Freeze point and order of work.** The specification is frozen when Tommaso decides the six author-choice rows and signs off the three sign-off rows, and the target contract with its scope is fingerprinted. Then, in order: implementation and measurement only, meaning the re-estimated earnings process, the national ACS rows and the new renter-share row, the external restrictions, the wedge wired into a newly pinned native source with bitwise-off tests, the accounting note and the first-birth-age reconciliation; joint calibration on the cluster with full fit and parameter tables; lifecycle validation against the untargeted profiles; and only then the same 48-cell mechanism grid on the refit with both outcome clocks and debt-by-age diagnostics, followed by the property-tax experiment under an explicit closure.

Unresolved after this pass: the cause of the credit cohort reversal and the plateau; the frozen source's debt taper; the national target values; whether the transition solver handles deterministic tenure kinks; the two unverified Kaplan, Mitman and Violante claims in the tenure memo; and the memory files I could not read.

Sources: [NBER w23345](https://www.nber.org/papers/w23345), [Review of Economic Studies version](https://academic.oup.com/restud/article-abstract/89/3/1120/6372706).


---
## Focused clarification following lead review

Provenance: Claude Max / fable (claude-fable-5-1); resumed session fbfab0f8-fd81-4d5e-9453-23b6912efaf9; completed 2026-09-20T15:31:41.860337+00:00; reviewed-not-adopted.

**Corrections**

1. **Earnings evidence, separated.** The original process pairs a Floden–Lindé AR(1) in after-tax units with a PSID gross-earnings fixed effect: a coherence problem, not proven double counting. The no-type candidate over-produces the old-age tail and the wealth level because its persistent variance is large, so "types are needed for the tail" no longer binds; that is a consequence of inflated variance, not a validation. The process is identified outside the model by the 13-lag autocovariances: the fixed effect absorbs the long-lag floor, and the no-type fit buys that floor with persistence near 0.97 and more variance, hence 39.5 against 14.3. My cell "types over-produce dispersion in the candidate" was wrong.

2. **Counts.** Now: nine structural coordinates plus normalized ψ, against twelve scored rows plus one exact normalization. Wedge variant: nine minus χ plus two wedge terms gives ten structural plus ψ, eleven free objects, against twelve scored rows plus normalization, thirteen if a renter-share row were later scored. The order condition holds; identification does not follow, since χ and the intercept enter ownership as near substitutes and ψ, ξ, both κ and β all move timing and childlessness. I withdraw "β identified only by wealth" and the inference that the bound on the first-child requirement proves two rooms rows conflict. Bound hits are search or misfit diagnostics; a sensitivity matrix at the selected point is the check.

3. **Nesting.** The cap is the infinite-slope limit of the wedge only with zero intercept and χ retained; at finite slope the renter choice set changes, so nesting is a limit, not an equality on the grid. The sandbox wedge run changed four things at once: cap removed, χ set to one, intercept 0.02, slope 0.05, so its loss of 402 does not isolate the slope. Low ownership with high rooms is consistent with the wedge hypothesis but also with preferences, the equilibrium price, or the state distribution; the scarce-space test speaks against the price-level version only.

4. **Not sign-offs.** Net estate valuation presumes liquidation at death; with no receiver it changes the bequest argument and the incentive to die in the house, so it needs a contract line. Deterministic tenure makes aggregate demand a step function on the finite grid, which is why the kink-state tools exist; clearing then needs a distribution-level rule. Both are author choices with a contract, and neither enters the variants below.

5. **Doses and curvature.** The plateau is between λ of 1 and 5, four and twenty annual after-tax earnings amounts; I wrote one and four. Curvature: the active solver at line 2456 applies $e(m)^{\sigma-1}$ with $e(m)=((2+0.7m)/2)^{0.7}$ to the CRRA composite, so at $\sigma=2$ flow utility is $-e(m)/X+\psi m$ and the scale increments at fixed $X$ are 0.234, 0.216, 0.203. Verified in active code; the frozen source line was not re-read. It is a fixed-composite statement, not the full cost including the floor and ξ. The coarse four-year mapping stays in both variants; a period-average process is not first-order Markov in the sampled state, so exactness would need an extra state or a fitted-chain approximation.

**Prerequisite, zero evaluations.** Read the saved refit arms: the lowest liquid asset with positive mass at λ of 1 and 5 against the cap and grid floor; debt and wealth by age and first-birth flow by wealth in the credit arms against baseline. This precedes Variant B's credit arms, not Variant A.

**Two isolated comparisons on one pinned switch-bearing source.** Common controls: nine structural parameters at retained values, ψ re-normalized to 2.1, unchanged target and weight fingerprint, same entry rule, same solver settings, no scored row added. The first stationary evaluation is a switches-off control that must reproduce 179.298.

| | Variant A: wedge slope nested against the cap | Variant B: no-type process, measurement-error adjusted |
|---|---|---|
| Hypothesis | A finite slope raises ownership and lowers rooms at unchanged timing rows; the deposit effect on births shrinks further | Housing misfit and the cohort sign reversal come from inflated transitory risk, not from removing the type |
| Fixed | χ retained, intercept zero, cap removed | Same coarse mapping; entry wealth distribution reported, since it varies with the process |
| Arms | Cap control; slope 0.05, 0.2, 1.0 per room above six | Retained and as-run exist; transitory variance scaled by BGM's 0.45; no-type refit to autocovariances net of that variance |
| Outcomes | Full 13 rows; renter share above six rooms against ACS; price; then φ=1 plus cohort on the best arm and its control | Full 13 rows; then λ of 0.25 and 1 plus cohort on the adjusted arm |
| Stop | After three slopes; slope-only insufficient if ownership stays below 0.55; cap preferred if renters above six rooms exceed ACS at every slope | Form stays in play if ownership and rooms recover half the gap to 179.298; otherwise retain the fixed effect provisionally; negative cohort sign again means not a transitory-risk artifact |
| Cost | 4 stationary plus 1 repetition; 4 fixed-price | 2 stationary plus 1 repetition; 4 fixed-price |

Total: 9 of 32 stationary, 8 of 16 fixed-price; the remainder is held back.
