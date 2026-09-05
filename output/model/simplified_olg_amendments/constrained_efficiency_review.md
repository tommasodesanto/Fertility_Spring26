# Cash transfers followed by housing-market clearing

September 5, 2026. Bounded independent review of the constrained-planner benchmark. No manuscript, model, calibration, or policy file was changed.

**Lead scope clarification:** references below to a conditional Pareto
improvement cover the four specified household groups only. Initial intermediary
title owners and future market clearing are omitted, so this is not a proved
Pareto improvement for the complete economy. The lead independently reproduced
the finite construction at three price changes and the obstruction's budgets
and coefficient signs; see `planner_benchmark_checks.json`.

**Main finding.** The financing wedge does not establish constrained Pareto inefficiency. With the proposed cash-only instrument, I can prove an improvement for the four current household groups in the saved illustration **conditional on fixed continuation prices and rebates**. It taxes current young households and compensates them exactly through lower housing costs; current old households gain. Housing moves **toward the old**. Initial intermediary title owners are not covered, so this is not an established Pareto improvement for a complete economy. I also construct a conditional economy satisfying the same current household equations, with a strictly binding young down-payment constraint, in which **no nearby balanced cash-transfer improvement for all four household groups exists**. Neither statement is a theorem about a complete OLG equilibrium path.

| Result | Status |
|---|---|
| Cash-only improvement for four household groups in the saved illustration, with continuation prices fixed and initial intermediary title owners outside the comparison | PROVED UNDER STATED ASSUMPTIONS; finite household/market calculation checked; not a complete-economy Pareto result |
| Any nearby conditional Pareto improvement in that illustration increases old housing | PROVED UNDER STATED ASSUMPTIONS |
| Explicit nearby-inefficiency obstruction with positive property tax and a binding young borrowing constraint | PROVED UNDER STATED ASSUMPTIONS in a conditional economy |
| Pareto improvement for every affected household along a genuine OLG path | NOT ESTABLISHED |
| Transfer eligibility, information, coverage of initial title owners, and future fiscal powers | NEEDS AUTHOR CHOICE |

## 1. The instrument and the welfare comparison

Let \(\ell_{it}\) be a household's net transfer in date-\(t\) goods, delivered **before housing closing**. A negative value is an immediately collected cash tax. The government neither borrows nor imposes future repayment. Transfers are lump sums conditional on predetermined household characteristics and, if permitted, its observed ownership taste; they are not payments conditional on subsequently buying housing. The mortgage share \(\phi_t\) stays fixed. Taxes must be payable from the donor's available cash. The examples below tax only young households with ample positive cash.

The young renter's budget becomes

\[
c+a'+q r_t h=y_i+b_i+T_t+\ell_{it}.
\]

The young owner's original budget and closing constraint become

\[
\begin{aligned}
c+a'+[(1-\phi_t)+q\tau^p]P_th&=y_i+b_i+T_t+\ell_{it},\\
a_{t+1}&=R_fa'-R_f\phi_tP_th,\\
(1-\phi_t)P_th&\le b_i+\ell_{it}.
\end{aligned}
\]

The grant enters total resources once. An old household's budget receives \(\ell_{it}\) on its right side. All original housing-size, nonnegative-saving, retained-title, and estate constraints remain. There is still one stock \(\bar H\), with \(h_R^{\max}<h_O^{\max}\); no separate tenure stocks are introduced.

The transfer account and existing property-tax account are separate:

\[
\int_{Y_t}\ell_{it}\,di+\int_{O_t}\ell_{jt}\,dj=0,
\qquad
(Y_t+O_t)T_t=q\tau^pP_t\bar H.
\]

Households then optimize and housing clears. The intermediary condition remains

\[
u_t\equiv q r_t=(1+q\tau^p)P_t-qP_{t+1}.
\]

No real moving cost is present in these calculations. Adding one requires its own resource account and would reduce the set of possible improvements. The source's proposed nonnegative reallocation cost is still provisional.

The first welfare test holds each household's fertility at its baseline value and holds cohort masses fixed. It compares realized lifetime utility, including the ownership taste, for every current young household, and remaining utility for every current old household. A higher average or weighted utility is insufficient. For a path comparison, the same requirement must also cover subsequent cohorts and any other agents whose claims change.

**Later payments are different.** A transfer delivered after closing changes the budget but not \(b_i\) in the down-payment condition. A loan with compulsory future repayment is an additional credit instrument. Neither is silently included in \(\ell\).

## 2. Exact local accounting with continuation prices fixed

This section is a **conditional date-\(t\) economy**: fix \(P_{t+1},P_{t+2},T_{t+1}\), the initial old states, fertility, and masses. Young households still optimize saving and their old-age allocations. Their altered future housing choices are **not** required to clear a future housing market in this calculation.

For compactness define

\[
\delta=1-\phi,\quad a_P=1+q\tau^p,\quad
t_P=\frac{q\tau^p\bar H}{Y+O},\quad A=1+\gamma+\omega_B.
\]

Thus \(dT=t_P\,dP\) and \(du=a_P\,dP\). Here \(a_P\) is a price coefficient, not household saving. Define \(z=(MV^Y-u)/(\delta P)>0\) for a young owner with a strictly binding cash constraint. It is the ratio of the utility multiplier on purchase cash to the utility multiplier on the budget. Positive saving and interior old-owner choices hold in both examples.

Let \(c_i\) denote the derivative of the exact cash transfer that keeps household \(i\)'s utility constant when \(P\) changes. This is a compensation coefficient, not consumption. Let \(d_i=\partial h_i/\partial\ell_i\) at fixed prices and \(k_i\) be its housing derivative along that exact compensation. On the relevant branches:

| Household | \(c_i\) | \(d_i\) | \(k_i\) |
|---|---|---|---|
| Young owner, cash constraint binding | \([a_Ph-t_P+z\delta h]/(1+z)\) | \(1/(\delta P)\) | \(c_i/(\delta P)-h/P\) |
| Young renter, size cap binding | \(a_Ph_R^{\max}-t_P\) | \(0\) | \(0\) |
| Old renter, size cap binding | \(a_Ph_R^{\max}-t_P\) | \(0\) | \(0\) |
| Old owner, retention and estate constraints slack | \(a_Ph^2-H-t_P\) | \(\gamma/(Au)\) | \(-(1-\gamma/A)h^2a_P/u\) |

**Derivation.** For a young owner, let \(\lambda=1/x\). Applying the envelope theorem to the original two budgets, and using the positive-saving Euler condition, gives

\[
\frac{dW^O}{\lambda}
=(1+z)d\ell+[t_P-a_Ph-z\delta h]dP.
\]

The house-price effect includes both the down payment and the mortgage repayment: their budget coefficients add to \(1+q\tau^p\), not \(1-\phi+q\tau^p\). Differentiating the binding closing condition gives

\[
dh=\frac{d\ell}{\delta P}-\frac hP dP.
\]

For an old owner, \(dV^O/(1/c^2)=d\ell+[H+t_P-a_Ph^2]dP\). Its original interior demands are \(h^2=\gamma(a+PH+T+\ell)/(Au)\). These expressions give the table. The capped-renter results follow directly from their original budgets. The ordinary rebate is included in every coefficient.

Let \(m_i\) denote group masses and define

\[
C=\sum_i m_i c_i,\qquad K=\sum_i m_i k_i.
\]

For a utility change, write \(g_i=dU_i/U_{i,\ell}\). A Pareto direction has \(g_i\ge0\). The fiscal and housing equations imply

\[
\sum_i m_i g_i+C\,dP=0,
\qquad
\sum_i m_i d_i g_i+K\,dP=0. \tag{1}
\]

These are accounting restrictions, not a presumed positive theorem. In particular, a positive financing wedge does not determine the signs of \(C\) or \(K\).

## 3. An explicit obstruction with the original positive property tax

**Local impossibility proposition.** Suppose all branches in the table remain valid in a neighborhood, every housing demand is nondecreasing in its own transfer, \(C>0\), and \(K<0\). Then there is no different sufficiently close equilibrium produced by balanced cash transfers that weakly improves every household and strictly improves anyone.

**Proof, including finite changes.** Let \(\ell_i^c(P)\) be the exact transfer maintaining baseline utility, with \(\ell_i^c(P_0)=0\). Monotonicity of utility in cash implies that any Pareto allocation has \(\ell_i\ge\ell_i^c(P)\). Since the sum of compensating transfers has positive derivative \(C\), a balanced Pareto reform cannot have \(P>P_0\) nearby. If \(P<P_0\), its compensated total housing demand exceeds \(\bar H\), because its derivative is \(K<0\). Giving households additional transfers to improve their utilities cannot reduce that demand. Thus housing cannot clear. At \(P=P_0\), balanced transfers and monotonic utility force every transfer to zero. Continuity preserves the strict signs locally. This proves a finite neighborhood result, not merely failure of one differentiable direction.

The following complete current-date example satisfies these hypotheses. Continuation prices and the next rebate are fixed as explicitly stated in Section 2; it is not presented as a stationary OLG equilibrium.

| Primitive/state | Value |
|---|---:|
| \(q,\tau^p,\phi\) | \(0.5,0.05,0.8\) |
| \(P_t,P_{t+1},P_{t+2}\) | \(1,1,1\) |
| \(Y,O\), owner shares in both cohorts | \(1,1\), both \(0.5\) |
| \(\bar H,h_R^{\max},h_O^{\max}\) | \(1.4,0.5,2\) |
| \(\alpha,\beta,\gamma,\omega_B,\chi,\kappa\) | \(0.4,0.4,0.3,0.4,0.15,0.5\) |
| \(y,b,\vartheta\) | \(3.73375,0.2,41/240\) |
| \(T_t,T_{t+1},u_t,u_{t+1}\) | \(0.0175,0.0175,0.525,0.525\) |
| Initial old-owner \((a,H)\); old-renter \(a\) | \((1.3625,1)\); \(1.9825\) |
| Logistic taste location and scale | \(-0.2749785880931923,0.12\) |

The young owner chooses \((n,h,x,a')=(0.5,1,2,1.65125)\). Its future interior choices are \((c^2,h^2,e)=(1.6,0.9142857143,1.28)\). Both retention and estate bounds are strictly slack; the down payment is exactly \(0.2\), and \(MV^Y-u=0.5416666667>0\). The young renter chooses \(n=0.2760297046\), \(h=0.5\), \(x=2.2595163746\), and \(a'=1.3878291698\); its old rental cap binds. The original fertility first-order conditions hold before fertility is frozen for the comparison.

Current old owners choose \((c^2,h^2,e)=(1.4,0.8,1.12)\). Current old renters have capped housing \(0.5\), consumption \(1.2410714286\), and estate \(0.9928571429\). Each of the four groups has mass \(0.5\), so housing sums to \(1.4\). Both tenure shares are strictly between zero and one: the reported taste location sets the optimized young ownership share to one half. The old distribution is an explicitly specified inherited state; it is not asserted to be generated by a stationary prior cohort.

The coefficients are

| Group | \(c_i\) | \(d_i\) | \(k_i\) |
|---|---:|---:|---:|
| Young owner | 0.4177528090 | 5.0000000000 | 1.0887640449 |
| Young renter | 0.4950000000 | 0 | 0 |
| Old owner | -0.1975000000 | 0.3361344538 | -1.2862745098 |
| Old renter | 0.4950000000 | 0 | 0 |

Hence \(C=0.6051264045>0\) and \(K=-0.0987552324<0\). The local impossibility proposition applies despite strict young financing rationing and a strictly larger young housing willingness to pay. This is a concrete obstruction to inferring a cash-transfer Pareto gain from that gap.

## 4. What happens in the saved illustration

The source is row 0 of `output/model/simplified_olg_amendments/transition_phi80.csv`, with the parameters in the proposal's final appendix. This is an illustration, not an estimated quantitative benchmark.

\(P_0=0.34268356297591945\), \(u_0=0.17990887056236\), \(T_0=0.00588696586437\), \(Y=O=1.4552639291234173\), \(q=0.5\), \(\tau^p=0.05\), \(\phi=0.8\), \(\bar H=2\). Entrants have \(y=1,b=0.06\), and \((\alpha,\beta,\gamma,\omega_B,\chi,\kappa)=(0.4,0.4,0.3,0.4,0.15,0.5)\). Rental and owner size limits are \(0.45,2\). Future prices and rebates are fixed at their baseline stationary values for the conditional exercise.

| Group | Mass | Housing | Baseline fertility / inherited state |
|---|---:|---:|---|
| Young owner | 1.0893612826199552 | 0.8754432147102463 | \(n=0.5494856788378992\), \(x=0.4933973791493271\) |
| Young renter | 0.3659026465034622 | 0.45 | \(n=0.352671788308486\), \(a'=0.3586354155314221\) |
| Old owner | 1.0893612826218146 | 0.6581964003570613 | \(a=0.3651334697777289,H=0.8754432147102463\) |
| Old renter | 0.3659026465016027 | 0.45 | \(a=0.7172708310628443\) |

The negligible differences between young/old owner masses are rounding/solution residuals from the saved row, not different intended stationary shares.

| Group | \(c_i\) | \(d_i\) | \(k_i\) |
|---|---:|---:|---:|
| Young owner | 0.3975930882 | 14.5907202452 | 3.2465001136 |
| Young renter | 0.4440709856 | 0 | 0 |
| Old owner | -0.2179709187 | 0.9808887560 | -3.0882034612 |
| Old renter | 0.4440709856 | 0 | 0 |

Thus \(C=0.5206469346>0\), \(K=0.1724422443>0\), and

\[
0<K/C=0.3312076435<d_{OO}=0.9808887560.
\]

**Exact conditional Pareto construction.** Choose \(P<P_0\) sufficiently close. Keep each young group's utility exactly constant with its compensating transfer. Since the old renter's cap remains binding, choose old-owner housing to clear the fixed stock. Infer its transfer from its original housing demand, and use the old-renter transfer to balance the cash budget. This construction solves both market and fiscal equations exactly.

For reproducibility, young-owner compensation is the root in \(\ell\) of

\[
(1+\beta A)\log x(P,\ell)+\alpha\log[h(P,\ell)-\kappa n]
=(1+\beta A)\log x_0+\alpha\log(h_0-\kappa n),
\]

where

\[
h(P,\ell)=\frac{b+\ell}{\delta P},\qquad
x(P,\ell)=\frac{y+b+\ell+t_PP+qT_{t+1}-\chi n-u(P)h(P,\ell)}{1+\beta A}.
\]

The young-renter transfer is \(\ell_{YR}=(a_Ph_R^{\max}-t_P)(P-P_0)\). Then set

\[
\begin{aligned}
h_{OO}&=\frac{\bar H-m_{YO}h_{YO}-(m_{YR}+m_{OR})h_R^{\max}}{m_{OO}},\\
\ell_{OO}&=Au(P)h_{OO}/\gamma-a_{OO}-PH_{OO}-t_PP,\\
\ell_{OR}&=-[m_{YO}\ell_{YO}+m_{YR}\ell_{YR}+m_{OO}\ell_{OO}]/m_{OR}.
\end{aligned}
\]

The implicit-function theorem gives the young-owner compensation root because utility is strictly increasing in \(\ell\). All formulas are continuously differentiable on these strictly maintained branches. Writing \(P=P_0-\epsilon\), the two old groups' first-order welfare surpluses in cash units are

\[
g_{OO}=\frac{K\epsilon}{m_{OO}d_{OO}}>0,\qquad
g_{OR}=\frac{(C-K/d_{OO})\epsilon}{m_{OR}}>0.
\]

Their strict signs persist for sufficiently small finite \(\epsilon\). Young utilities are constant by construction, not merely constant to first order. This establishes the conditional Pareto claim.

For \(\epsilon=10^{-4}\), the independent scalar calculation gives:

| Group | Cash transfer per household | Utility change | New housing |
|---|---:|---:|---:|
| Young owner | -0.00003973313023 | 0 to numerical precision | 0.87511885200968 |
| Young renter | -0.00004440709856 | 0 exactly by budget | 0.45 |
| Old owner | +0.00003745387111 | +0.00003970938445 | 0.65852076305763 |
| Old renter | +0.00005119288340 | +0.00020839343316 | 0.45 |

The mass-weighted transfer residual is zero and the housing residual is \(4.44\times10^{-16}\). All relevant cash, saving, retention, and estate inequalities remain feasible. Current old-owner estate slack is \(0.0900857009\); the young owner's future retention and estate slacks are \(0.2168378022\) and \(0.0902328382\). A second \(\epsilon=10^{-5}\) gives the same signs.

**The direction is toward old households.** More generally, since \(C>0\), equation (1) requires a weakly lower price for any nearby conditional Pareto improvement. Interior old owners increase housing when their utility is held constant and price falls; giving them a utility gain raises their housing further. Renters' housing remains at the cap. Housing clearing therefore forces total young-owner housing to fall. This result is an obstruction to the desired **old-to-young** cash-transfer implementation, not a proof of constrained efficiency: the explicit Pareto improvement above exists.

## 5. Tenure, claims, and the missing OLG proof

**Tenure.** Both examples retain mixed tenure tastes and the physical size segmentation. The proofs concern a neighborhood with the displayed branches maintained; they are not global statements across tenure switches or cap changes. The original model fixes old tenure and bars an old owner from buying additional housing, so neither old-tenure switching nor new old title is used. For young households, freezing each household's own optimized baseline fertility gives a useful local qualification: when the two tenure branches have distinct unique fertility optima, a household at the original taste threshold strictly prefers its own tenure when evaluated at its own frozen fertility. Strict concavity therefore supplies a positive local tenure margin. If fertility is instead fixed to one common value or tastes/transfers are defined differently, the tenure margin must be checked again; smooth logit shares alone are not that check. The finite receipts verify conditional branches and constraints; they do not constitute a global search over all alternative tenure policies.

**Initial titles and intermediaries.** The intermediary earns its prescribed return on a unit acquired at the new price because the user-cost relation is imposed. That does not insure the capital value of titles it already owned before an unexpected reform. In the saved illustration, current old households' inherited title totals about \(0.954\), below the stock \(2\). The residual stock is not an unassigned free asset. The source's external trading closure does not identify the ultimate welfare recipient of every initial intermediary/estate claim. The conditional construction lowers \(P_t\), and therefore can lower those initial claims' value. If their owners belong in the Pareto comparison, their compensation must be funded. In particular, the conditional construction cannot be advertised as leaving every initial title holder unaffected.

**Future cohorts and clearing.** Even though each current young household's lifetime utility is held constant, its saving, old consumption, estate, and old housing change. Consequently the next old distribution changes. The fixed continuation prices used above generally will not clear that changed next-date housing market. A complete OLG proof must solve \(P_{t+1},P_{t+2},\ldots\), update the property-tax rebates at each date, and check all subsequent cohorts. Fixed fertility keeps their masses fixed; it does not fix their endowments, demand, prices, or welfare. Estates remain warm-glow expenditure and must not be converted into descendant initial wealth.

Future price effects cannot be suppressed in a path derivative. For example, an interior-old young household has the continuation term

\[
q\,[dT_{t+1}-h^2_{t+1}\,du_{t+1}]
\]

in its welfare change measured in current consumption units, in addition to its current housing-user-cost effect; \(du_t=(1+q\tau^p)dP_t-q\,dP_{t+1}\). Binding estate restrictions add further terms involving their shadow values and future prices. These terms are absent only because Section 2 explicitly fixes the continuation.

**External resources.** Balanced household cash transfers do not themselves create government debt. Private changes in mortgages, bonds, and estates must still be valued at \(q\); existing creditors must receive their promised returns. A full path must state initial outside claims and the terminal/no-Ponzi condition. Neither a zero current transfer budget nor free entry in newly acquired rental housing replaces this accounting.

## 6. Verification and the remaining decisions

Three checks were performed:

1. Derived the envelopes, exact compensation functions, finite-neighborhood obstruction, and exact conditional Pareto construction.
2. Checked them against the proposal's young budgets, mortgage repayment, property-tax rebate, old retained-title/estate constraints, and mixed-tenure formulation. Both examples retain a positive property tax. No direct housing assignment or additional public lending appears in the constructed cash-transfer policy: housing quantities follow private choices at the clearing price.
3. Used independently written scalar arithmetic, without importing or running the model solver. The saved illustration's exact compensation was checked at two finite price changes and its original two-budget utility reconstructed. At price steps \(\pm10^{-6}\), the compensated young housing slopes were \(3.24647135\) and \(3.24652888\), bracketing the analytical \(3.24650011\); original utility differences were at most \(4.44\times10^{-16}\), and budget residuals at most \(2.3\times10^{-16}\). The explicit counterexample provides the separate adversarial check against a universal positive claim. Calculations took less than one second each; no equilibrium sweep or calibration ran.

The main source equations are in `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex`, especially the household, equilibrium, and compensated-sale sections. The live September 5 instruction is in `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/ACTIVE_DECISION_LEDGER.md` under “Both efficiency comparisons requested.”

Minimal decisions before attempting an OLG theorem are:

- Whether grants/taxes really occur before closing and are eligible cash; whether targeting may use the observed taste or predetermined household groups.
- Whether the admissible schedule is one unexpected transfer date or a committed sequence of separately balanced transfers, and whether all initial title owners and all subsequent cohorts must be compensated.
- The ownership and treatment of initial intermediary/estate claims and the external terminal condition.
- Whether public bridge credit, secured advances against sale proceeds, compulsory future repayment, or changes to mortgage eligibility are permitted as **additional** instruments. These may support a different theorem; they are absent here.

The direct-allocation compensated-sale proof uses additional financing and controls the transaction. It cannot be relabeled a proof for this transfer class. A sequence of unrestricted compensations should likewise not be introduced after observing losses without treating that sequence as a different fiscal commitment.

The applicable conceptual reference is Alberto Bisin's [2014 general-equilibrium lecture notes, Section 3.2.4](https://bpb-us-e1.wpmucdn.com/wp.nyu.edu/dist/c/16384/files/2019/11/Lecture-notes-Oct2014.pdf), printed pp. 43–48: the constrained planner anticipates the prices after transfers. Their generic incomplete-market discussion does not prove an inefficiency theorem for this deterministic OLG model. The results above are derived directly from the present household equations.

**Recommendation.** Preserve constrained inefficiency as an open OLG question. The proved conditional illustration supplies a useful warning: a cash-transfer Pareto gain can exist while moving housing in the opposite direction from the intended young-housing mechanism. The separate direct-allocation benchmark is needed to establish that mechanism under its own explicit permissions.
