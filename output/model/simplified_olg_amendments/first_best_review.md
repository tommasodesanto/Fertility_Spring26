# Direct allocation at fixed fertility

Status: **PROVED UNDER STATED ASSUMPTIONS** for the local Pareto result below. A complete first-best allocation for the whole economy remains **NEEDS AUTHOR CHOICE**. This report concerns the September 5 proposal, with no capital-gains tax, one aggregate housing stock, and the original household utility and value-function notation. It does not change the model or the author-owned manuscript.

The useful result is simple. If an eligible young owner values another unit of housing more than an old owner values keeping it, by enough to cover its real transfer cost, the planner can move a little housing to the young household and compensate the old household with current consumption. Every household can retain its fertility, tenure, estate, and future consumption and occupancy. No welfare weights or complete solution of a planner's problem are needed to prove this allocation inefficient.

The inequality is a comparison of housing values in consumption units. Young age, greater raw marginal utility of housing, and a binding borrowing constraint alone do not establish it.

## 1. Which allocation problem is well defined?

Fix the competitive allocation and its entire population path. Keep each household's fertility $n_i$, actual tenure-taste draw $\xi_i$, and identity fixed. For a young household born at $t$, welfare remains the original lifetime utility:

\[
u_t^Y(c_i,h_i,n_i)+\beta u^2(c_i^2,h_i^2,e_i)
+\xi_i\mathbf 1\{m_i=O\}.
\]

For the initial old cohort, the remaining welfare object is $u^2(c_j^2,h_j^2,e_j)$; past utility is sunk. Using actual taste draws makes the comparison household by household, rather than only a comparison of average logit values. The construction below preserves tenure, so the taste term is unchanged.

### Physical housing feasibility

There is one fixed stock, not separate fixed owner and rental stocks. At every date, direct allocations must satisfy:

\[
\int_{\mathcal Y_t}h_i\,di+\int_{\mathcal O_t}h_j^2\,dj\leq\bar H.
\]

The measures in these integrals have the fixed cohort masses $Y_t$ and $O_t$. Young renters obey $h_i\leq h_R^{\max}$; young owners obey $h_i\leq h_O^{\max}$; old renters obey $h_j^2\leq h_R^{\max}$. An old owner has $0<h_j^2\leq H_j$, where $H_j$ is the housing it bought when young. Owners cannot buy additional housing when old or rent their retained housing to another household. Every young bundle also obeys $c_i>\chi n_i$ and $h_i>\kappa n_i$.

These restrictions give a precise feasible set when lifetime tenure is fixed at the competitive assignment. Permitting tenure changes is a separate enlargement: one must specify when tenure may change, what happens to title, and whether a conversion uses real resources. Nothing below requires that enlargement. The transfer is owner to owner and leaves rental occupancy and the rental intermediary's holdings unchanged.

### Dated goods and external claims

A planner need not obey each household's competitive-price budget. It must obey the economy's real dated resource budget, including the world bond terms. One consistent accounting convention is to let $A_t$ be the face value of the economy's net external financial claims delivered at $t$, after consolidating households, government, and rental intermediaries. A positive $A_t$ is an external asset. Let $\mathcal E_t$ be genuinely exogenous net goods inflows at $t$, and let $\mathcal C_t$ be real reallocation costs. The dated resource restriction is:

\[
\int_{\mathcal Y_t}c_i\,di+
\int_{\mathcal O_t}c_j^2\,dj+
\int_{\mathcal O_{t-1}}e_j\,dj+
\mathcal C_t+qA_{t+1}\leq\mathcal E_t+A_t.
\]

An estate chosen by an old household at $t-1$ is a goods expenditure at $t$, as the proposal specifies. Its housing sale is not an additional source of aggregate goods. The good transferred by the buyer and received by the estate is counted once. Entrant wealth $b_i$ is not financed by the modeled estates. If $b_i$ is an external entrant endowment, its entry flow belongs in $\mathcal E_t$; if it is a domestic claim or transfer, it belongs in consolidated initial claims or transfers instead. It cannot be counted twice.

Initial external claims, already promised external payments, and initial outstanding estates must be honored. Over a finite comparison, retain the baseline terminal external position; over an infinite comparison, impose the appropriate external solvency/no-Ponzi restriction. For example, a sufficient standard boundary is $\liminf_{T\to\infty}q^{T-t}A_T\geq0$, when the corresponding present values exist. Removing a household borrowing restriction is not permission to obtain unrepayable foreign resources.

All property-tax payments and rebates are internal transfers after consolidation. The proposal's identity $(Y_t+O_t)T_t=q\tau^pP_t\bar H$ does not create or destroy goods. Domestic house prices, rents, mortgage advances, and title-sale receipts must not be added as real resource costs in the planner's aggregate budget. Real transaction costs belong there separately. The world bond price $q$, in contrast, continues to price external dated trade.

For a **complete** planner program, the origins of initial claims and entrant wealth and the settlement of estates must be pinned down explicitly. The local proof is stronger on this accounting dimension: its current goods budget balances exactly, every estate expenditure stays fixed, and its net external position is unchanged at every date. It requires no increase in national foreign borrowing.

### Estate treatment is not automatically a physical restriction

The original old-owner problem contains the condition:

\[
e_j\geq P_{t+1}h_j^2.
\]

It says that the estate contains the retained house valued at its future sale price, with a nonnegative financial remainder. The retained title is physical; the lower bound stated in goods uses a market valuation. A global direct-allocation planner has not yet been assigned a new competitive housing price path with which to value it.

Two different completions should not be conflated. A planner could deliver a real goods estate and arrange title settlement collectively; that requires saying whether it may replace or tax the title proceeds and relax the original estate-composition rule. Alternatively, it could retain the institutional composition rule at a stated settlement price and specify the settlement counterparties. Treating the planner's scarcity multiplier as that settlement price is another substantive choice, not an identity. None has been selected by the author.

The local proof avoids choosing among these global completions. It fixes estate expenditure and uses the original price path only to verify an admissible settlement of the changed titles. The seller's lower retained housing relaxes its original estate bound; the buyer's future occupied housing and estate are unchanged. No estate-composition permission is silently removed.

### A precise Pareto problem that is sufficient here

Choose one young owner $i$ and one old owner $j$. For clarity they have equal mass; transfers between unequal-mass groups are scaled to equal aggregate housing units. Hold every other current allocation, all subsequent consumption and occupancy, all estates, fertility, and tenure fixed. Choose current $c_i',c_j^{2\prime}$ and a transfer $\epsilon\geq0$ to maximize the young household's lifetime utility subject to the old household retaining its baseline utility. The nontrivial constraints are:

\[
\begin{aligned}
h_i'&=h_i+\epsilon\leq h_O^{\max},&
h_j^{2\prime}&=h_j^2-\epsilon>0,\\
c_i'+c_j^{2\prime}+C(\epsilon)&\leq c_i+c_j^2,&
u^2(c_j^{2\prime},h_j^2-\epsilon,e_j)&\geq u^2(c_j^2,h_j^2,e_j),\\
c_i'&>\chi n_i,&
h_i^2&\leq h_i+\epsilon.
\end{aligned}
\]

The last condition is the young owner's unchanged old-age occupancy constraint. The original donor conditions $h_j^{2\prime}\leq H_j$ and $e_j\geq P_{t+1}h_j^{2\prime}$ are retained and hold automatically when housing is reduced. The buyer's old estate-composition condition is unchanged. Section 3 verifies title and financial settlement against the original equations.

This is a direct allocation problem with individual utility protection, physical constraints, and a goods budget. There is no household down-payment constraint in this benchmark. It is a fully specified restricted Pareto problem; failure of efficiency in this smaller set proves failure in every broader direct-allocation set that admits these changes. It is not a characterization of the entire economy's first-best allocation.

## 2. A local Pareto proposition

The young household's adult consumption and space are $x_i=c_i-\chi n_i>0$ and $s_i=h_i-\kappa n_i>0$. The two consumption-unit housing values computed from the original preferences are:

\[
MV_i^Y=\frac{\alpha x_i}{s_i},\qquad
MV_j^O=\frac{\gamma c_j^2}{h_j^2}.
\]

Each is the amount of its own current consumption a household would exchange for a marginal housing unit, holding its other arguments fixed. They are invariant to rescaling a household's utility representation by a positive constant.

**Proposition — PROVED UNDER STATED ASSUMPTIONS.** Suppose the young owner has $h_i<h_O^{\max}$, the donor is an old owner with $h_j^2>0$, and both start from feasible bundles. The planner may redistribute current goods and settle titles and financial claims as in Section 3 while relaxing private borrowing limits. Let $C(0)=0$, $C\geq0$, and $C'(0+)=K\geq0$ exist and be finite. If:

\[
MV_i^Y>MV_j^O+K,
\]

then sufficiently small positive transfers of housing from $j$ to $i$ are Pareto improvements. The old seller is exactly indifferent, the young household is strictly better off, and every other household, estate, and outside creditor can retain its allocation or contractual payment. Fertility, cohort masses, and tenure are unchanged.

**Proof.** With its estate fixed, exactly compensating the old household for losing $\epsilon$ units requires extra current consumption:

\[
D_j(\epsilon)=c_j^2\left[\left(\frac{h_j^2}{h_j^2-\epsilon}\right)^\gamma-1\right].
\]

Give it $D_j(\epsilon)$, and reduce the young household's current consumption by $D_j(\epsilon)+C(\epsilon)$. Current goods and housing balance exactly. Small enough $\epsilon$ preserves positive consumption, space, and donor housing and respects the buyer's physical upper bound. The old household's consumption gain exactly cancels its housing utility loss; all future utility arguments are unchanged. The young household's lifetime gain is:

\[
\Delta U_i=
\log\left(1-\frac{D_j(\epsilon)+C(\epsilon)}{x_i}\right)
+\alpha\log\left(1+\frac{\epsilon}{s_i}\right).
\]

Since $D_j'(0)=\gamma c_j^2/h_j^2$, its right derivative is:

\[
\Delta U_i'(0+)=\frac{MV_i^Y-MV_j^O-K}{x_i}>0.
\]

Continuity and the positive derivative give the claimed strict gain for sufficiently small positive transfers. The original title and claim constraints are checked below. This completes the proof without imposing an interior-saving Euler equation.

**Competitive-equilibrium corollary — PROVED UNDER STATED ASSUMPTIONS.** An old owner with slack retention and estate-composition constraints satisfies $MV_j^O=q r_t$. For a young owner with positive saving, slack old-age retention and estate-composition constraints, a slack physical owner cap, and a strictly positive down-payment multiplier, its first-order conditions imply $MV_i^Y>q r_t$. Therefore such a pair implies a local Pareto improvement at zero real reallocation cost, or whenever $K<MV_i^Y-q r_t$.

To see exactly what the financing statement uses, let $\lambda_i=1/x_i$ be the multiplier on the young goods budget and $\eta_i>0$ the multiplier on $(1-\phi_t)P_th_i\leq b_i$. The saving and housing conditions, under the stated slackness assumptions, give:

\[
c_i^2=\frac{\beta x_i}{q},\qquad
MV_i^Y=q r_t+\frac{\eta_i}{\lambda_i}(1-\phi_t)P_t.
\]

The term on the right is the consumption value of relaxing the financing limit. A numerically binding inequality need not have a strictly positive multiplier. Nor does the corollary prove that a positive measure of eligible pairs exists for every parameterization. That is a condition to establish for the equilibrium being discussed, rather than an implication of age alone.

The proposition itself is less restrictive than this corollary: it can use an old donor at a retention or estate boundary, or a young buyer without interior saving, if the actual consumption-unit value gap and the permitted financial settlement satisfy its assumptions.

A necessary efficiency condition follows in the same way. At a zero-cost direct Pareto optimum, if a current housing unit and current goods can be transferred in either direction while keeping the relevant future retention and estate constraints feasible, the two current consumption-unit housing values must agree. Otherwise the direction with the higher recipient value gives the improvement just proved. With physical bounds, this is a one-sided inequality instead. Equality on one pair does not establish whole-economy efficiency: other pairs, intertemporal trades, estate choices, or tenure changes can still offer gains. This condition characterizes the locally relevant margin without assigning social weights.

## 3. Original-constraint and external-claim check

Write the unchanged equilibrium user cost as $u_t=q r_t=(1+q\tau^p)P_t-qP_{t+1}$, and define a reference-price transfer $L=D_j(\epsilon)-u_t\epsilon$. Prices here document how to settle the direct allocation; they do not constrain its real goods budget or claim to be prices of a new equilibrium.

| Original object | Exact check under the perturbation |
|---|---|
| Current housing stock | The young owner's housing rises by $\epsilon$, the old owner's occupancy falls by $\epsilon$; the single aggregate stock is unchanged. |
| Rental segmentation | No rental household or rental holding changes. The young owner remains below $h_O^{\max}$; there is no conversion across tenure or new rental stock. |
| Old donor retention | $h_j^2-\epsilon\leq H_j$, including when the original bound binds. |
| Old donor estate | Keep $e_j$ fixed. Its original lower bound falls from $P_{t+1}h_j^2$ to $P_{t+1}(h_j^2-\epsilon)$. The seller replaces the missing title proceeds with a bond delivering $P_{t+1}\epsilon$. |
| Young buyer's subsequent occupancy | Its inherited title next period becomes $H_i'=h_i+\epsilon$. Keep its old housing $h_i^2$ fixed. The retention bound weakens, and it can sell the additional title when it becomes old. |
| All subsequent title availability | At $t+1$, the buyer sells the extra unit when the donor's estate would have sold that same unit. The unit remains available to the same later occupants. |
| Young buyer's subsequent estate | Both $h_i^2$ and $e_i$ stay fixed, so the original estate-composition inequality is unchanged. |
| Current goods | $\Delta c_i+\Delta c_j^2+C=-D_j-C+D_j+C=0$. All later goods use, including estates, is unchanged. |
| Property taxes and rebates | The taxable stock and price path are unchanged. The extra owner tax reserve $q\tau^pP_t\epsilon$ replaces the donor's former reserve; total receipts and common rebates are unchanged. |
| Fertility, population, and tastes | Each $n_i$, lifetime tenure, and actual $\xi_i$ is fixed. All cohort sizes and the demographic recursion remain unchanged. |

The original young-owner equations also give an exact financial representation. Change its saving account by:

\[
\Delta a_i'=\phi_tP_t\epsilon-qP_{t+1}\epsilon.
\]

Its future financial wealth plus housing liquidation value then has zero change:

\[
\Delta(a_{i,t+1}+P_{t+1}H_i)
=\frac{\Delta a_i'}q-\frac{\phi_tP_t\epsilon}q+P_{t+1}\epsilon=0.
\]

Thus its unchanged future bundle satisfies the original old budget; allowing more inherited title does not invalidate any of its physical or estate bounds. The change in its current budget's left-hand side is:

\[
-D_j-C+\Delta a_i'+[(1-\phi_t)+q\tau^p]P_t\epsilon
=-L-C.
\]

A current transfer $-L-C$ to the buyer and $L$ to the seller, with the remaining $C$ used for the real transfer activity, reproduces the direct goods allocation. For the seller, the old-budget left-hand side changes by $D_j-u_t\epsilon=L$. This algebra is valid whether $L$ is positive or negative; targeted redistribution belongs to the direct-allocation permission.

There is no increase in net foreign financing. The seller purchases a replacement estate claim costing $qP_{t+1}\epsilon$; the buyer's saving change and additional mortgage have net cost $-qP_{t+1}\epsilon$. Their consolidated change in external asset purchases is exactly:

\[
\Delta a_i'-\phi_tP_t\epsilon+qP_{t+1}\epsilon=0.
\]

Tax reserves offset separately. Existing foreign claims are honored, and any new gross lending receives the original market return. If a particular lender cannot provide the bridge or a title is itself subject to an additional unmodeled lien, that restriction must be specified; it is not established by the proposal.

**Private financial restriction exposed.** If $a_i'=0$ and $\phi_tP_t<qP_{t+1}$, this representation requires negative $a_i'+\Delta a_i'$. The planner must then replace it with permitted public/unsecured finance or otherwise relax that private financial bound. The aggregate resource proof does not claim the original private nonnegative-saving constraint still holds. If one insists on that private constraint, require sufficient initial saving or $\phi_tP_t\geq qP_{t+1}$. The original down-payment constraint is intentionally relaxed: a strictly constrained buyer could not acquire the extra unit under its unchanged cash limit. These are precisely why this is not a transfers-followed-by-markets theorem.

## 4. What the existing compensated-credit proof establishes

The proposal's existing proof sets $L=D_j-u_t\epsilon$, extends purchase credit, keeps the young household's current consumption fixed, and reduces its old consumption by $[D_j+C]/q$. It checks an exact transaction, seller compensation, future resale, estate replacement, taxation, and creditor repayment. Under its assumptions, its young utility gain is:

\[
\alpha\log(1+\epsilon/s_i)+
\beta\log\left(1-\frac{D_j(\epsilon)+C(\epsilon)}{q c_i^2}\right).
\]

Its derivative equals the consumption-unit housing gap only when the interior-saving equation $q c_i^2=\beta x_i$ applies. Under that equation it gives exactly the same first-order value gap as the current-consumption construction. The two proofs use different dated compensation paths; neither is a complete first-best solution.

It already suffices to refute first-best efficiency if that transaction belongs to the chosen direct-allocation feasible set. It does not establish that every competitive equilibrium is inefficient, that every old owner should sell, that any arbitrary weighted planner gives more housing to young households, or that cash transfers followed by market clearing improve every household. The current-goods proof removes unnecessary Euler-equation dependence from the basic allocation test while keeping the stronger implementation question separate.

Unequal raw housing marginal utilities do not prove inefficiency. For example, with $\alpha=\gamma=1$, let the young household have $x_i=s_i=1$ and the old household have $c_j^2=h_j^2=2$. The raw housing marginal utilities are $1$ and $1/2$, yet both consumption-unit housing values equal $1$. This allocation maximizes the young current log utility plus twice the old current log utility with totals $x_i+c_j^2=3$ and $s_i+h_j^2=3$; it is Pareto efficient in that fixed-continuation specialization. Raw utility comparisons also change under an irrelevant rescaling of one person's utility.

## 5. Precisely what “more housing for young” means

The proposition constructs an improving allocation with $\Delta H_t^Y=\epsilon>0$ and $\Delta H_t^O=-\epsilon$, for a pair of positive-mass households or appropriately scaled groups. With a continuum, a single type of zero mass does not establish a positive aggregate transfer. There must be a positive-measure set of eligible pairs. A uniform positive gap, physical slack, and bounded curvature on such a set allow a common small transfer; continuity and positive type density around an eligible interior pair are one sufficient way to obtain such a set.

The direction follows the value comparison, rather than a chosen social preference for young people. It proves that **there exists a Pareto-improving reallocation toward the young**. It does not prove that every Pareto improvement, or every entire-economy efficient allocation, gives young households more aggregate housing than the competitive allocation. Other pairs, dates, estates, and consumption allocations may move differently in a complete optimum.

The dependence of a weighted optimum on weights is visible even in the original preferences with continuations fixed. With available adult housing $S$, zero transfer cost, slack caps, and positive weights $w_Y,w_O$, the weighted optimum allocates:

\[
s_i^*=\frac{w_Y\alpha}{w_Y\alpha+w_O\gamma}S,
\qquad h_j^{2*}=S-s_i^*.
\]

As $w_Y/w_O$ falls, young adult space falls. The young total remains $h_i^*=s_i^*+\kappa n_i$, so child space is retained. The model does not select those weights. A globally efficient weighted allocation can therefore have less young housing than a given competitive allocation, even though that competitive allocation also admits a small Pareto-improving transfer toward an eligible young group.

## 6. A diagram derived from the actual preferences

**PROVED UNDER STATED ASSUMPTIONS:** an exact diagram is available for the restricted two-household Pareto problem in Section 1. Put housing transferred to the young household, $\epsilon$, on the horizontal axis. Put current consumption goods per unit of housing on the vertical axis. At each $\epsilon$, keep the old household exactly indifferent and pay the real cost from the young household. The plotted values must be:

\[
\begin{aligned}
c_j^2(\epsilon)&=c_j^2\left(\frac{h_j^2}{h_j^2-\epsilon}\right)^\gamma,\\
x_i(\epsilon)&=x_i-D_j(\epsilon)-C(\epsilon),\\
MV_i^Y(\epsilon)&=\frac{\alpha x_i(\epsilon)}{s_i+\epsilon},\\
MV_j^O(\epsilon)&=\frac{\gamma c_j^2(\epsilon)}{h_j^2-\epsilon}.
\end{aligned}
\]

Assume, for this figure only, a differentiable, nondecreasing convex real cost $C$. The young value decreases because it has more space and pays compensation; the old value increases because it has less housing and more compensating consumption. Hence $MV_j^O(\epsilon)+C'(\epsilon)$ rises. The young utility derivative along this exact compensation path is:

\[
\frac{dU_i}{d\epsilon}=
\frac{MV_i^Y(\epsilon)-MV_j^O(\epsilon)-C'(\epsilon)}{x_i(\epsilon)}.
\]

An initial positive gap therefore gives the rightward improving arrow. If the curves cross before a physical bound is reached, their unique crossing is the optimum of this restricted Pareto problem. Otherwise its optimum is at the feasible physical boundary. In the zero-cost case the intersection equates the two consumption-unit housing values. With positive marginal cost it equates the young value to the old value plus cost.

The correct label is “best allocation holding other allocations and the old household's utility fixed,” or a shorter clearly defined “conditional efficient allocation.” Calling that crossing the economy's full first best would add an unproved claim. A point labelled “CE” at zero is conditional on the plotted baseline actually being an equilibrium; arbitrary parameter choices supply only an illustrative feasible baseline. No assertion that one representative pair describes every age group's aggregate demand is justified without further aggregation assumptions.

The area between these consumption-unit value curves is not, by itself, an exact measure of utility or a global welfare loss. The exact young utility gain is the integral of the displayed gap divided by $x_i(\epsilon)$. This division varies along the path. The original logarithmic expression gives the exact finite gain without an area approximation.

## 7. Boundary and counterexample review

| Claim or possible failure | Status and implication |
|---|---|
| Positive gap with an admissible owner-to-owner direction | **PROVED UNDER STATED ASSUMPTIONS.** The exact goods compensation proves a local Pareto improvement. |
| Strictly binding young credit alone | **OBSTRUCTION.** It gives the simple value wedge only under the relevant saving/continuation conditions, and still needs an eligible lower-value old donor. |
| Old retention or estate bound binds | **OBSTRUCTION to the simple CE corollary, not automatically the proposition.** The old value can exceed $q r_t$. Compare its actual $\gamma c_j^2/h_j^2$; reducing housing still respects its bounds. |
| Young owner at the physical size cap | **OBSTRUCTION.** The proposed direction is not physically feasible, however high the value gap. A renter at its cap cannot silently be treated as an unconstrained owner. |
| Zero gap, zero cost | **OBSTRUCTION to a strict local gain.** In the convex-cost compensation diagram, the young value falls and the old value rises, so a positive transfer from a zero-gap point reduces young utility. |
| Large marginal cost | **OBSTRUCTION to this direction.** $K\geq MV_i^Y-MV_j^O$ eliminates the strict first-order conclusion. It does not prove global efficiency. |
| Positive fixed moving cost | **OBSTRUCTION to the marginal theorem.** If $C(0+)>0$, test a finite feasible transfer instead. Exact sufficiency is $D_j(\epsilon)+C(\epsilon)<x_i[1-(1+\epsilon/s_i)^{-\alpha}]$. |
| No eligible positive-measure pairs | **OBSTRUCTION to an aggregate young-housing claim.** The model needs enough physical slack and a strict value ordering on a non-negligible set. |
| Relaxing external solvency | **Not permitted.** The proof preserves the net external asset path and every existing payment. A global planner must separately respect its specified external boundary. |
| Altering fertility or population to generate gains | **Outside the result.** Every fertility choice and every cohort mass is fixed. |
| A re-clearing competitive equilibrium after transfers | **Not proved here.** The direct allocation may violate individual optimality at reference prices and relaxes private financing. The sibling benchmark must address market responses. |

The constrained old-owner first-order condition makes the third row explicit. If $\rho_j\geq0$ is its retention multiplier and $\mu_j\geq0$ its estate-composition multiplier, its housing condition implies:

\[
MV_j^O=q r_t+c_j^2\rho_j+c_j^2P_{t+1}\mu_j.
\]

Thus neither old age nor a lower housing-expenditure share guarantees a willing donor at value $q r_t$.

## 8. Three verification passes and minimal choices

**Pass 1 — direct derivation.** The old compensation was solved exactly from the original logarithmic utility. Differentiating it gives $D_j'(0)=\gamma c_j^2/h_j^2$. Substitution into the young original utility gives the exact finite gain and the stated derivative. No reduced fertility problem, welfare weights, or old-saving formula was used in this proof.

**Pass 2 — original constraints.** Section 3 checks original occupancy, retention, estate composition, current and future household accounting, taxes, the external account, tenure tastes, and all future cohorts. It explicitly identifies the relaxed down-payment/no-unsecured-borrowing restrictions. The title accounting uses the original model's permitted partial old-owner sale and subsequent estate sale, with no additional housing stock or rental conversion.

**Pass 3 — adversarial and independent algebra.** The boundary cases above include the raw-MU counterexample, physical-cap obstruction, strict-multiplier qualification, positive fixed costs, and the possible negative saving-account adjustment. An independent in-memory calculation with $q=0.5$, $P_t=P_{t+1}=1$, $\tau^p=0.1$, $\alpha=0.4$, $\gamma=0.3$, $x_i=1$, $s_i=0.5$, $h_j^2=1$, $c_j^2=0.55/0.3$, and $C(\epsilon)=0.05\epsilon+0.1\epsilon^2$ reproduced exact seller indifference, all financial identities at both $\phi=0.4$ and $0.8$, and a young utility slope converging to $0.2$. At $\epsilon=10^{-7}$, the computed slope was $0.1999998546$. The exact compensation diagram crossed at $\epsilon=0.07082663$, where both cost-adjusted values were $0.66927892$. These are algebra checks on constructed bundles, not an equilibrium or calibration result.

The lead's independently generated [planner_benchmark_checks.json](planner_benchmark_checks.json) was inspected as a further check. It reports 27 finite perturbations of constructed conditional households, including property-tax rates $0,0.02,0.05$, 18 improving cases, and 9 cases where cost exceeds the gap. The largest accounting residual is $1.78\times10^{-15}$. These checks also do not establish a full competitive equilibrium or the existence of eligible pairs in every equilibrium.

The remaining author choices can be kept small:

1. **For the local Pareto claim:** admit direct redistribution of current goods and financial/title settlement while relaxing private borrowing limits. Fixed tenure and fixed goods expenditure by estates suffice; no tenure conversion or social weights are needed.
2. **For costs:** choose whether the baseline has zero cost or a smooth nonnegative cost. A positive fixed cost needs the finite test and cannot be hidden inside $K$.
3. **For a complete first-best program:** specify estate settlement and initial/exogenous external resources; decide whether young tenure may be reassigned and which original title restrictions remain. Then select a point on the Pareto frontier, for example by protecting all other households at competitive utility, or supply explicit welfare weights. Infinite-horizon existence also requires its own regularity and resource conditions.
4. **For the figure:** the exact conditional diagram above is ready analytically. A figure claiming the entire economy's first best remains conditional on those additional decisions and a fuller derivation.

The strongest defensible conclusion is therefore a conditional equilibrium-inefficiency theorem with an explicitly feasible transfer toward young owners. The model does not yet deliver an unconditional age ordering of housing values or a uniquely defined whole-economy first-best allocation.

### Sources inspected

- `memory/AGENT_MEMORY.md`, current theory instructions; `memory/daily/2026-09-05.md`; live `CALIBRATION_STATUS.md` for separation from the quantitative task.
- `AGENTS.md` and `docs/style/econ_writing_style_guide.md`.
- `docs/model/ACTIVE_DECISION_LEDGER.md`, especially the September 5 planner entries.
- `latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex`, original preferences, household constraints, equilibrium accounting, Proposition 1, and its appendix proof.
- `output/model/simplified_olg_amendments/planner_benchmark_checks.json`, independently written lead evidence.

Only this report was written by this worker. No paper source, protected manuscript, calibration file, code, memory, or decision ledger was edited.
