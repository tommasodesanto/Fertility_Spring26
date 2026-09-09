# Consolidated context: current model, author decisions, and checked work

Prepared for a new Pro review, September 8, 2026. This digest is research context, not a claim that the theory is complete. The accompanying consolidated prompt has priority. The exact source excerpts at the end describe the maintained household model and an existing sufficient existence argument; they do not reinstate any superseded planner.

## A. Decisions that must survive the handoff

The author wants a simple illustrative theory of housing allocation, ownership and fertility. The intended message is that finance and limited rental sizes can prevent young families from occupying the housing a planner would assign them. The final exposition should resemble a transparent problem set with a few extensions. It should retain the notation and substance the author chose, and avoid technical language that adds no content.

Earlier work pursued compensated/Pareto reallocations and constrained inefficiency. Those results are now background or possible appendix material. The present main route is equally weighted utilitarian welfare. A stationary lifetime comparison that improves later cohorts while sacrificing initial old households was explicitly rejected as answering a different question. The current comparison is a full allocation at one date, optimizing current consumption as well as housing. Do not switch back to a housing-only objective or assume future-cohort compensation solves the current question.

The planner's accepted direct powers are broad enough to relax private finance, including the old financial-saving floor, while honoring obligations and preserving fixed continuation opportunities and net estates. Old competitive households themselves still cannot borrow. They may sell and resize their houses. The author was repeatedly confused by descriptions that made old households mechanically unable to sell; that is not the model. The retained-home size is not a physical constraint.

Tenure is fixed for the dated comparison and persists across the two ages. “Same tenure distribution” means the matched stationary cohorts have the same distribution of owner/renter status, not that young and old occupy the same amount of housing. A statement about a young individual, a matched pair, a positive-mass group, and the entire young cohort are different statements. We would like both an individual characterization and a useful aggregate theorem.

The author initially expected old households to occupy larger total homes. That remains a prediction to investigate. It is acceptable to use a weaker comparison with the young household's adult space if that is enough for the housing-allocation result. Do not assume the desired equilibrium ordering or a marginal-utility gap as the entire theorem. Derive conditions from preferences, endowments, credit and physical constraints. Multiplier conditions and numerical reference-point/continuity arguments were rejected as the main final proposition.

The author is concerned about a standalone restriction on \(\beta R_f\) in either direction. We have not adopted one. Recent inequalities using \(\beta\ge q\) are partial results to improve, not instructions to impose that condition. High annual discount factors in the richer quantitative model do not automatically map into the two-period parameter without matching horizons.

The demographic illustration must start with an exogenous fertility-taste decline that initiates a transition from one steady state to another. Housing policy is introduced along this transition and may lead to a different eventual population level. We are not trying to explain the initial fertility shock. Both new positive stationary endpoints have replacement fertility; the difference in population comes from fertility along the path. Do not replace this with a static policy comparison or two arbitrarily selected fertility paths.

Transaction costs, larger rental homes, another lifecycle stage, and an alternative policy instrument are not adopted. They may be useful if a precise obstruction is established. The current baseline already includes old-age income and current earnings available for the down payment. Entry wealth is exogenous and the estate is warm glow; do not accidentally impose an endogenous dynastic inheritance law.

The separate general-preferences branch has just been corrected at the author's request. Begin with gross \(U^y(c,h,n)\), optionally \(a(c,n)+b(h,n)+v(n)\), without linear offsets. The shifted-log model is a specialization to test afterward. Both the main model and this exploratory branch must remain clearly identifiable.

## B. Useful formulas and their status

All identities below should be independently verified against the original budgets. They summarize existing progress, not an instruction to accept every earlier auxiliary bound.

### Stationary reduction

Write \(d_p=1-q+q\tau^p\), \(p=d_pP\), \(L_R=p\), and \(L_O=(1-\phi+q\tau^p)P\). Current cash is \(w=y^y+b+T\); old income including its rebate is \(v=y^o+T\). Conditional on tenure, use total resources on entering old age \(z\). The young restrictions reduce to
\[
c+ph+qz=w+qv,\qquad c+L_dh\le w.
\]
Let \(K=1+\gamma+\omega_B\), \(E=1+\alpha+\vartheta\), and \(D=E+\beta K\). Define
\[
\Gamma_R=\gamma,\qquad
\Gamma_O=\Gamma=\min\left\{\gamma,
\frac{(\gamma+\omega_B)(1-q+q\tau^p)}{1+q\tau^p}\right\}.
\]
When old housing is uncapped, \(c^o=z/K\) and \(ph^o=\Gamma_dc^o\), including either old-owner estate regime. A sufficient primitive restriction making the financial-estate floor slack is
\[
\omega_B(1-q+q\tau^p)>q\gamma.
\]
Its economic role must be explained; it is not a condition for the old to be able to sell.

### Full dated planner

Let \(x_i=c_i^y-\chi n_i\), \(s_i=h_i^y-\kappa n_i\). An overbar is a mean under the common stationary law \(Q\), not an expectation about aggregate uncertainty. At fixed fertility the planner's adult consumption is common:
\[
x^F=c^{o,F}=(\bar x+\bar c^o)/2.
\]
Its housing choices are
\[
h_i^{y,F}=\min\{H_{d_i},\kappa n_i+\alpha/\lambda\},
\qquad h_i^{o,F}=\min\{H_{d_i},\gamma/\lambda\},
\]
where \(\lambda\) clears housing. This solves the full consumption–housing problem because the fixed-fertility log objective separates. That separation need not hold for arbitrary utility.

For owners the proposed settlement uses
\[
\Delta a'_i=-P_{t+1}\Delta h_i^y,\qquad
\Delta a_j^e=-qP_{t+1}\Delta h_j^o,
\qquad t_k=\Delta c_k+u_t\Delta h_k.
\]
It preserves young future net resources and old net estates. Rental allocations require the corresponding intermediary accounting. With total goods and housing fixed, transfers sum to zero. These are dated settlement identities; they are not an equilibrium tax-policy implementation or permission for competitive old households to borrow.

### Partial equilibrium-to-planner comparisons

The latest attached Pro response used the strong condition
\[
\alpha\ge\gamma,\qquad \beta\Gamma/q\ge\alpha+\vartheta
\]
to put old total housing above young total housing and derive further results. We then investigated weaker adult-space comparisons. Under appropriate slack caps, \(\beta\Gamma/q>\gamma\) is sufficient for the resource comparisons used by the aggregate housing and joint-planner fertility argument. With \(\alpha\ge\gamma\), binding competitive caps can be handled in some comparisons; binding **planner** caps remain an important unresolved generalization. Do not collapse these different cap requirements.

If the old financial-estate floor is everywhere slack, \(\Gamma=\gamma\). At \(\beta=q\), the uncapped-old first-order conditions imply
\[
\alpha h_i^o-\gamma s_i
=(\mu_iL_{d_i}+\eta_i^y)s_i h_i^o\ge0,
\qquad c_i^o=x_i(1+\mu_i/\Lambda_i).
\]
Here \(\Lambda_i\) is the young lifetime-budget multiplier, \(\mu_i\) the young financing multiplier, and \(\eta_i^y\) the young housing-cap multiplier. With \(\alpha\ge\gamma\) and the relevant planner caps slack, positive mass with \(\mu_i>0\) makes the aggregate housing and joint-fertility comparisons strict. These are useful partial results; the goal is to avoid using patience relative to \(q\) as a standalone final assumption if possible.

In the explicit zero-tax, \(\phi=q\), market-cap-slack, positive-financial-estate subcase,
\[
x_i=\min\{w_i/E,(w_i+qv_i)/D\},\qquad
p=(\nu\vartheta\bar x-\chi)/\kappa,
\qquad n_i=\vartheta x_i/(\chi+p\kappa).
\]
At \(\beta=q\), \(F\{v>(K/E)w\}>0\) gives positive strictly constrained mass. This demonstrates a tractable primitive income–cash condition, but \(\phi=q\) is a special benchmark, not an adopted restriction on the main model. If all financial and physical constraints are slack at \(\beta=q\), aggregate housing and average fertility can be unchanged at the planner optimum even with heterogeneous endowments. Individual redistribution may nevertheless improve utilitarian welfare.

The earlier dated Pro response contains a counterfamily with binding young finance but a reversed aggregate housing direction when the old estate regime differs. Preserve and check it before asserting that young borrowing constraints alone suffice. It also contains conservative primitive inequalities; the author found those too elaborate for the main statement. They may inform a proof without becoming the presentation.

### Fertility identities and the corrected aggregation step

For the maintained log utility, a private interior fertility optimum at a given gross bundle satisfies
\[
\frac{\vartheta}{n}=\frac{\chi}{x}+\frac{\alpha\kappa}{s}.
\]
Its conditional response is
\[
dn=\frac{\chi\,dc/x^2+\alpha\kappa\,dh/s^2}
{\vartheta/n^2+\chi^2/x^2+\alpha\kappa^2/s^2}.
\]
For housing rising by \(\Delta h\ge0\) and consumption falling by \(\delta\ge0\), the finite bundle comparison, evaluated at the initial fertility and positive adult bundle, gives
\[
n_1\ge n_0\quad\Longleftrightarrow\quad
\delta\le\frac{\alpha\kappa x^2\Delta h}
{\chi s(s+\Delta h)+\alpha\kappa x\Delta h},
\]
with the usual feasibility/interiority conditions. This is not a funded policy theorem.

If the joint dated planner is uncapped, it chooses common fertility \(n^J\) satisfying
\[
\frac{\vartheta}{n^J}=
\frac{2\chi}{\bar c-\chi n^J}
+\frac{\kappa(\alpha+\gamma)}{\bar h-\kappa n^J},
\qquad \bar c=C^{eq}/N,\quad \bar h=\bar H/N.
\]
The heterogeneity step uses Cauchy–Schwarz after multiplying the reference household fertility equation by **\(n_i^2\)**:
\[
\vartheta\bar n=
\chi\int n_i^2/x_i\,dQ+\alpha\kappa\int n_i^2/s_i\,dQ
\ge\bar n^2(\chi/\bar x+\alpha\kappa/\bar s).
\]
The earlier response said multiply by \(n_i\); the squared version is the corrected argument. Additional comparisons of old versus young mean resources are still needed to sign the joint optimum's fertility gain. The relevant cap assumptions must be checked for the **joint** optimum, not just the reference or fixed-fertility optimum.

### General-preferences identity

For the same stationary reduced budgets and regular interior goods choices, write \(m_i=V_d'(z_i)=U_c^o\), and put \(\rho_i\ge0\) on the old-owner estate floor \(e-Ph^o\ge0\), zero for renters. Then
\[
U_h^y-U_h^o=
(\beta/q-1)p m_i+L_d\mu_i+\eta_i^y-P\rho_i-\eta_i^o.
\]
No logarithm is required. The identity does not prove the full planner's aggregate direction, especially with nonseparable consumption–housing utility. It also shows why eliminating logarithms alone does not eliminate every appearance of \(\beta/q\).

### Existence and remaining transition work

The existing sufficient stationary existence proof is reproduced below. Its low-price fertility condition contains \(D=E+\beta K\), so it has a finite upper bound on \(\beta\) at fixed other primitives. This is a conservative certificate, not a necessary existence condition and not the infinite-lived precautionary-saving restriction.

Strict concavity gives conditional household real-choice uniqueness and a unique dated planner real allocation where an optimum exists. The explicit \(\phi=q\), zero-tax, cap-slack benchmark gives one positive stationary price within that regime when \(\nu\vartheta\bar x>\chi\). This does not exclude additional equilibria involving capped households. General stationary uniqueness and deterministic transition existence, uniqueness and convergence have not been established. None of the recent static results proves a funded general-equilibrium property-tax transition. Older transition arguments used earlier household specifications and should not be silently imported.

## C. Literature starting points and visual structure

These links were checked in the preceding work; verify the actual statements if using them in the new result.

- Coven, Golder, Gupta and Ndiaye, August 1, 2026 version: <https://abdouecon.github.io/research/papers/Property_Tax.pdf>. The simple two-period section establishes capitalization and intergenerational redistribution with consumption utility. Its quantitative lifecycle model includes housing, finance, transaction costs and estates. It does not supply a theorem that our equally weighted planner gives more housing to young households. Welfare and allocation effects of the reform need separate assessment.
- van Doornik, Fazio, Ramadorai and Skrastins, Housing and Fertility, December 2024 version: <https://bcb.gov.br/content/publicacoes/WorkingPaperSeries/WP612.pdf>. Goods and space are useful fertility ingredients. The simple model does not provide our endogenous tenure/mortgage architecture. Verify exact formulations before borrowing a result.
- Aiyagari (1994), especially pp. 668–670: <https://www.liuyanecon.com/wp-content/uploads/Aiyagari-1994.pdf>. Its stationary precautionary-saving condition concerns an infinite-lived income-risk environment. It should not be applied as a theorem about this finite-lived model with exogenous entry wealth.

The author wants two familiar illustrations and only 5–7 theory slides. First show housing moving from an old household to a young household, with labeled market and planner allocations and housing marginal-utility curves that match the welfare criterion. An older compensated-Pareto figure is a visual reference, not the current welfare theorem; do not carry over compensation claims automatically.

The second illustration combines a fertility decline and an intervention during the resulting transition. It should show an initial steady state \(S_-\), a baseline destination \(S_0\), and a reform introduced at a common inherited state that leads toward \(S_1\). Two panels should make the fertility and population/equilibrium movements readable together. Curves, axes and arrows must correspond to the proven mechanism; do not represent all transition points as lying on a stationary schedule. The existing deck's recent curves were acknowledged to assume fertility paths. The author explicitly rejected treating those as a solved equilibrium transition. Numerical illustrations belong mainly in the quantitative model, so an analytical or clearly labeled schematic construction is preferable here.

## D. Maintained household model: source excerpt

The following is extracted from `latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex`, from Environment through the stationary household reduction. Only old-age superscripts have been harmonized from the stale `2` to the author's requested `o`. This is the household model, **not** the new full dated planner definition, which is supplied by the consolidated prompt above. Allowing the existing property-tax parameter to be date-dependent for a permanent reform is the explicit policy comparison requested in that prompt.

```latex
\section{Environment}

\textbf{Households.} At date $t$ there are $Y_t$ young and $O_t$ old households.
A young household has liquid wealth $b_i>0$, income $y_i^y>0$ when young, and
income $y_i^o\ge0$ when old. The triple $(y_i^y,b_i,y_i^o)$ has distribution $F$
in each entering cohort. Future income is known. Households choose completed
fertility once and retain their tenure when old.

\textbf{Preferences.} Each child uses $\chi>0$ units of goods and $\kappa>0$
units of space. Total nondurable expenditure is $c$ and housing is $h$; the
bundles left for adults are:
\begin{equation}
x=c-\chi n>0,\qquad s=h-\kappa n>0.
\end{equation}
Young and old utility are:
\begin{align}
u_t^y(c,h,n)&=\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta_t\log n,
\label{eq:new_uy}\\
u^o(c^o,h^o,e)&=\log c^o+\gamma\log h^o+\omega_B\log e.
\label{eq:new_uo}
\end{align}
All preference weights are positive. Households discount future utility by
$\beta>0$. Fertility $n>0$ is continuous; the superscript $o$ labels old-age
quantities. The estate $e$ consists of financial assets and the proceeds from
selling housing at death. It enters the parent's utility but does not determine
an entering household's wealth $b_i$.

A household draws an ownership preference $\xi_i$ before making its choices.
The draw is independent of its endowments and logistic with location $\bar\xi$
and scale $\sigma_\xi>0$.

\textbf{Housing and financial markets.} The housing stock is $\bar H>0$.
Rental units satisfy $h\le h_R^{\max}$ and owner units satisfy
$h\le h_O^{\max}$, where $0<h_R^{\max}<h_O^{\max}$. Goods and bonds trade with
the rest of the world. The bond price is $q\in(0,1)$ and its gross return is
$R_f=1/q$. House prices are $P_t>0$. Owners pay property tax at rate $\tau^p$.
Competition among rental intermediaries gives:
\begin{equation}
u_t\equiv qr_t=(1+q\tau^p)P_t-qP_{t+1}>0.
\label{eq:new_usercost}
\end{equation}
Here $r_t$ is rent paid at the end of the period, and $u_t$ is its value at
the beginning. Property-tax revenue is rebated equally to young and old
households; $T_t$ denotes the rebate at the beginning of the period.

\textbf{Home finance.} Current income and liquid wealth are available at
purchase. A mortgage finances at most the share $\phi_t\in(0,1)$ of the house
price. Principal and accumulated interest are repaid on entering old age.
Unsecured borrowing is unavailable. Old owners can buy or sell housing within
the owner size limit, using their income and wealth without new borrowing.
Thus an old owner's choice is not bounded by the size of its previous home.

\textbf{Demography.} A child produces $\nu>0$ young households next period,
after survival and household formation. Cohort masses satisfy:
\begin{equation}
Y_{t+1}=\nu\bar n_tY_t,\qquad O_{t+1}=Y_t,
\label{eq:new_demography}
\end{equation}
where $\bar n_t$ is average fertility among the young. Population here counts
adult households.

\section{Household choices}

\textbf{Young renters.} Financial wealth on entering old age is $a'$.
Conditional on renting, the young household solves:
\begin{align}
W_t^R(i)=\max_{c,h,n,a'}\;&u_t^y(c,h,n)+\beta V_{t+1}^R(a';i),\nonumber\\
&c+qa'+u_th=y_i^y+b_i+T_t,\qquad a'\ge0,\quad h\le h_R^{\max}.
\label{eq:new_young_renter}
\end{align}
The renter pays for current housing services and saves for old age.

\textbf{Young owners.} For an owner, $a'$ is financial wealth net of mortgage
repayment. Its problem is:
\begin{align}
W_t^O(i)=\max_{c,h,n,a'}\;&u_t^y(c,h,n)+\beta V_{t+1}^O(a',h;i),\nonumber\\
&c+qa'+(1+q\tau^p)P_th=y_i^y+b_i+T_t,\nonumber\\
&qa'+\phi_tP_th\ge0,\qquad h\le h_O^{\max}.
\label{eq:new_young_owner}
\end{align}
To see the mortgage directly, let $d$ be principal borrowed at purchase and
$k$ the amount invested in bonds. Then $0\le d\le\phi_tP_th$, $k\ge0$,
$a'=(k-d)/q$, and the budget is
$c+k+(1+q\tau^p)P_th=y_i^y+b_i+T_t+d$.
Eliminating $k,d$ gives \eqref{eq:new_young_owner}.

\textbf{Old households.} An old owner has net financial wealth $a$ and title
to $H$ units of housing. Let $a^e\ge0$ be financial saving during old age.
The estate is:
\begin{equation}
e=\begin{cases}
q^{-1}a^e,&\text{renter},\\
q^{-1}a^e+P_{t+1}h^o,&\text{owner}.
\end{cases}
\label{eq:new_estate}
\end{equation}
The retained house is sold at death, at the end of old age. The two old-age
problems are:
\begin{align}
V_t^R(a;i)=\max_{c^o,h^o,e}\;&u^o(c^o,h^o,e),\nonumber\\
&c^o+qe+u_th^o=a+y_i^o+T_t,\quad h^o\le h_R^{\max},
\label{eq:new_old_renter}\\
V_t^O(a,H;i)=\max_{c^o,h^o,e}\;&u^o(c^o,h^o,e),\nonumber\\
&c^o+qe+u_th^o=a+P_tH+y_i^o+T_t,\nonumber\\
&h^o\le h_O^{\max},\qquad e\ge P_{t+1}h^o.
\label{eq:new_old_owner}
\end{align}
The owner's estate restriction is equivalent to $a^e\ge0$. Its cash budget
before substituting for the estate is
$c^o+a^e+(1+q\tau^p)P_th^o=a+P_tH+y_i^o+T_t$.
Income, liquid wealth, and sale proceeds therefore enter the same budget.

\textbf{Tenure.} The household owns if $W_t^O(i)+\xi_i\ge W_t^R(i)$.
Its ownership probability is:
\begin{equation}
\pi_t^O(i)=
\frac{\exp\{[W_t^O(i)+\bar\xi]/\sigma_\xi\}}
{\exp\{W_t^R(i)/\sigma_\xi\}+\exp\{[W_t^O(i)+\bar\xi]/\sigma_\xi\}}.
\label{eq:new_tenure}
\end{equation}
The taste draw affects tenure but not choices conditional on tenure.

\Needspace{10\baselineskip}
\section{Equilibrium and the welfare comparison}

\begin{definition}
An equilibrium consists of household choices, tenure probabilities, prices,
rebates, and cohort masses satisfying the household problems, rental pricing,
and demographic equations above. Housing and the property-tax budget clear:
\begin{equation}
Y_t\bar h_t^y+O_t\bar h_t^o=\bar H,
\qquad (Y_t+O_t)T_t=q\tau^pP_t\bar H.
\label{eq:clearing}
\end{equation}
Here $\bar h_t^y$ and $\bar h_t^o$ are average housing occupied by young and old households.
The old distribution is generated by the preceding cohort's choices.
\end{definition}
At a positive stationary equilibrium, $Y=O=N$, $\bar n=1/\nu$, and
$N=\bar H/(\bar h^y+\bar h^o)$.

For the stationary results assume $0\le\tau^p<2$. At stationarity, let
$w=y^y+b+T$ be current cash and $v=y^o+T$ be old income
including the rebate. Write:
\begin{equation}
\begin{gathered}
p=(1-q+q\tau^p)P,\qquad L=(1-\phi+q\tau^p)P,\\
z=a'+Ph+v,\qquad K=1+\gamma+\omega_B,\qquad
D=1+\alpha+\vartheta+\beta K.
\end{gathered}
\label{eq:reduction_objects}
\end{equation}
Here $p$ is the cost of housing services, $L$ is the cash required per unit
of owner housing, and $z$ is the owner's resources on entering old age.
The young owner's two financial restrictions become:
\begin{equation}
c+ph+qz=w+qv,\qquad c+Lh\le w.
\label{eq:reduced}
\end{equation}
The second inequality is the original mortgage limit. If old housing and the
estate restriction are slack, $V(z)=K\log z+C$, with
$c^o=z/K$, $h^o=\gamma z/(Kp)$, and $e=\omega_Bz/(Kq)$.
```

## E. Existing sufficient stationary existence certificate: source excerpt

This excerpt is a proof to independently audit. Its reference to the original numbered Proposition does not adopt that older proposition. No claim of general uniqueness or transition stability accompanies it.

```latex
The following sufficient conditions produce the equilibrium and constrained
owner group in Proposition~\ref{prop:direct}. They are conservative analytical
bounds, with no equilibrium price or multiplier as an input. Assume bounded
endowments, $w_{0i}=y_i^y+b_i\ge\underline w>0$, $v_{0i}=y_i^o\ge0$, and
$0\le\tau^p<2$. Define:
\begin{equation}
\begin{gathered}
d_p=1-q+q\tau^p,\quad d_L=1-\phi+q\tau^p,\quad
d_{\max}=\max\{d_p,d_L\},\\
\bar M_0=\int(w_0+qv_0)\,\mathrm dF,\quad
\bar T=\frac{\tau^p\bar M_0}{(1-q)(2-\tau^p)},\quad
\bar M=\bar M_0+(1+q)\bar T.
\end{gathered}
\end{equation}
Require enough potential fertility when housing is inexpensive:
\begin{equation}
h_R^{\max}>\kappa/\nu,\qquad
\vartheta\nu>\chi D/\underline w+
\frac{\alpha\kappa}{h_R^{\max}-\kappa/\nu}.
\label{eq:existence_conditions}
\end{equation}
Then a positive stationary equilibrium exists, and every such equilibrium has
$0\le T\le\bar T$ and $P_-<P<P_+$, where:
\begin{equation}
P_-\equiv\frac{\alpha\underline w}{D d_{\max}h_R^{\max}},
\qquad P_+\equiv\frac{\nu\bar M}{d_p\kappa}.
\label{eq:price_bounds}
\end{equation}

To verify existence, combine the original budgets in either tenure:
\begin{equation}
x+ps+(\chi+p\kappa)n+qc^2+qp h^2+q^2e=w+qv.
\label{eq:lifetime}
\end{equation}
Scaling the first-order conditions gives
$\lambda(w+qv)+\mu w\le D$, while $1/x=\lambda+\mu$; hence
$x\ge w/D\ge\underline w/D$. Below $P_-$, the housing condition
$\alpha x/s\le d_{\max}P$ when uncapped, or the binding cap itself,
implies $h\ge h_R^{\max}$. Equation \eqref{eq:new_fertility} and
\eqref{eq:existence_conditions} then give fertility above replacement for
all conditional choices. Above $P_+$, \eqref{eq:lifetime} gives mean
fertility below replacement for every $T\le\bar T$.

The rebate map, writing $\bar h=\bar h^y+\bar h^o$, satisfies:
\begin{equation}
\mathcal T(P,T)=q\tau^pP\bar h/2
\le\frac{\tau^p}{2d_p}\{\bar M_0+(1+q)T\}\le\bar T.
\end{equation}
The final inequality uses
$2d_p-\tau^p(1+q)=(1-q)(2-\tau^p)>0$.
Conditional allocations and tenure probabilities are continuous. On the
compact rectangle of prices and rebates, move the price in the direction of
$\nu\bar n-1$, clipping it to $[P_-,P_+]$, and update the rebate by
$\mathcal T$. A Brouwer fixed point has an interior price by the strict
boundary signs, and therefore replacement fertility. Setting
$N=\bar H/\bar h$ completes the stationary equilibrium. This argument
establishes existence, not uniqueness or transition stability.
```
