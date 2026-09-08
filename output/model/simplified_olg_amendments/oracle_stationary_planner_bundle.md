# Full planner with fixed fertility: settle the benchmark before the theorem

You are reviewing a simple, illustrative two-age OLG model of housing and
fertility. We need a precise analytical foundation for a paper and a few theory
slides. Previous iterations have produced conditional housing variations but
have not specified the full planner satisfactorily. Resolve that problem first.
Do not simply improve the wording of the attached proposition.

The author wants to compare a competitive steady state with a planner's
stationary allocation, initially holding fertility fixed. The planner chooses
consumption as well as housing and the other components of the allocation.
The hoped-for economic conclusion is that borrowing limits and differences in
the housing available to renters and owners can leave too little housing with
young households. This is a hypothesis to establish or qualify, not an answer
you must obtain. We have spent too many iterations selecting assumptions to
rescue a predetermined conclusion.

Please give one recommended formulation and at most one substantive welfare
alternative. The author is willing to compare welfare weights and generally
prefers conventions established in dynastic OLG. Explain what that preference
does and does not imply for this particular household model. Preserve the model
unless a change is necessary; identify every proposed change explicitly.

## 1. Author instructions and scope

- First compare stationary allocations with the SAME positive masses
  \(Y=O=N\), housing stock, income distribution, and aggregate resources.
  Fix the fertility of each household type at its competitive choice, not just
  average fertility. A type includes endowments and the ownership taste draw.
  If the planner changes tenure, its fertility remains the frozen individual
  choice. This is a conditional welfare exercise, not a claim that fertility
  would still be privately optimal after reallocation.
- Stationarity requires \(\nu\bar n=1\), but this does not determine \(N\).
  Take \(N\) from the reference competitive steady state. Do not let the
  planner get more housing per household by choosing a smaller population.
- The direct planner chooses the FULL feasible allocation: young and old
  consumption, housing, permitted tenure assignments, saving and estates (or
  their correctly defined real counterparts). It can redistribute resources
  and relax individual financing restrictions. Physical housing limits and
  tenure persistence across a household's two ages remain. Assume full
  information for this direct benchmark, including realized ownership tastes.
  Fixing tenure or consumption in a proof variation is permissible; defining
  the entire planner by that restricted variation is not.
- Welfare weights must be justified independently of the desired housing
  direction. Distinguish failure to maximize the chosen utilitarian objective
  from Pareto inefficiency and from constrained inefficiency.
- The two-period structure, heterogeneous young wealth and young/old income,
  housing size limits by tenure, current income available for purchase, standard
  mortgage repayment, positive child goods and space costs, and warm-glow
  estate preferences are the maintained starting point. Retain the author's
  notation. Do not remove heterogeneity, assume a very low household discount
  factor, add liquidation costs, or replace the mortgage to manufacture a sign.
- THIS RUN DOES NOT need a transfer-policy theorem, endogenous-fertility
  welfare, a demographic transition proof, numerical examples, calibration,
  plots, or slides. Those come after this benchmark is settled. Explain how a
  stationary comparison differs from welfare along an attainable transition,
  but do not spend the run solving the transition.

## 2. What the attachment establishes, and what it does not

The attached `simplified_olg_utilitarian.tex` is the exact current discussion
note. Its environment and household equations are the maintained specification.
Its welfare section and propositions use an earlier, dated remaining-lifetime
criterion; they are material to reassess, not an authoritative definition of
the new full planner. Later transfer and fertility results are background only.

Each household lives young and old. Lifetime utility, including the tenure
taste once, is
\[
u^y(c,h,n)+\xi\mathbf1\{\mathrm{owner}\}+\beta u^o(c^2,h^2,e),
\]
\[
u^y=\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta\log n,
\qquad
u^o=\log c^2+\gamma\log h^2+\omega_B\log e.
\]
The note's superscript \(2\) denotes old-age quantities; \(y,o\) label utilities.
The endowment triple \((y^y,b,y^o)\) has distribution \(F\), and the independent
ownership taste \(\xi\) is logistic. Future income is known. Use integrals or
bars for aggregation, not expectations suggesting aggregate uncertainty.

The estate enters parental utility directly. There is NO continuation utility
of children and NO equation making entrant wealth \(b\) inherit the estate.
Thus this is currently warm-glow OLG, not a Barro–Becker dynasty. If genuine
dynastic altruism is essential to a suggested benchmark, show the extra
preference and inheritance equations and label this a different model. Do not
silently replace \(\omega_B\log e\) by descendants' utility.

Goods and bonds trade with the outside world at bond price \(q\in(0,1)\).
At stationarity let
\[
p=(1-q+q\tau^p)P,\quad L=(1-\phi+q\tau^p)P,
\quad w=y^y+b+T,\quad v=y^o+T,
\]
\[
z=a'+Ph+v,\qquad K=1+\gamma+\omega_B.
\]
Here \(p\) is the housing service cost, \(L\) the cash requirement per owner
unit, and \(z\) resources on entering old age. The young owner's constraints
reduce exactly to
\[
c+ph+qz=w+qv,\qquad c+Lh\le w.
\]
Thus current income CAN finance the down payment. The original constraints,
renter problems, estate equation, rebates and housing clearing are attached.

## 3. Resolve welfare weights first

Use primary literature to distinguish the objects below; do not call all of
them interchangeable versions of utilitarian welfare.

1. A stationary cohort's lifetime welfare, integrating
   \(u^y+\xi\mathbf1\{\mathrm{owner}\}+\beta u^o\) across the cohort.
2. Welfare over an infinite sequence of cohort lifetimes, with explicit social
   generation weights, a social discount factor distinct from \(\beta\), and
   an explicit treatment of the initial old. Derive the relative age weights
   in its dated first-order conditions. Explain why maximizing welfare only
   over stationary allocations need not give the steady state approached by
   this dynamic planner.
3. Welfare of initial dynasties when parental utility recursively includes
   descendants. Explain separately the social weights across dynasties and
   the private altruism weights within each dynasty.

Recommend ONE primary criterion for our immediate stationary comparison and
ONE useful comparator if needed. Tell us whether the recommendation preserves
the household model. Do not introduce a new dynasty merely because the author
asked for literature-consistent weighting. With fixed population, state when
total and average welfare differ only by a constant normalization.

The current note gives weight one to each living household's remaining utility:
young lifetime utility plus current old utility. At stationarity its old flow
terms consequently have combined coefficient \(1+\beta\), whereas a single
cohort lifetime criterion has coefficient \(\beta\). Neither criterion should
be conflated with equal weights on every cohort along a whole transition.

Useful PRIMARY starting points, to verify rather than invoke as authority:

- Becker and Barro, *A Reformulation of the Economic Theory of Fertility*,
  QJE (1988); working-paper version: https://www.nber.org/papers/w1793 and
  https://www.nber.org/system/files/working_papers/w1793/w1793.pdf.
  Parental altruism generates a dynastic utility formulation.
- Golosov, Jones and Tertilt, *Efficiency with Endogenous Population Growth*,
  Econometrica (2007):
  https://tertilt.vwl.uni-mannheim.de/research/optimality_Econometrica.pdf.
  Section 3.4, especially Results 1–2 around printed pp. 1055–1056, distinguishes
  planning problems with weights on potential people versus initial agents.
  These results do not select a unique equal-weight social objective for us.
- Farhi and Werning, *Inequality and Social Discounting*, JPE (2007):
  https://www.journals.uchicago.edu/doi/10.1086/518741;
  author manuscript: https://web.mit.edu/iwerning/Public/inequality_social_screen_old.pdf.
  Check its separation of parental altruism and social concern for descendants.

Give precise equation/section references for the convention you actually use.
Report a lack of a unique convention honestly. We need a small, relevant
literature foundation, not a literature review.

## 4. Specify full feasibility instead of assuming it

An independent accounting check found that the current note's household
budgets alone do not uniquely specify the price-free planner. Please resolve:

- Who receives estates and who supplies the exogenous entrant wealth \(b\)?
  Internal bequests are transfers; they are not lost goods. Outside recipients
  imply external outflows. Identify which closure preserves the note.
- Who owns rental intermediaries and housing titles outside living owners?
  Their pricing condition alone does not specify ownership or their resource
  contribution. Preserve the aggregate housing stock; do not count purchases
  of existing housing as newly produced goods.
- What is the estate object in a direct allocation? The owner's definition
  \(e=q^{-1}a^e+P_{t+1}h^2\) is wealth including house resale value; the
  renter's estate is \(e=q^{-1}a^e\). A planner
  cannot choose arbitrary nominal house prices to create utility/resources.
  Explain a consistent real delivery/valuation rule. If reference-price
  valuation is retained, say exactly what that means for the benchmark.
- Which individual financial constraints does the planner relax? The owner's
  \(e\ge P_{t+1}h^2\) is equivalent to nonnegative financial saving when old;
  it is not itself a physical housing limit. Do not retain it as physical
  feasibility without explanation.
- What fixes aggregate external wealth? Do not give the planner an arbitrary
  asset endowment. Individual saving and wealth distribution may change while
  the aggregate endowment is held fixed. Distinguish that restriction from
  requiring an attainable path from a common inherited state.

For orientation ONLY, if all domestic claims can be consolidated, a candidate
dated ledger is
\[
C_t+qB_{t+1}^{\mathrm{ext}}+X_t
=Y_t^g+B_t^{\mathrm{ext}}+I_t,
\]
where \(C_t\) is total nondurable expenditure (already including child goods),
\(Y_t^g\) goods income, \(B_t^{\mathrm{ext}}\) external bond payoffs available
at \(t\), and \(I_t,X_t\) actual outside inflows and outflows. This is NOT a
closed specification until recipients, ownership and timing are fixed. At
stationarity it implies
\(C+X=Y^g+I+(1-q)B^{\mathrm{ext}}\).
Freely choosing endowed \(B^{\mathrm{ext}}\) can make welfare unbounded.
With domestic taxpayers and rebate recipients, fully rebated property taxes
are internal transfers, not resource losses. Any foreign titleholders' tax
payments must be recorded consistently with outside rental/title cash flows.

Recommend the smallest explicit completion consistent with the household
model, label its added assumptions, and then write the full planner with all
its controls, objective, feasibility constraints, and fixed objects. If no
such completion preserves the claimed economy, identify the precise conflict
and the smallest proposed revision. Do not hide it in a theorem's hypotheses.

## 5. Analytical result and a diagnostic that needs interpretation

On the uncapped old-owner branch with positive financial estate,
\[
c^2=z/K,\quad h^2=\gamma z/(Kp),\quad e=\omega_Bz/(Kq).
\]
This branch requires
\[
\omega_B(1-q+q\tau^p)>q\gamma.
\]
It is a substantive restriction, not a harmless regularity condition.
With \(m=1/c^2\), \(s=h-\kappa n\), an uncapped young owner satisfies
\[
\alpha/s=(\beta/q)pm+\mu L,\qquad \gamma/h^2=pm,
\]
where \(\mu\ge0\) is the multiplier on \(c+Lh\le w\).

For the stationary cohort lifetime criterion, a candidate permanent housing
variation across otherwise identical young and old owner types gives
\[
\Delta h=\varepsilon,\quad\Delta h^2=-\varepsilon,
\quad\Delta a'=-P\varepsilon,\quad\Delta a^e=qP\varepsilon,
\]
with internal age transfers \(+p\varepsilon\) to young and
\(-p\varepsilon\) to old. Each entering old household then inherits an extra
\(\varepsilon\) housing and \(-P\varepsilon\) financial wealth. At reference
prices, entering-old total resources, estates and both ages' consumption stay
unchanged; \(q\Delta a'+\Delta a^e=0\). Repeated cohort by cohort, this is
a stationary candidate variation. It relaxes private financing when needed.

The resulting derivative appears to be
\[
\frac{\alpha}{s}-\beta\frac{\gamma}{h^2}
=\beta(1/q-1)pm+\mu L.
\]
Independently verify the algebra AND feasibility under your proposed complete
closure. It is positive for \(q<1\) even if \(\mu=0\). Interpret that fact:
it cannot by itself attribute the welfare difference to a binding mortgage.
Nor does it account for the initial old's transition welfare.

We need more than another sufficient condition saying that marginal utilities
are unequal. Seek a simple statement linking the model's parameters and
competitive choices to failure of the chosen FULL planner optimality
conditions. Wherever possible give interpretable primitive inequalities;
the attached note has existing income/capacity restrictions to inspect, not
automatically reuse under a new benchmark. Check the regular regime with zero
financial estate as well. If a positive-estate or cap restriction is necessary,
say how much of the result it limits and why.

Most importantly distinguish:
(a) an improving feasible housing variation; (b) the competitive allocation
not solving the full planner problem; (c) the full optimum allocating MORE
aggregate housing to young households after consumption, estates and tenure
are jointly optimized. Prove (c) if possible, but do not infer it from (a).
Separate ordinary redistribution under the welfare weights from the extra
effect of financing and tenure restrictions. If equal utility weights already
generate a gap without those restrictions, say so clearly. A frictionless
comparison with the SAME objective and resources is a useful diagnostic, not
an invitation to choose convenient weights after observing the answer.

## 6. Required output

Give a short verdict, the recommended welfare criterion and one comparator,
the fully specified planner, and one main analytical proposition with proof.
Use a compact appendix only for essential algebra or accounting. End with the
few author decisions that actually remain, separating normative choices from
missing accounting assumptions and changes to household behavior. Keep the
main answer roughly within 2,500 words, with essential proof detail beyond
that if necessary. Use plain economics prose and explicit equations. Take the
reasoning time needed to check the accounting and theorem; verbosity is not
the objective. No numerical-point-plus-open-neighborhood proof, no claims of
universal inefficiency without proof, and no promise of an unproved fertility
effect. A clear limitation is preferable to another superficially complete
formulation that requires restarting tomorrow.

### File: latex/JMP_DS_suggestions/simplified_olg_utilitarian.tex
```
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern,microtype,amsmath,amssymb,amsthm,booktabs,tabularx,enumitem,needspace}
\usepackage[colorlinks=true,linkcolor=blue,urlcolor=blue]{hyperref}
\newtheorem{proposition}{Proposition}
\newtheorem{definition}{Definition}
\newcommand{\dd}{\mathrm d}
\begin{document}
\begin{center}
{\Large Housing Allocation, Transfers, and Fertility}\par
\smallskip
{\large An illustrative model}\par
\smallskip
{\small Discussion draft for Tommaso De Santo \quad September 8, 2026}
\end{center}
\medskip

Borrowing limits can keep housing with old owners even when young households
would gain more from the space. We compare the market allocation with an
equally weighted utilitarian benchmark. We then ask whether taxes and transfers
can improve that allocation while households continue to choose through markets,
and how the resulting changes in housing and consumption affect fertility.

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
u^o(c^2,h^2,e)&=\log c^2+\gamma\log h^2+\omega_B\log e.
\label{eq:new_uo}
\end{align}
All preference weights are positive. Households discount future utility by
$\beta>0$. Fertility $n>0$ is continuous; the superscript 2 labels old-age
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
q^{-1}a^e+P_{t+1}h^2,&\text{owner}.
\end{cases}
\label{eq:new_estate}
\end{equation}
The retained house is sold at death, at the end of old age. The two old-age
problems are:
\begin{align}
V_t^R(a;i)=\max_{c^2,h^2,e}\;&u^o(c^2,h^2,e),\nonumber\\
&c^2+qe+u_th^2=a+y_i^o+T_t,\quad h^2\le h_R^{\max},
\label{eq:new_old_renter}\\
V_t^O(a,H;i)=\max_{c^2,h^2,e}\;&u^o(c^2,h^2,e),\nonumber\\
&c^2+qe+u_th^2=a+P_tH+y_i^o+T_t,\nonumber\\
&h^2\le h_O^{\max},\qquad e\ge P_{t+1}h^2.
\label{eq:new_old_owner}
\end{align}
The owner's estate restriction is equivalent to $a^e\ge0$. Its cash budget
before substituting for the estate is
$c^2+a^e+(1+q\tau^p)P_th^2=a+P_tH+y_i^o+T_t$.
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
$c^2=z/K$, $h^2=\gamma z/(Kp)$, and $e=\omega_Bz/(Kq)$.

For the direct housing comparison, hold fertility and tenure fixed at their
market choices. Give each living household weight one on its remaining
utility:
\begin{equation}
\mathcal W_t=\int\bigl[u_t^y+\beta V_{t+1}
 +\xi_i\mathbf 1\{\text{owner}\}\bigr]\,\dd M_t^Y
 +\int u_t^o\,\dd M_t^O.
\label{eq:welfare}
\end{equation}
The measures count households. Thus equal weights do not mean equal weights
on the two age-group averages. This is the proposed welfare convention for
this note. Equal weights on utility as measured at birth would instead put
$\beta$ on the initial old; Appendix~\ref{app:weights} explains the difference.

The direct planner can assign housing and waive individual borrowing limits,
while respecting the housing stock, tenure-specific size limits, estates, and
aggregate resources. A transfer authority instead announces lump-sum payments
and lets households choose through markets under the original restrictions.
We first compare allocations that leave all later cohorts unchanged, so their
welfare cancels. Neither comparison attaches a social value to additional births.

\section{Housing allocation}

The relevant young owners have a strictly binding mortgage limit. They are
below their housing cap, and their old-age housing cap and estate floor are
slack. These properties follow from explicit
income and capacity inequalities in Appendix~\ref{app:primitives}; they need
not be assumed from a computed equilibrium.

\begin{proposition}[Direct utilitarian improvement]
\label{prop:direct}
Consider a positive stationary equilibrium with a positive share of these
owners. If $\beta\ge q$, a small reallocation of housing from their current
old counterparts to the young raises $\mathcal W_t$. Consumption, fertility,
tenure, estates, all future real allocations, and net external assets can
remain unchanged.
\end{proposition}

\begin{proof}
Let $m=1/c^2$ be the marginal utility of cash when old, and let $\mu>0$
be the multiplier on the young owner's cash restriction in
\eqref{eq:reduced}. The young and old first-order conditions imply:
\begin{equation}
\frac\alpha s=\frac\beta q\,pm+\mu L,
\qquad \frac\gamma{h^2}=pm.
\end{equation}
Stationarity supplies an old owner with the same endowments and old resources.
The difference in housing marginal utilities is therefore:
\begin{equation}
\boxed{\quad
\frac\alpha s-\frac\gamma{h^2}
=\left(\frac\beta q-1\right)pm+\mu L>0.
\quad}
\label{eq:direct_gap}
\end{equation}
Giving the young $\varepsilon$ units and taking the same amount from the old
changes welfare by
$\alpha\log(1+\varepsilon/s)+\gamma\log(1-\varepsilon/h^2)$, which is
positive for small $\varepsilon>0$. Appendix~\ref{app:direct_settlement}
verifies the financial settlement. Stationary young and old owner measures
coincide, so the comparison integrates over heterogeneous endowments.
\end{proof}

Young households value the same space more because current resources are
scarce relative to their resources when old. At $\beta=q$, the entire gap in
\eqref{eq:direct_gap} is due to the mortgage constraint. At $\beta>q$, the
household's relative valuation of the two ages also contributes under the
chosen social weights. The old lose housing utility; compensation is not
required by this utilitarian criterion.

This proves that the market allocation does not maximize the planner's
objective. It also identifies an improving direction toward young housing.
It need not identify the direction of the global optimum when other
households and binding physical caps are included. For the conditional
housing optimum, hold consumption, children, and estates fixed, and write
$C_n=\kappa\int n_i\,\dd M^Y$ for space used by children. If the unconstrained
allocation respects all retained caps and estate floors, young aggregate
housing is:
\begin{equation}
H_Y^*=C_n+\frac{\alpha Y}{\alpha Y+\gamma O}(\bar H-C_n).
\label{eq:planner_housing}
\end{equation}
At stationarity, this exceeds equilibrium young housing if all old housing
caps and estate floors are slack and the hypotheses of
Proposition~\ref{prop:direct} hold. Appendix~\ref{app:direct_settlement}
gives the aggregate comparison and the solution with active caps.

Along a transition, stationarity can no longer identify current old marginal
utility with the young household's future old marginal utility. The exact
local test becomes:
\begin{equation}
p_t\left(\frac\beta q m_{i,t+1}-m_{j,t}\right)+\mu_{i,t}L_t>0,
\label{eq:dated_gap}
\end{equation}
where $p_t=u_t$ and $L_t=(1-\phi_t+q\tau^p)P_t$. This applies at a date on
any existing path; it requires no convergence assumption. The stationary
parameter condition alone does not sign it along every path.

\section{Taxes and transfers through markets}

A grant can relax the young household's cash shortage. A tax on an old owner
can release housing. Their demands must still clear the market, and the
program must finance both current and promised payments.

There is a simple constrained utilitarian comparison that preserves the
original timing of all young choices.

\begin{proposition}[Redistribution among old owners]
\label{prop:old_redistribution}
Suppose a positive stationary equilibrium has two positive groups of old
owners with slack housing caps and estate floors. If every household in one
group has more old-age resources than every household in the other, a small,
unanticipated one-time balanced cash transfer from the richer group to the
poorer raises utilitarian welfare. Every household chooses optimally under the original constraints,
and the original aggregate housing, fertility, and price path are preserved.
\end{proposition}

\begin{proof}
At the same prices, old demands are proportional to resources:
\begin{equation}
(c^2,h^2,e)=z(1/K,\gamma/(Kp),\omega_B/(Kq)).
\end{equation}
A balanced transfer therefore preserves all three aggregates. Select equal
positive submeasures from the two groups with uniform strict margins. For a
matched unit of mass with $z_L<z_H$, transferring $b$ raises welfare at the rate:
\begin{equation}
K\left[\frac1{z_L+b}-\frac1{z_H-b}\right]>0
\qquad\text{for }0<b<(z_H-z_L)/2.
\end{equation}
Choose a smaller transfer if needed to preserve the slack caps. It is a
one-time intervention for current old households, so current young face the
same prices and lifetime resources and retain all their choices. Estates do
not determine entrant wealth, so later cohorts are also unchanged.
\end{proof}

The gain comes from redistributing wealth within old age. Young housing is
unchanged. The condition on resource dispersion also has a simple primitive
case: when $\phi=q$, owners with binding mortgages enter old age with
$z=y^o+T$. Two separated old-income groups satisfying
Appendix~\ref{app:primitives} therefore suffice. This result needs no
restriction on $\beta/q$.

Moving housing toward young households through transfers requires more.
The following construction is conditional on fertility and tenure having already
been chosen. One interpretation is an unexpected intervention after those
commitments and before housing purchase. This additional timing restriction
is explicit: it is not yet a theorem for the original simultaneous fertility,
tenure, and housing choices.

\begin{proposition}[A transfer improvement with committed fertility]
\label{prop:transfers}
Start from the stationary equilibrium in Proposition~\ref{prop:direct}.
Suppose:
\begin{equation}
\boxed{\quad
\phi\ge q,\qquad \alpha(1+\omega_B)\le
\gamma\frac{1-\phi+q\tau^p}{1-q+q\tau^p}.
\quad}
\label{eq:transfer_conditions}
\end{equation}
Suppose also that a positive group of other old owners is strictly at its
housing cap. A small program of current and next-period grants to selected young
owners, financed by taxes on current old owners, raises utilitarian welfare
and moves housing toward the young. Households optimize conditional on their
committed fertility and tenure; all private financial restrictions remain.
The government saves current revenue to pay its next-period grants.
\end{proposition}

The preference restriction makes the grants' remaining funding requirement
nonnegative after taxing the old owners who release housing. The additional
wealthy old owners cover that requirement without changing their housing,
because their cap remains binding. Their marginal utility of cash is
automatically lower than that of the uncapped old donors:
\begin{equation}
m_F<\frac\gamma{p h_O^{\max}}<m_i.
\label{eq:cap_marginal_order}
\end{equation}
Appendix~\ref{app:transfers} gives an income bound producing this group and proves that the restrictions can hold
jointly. This construction includes redistribution from households with low
marginal utility of consumption; its gain cannot all be attributed to housing.

Here is the policy for one selected young owner. Put $\lambda=\beta K/(qz)$.
For a small increase $\varepsilon$ in adult space, choose:
\begin{equation}
s_\varepsilon=s+\varepsilon,\qquad
x_\varepsilon=\frac L{\alpha/(s+\varepsilon)-(p-L)\lambda}.
\label{eq:transfer_allocation}
\end{equation}
The current grant $G$ and next-period grant $J$ are:
\begin{equation}
G=x_\varepsilon-x+L\varepsilon,
\qquad qJ=(p-L)\varepsilon.
\label{eq:transfer_grants}
\end{equation}
Both are nonnegative. These are fixed payments, calculated from the baseline,
not subsidies that vary with the household's subsequent choice. The household
voluntarily chooses \eqref{eq:transfer_allocation}, remains at its mortgage
limit, and enters old age with its original resources $z$.

An uncapped old owner's housing demand is $gz$, where $g=\gamma/(Kp)$.
Tax its matched old counterpart $\varepsilon/g$, so that housing falls by
exactly $\varepsilon$. The wealthy capped group pays the residual:
\begin{equation}
R=G+qJ-\varepsilon/g=x_\varepsilon-x+p\varepsilon-\varepsilon/g\ge0.
\label{eq:residual}
\end{equation}
Thus housing clears at its original price. Current net revenue is $qJ$;
saving it pays $J$ next period. The original old-age resources and choices
of the young are preserved, and later cohorts face the original economy.
The initial old leave smaller estates, which are included in their welfare.

This establishes one feasible transfer improvement. It does not characterize
a uniform age-transfer rule or the optimal transfer schedule. Allowing
fertility to respond changes housing demand and future cohort masses.
Appendix~\ref{app:free_choice} gives a different construction that lets
individual fertility and tenure remain private choices while preserving total
births. The next section distinguishes these welfare comparisons from a
policy that raises aggregate fertility.

\section{Fertility on the adjustment path}

Children use both goods and space. The fertility condition is:
\begin{equation}
\frac\vartheta n=\frac\chi{c-\chi n}+\frac{\alpha\kappa}{h-\kappa n}.
\label{eq:new_fertility}
\end{equation}
Additional housing lowers the space cost of another child; a fall in
consumption raises its goods cost. This distinction matters even when the
household gains utility.

\begin{proposition}[A gift and a dated fertility comparison]
\label{prop:fertility}
At fixed tenure, prices, and old-age transfers, a current cash gift strictly
raises fertility. Housing rises until its physical cap binds. The result
allows either mortgage regime and all old-age housing and estate constraints.
The original pair of current and old-age grants in
\eqref{eq:transfer_grants} also strictly raises a treated owner's fertility
at the fixed price-and-rebate path, after all its choices are reoptimized.

More generally, compare two privately chosen bundles at the same fertility
preference $\vartheta$. If the new bundle can accommodate baseline fertility
$n_0$, then:
\begin{equation}
\boxed{\quad
n_1\ge n_0\quad\Longleftrightarrow\quad
\frac\chi{c_1-\chi n_0}+\frac{\alpha\kappa}{h_1-\kappa n_0}
\le \frac\chi{c_0-\chi n_0}+\frac{\alpha\kappa}{h_0-\kappa n_0}.
\quad}
\label{eq:finite_fertility}
\end{equation}
This test applies at a date on two existing equilibrium paths, without
requiring either path to converge.
\end{proposition}

A rise in both consumption and housing is sufficient for higher fertility.
If one falls, \eqref{eq:finite_fertility} states exactly how much the other
must compensate. Appendix~\ref{app:fertility} proves the gift result across
constraint regimes. A gift paid for by taxing today's old is different from
a loan repaid by the young recipient: anticipated repayment can reverse its
fertility response.

\textbf{Tenure changes.} A common gift raises conditional fertility in both
tenures, but households may switch between them. For each endowment type,
mean fertility obeys the exact identity:
\begin{equation}
\Delta\bar n_i=\pi_1^O\Delta n_i^O+(1-\pi_1^O)\Delta n_i^R
 +(\pi_1^O-\pi_0^O)(n_{i0}^O-n_{i0}^R).
\label{eq:tenure_fertility}
\end{equation}
The last term can be negative. A grant available only to buyers raises
ownership at fixed prices; if, for every affected endowment type, owners initially choose at least
as many children as renters, both contributions are nonnegative. This is an
additional tenure comparison, not a consequence of utilitarian improvement.

\textbf{The same grant pair raises fertility.} Give the household exactly
the grants in \eqref{eq:transfer_grants}. Let $W_\varepsilon(n)$ be its
value after optimizing all other choices conditional on fertility $n$.
At baseline fertility $n_0$, the conditional optimum has
$x_\varepsilon>x_0$ and $s_\varepsilon>s_0$. The envelope theorem gives:
\begin{equation}
W_\varepsilon'(n_0)=
\chi\left(\frac1{x_0}-\frac1{x_\varepsilon}\right)
+\alpha\kappa\left(\frac1{s_0}-\frac1{s_\varepsilon}\right)>0.
\label{eq:grant_pair_fertility}
\end{equation}
The profiled value is concave, so the new optimal fertility is strictly
higher. This includes the promised old-age grant and all household
reoptimization. It does not require consumption and housing both to rise
at the final, freely chosen fertility.

This positive treatment effect does not preserve the market implementation.
For example, when $\phi=q$ and $\gamma=\alpha(1+\omega_B)$, the balanced
transfer in Appendix~\ref{app:transfers} clears housing at fixed fertility.
With fertility freely chosen, the difference between young and old housing
responses to cash instead is:
\begin{equation}
\frac{\vartheta(\kappa p-\alpha\chi)}
{p(\chi+\kappa p)(1+\alpha)(1+\alpha+\vartheta)}.
\end{equation}
It can have either sign. Even if it is zero, higher fertility changes the
next cohort's housing demand. A policy that raises aggregate fertility
therefore needs a new demographic equilibrium path. The separate construction
in Appendix~\ref{app:free_choice} avoids fixing individual fertility by
offsetting births within the young cohort; total fertility remains unchanged.

\textbf{The intended transition.} Start from a stationary economy. An
exogenous decline in $\vartheta_t$ initiates demographic adjustment. At a
later date, compare a housing policy with continued adjustment without it,
from the same inherited state and with the same preference shocks. The
initial fertility decline is part of the experiment, rather than something
the housing mechanism must explain.

For two existing paths sharing $Y_{t_0}$, cohort accounting gives:
\begin{equation}
\frac{Y_t^{\rm policy}}{Y_t^{\rm baseline}}
=\prod_{s=t_0}^{t-1}
\frac{\bar n_s^{\rm policy}}{\bar n_s^{\rm baseline}},
\qquad O_{t+1}^j=Y_t^j.
\label{eq:population_path}
\end{equation}
This orders populations at finite dates if fertility is ordered along the
path. Convergence is unnecessary for either this identity or the dated
fertility test. Determining equilibrium consumption, housing, and prices
still requires an equilibrium path, including the policy's fiscal response.

If both paths converge to positive stationary equilibria, both fertility
limits equal $1/\nu$. Their limiting adult-household populations satisfy:
\begin{equation}
\frac{(Y+O)^{\rm policy}}{(Y+O)^{\rm baseline}}
=\frac{(\bar h^y+\bar h^o)^{\rm baseline}}
{(\bar h^y+\bar h^o)^{\rm policy}}.
\label{eq:population_limit}
\end{equation}
\clearpage
\appendix
\section{Primitive conditions for the stationary result}
\label{app:primitives}

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

For the constrained owner group, put:
\begin{equation}
\ell=d_L/d_p,\qquad C_{\min}=1+\alpha\ell+
\vartheta\min\{1,\ell\},\qquad k=D/C_{\min}-1.
\end{equation}
On a positive-$F$-mass set $\mathcal S$, require:
\begin{equation}
qv_{0i}>kw_{0i}+\max\{k-q,0\}\bar T.
\label{eq:income_bound}
\end{equation}
For $M_{\mathcal S}\ge\sup_{i\in\mathcal S}
\{w_{0i}+qv_{0i}+(1+q)\bar T\}$, also require:
\begin{equation}
\omega_Bd_p>q\gamma,\qquad
h_O^{\max}>\frac{M_{\mathcal S}}{d_pP_-}
\max\{1,\gamma/(qK)\}.
\label{eq:cap_bound}
\end{equation}
The latter bounds make the old estate floor and both owner caps slack in
actual choices and in the comparison without the cash restriction.
In that comparison, optimal cash expenditure is:
\begin{equation}
\frac{w+qv}{D}\left[1+\alpha\frac Lp+
\vartheta\frac{\chi+L\kappa}{\chi+p\kappa}\right].
\end{equation}
The bracket is at least $C_{\min}$. Condition \eqref{eq:income_bound} makes
this expenditure strictly greater than $w$ for every admissible rebate.
Thus the original mortgage limit is strictly restrictive. The income
condition says that households cannot bring enough of their future resources
forward to finance their preferred current bundle. Logistic tastes give
these endowments positive owner mass; their stationary predecessors supply
the matched current old. Together with $\beta\ge q$, these primitive
inequalities imply Proposition~\ref{prop:direct}.

\section{Direct allocation: settlement and aggregate housing}
\label{app:direct_settlement}

At any date, increase a selected young owner's housing by $\varepsilon$ and
reduce the old donor's housing by the same amount. Keep consumption and
estates fixed. Set:
\begin{equation}
\Delta a'=-P_{t+1}\varepsilon,\qquad
\Delta a^e=qP_{t+1}\varepsilon,
\end{equation}
and transfer $p_t\varepsilon$ from old to young. The young budget's extra
cost is $q\Delta a'+(1+q\tau^p)P_t\varepsilon=p_t\varepsilon$.
The old budget releases the same amount. The old estate is unchanged since
$q^{-1}\Delta a^e-P_{t+1}\varepsilon=0$; the young household's future
resources are unchanged since $\Delta a'+P_{t+1}\varepsilon=0$.
Old saving finances young borrowing, so net external assets and existing
creditor payments are unchanged. The young mortgage restriction may be
violated, as permitted for the direct planner.

At stationarity, the same measure $N\pi^O(i)\,\dd F(i)$ describes both
owner cohorts. A positive eligible set has a positive subset with uniform
strict margins on the mortgage, caps, and log arguments. The finite log
welfare change therefore permits a common sufficiently small housing
transfer on this subset. No equal-mass assumption on unrelated income
groups is needed.

For the conditional housing optimum, let $\lambda_H>0$ be the stock
constraint's multiplier. The exact allocation with caps is:
\begin{equation}
s_i^*=\min\{h_{m_i}^{\max}-\kappa n_i,\alpha/\lambda_H\},
\qquad h_j^{2*}=\min\{\bar h_j^2,\gamma/\lambda_H\},
\end{equation}
where $\bar h_j^2=h_{m_j}^{\max}$ for renters and
$\bar h_j^2=\min\{h_O^{\max},e_j/P_{t+1}\}$ for owners with retained
estates. The multiplier makes total housing equal $\bar H$.
If no cap binds, this yields \eqref{eq:planner_housing}.

For the stationary aggregate comparison, let $Q$ be the common distribution
of endowments and realized tenure across households in either cohort. When every old
cap and estate floor is slack, the first-order conditions imply:
\begin{equation}
\alpha/s_i=\frac\beta q pm_i+\mu_iL_{m_i}+\eta_i^y
\ge pm_i=\gamma/h_i^o,
\end{equation}
where $L_O=L$, $L_R=p$, and $\eta_i^y\ge0$ is the young housing-cap
multiplier. If the uncapped planner allocation is feasible, then:
\begin{equation}
H_Y^*-H_Y^{\rm eq}=
\frac N{\alpha+\gamma}\int(\alpha h_i^o-\gamma s_i)\,\mathrm dQ>0.
\end{equation}
Strictness follows from the positive constrained-owner mass. At $\beta=q$
and with slack young caps, the integrand equals $\mu_iL_{m_i}s_i h_i^o$.
This is an aggregate statement for the stated regime. With active planner
caps, the preceding minimum formulas apply and this aggregate ordering has
not been established generally.

\section{Transfer implementation and funding}
\label{app:transfers}

Keep fertility and tenure committed. At the baseline let
$\Lambda=1/x$, $m=K/z$, and $\lambda=\beta m/q$. The first-order
conditions are $\Lambda=\lambda+\mu$ and
$\alpha/s=\lambda p+\mu L$, with $\mu>0$.
For \eqref{eq:transfer_allocation}--\eqref{eq:transfer_grants}, choose
$\Delta a'=-\phi P\varepsilon/q$. The original budget and mortgage give:
\begin{equation}
\Delta c+q\Delta a'+(1+q\tau^p)P\varepsilon=G,
\quad q\Delta a'+\phi P\varepsilon=0,
\quad \Delta a'+P\varepsilon+J=0.
\end{equation}
The updated choices satisfy every first-order condition. For sufficiently
small $\varepsilon>0$, all strict margins persist. Strict concavity of the
original log problem with affine constraints establishes the conditional
optimum, including its original estate restriction.

Along this direction, write $A=\dd x_\varepsilon/\dd\varepsilon$ and
$\rho_\varepsilon=\alpha x_\varepsilon/s_\varepsilon$. Then:
\begin{equation}
A=\frac{\rho_\varepsilon^2}{\alpha L}\ge\frac L\alpha
\ge\frac{(1+\omega_B)p}{\gamma}=1/g-p.
\end{equation}
The first inequality uses $L\le\rho_\varepsilon\le p$, since $\phi\ge q$.
The second is \eqref{eq:transfer_conditions}. Thus $R\ge0$.
Let $m_F$ be the funders' mean marginal utility of cash. For one unit of
recipient mass, the welfare derivative is:
\begin{equation}
\begin{aligned}
\mathcal W'(0)
&=\Lambda A+\alpha/s-m/g-m_F(A+p-1/g)\\
&=(\Lambda-m)A+(\alpha/s-pm)+(m-m_F)(A+p-1/g)>0.
\end{aligned}
\end{equation}
The first two terms are strictly positive when $\beta\ge q$ and $\mu>0$;
the third is nonnegative. Integrate over recipients and collect their total
residual equally across a positive funder group. Uniform strict margins
and continuity give a finite welfare gain, even when the two groups have
different masses.

All original mortgages are repaid. At date one the extra title held by the
selected young is sold when they choose their unchanged old housing; these
sales offset the smaller death sales of the initial old. The government's
saving pays $J$, and its balance is zero. From date two the inherited
household state is the baseline state. For each matched transfer, aggregate
consumption and the initial old's estate payments satisfy:
\begin{equation}
\Delta C_0=\frac{\omega_B\Delta x}{1+\omega_B},\qquad
\Delta E_1=-\frac{\omega_B\Delta x}{q(1+\omega_B)},\qquad
\Delta C_0+q\Delta E_1=0.
\end{equation}
Thus the extra consumption is financed by smaller estates, with the old's
loss included in welfare.

\textbf{A primitive funder bound.} Under Appendix~\ref{app:primitives}, put
$\delta=\phi/q-1\ge0$. A positive endowment group satisfying:
\begin{equation}
y_f^o>\delta P_+h_O^{\max}+\frac{Kd_pP_+h_O^{\max}}\gamma
\label{eq:funder_bound}
\end{equation}
has strictly capped old owner housing. Since its cap threshold is
$z_*=Kph_O^{\max}/\gamma$, its marginal utility is
$m_f=(1+\omega_B)/(z_f-ph_O^{\max})<\gamma/(ph_O^{\max})$.
An uncapped donor has $z_i<z_*$ and hence $m_i=K/z_i$ above that threshold,
which proves \eqref{eq:cap_marginal_order}. This follows from
$z_f\ge y_f^o-\delta P_+h_O^{\max}$ and the old log demands. The
price and rebate bounds use the full distribution, including the funders.

These restrictions are jointly possible. Choose $\phi>q$ close enough to
$q$ that $\ell>1/(1+q)$, and choose
$0<\alpha<\gamma\ell/(1+\omega_B)$ with
$\omega_Bd_p>q\gamma$. For large $\vartheta$,
$k\to\ell^{-1}-1<q$. Choose $h_R^{\max}>\kappa/\nu$, large enough
$\vartheta$, and then small positive $\chi$ to satisfy
\eqref{eq:existence_conditions}. Let ordinary types have $w_0=W$ and
$v_0=V_S$ with $qV_S>kW$. For fixed small $c_*>0$, give a funder group
old income $V>V_S+c_*$ and mass $c_* /(V-V_S)$, replacing that mass of
ordinary types. Then $\bar M_0=W+qV_S+qc_*$ is independent of $V$.
Choose the finite owner cap from \eqref{eq:cap_bound}. Every bound on the
right of \eqref{eq:funder_bound} is now independent of $V$, so sufficiently
large finite $V$ satisfies it with positive funder mass. The restrictions
are analytical sufficient conditions, not a claim of empirical mildness.

\textbf{A simpler special case.} If $\phi=q$ and
$\alpha(1+\omega_B)=\gamma$, then $L=p$ and $J=R=0$.
A current grant $b>0$ and an equal tax on the matched old suffice, without
additional funders. Put $w_n=w-(\chi+p\kappa)n$; young adult choices are
$x=w_n/(1+\alpha)$ and $s=\alpha w_n/[p(1+\alpha)]$. Both ages have
housing response $g=\alpha/[p(1+\alpha)]=\gamma/(Kp)$ to cash. Welfare
changes by:
\begin{equation}
\Delta\mathcal W(b)=(1+\alpha)\log(1+b/w_n)+K\log(1-b/v)>0
\end{equation}
for sufficiently small $b>0$, since strict finance and $\beta\ge q$ give
$(1+\alpha)/w_n>K/v$. This is a useful exact illustration, but it requires
two parameter equalities.

\section{Private fertility and the limits of the policy comparison}
\label{app:fertility}

At fixed prices, let $\mathcal V(z)$ retain every old-age constraint. The
reduced young problem is \eqref{eq:reduced}, with $L=p$ for renters. This
old log problem is strictly concave, continuous, and piecewise smooth.
When finance binds and young housing is uncapped, put $\Delta=L-p$ and
$\Gamma=-\beta\mathcal V''(z)/q^2>0$. The curvature terms for choices
$(h,n)$ are:
\begin{equation}
\begin{gathered}
A=L^2/x^2+\alpha/s^2+\Gamma\Delta^2,\quad
B_c=-L\chi/x^2+\alpha\kappa/s^2,\\
D_n=\chi^2/x^2+\alpha\kappa^2/s^2+\vartheta/n^2,
\qquad J_n=AD_n-B_c^2>0.
\end{gathered}
\end{equation}
Implicit differentiation with respect to current cash yields:
\begin{equation}
h_w=\frac{LD_n+\chi B_c}{x^2J_n}>0,
\qquad n_w=\frac{LB_c+\chi A}{x^2J_n}>0,
\end{equation}
because
$LD_n+\chi B_c=L\vartheta/n^2+\alpha\kappa(\chi+L\kappa)/s^2>0$
and $LB_c+\chi A=\alpha(\chi+L\kappa)/s^2+\chi\Gamma\Delta^2>0$.
If young housing is capped and finance binds,
$n_w=\chi/(x^2D_n)>0$.

If both young constraints are slack, expenditures on adult goods, adult
space, and children have fixed positive shares of total current expenditure;
that total rises with wealth by strict concavity of old utility. If only the
young housing cap binds, put $k_c=\chi/(x^2D_n)$. Then:
\begin{equation}
c_w=\frac\Gamma{\Gamma+(1-\chi k_c)/x^2}>0,
\qquad n_w=k_c c_w>0.
\end{equation}
Here $0<\chi k_c<1$. At changes in old constraint regimes, $\mathcal V'$
is continuous and the one-sided responses remain positive. Unique choices
join continuously at all regime changes, extending monotonicity to finite
gifts. Housing is strictly increasing whenever its young cap is slack.

For the finite test, define
$f(n;c,h)=\vartheta/n-\chi/(c-\chi n)-\alpha\kappa/(h-\kappa n)$.
This function strictly decreases in $n$ over its feasible interval, and the
chosen fertility is its unique zero. Evaluating the new function at $n_0$
gives \eqref{eq:finite_fertility}. If $c_1\le\chi n_0$ or
$h_1\le\kappa n_0$, the new bundle cannot accommodate $n_0$, so $n_1<n_0$.

For a current transfer $\dd b$ with future repayment of present value
$r\,\dd b$, with finance binding and young housing uncapped, the calculation gives:
\begin{equation}
\frac{\dd n}{\dd b}=
\frac{(LB_c+\chi A)/x^2+r\Gamma\Delta B_c}{J_n}.
\end{equation}
The gift has $r=0$; the term involving repayment can be negative. Similarly,
for a common gift before tenure choice, its derivative includes the term
$\pi^O(1-\pi^O)(n^O-n^R)(1/x_O-1/x_R)/\sigma_\xi$, whose sign is not
fixed. These are reasons to specify the policy before inferring aggregate
fertility from the individual gift result.

One could instead redesign the two grants to hold old resources fixed
while fertility adjusts. Along the direction in
\eqref{eq:transfer_allocation}, both adult goods and adult space rise, and:
\begin{equation}
\dd n=\frac{n^2}{\vartheta}
\left(\frac\chi{x^2}\dd x+\frac{\alpha\kappa}{s^2}\dd s\right)>0.
\end{equation}
The redesigned payments and residual funding would be:
\begin{equation}
\begin{aligned}
G&=\Delta x+L\varepsilon+(\chi+L\kappa)\Delta n,\\
qJ&=(p-L)(\varepsilon+\kappa\Delta n),\\
R&=\Delta x+(p-1/g)\varepsilon+
 [\chi+\kappa(p-1/g)]\Delta n.
\end{aligned}
\label{eq:fertile_residual}
\end{equation}
This is a different policy from holding the original grant pair fixed.
Its corrected funding term matters even at unchanged prices. In the simple case $L=p$ and
$\alpha(1+\omega_B)=\gamma$, its derivative along an adult-space increase
is $R'=(\chi-\kappa p/\alpha)n'$. It is negative if
$\chi/\kappa<p/\alpha$, although the fixed-fertility funding condition
holds at equality. Reusing that proof with free fertility would therefore
miss a current funding term as well as the change in future population.

\section{Voluntary fertility with unchanged cohort size}
\label{app:free_choice}

This construction removes the commitment restriction in
Proposition~\ref{prop:transfers}, at the cost of a more specific policy and
$\phi=q$. The intervention is unexpected, announced after the ownership
taste draw and before the young make any choices. The authority can select
households using predetermined identities,
endowments, and their realized ownership tastes. Selected young households
strictly prefer ownership at baseline. Transfers apply regardless of their
subsequent tenure choice, so both tenure and fertility remain freely chosen.
This targeting requires more information than an anonymous age-transfer rule.

Put $H_O=h_O^{\max}$ and $E=1+\alpha+\vartheta$. Since $\phi=q$, we have
$L=p$, and a strictly constrained young owner's old resources equal $v$.
For an uncapped young owner, a unit of cash raises fertility by $a_n$ and
housing by $a_h$:
\begin{equation}
a_n=\frac{\vartheta}{E(\chi+p\kappa)},\qquad
a_h=\frac{\alpha/p+\vartheta\kappa/(\chi+p\kappa)}E,
\qquad g=\frac\gamma{Kp}.
\end{equation}
Select equal positive submeasures of four groups: an uncapped, constrained
young owner A; a constrained young owner B strictly at $H_O$; an uncapped
old counterpart of A, with resources $v_A$ and slack estate floor; and an
old owner strictly at $H_O$ with slack estate floor. The latter can be B's
stationary counterpart. Uniform strict ownership margins for the young
preserve their voluntary tenure choices under small transfers.

At B's cap, let $x_B=w_B-pH_O-\chi n_B$ and $s_B=H_O-\kappa n_B$.
Its cash response of fertility is $b_B>0$. Define the tax per unit grant
needed to offset births, $q_B$, and the old tax needed to release housing, $d$:
\begin{equation}
b_B=\frac{\chi/x_B^2}
 {\vartheta/n_B^2+\chi^2/x_B^2+\alpha\kappa^2/s_B^2},
\qquad q_B=\frac{a_n}{b_B},\qquad d=\frac{a_h}{g}.
\end{equation}
A sufficient condition for the construction is:
\begin{equation}
0<q_B<1,\qquad q_B+d>1,\qquad
\frac E{w_A}>\frac{q_B}{x_B}+\frac{Kd}{v_A}.
\label{eq:free_choice_conditions}
\end{equation}
The last inequality compares the young recipient's gain per unit of cash
with the losses of the two taxpayers. The primitive construction below
makes these restrictions explicit and shows that they can hold jointly.
There is no additional restriction on $\beta/q$.

Give young A a small grant $G$. Tax young B by the amount $Q(G)$ satisfying:
\begin{equation}
n_B(w_B-Q(G))=n_B(w_B)-a_nG.
\end{equation}
This uses B's privately optimal fertility after the tax. Since $b_B>0$, the
small tax is unique and $Q'(0)=q_B$. Tax old A by $dG$ and give the capped
old group $S(G)=Q(G)+dG-G>0$. The program balances immediately:
$G+S=Q+dG$. A's housing rises by $a_hG$, B stays at its cap, old A's
housing falls by $gdG=a_hG$, and the capped old recipients retain their
housing. Thus the original price clears housing. Births also remain unchanged.

For both selected young owners, the mortgage satisfies $a'=-Ph$, so
$a'+Ph+v=v$ is unchanged. They therefore make their original old-age
choices. Current young of other types and all future entrants face the
same prices and resources, while their cohort masses are unchanged.
A's larger inherited title offsets its larger mortgage, including at resale.
Initial old consumption and estate payments change, with
$\Delta C_0+q\Delta E_1=0$ by the original budgets and the balanced transfer
account. No new government debt or unfinanced external resources are needed.

Let $m_F>0$ be the capped old recipient's marginal utility of cash. The
welfare derivative, including both young households' fertility choices, is:
\begin{equation}
\mathcal W'(0)=\frac E{w_A}-\frac{q_B}{x_B}-\frac{Kd}{v_A}
 +m_F(q_B+d-1)>0.
\end{equation}
The first three terms are positive in sum by
\eqref{eq:free_choice_conditions}, and the last is positive. Strict regime
and ownership margins and continuity give a finite positive transfer.
Thus this is an equilibrium utilitarian improvement that moves housing
toward the young with the original choice timing. Aggregate fertility is
held fixed through voluntary offsetting responses, rather than by fixing
individual fertility. Tenure probabilities must be computed from individual choice inequalities
and the original taste distribution. The endowment-only logit formula does
not apply to taste-dependent transfers. Here all tenure choices remain at
baseline.

\textbf{An explicit cap interval.} Put $r=\chi/(p\kappa)$,
$a=\alpha+\vartheta$, and define:
\begin{equation}
f(t)=t+\frac{t(1-t)}{\vartheta-at},\quad
 t_* =\frac{\vartheta}{\alpha(r+1)+\vartheta},\quad
 t_\dagger=\frac{\vartheta}{a}
 \left[1-\sqrt{\frac{\alpha r}{E(\alpha r+a)}}\right].
\end{equation}
If $\alpha r>1$, then $t_*<t_\dagger<\vartheta/a$.
Choose $t_B\in(t_*,t_\dagger)$ and set
$w_B=pH_O+(\chi H_O/\kappa)f(t_B)$ and $n_B=H_Ot_B/\kappa$.
These solve B's fertility condition with a strict cap. Here
$f'>0$, $f''>0$, and $q_B=a_n\chi f'(t_B)<1$. At the lower endpoint:
\begin{equation}
q_*=1-pa_h\left(1-\frac1{\alpha r}\right),\qquad
q_*+d-1=a_h\left(\frac1g-p+\frac p{\alpha r}\right)>0.
\end{equation}
Consequently the first two conditions in \eqref{eq:free_choice_conditions}
hold throughout this explicit interval. The finite tax is:
\begin{equation}
Q(G)=\frac{\chi H_O}{\kappa}
\left[f(t_B)-f\left(t_B-\frac{\kappa a_nG}{H_O}\right)\right].
\end{equation}
For $G<H_O(t_B-t_*)/(\kappa a_n)$, B remains strictly capped and
$0<Q'(G)<1$, $S'(G)>0$. Further restrict $G$ by A's housing cap and
mortgage margin, old A's positive resources, and the ownership margins.

\textbf{A nonempty family of primitive economies.} Fix any finite
$\beta>0$, $q\in(0,1)$, and positive $p,H_O$ and preference parameters
with $\alpha\chi>p\kappa$ and $\omega_B(1-q)>q\gamma$.
Choose $\eta\in(0,1)$ and set $v_A=\eta H_O/g$; choose B's cash from
the interval just given. The following inequalities are sufficient:
\begin{equation}
\begin{gathered}
v_B>\max\{H_O/g,\ pH_O+\beta(1+\omega_B)x_B/q\},\\
0<w_A<\min\left\{H_O/a_h,\ qEv_A/(\beta K),\
 E/(q_B/x_B+Kd/v_A)\right\}.
\end{gathered}
\end{equation}
They make old B capped, both young mortgages strict, young A uncapped, and
the welfare inequality strict. Every upper bound on $w_A$ is positive.
Give the two endowment types positive mass, choose any
$0<h_R^{\max}<H_O$, and retain logistic tastes with positive scale.
Both tenure menus have finite values, so both tenures occur and positive
subsets of strict owners can be selected.

To close the stationary economy, evaluate these menus at the auxiliary total
resources $(w_i,v_i)$ and service price $p$. Let $\bar h_{\rm life}$ be
mean lifetime housing. Choose a small positive property tax, set
$P=p/(1-q+q\tau^p)$, and set:
\begin{equation}
T=q\tau^pP\bar h_{\rm life}/2,\qquad
N=\bar H/\bar h_{\rm life}.
\end{equation}
The bound
$\tau^p<(1-q)\min_i\{w_i,v_i\}/(2qpH_O)$ ensures that $T$ is smaller
than every chosen resource level. Set primitive current income plus wealth
to $w_i-T$, split into positive components, and old income to $v_i-T$.
All real menus remain unchanged: $p,L,w_i,v_i$ are unchanged, and the
old estate floor remains slack.

For any given $\nu$, multiply both provisional child costs by
$\nu\bar n_0$, where $\bar n_0$ is mean fertility before this rescaling.
Fertility is divided by the same factor, leaving goods and space used by
children, all other real choices, and tenure-value differences unchanged.
Mean fertility becomes $1/\nu$. The young choices generate the stationary
old distribution, completing equilibrium. This is an analytical construction
of primitive economies satisfying the theorem, not a claim about already
fixed calibration parameters or a general mortgage share.

\section{Welfare weights and the retained Pareto results}
\label{app:weights}

The consumption-log coefficients are one at both ages. Under the benchmark
in \eqref{eq:welfare}, each living household receives equal weight on
remaining utility. If instead the initial old receive weight $\rho>0$,
the stationary direct gap is:
\begin{equation}
\alpha/s-\rho\gamma/h^2=(\beta/q-\rho)pm+\mu L.
\end{equation}
Equal weights on lifetime utility measured at birth imply $\rho=\beta$,
since old utility carried that factor in the original lifetime objective.
The condition then differs from the remaining-utility convention $\rho=1$.
Choosing between them is a welfare judgment, not a harmless change in notation.

The earlier note, \emph{Housing and Fertility with Conventional Mortgage
Finance}, retains the separate compensated Pareto constructions and their
stronger assumptions. Those results can serve as appendix material. The
present utilitarian argument does not use them and does not relabel its
uncompensated gains as Pareto improvements. Comparing policies that change
who is born would additionally require a population-welfare criterion;
the fertility and population statements here make no such ranking.

\end{document}
```
