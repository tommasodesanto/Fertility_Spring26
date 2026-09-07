[SYSTEM]
You are Oracle, a focused one-shot problem solver. Emphasize direct answers and cite any files referenced.

[USER]
# A finite demographic transition with renters and owners

You are an independent mathematical economic theorist. Please help settle one
specific proof question in a simple two-period overlapping-generations model
of housing and fertility. Work on the mathematics; do not rewrite the paper,
invent a different household problem, or give a literature survey. A short,
transparent result is more useful than an elaborate theorem with no clear
economic content. Treat all supplied claims as claims to check.

## The question

Can the finite transition argument in `transition_extensions.md`, section 4,
be extended to an economy with a substantial share of renters and strictly
positive child goods costs? Ideally give simple sufficient conditions that
cover a class of economies. A proved, explicit finite neighborhood around
the genuine mixed-tenure example in `mixed_transition_proof.md` would also
be useful. An existence statement for an unspecified sufficiently small
shock around that example is already available; repeating it would not
resolve this question.

The economic experiment has two steps. Start at a positive stationary
equilibrium with fertility weight \(\vartheta_0\) and financed share
\(\phi_0\). At date zero an unexpected permanent decline to
\(\vartheta_1<\vartheta_0\) starts a demographic transition. At a later date
\(t_p\geq1\), compare continuation of that baseline with an unexpected
permanent credit relaxation \(\phi_1>\phi_0\). Both continuations start
with the same actual inherited cohorts, financial claims, housing titles,
and tenure distribution at the intervention date. Household expectations
are correct following each surprise; previously chosen financial claims
remain those chosen under the earlier forecast.

The desired conclusions are:

1. Both infinite equilibrium paths exist and converge to positive stationary
   endpoints, with the original household constraints satisfied throughout.
2. The original preference decline lowers initial fertility and the baseline
   terminal population relative to the original stationary equilibrium.
3. Credit relaxation during that baseline raises fertility at the intervention
   date and leads to a larger terminal population than continuation without it.

Try first for a finite parameter range and a result valid at any later baseline
date. If that is too much, state precisely which of these conclusions you can
prove and why the rest fails or remains open. The priority is a useful theorem,
not defending every desired sign. If the proposed generality is false, give a
counterexample or an explicit additional restriction. Separate failures of a
proof method from failures of the economic claim.

## The model must stay fixed

The attached TeX gives the full model, budgets, value functions and equilibrium.
The following conventions matter:

- Young utility is
  \(\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta_t\log n\).
  Old utility is \(\log c^2+\gamma\log h^2+\omega_B\log e\), discounted
  by \(\beta\) when young. Preserve the author's \(W^m,V^m\) functions.
- Completed fertility, housing, goods, saving and tenure are joint household
  choices. A single additive ownership-taste draw is observed before those
  choices. Both parameters of its logistic distribution stay fixed within
  each economy across both shocks. Do not reset them to keep a tenure share
  constant after a reform.
- The entrant distribution of income and liquid wealth is fixed. Homogeneous
  entrants are acceptable for this proof. Estates are warm-glow expenditure
  and are not automatically inherited as entrants' liquid wealth.
- Only pre-existing liquid wealth finances the down payment:
  \((1-\phi_t)P_t h\leq b\). Current income and rebates arrive after purchase.
  Rental and owner size limits are physical constraints. Old owners may
  downsize but cannot buy additional housing or rent out retained housing.
- The consumption good is the numeraire; goods and bonds trade with the rest
  of the world at a fixed bond price \(q\). Do not introduce domestic goods
  or bond clearing. Housing clears against a fixed \(\bar H\).
- The date-\(t\) housing service cost is
  \(u_t=q r_t=(1+q\tau^p)P_t-qP_{t+1}>0\). Property-tax revenue is rebated
  equally to young and old, with
  \(T_t=q\tau^pP_t\bar H/(Y_t+O_t)\). Zero tax is acceptable as a clearly
  stated intermediate theorem; positive child goods costs and material
  renting are the main extension being sought.
- Inherited old financial assets include repayment of the mortgage originally
  contracted. Do not recalculate that mortgage at the new financed share.
  Existing titles revalue at the new price. Consequently the initial old-owner
  housing coefficient \(M_0\) in the six-variable system is not fixed when
  the surprise changes \(P_0\). Retain the actual initial-old boundary.
- Demography is \(Y_{t+1}=\nu\bar n_tY_t\), \(O_{t+1}=Y_t\). Every positive
  stationary endpoint has \(\bar n^*=1/\nu\). The long-run comparison is
  population, not permanently above-replacement fertility. With a common
  inherited young cohort, the limiting young-cohort ratio equals the limit
  of the products of the relative fertility rates.

It is acceptable to retain the strict household branches used in the proofs:
young owners have restrictive down payments, positive saving and a slack
physical cap; young and old renters are at their rental caps; old owners have
slack retention and estate bounds. Do not assume these inequalities without
checking that your proposed path region implies them. Proving transitions
across changes in binding constraints is not required for a useful answer.

## What has already been obtained

Read the attachments in this order:

1. `simplified_olg_amendment_proposal.tex`: the original household and
   equilibrium specification and the intended two-stage comparison. The
   allocation proof is background; planner design is outside this question.
2. `mixed_transition_proof.md`: exact six-variable equilibrium map, initial-old
   boundary, a one-sided sequence argument, and a genuine mixed economy with
   owner share \(11/21\), \(\chi=3/20\), and \(\tau^p=467/9250\). The
   forward derivative has four stable and two unstable roots, including one
   zero stable root. Exact algebra and interval calculations certify the local
   credit-impact and terminal-population signs for a family of taste scales.
3. `transition_extensions.md`: the newest additions. Section 2 proves wider
   mixed-tenure stationary signs when \(\chi=\tau^p=0\). Section 3 supplies
   a simpler sufficient condition for convergence in the all-owner limit.
   Section 4 gives an explicit finite preference decline, later credit reform,
   and an infinite-sequence contraction with inherited claims handled exactly.
4. `local_transition_proof.md`: earlier all-owner derivations, the exact
   parameter example used by section 4 above, and limitations/counterexamples.
   Its older statements about unspecified neighborhoods should be read together
   with the newer finite construction, not as negating that construction.
5. `verify_simplified_olg_transition_extensions.py`: the latest symbolic and
   rational-interval verification calculations. It imports earlier repository
   helpers which are not included. You are not being given a self-contained
   runnable repository; this file supplies inspectable formulas and exact
   bounds, not independent authority for any claim.

Two limitations are already established. Positive child goods costs can make
stationary credit/population signs fail outside suitable restrictions. Also,
mixed equilibrium paths can oscillate around their endpoints. We therefore
do not ask for fertility to be higher at every future date or for monotone
convergence. A welfare improvement does not imply higher fertility, and the
competitive credit reform is not asserted to be a Pareto improvement.

## What would count as progress

Please begin by identifying any substantive error in the supplied equations
or arguments that changes the target. Otherwise spend most of the answer
trying to prove the extension. A better choice of variables or an
infinite-sequence boundary-value operator may help. The six-dimensional
forward map has unstable roots, so simply assuming its whole derivative is
a contraction cannot justify the desired equilibrium transition.

Give a precise proposition, define its assumptions, and provide the proof
with the estimates on which it depends. If using validated numerical bounds,
state the finite box, the uniform residual and derivative bounds, how the
infinite tail is controlled, and how the original household inequalities and
initial-old boundary are verified. A finite-horizon simulation with an imposed
terminal steady state, pointwise eigenvalues, or a plotted path is not an
infinite-horizon existence or convergence proof. Distinguish analytical
conditions, proved interval bounds, and numerical evidence.

Look for conditions that can ultimately be explained in a few lines to an
economist. Do not hide the result in an assumption that essentially restates
the desired signs, and do not insert an unexplained borrowing multiplier into
a purported primitive condition. A modest transparent theorem is preferable
to an unsupported global claim. Aim for a compact mathematical answer, with
one clear recommendation about what belongs in the illustrative model and
what is better left as an extension.

### File: latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex
```
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern,microtype,amsmath,amssymb,amsthm,booktabs,tabularx,enumitem,graphicx,needspace,placeins}
\usepackage[colorlinks=true,linkcolor=blue,urlcolor=blue]{hyperref}
\newtheorem{proposition}{Proposition}
\newtheorem{corollary}[proposition]{Corollary}
\newtheorem{definition}{Definition}
\newcommand{\dd}{\mathrm d}
\newcommand{\Rnet}{\mathcal R}
\title{Housing Allocation, Fertility, and Population\\\large An illustrative model}
\author{Working note for Tommaso De Santo}
\date{September 6, 2026}
\begin{document}
\begin{center}
{\Large Housing Allocation, Fertility, and Population}\par
\smallskip
{\large An illustrative model}\par
\smallskip
{\small Working note for Tommaso De Santo \quad September 6, 2026}
\end{center}
\medskip

Borrowing limits can keep housing with old owners even when young households
would give up more consumption to use it. We show that a small reallocation
can then improve welfare while holding fertility fixed. We also study housing
policy during a demographic transition: an earlier decline in the desire for
children starts the adjustment, and a later policy changes its course.

\section{Environment}

\textbf{Households.}
At date $t$, the economy contains $Y_t$ young households and $O_t$ old
households. A young household has income $y_i>0$, liquid wealth $b_i>0$, and
chooses total nondurable expenditure $c$, total housing $h$, completed
fertility $n>0$, saving $a'\geq0$, and tenure $m\in\{R,O\}$. Fertility is a
continuous choice made once. The pair $(y_i,b_i)$ is drawn from the same
distribution $F$ in each cohort.

\textbf{Preferences.}
Each child uses $\chi>0$ units of goods and $\kappa>0$ units of housing. Define
the parts of the household bundle available to the adults as
\begin{equation}
 x=c-\chi n>0,\qquad s=h-\kappa n>0.
 \label{eq:amend_usable_bundles}
\end{equation}
Young and old flow utilities are
\begin{align}
 u_t^Y(c,h,n)&=\log(c-\chi n)+\alpha\log(h-\kappa n)
 +\vartheta_t\log n,
 \label{eq:amend_preferences}\\
 u^2(c^2,h^2,e)&=\log c^2+\gamma\log h^2+\omega_B\log e.
 \label{eq:amend_old_preferences}
\end{align}
Both child requirements reduce the bundle left for adults. All log arguments
and preference weights are positive. The weights $\alpha$ and $\gamma$ value
young and old housing, $\vartheta_t$ values children, and $\omega_B$ values
estates. The factor $\beta>0$ discounts old-age utility. The estate $e$ is valued in goods at
the date after the old period. Parents value what they leave, but their
children do not inherit it as entrant wealth in this model.

Each household draws once an idiosyncratic preference $\xi_i$ for ownership,
independently of $(y_i,b_i)$, and observes it before making its choices.
The draw is logistic with location $\bar\xi$ and scale $\sigma_\xi>0$.

\textbf{Demography.}
A child born at $t$ produces $\nu>0$ young households at $t+1$, after survival
and household formation. Thus
\begin{equation}
 Y_{t+1}=\nu\bar n_tY_t,\qquad O_{t+1}=Y_t,
 \label{eq:amend_demography}
\end{equation}
where $\bar n_t$ is average fertility among the young. Population counts adult households; Appendix~\ref{app:amend_units} gives
the conversion to resident persons.

\textbf{Housing and asset markets.}
The fixed housing stock is $\bar H$. Rental units satisfy
$h\leq h_R^{\max}$, owner units satisfy $h\leq h_O^{\max}$, and
$h_R^{\max}<h_O^{\max}$. The consumption good is the numeraire. Goods and one-period bonds trade
with the rest of the world. The bond price is $q\in(0,1)$ and its gross
return is $R_f=q^{-1}$; housing clears against the domestic stock. The owner
of record pays property tax at rate $\tau^p$. A competitive rental intermediary buys a unit for $P_t>0$ at the
beginning of date $t$, then collects $r_t$, pays $\tau^pP_t$, and holds a unit
worth $P_{t+1}$ at the end of the period. Free entry therefore imposes
\begin{equation}
 r_t+P_{t+1}=R_fP_t+\tau^pP_t
 \quad\Longleftrightarrow\quad
 r_t=(R_f+\tau^p)P_t-P_{t+1}.
 \label{eq:amend_usercost}
\end{equation}
Rent is paid in date-$(t+1)$ goods, so its date-$t$ value is $q r_t>0$;
at a steady state, $r=(R_f-1+\tau^p)P$. There is no capital-gains tax. Property-tax revenue is rebated equally to
all young and old households. The date-$t$ value of this rebate is $T_t$.

\textbf{Home finance.}
A young buyer finances the share $\phi_t\in(0,1)$ and posts
$(1-\phi_t)P_th$ at purchase. Closing occurs before current income and rebates
are available, and unsecured borrowing is unavailable. Thus only liquid
wealth $b_i$ can satisfy the down payment. The owner cannot buy more housing
when old or rent its retained housing to another household. At death the
estate sells it back into the fixed stock, and the proceeds enter warm-glow
estate expenditure.

\section{Household problems}

Households take prices and transfers as given. Renters carry financial wealth
into old age; owners carry the house and its mortgage.

\textbf{Young renters.}
Conditional on renting, the young household solves
\begin{align}
 W_t^R(i)=\max_{c,h,n,a'}\;&u_t^Y(c,h,n)+\beta V_{t+1}^R(R_fa') \nonumber\\
 \text{s.t.}\quad&
 c+a'+q\,r_th=y_i+b_i+T_t,\nonumber\\[-2pt]
 &h\leq h_R^{\max},\qquad a'\geq0.
 \label{eq:amend_young_renter}
\end{align}
The renter enters old age with financial wealth $R_fa'$.

\textbf{Young owners.}
Conditional on owning, the young household solves
\begin{align}
 W_t^O(i)=\max_{c,h,n,a'}\;&u_t^Y(c,h,n)+\beta V_{t+1}^O(a_{t+1},h) \nonumber\\
 \text{s.t.}\quad&
 c+a'+[(1-\phi_t)+q\tau^p]P_th=y_i+b_i+T_t,\nonumber\\[-2pt]
 &a_{t+1}=R_fa'-R_f\phi_tP_th,\qquad
 (1-\phi_t)P_th\leq b_i,\nonumber\\[-2pt]
 &h\leq h_O^{\max},\qquad a'\geq0.
 \label{eq:amend_young_owner}
\end{align}
The property tax is set aside in the bond. The mortgage is repaid when the
household becomes old, leaving financial wealth $a_{t+1}$ and title $H=h$.

\textbf{Old households.}
An old household has financial wealth $a$. An owner also carries housing $H$
from young age. Renters and owners solve, respectively,
\begin{align}
 V_t^R(a)&=\max_{c^2,h^2,e}u^2(c^2,h^2,e),\nonumber\\
 &c^2+qe+q\,r_th^2=a+T_t,\quad 0<h^2\leq h_R^{\max},
 \label{eq:amend_old_renter}\\
 V_t^O(a,H)&=\max_{c^2,h^2,e}u^2(c^2,h^2,e),\nonumber\\
 &c^2+qe+q\,r_th^2=a+P_tH+T_t,\quad
 0<h^2\leq H,\quad e\geq P_{t+1}h^2.
 \label{eq:amend_old_budget}
\end{align}
The owner may sell part of its house, and retained housing must be included
in its estate. The budget includes sale receipts, property tax, and the value
of housing eventually sold by the estate.

\textbf{Tenure choice.}
A household owns when $W_t^O(i)+\xi_i\geq W_t^R(i)$. The logistic ownership
taste gives the owner share
\begin{equation}
 \pi_t^O(i)=
 \frac{\exp\{[W_t^O(i)+\bar\xi]/\sigma_\xi\}}
 {\exp\{W_t^R(i)/\sigma_\xi\}+\exp\{[W_t^O(i)+\bar\xi]/\sigma_\xi\}}.
 \label{eq:amend_tenure}
\end{equation}
The taste draw smooths aggregate tenure shares without changing conditional
renter or owner policies.

\section{Equilibrium}

Prices and transfers equate housing demand to the fixed stock and balance the
government budget. Let $\bar h_t^Y$ and $\bar h_t^O$ denote average housing per
young and old household, averaging over types and tenure choices. Let $G_t$
be the old cohort's normalized distribution over net financial wealth,
inherited housing titles, and tenure, with separate mass $O_t$.
Housing clearing and the government budget are
\begin{equation}
 Y_t\bar h_t^Y+O_t\bar h_t^O=\bar H,\qquad
 (Y_t+O_t)T_t=q\tau^p P_t\bar H.
 \label{eq:amend_clearing}
\end{equation}
The second condition rebates all property-tax receipts in current value.
Goods and bonds trade externally; their markets do not impose an additional
domestic clearing equation.

\begin{definition}[Equilibrium]
Given policies, fertility preferences, the initial state $(Y_0,O_0,G_0)$,
and the fixed entrant distribution $F$, an equilibrium consists of price and
transfer sequences $(P_t,r_t,T_t)$, household choices, tenure probabilities,
and cohort states $(Y_t,O_t,G_t)$ satisfying:
\begin{enumerate}[label=(\roman*),nosep,leftmargin=*]
 \item Young and old households solve their stated problems, and tenure
 follows \eqref{eq:amend_tenure}.
 \item Rental entry satisfies \eqref{eq:amend_usercost}.
 \item Housing clears and the government budget balances as in
 \eqref{eq:amend_clearing}.
 \item Cohort sizes follow \eqref{eq:amend_demography}; $G_{t+1}$ is the
 distribution of assets, titles, and tenure chosen by the young at $t$.
\end{enumerate}
After an unexpected reform, old households keep their inherited assets and
titles and reoptimize using the new price path. Existing contractual
obligations remain fixed.
\end{definition}

\section{Housing allocation}

We ask whether the planner can improve the use of the existing housing stock
by moving space from old owners to young households. The planner chooses
housing and consumption directly. It respects the physical limits on each
tenure and can relax private financing restrictions. Fertility, tenure,
estates, and all future real allocations remain fixed.

Let $i$ be a young owner and $j$ an old owner. Their marginal values of
housing, measured in units of current consumption, are:
\begin{equation}
 MV_i^Y=\frac{\alpha(c_i-\chi n_i)}{h_i-\kappa n_i}
       =\frac{\alpha x_i}{s_i},\qquad
 MV_j^O=\frac{\gamma c_j^2}{h_j^2}.
 \label{eq:amend_marginal_values}
\end{equation}
These values measure the consumption each household would give up for a
little more housing. Comparing them requires no interpersonal comparison of
utility levels.

\begin{proposition}[Housing misallocation]
\label{prop:amend_welfare}
Suppose matched groups of young and old owners have equal positive mass
and room for a small transfer of housing. Assume uniformly positive young
saving, adult consumption, and room below the young physical housing cap.
If $MV_i^Y-MV_j^O$ is uniformly positive on matched pairs, a sufficiently small reallocation toward the young can
make them better off and leave every other household equally well off, in
the absence of real reallocation costs.
\end{proposition}

\begin{proof}
Give the young owner $\epsilon$ units of the old owner's housing. Keeping
the old owner's utility fixed requires additional consumption:
\begin{equation}
 D_j(\epsilon)=c_j^2\left[
 \left(\frac{h_j^2}{h_j^2-\epsilon}\right)^\gamma-1\right].
 \label{eq:amend_direct_compensation}
\end{equation}
Take this amount from the young owner. Goods and housing balance. Since
$D_j'(0)=MV_j^O$, the young owner's utility gain has derivative
$(MV_i^Y-MV_j^O)/x_i>0$ at zero. A small change preserves feasibility.
The seller replaces the lost estate title with a bond, and the buyer's
later sale of the extra title repays its financing. These claims offset,
so future allocations and existing creditors are unchanged.
Appendix~\ref{app:amend_settlement} gives the payments.
\end{proof}

\textbf{Why the values can differ.}
With positive saving and interior old-age choices, a young owner who would
buy more housing if its down-payment limit were relaxed has $MV_i^Y>q r_t$,
provided the physical owner cap is slack. An old owner with slack retention and estate
bounds instead satisfies:
\begin{equation}
 MV_j^O=\frac{\gamma c_j^2}{h_j^2}=q r_t.
 \label{eq:amend_oldmrs}
\end{equation}
Young buyers need cash at purchase; old owners already hold their homes.
The resulting value gap gives the reallocation in Proposition~\ref{prop:amend_welfare}.
The argument applies at any date when these conditions hold, including
during a demographic transition.

The rental size limit restricts an alternative way to obtain space. A
household can be unable to rent a larger home and unable to fund its
purchase. The owner-to-owner comparison above holds tenure fixed; a renter
already at its physical cap cannot receive more housing in that tenure.
An ownership taste can make constrained ownership preferable even when
larger rentals are available.

If moving $\epsilon$ housing uses $C(\epsilon)$ goods, with $C(0)=0$ and
$C'(0)=K_{ijt}\geq0$, the condition becomes $MV_i^Y>MV_j^O+K_{ijt}$.
The young pays this real cost as well as the old owner's compensation.
A fixed moving cost requires a finite gain large enough to cover it.

The proposition establishes inefficiency relative to a planner who can
relax private finance. It identifies a direction of improvement without
solving the entire first-best allocation. A market-based reform with
households choosing freely requires a separate welfare comparison.

\begin{figure}[!htbp]
 \centering
 \includegraphics[width=.73\textwidth]{output/model/simplified_olg_amendments/theory_slides_misallocation.pdf}
 \caption{Housing moves from the old owner to the young buyer, with the old
 owner compensated in consumption. Each curve uses the household's own
 housing on the horizontal axis. The arrows describe the same transfer of
 space. Fertility and future real allocations are fixed.}
 \label{fig:amend_allocation}
\end{figure}
\FloatBarrier

\section{Housing access and fertility}

More space makes children easier to accommodate. Paying for that space
reduces the goods available to the family. Fertility rises when the first
effect dominates the second.

Consider the young owner's conditional problem at a given date. Assume
positive saving, a binding down payment, a slack physical housing cap, and
interior old-age choices. Suppress household and date indices. Let
$k=\beta(1+\gamma+\omega_B)$ and
$w=y+b+T_t+qT_{t+1}$, lifetime resources before housing. Eliminating saving
from the original problem gives:
\begin{equation}
 (1+k)x+\chi n+u h=w,\qquad u=q r_t.
 \label{eq:amend_owner_resources}
\end{equation}
The net cost of housing is $u$ per unit, after accounting for its eventual
sale. The fertility first-order condition is:
\begin{equation}
 \frac{\vartheta}{n}=\frac{\chi}{x}+\frac{\alpha\kappa}{s},
 \qquad s=h-\kappa n.
 \label{eq:amend_fertility_foc}
\end{equation}
The marginal benefit of a child equals its goods and space requirements,
valued by the household.

Let $p$ be the net lifetime payment for a marginal unit of extra current
housing, including any accompanying transfers. With prices and other
resources fixed, $p=u$; compensation or a subsidy changes that payment.
The budget response is $(1+k)\dd x=-\chi\dd n-p\dd h$.
The implied fertility response is:

\begin{proposition}[A fertility condition]
\label{prop:amend_fertility}
Under these household conditions, more housing increases fertility locally
if and only if:
\begin{equation}
 \frac{\alpha\kappa}{(h-\kappa n)^2}
 >\frac{\chi p}{[1+\beta(1+\gamma+\omega_B)](c-\chi n)^2}.
 \label{eq:amend_fertility_condition}
\end{equation}
The comparison holds at any date, while tenure and the binding constraints
remain unchanged and saving and old-age choices are reoptimized.
\end{proposition}

The condition uses household allocations and the net payment, with no
borrowing multiplier. It is a conditional prediction of the household problem;
a policy's equilibrium effect also includes changes in prices, transfers,
and tenure. Appendix~\ref{app:amend_fertility} gives a sufficient restriction
using only parameters in a stationary specialization, and its extension to
dates with nonincreasing expected house prices. The fixed-fertility welfare
comparison in Proposition~\ref{prop:amend_welfare} does not require this sign.

\section{Fertility decline and policy during the transition}
\label{sec:amend_transition}

The economy need not be in a steady state when policy changes. Consider an
initial steady state with fertility weight $\vartheta_0$ and financed share
$\phi_0$. At date 0 an unexpected permanent fall to
$\vartheta_1=\vartheta_0-\epsilon$, with $\epsilon>0$, starts a baseline
transition. Explaining this initial change in preferences is outside the
exercise.

At a later date $t_p\geq1$, compare continued policy $\phi_0$ with a
permanent increase to $\phi_1=\phi_0+\delta$, where $\delta>0$.
The increase is announced and implemented at $t_p$. The two continuations
therefore start from the same young and old cohorts and the same inherited
assets and housing titles. They face the same $\vartheta_1$ and all other
non-policy primitives. This construction represents a housing policy
introduced during an adjustment already in progress.

The unanticipated policy may cause prices to jump at $t_p$, while inherited
assets and titles remain fixed. With earlier announcement, the two paths
would begin to differ at that date.

\subsection{Steady states and population}
\label{sec:amend_population}

For constant prices and policies, the household problems determine mean
fertility and housing at any candidate $(P,T)$. The old distribution comes
from repeating the young household choices. A positive steady state solves:
\begin{equation}
 \bar n(P,T;\vartheta,\phi)=\frac1\nu,\qquad
 T=\frac{q\tau^pP}{2}\,[\bar h^Y(P,T)+\bar h^O(P,T)].
 \label{eq:amend_stationary_system}
\end{equation}
The first condition is replacement; the second rebates property-tax revenue.
Once a positive root exists, housing clearing determines population:
\begin{equation}
 N_{\rm hh}^*=Y^*+O^*=2Y^*
 =\frac{2\bar H}{\bar h^{Y*}+\bar h^{O*}}.
 \label{eq:amend_population_scale}
\end{equation}
Thus two steady states can have the same replacement fertility and different
population levels. With the fixed stock, a larger population uses less
housing per household on average.

In the mixed-tenure economies used below, a small fall in $\vartheta$
lowers stationary price and population. Housing per young and per old
household rises. A small increase in $\phi$ raises stationary price,
population, and young housing; old housing falls by more than young housing
rises. These are local comparisons around an economy with both renters and
owners.

\subsection{The two continuations}

Let superscripts 0 and 1 denote the baseline and policy paths from $t_p$.
Because their initial young cohorts are identical, the demographic law gives:
\begin{equation}
 \frac{Y_{t_p+k}^1}{Y_{t_p+k}^0}
 =\prod_{j=0}^{k-1}
       \frac{\bar n_{t_p+j}^1}{\bar n_{t_p+j}^0},\qquad
 O_{t+1}^a=Y_t^a\quad(a=0,1).
 \label{eq:amend_transition_population}
\end{equation}
Relative cohort size records the accumulated fertility differences.
This identity is useful along the entire transition, without assuming
that either path is stationary.

\begin{proposition}[A decline followed by a housing policy]
\label{prop:amend_combined}
There is an open set of economies with positive child goods costs,
positive rebated property tax, and positive masses of renters and owners
for which the following holds. A sufficiently small permanent fall in
$\vartheta$ produces a converging baseline transition with lower initial
fertility and a smaller terminal population. A sufficiently small permanent
increase in $\phi$, introduced at any later date on that baseline, has a
converging continuation from the same inherited state and satisfies:
\begin{equation}
 \bar n_{t_p}^1>\bar n_{t_p}^0,\qquad
 N_{{\rm hh},1}^*>N_{{\rm hh},0}^*.
 \label{eq:amend_policy_comparison}
\end{equation}
The same constraints continue to bind along these local transitions.
\end{proposition}

Appendix~\ref{app:amend_transition} gives the assumptions and proof. The
policy need not restore the original population or raise fertility above
replacement.

Both terminal fertility levels equal $1/\nu$. With exponential convergence,
\eqref{eq:amend_transition_population} also implies:
\begin{equation}
 \sum_{t=t_p}^{\infty}\log\frac{\bar n_t^1}{\bar n_t^0}
 =\log\frac{N_{{\rm hh},1}^*}{N_{{\rm hh},0}^*}>0.
 \label{eq:amend_cumulative_fertility}
\end{equation}
The cumulative relative fertility effect is positive. Individual dates can
have effects of either sign: the example admits damped oscillations as it
approaches the new steady state.

\begin{figure}[htbp]
 \centering
 \includegraphics[width=\textwidth]{output/model/simplified_olg_amendments/combined_transition_figure.pdf}
 \caption{A fertility decline followed by a housing policy. From the old
 steady state $S_-$, the preference decline at $t=0$ starts the move toward
 $S_0$. The policy changes
 the intervention-date equilibrium from $I^0$ to $I^1$, with the same inherited
 cohorts and claims, and leads to $S_1$. Dashed lines are constant-price
 schedules. Arrows connect selected first-order equilibrium responses and
 do not imply monotone adjustment. Appendix~\ref{app:amend_figure} defines
 the comparison.}
 \label{fig:amend_credit}
\end{figure}

The permanent policy changes the terminal household problem. A temporary
reform that restores all final primitives cannot change a unique stationary
endpoint if both paths converge to it. The population comparison also does
not establish a welfare gain from the credit reform. The compensated
allocation and the competitive policy answer different questions.

\clearpage
\appendix
\section{Household arguments}
\label{app:amend_households}

\subsection{Financial settlement of the reallocation}
\label{app:amend_settlement}

The direct-allocation proof uses consumption and housing as resources.
Prices enter only to check that inherited claims and future allocations can
be respected. Set $u_t=q r_t$ and write the payment to the seller in excess
of the title price as:
\begin{equation}
 L(\epsilon)=D_j(\epsilon)-u_t\epsilon.
 \label{eq:amend_compensation}
\end{equation}
The seller receives $P_t\epsilon+L(\epsilon)$ and releases the tax reserve
$q\tau^pP_t\epsilon$. After consuming $D_j(\epsilon)$, it buys a bond
costing $qP_{t+1}\epsilon$, replacing the title lost from its estate.
The young buyer spends $D_j(\epsilon)$ less on current consumption and
finances the rest of the purchase and compensation. Its net borrowing is:
\begin{equation}
 (P_t+q\tau^pP_t)\epsilon+L(\epsilon)-D_j(\epsilon)
 =qP_{t+1}\epsilon.
 \label{eq:amend_settlement}
\end{equation}
The two new bond positions offset. At $t+1$ the buyer sells the extra title
for $P_{t+1}\epsilon$ and repays the bond. The old seller's estate receives
the same amount through its bond. The timing and quantity of future housing
sales, estate expenditures, total tax receipts, and common rebates are unchanged.
Existing loans retain their promised repayments. Positive saving provides
room to adjust the original financial portfolio; the planner relaxes the
private borrowing limit. Uniform strict margins on matched groups permit
a common sufficiently small transfer.

With a real cost $C(\epsilon)$, the young owner's exact gain is:
\begin{equation}
 \Delta U_i=
 \log\!\left(1-\frac{D_j(\epsilon)+C(\epsilon)}{x_i}\right)
 +\alpha\log\!\left(1+\frac{\epsilon}{s_i}\right).
 \label{eq:amend_exact_gain}
\end{equation}
Its derivative at zero is $(MV_i^Y-MV_j^O-K_{ijt})/x_i$.
If a fixed cost is present, the finite transfer must keep consumption
positive and make this exact gain positive.

\subsection{Fertility and a sufficient parameter restriction}
\label{app:amend_fertility}

For a fixed tenure, positive saving, and unchanged old-age constraints,
let $\widetilde w$ denote lifetime resources net of any fixed old-age rent.
Let $\widehat h$ be the current housing limit: $h_R^{\max}$ for a renter
and $\min\{h_O^{\max},b/[(1-\phi)P_t]\}$ for an owner. The conditional
problem, after eliminating saving and setting $s=h-\kappa n$, is:
\begin{equation}
 \begin{aligned}
 \max_{x,s,n}\quad &(1+k)\log x+\alpha\log s+\vartheta\log n\\
 \text{s.t.}\quad &(1+k)x+u s+(\chi+\kappa u)n=\widetilde w,
 \qquad s+\kappa n\leq\widehat h.
 \end{aligned}
 \label{eq:amend_reduced_problem}
\end{equation}
The goods and space used by each child cost $\chi+\kappa u$.
For interior old owners, $k=\beta(1+\gamma+\omega_B)$ and
$\widetilde w=y+b+T_t+qT_{t+1}$. For capped old renters,
$k=\beta(1+\omega_B)$ and the old rent must be subtracted:
\begin{equation}
 \widetilde w=y+b+T_t+qT_{t+1}-q^2r_{t+1}h_R^{\max}.
 \label{eq:amend_capped_old}
\end{equation}
This reduction requires different formulas if an old owner's retention or
estate bound binds.

Differentiating the fertility condition and using
$(1+k)\dd x=-\chi\dd n-p\dd h$ gives:
\begin{equation}
 \frac{\dd n}{\dd h}
 =\frac{\alpha\kappa/s^2-\chi p/((1+k)x^2)}{\Delta},\qquad
 \Delta=\frac{\vartheta}{n^2}+\frac{\chi^2}{(1+k)x^2}
              +\frac{\alpha\kappa^2}{s^2}>0.
 \label{eq:amend_fertility_sign}
\end{equation}
The numerator separates the benefit of additional space from the cost of
paying for it.

Differentiating the resource and fertility equations, allowing prices,
resources, and the preference weight to change, gives:
\begin{equation}
 \Delta\,\dd n=\frac{\dd\vartheta}{n}
 +\left(\frac{\alpha\kappa}{s^2}-\frac{\chi u}{(1+k)x^2}\right)\dd h
 +\frac{\chi}{(1+k)x^2}(\dd\widetilde w-h\,\dd u).
 \label{eq:amend_dated_fertility}
\end{equation}
This separates preference, space, and resource effects at any date.
Equation~\eqref{eq:amend_fertility_sign} sets the preference weight fixed and
collects the net resource change into the payment $p$.

For a sufficient condition without the borrowing multiplier, write
$MV^Y=\alpha x/s$. If $\kappa u/\chi>\alpha/(1+k)$ and the housing cap
binds, then $(1+k)\kappa MV^Y>\alpha\chi$. It follows that fertility rises
for any $p<MV^Y$. A condition entirely in parameters can imply this inequality.

\begin{proposition}[A sufficient parameter restriction]
\label{prop:amend_primitive}
Suppose property tax and rebates are zero, old-owner choices are interior,
saving is positive, the physical owner cap is slack, and the down payment
is strictly restrictive.
At stationary prices, the following restriction is sufficient for fertility
to rise after a local housing increase with net payment $p<MV^Y$:
\begin{equation}
 \frac{(1-q)b}{(1-\phi)(y+b)}
 \geq\frac{\alpha}{1+\beta(1+\gamma+\omega_B)+\alpha}.
 \label{eq:amend_fully_primitive}
\end{equation}
It is also sufficient at a transition date if the household expects
$P_{t+1}\leq P_t$.
\end{proposition}

\begin{proof}
Let $v=\kappa u/\chi$ and let $a_h=u\widehat h/\widetilde w$ be the
housing expenditure share at the cap. Within the reduced problem, removing
the current cap gives the optimal housing expenditure share:
\begin{equation}
 g(v)=\frac{\alpha+\vartheta v/(1+v)}{1+k+\alpha+\vartheta}.
 \label{eq:amend_unconstrained_share}
\end{equation}
Strict binding gives $a_h<g(v)$. The function is increasing, and
$g(\alpha/(1+k))=\alpha/(1+k+\alpha)$. At stationary prices,
$a_h=(1-q)b/[(1-\phi)(y+b)]$. The restriction therefore implies
$v>\alpha/(1+k)$. If $P_{t+1}\leq P_t$, then
$u/P_t=1-qP_{t+1}/P_t\geq1-q$, so the same lower bound on $a_h$
holds at that date. The preceding sufficient condition gives the result.
\end{proof}

This parameter restriction excludes very low liquid wealth relative to
lifetime resources. More generally, welfare can increase while fertility
falls: the exact sign still depends on both space and payment.
If a common rental cap increases at both ages, it also raises future rent;
under binding caps the appropriate net payment is $p=u_t+qu_{t+1}$.

A lower preference weight weakly lowers an individual household's fertility
when its full price and transfer sequence is fixed, even allowing tenure
choice. To see this, write its objective as $A(z)+\vartheta\log n(z)$ for a
feasible lifetime plan $z$. Adding the two optimality inequalities at
weights $\vartheta_0$ and $\vartheta_1$ gives:
\begin{equation}
 (\vartheta_1-\vartheta_0)(\log n_1-\log n_0)\geq0.
 \label{eq:amend_preference_fertility}
\end{equation}
The feasible set is unchanged. The general-equilibrium sign also requires
the price and tenure responses considered below.

\subsection{Population units}
\label{app:amend_units}

Let $n^L=a_n n$ count literal children. Consistent rescaling uses
$\nu^L=\nu/a_n$, $\chi^L=\chi/a_n$, and $\kappa^L=\kappa/a_n$.
Replacement literal fertility is $a_n/\nu$; the change in $\log n$ is a
choice-independent constant. If young and old households contain
$a_Y$ and $a_O$ adults and children reside with parents in their birth
period, resident population is:
\begin{equation}
 N_{{\rm persons},t}=a_YY_t+a_OO_t+a_n\bar n_tY_t.
 \label{eq:amend_persons}
\end{equation}
At a positive steady state,
$N_{\rm persons}^*/N_{\rm hh}^*=(a_Y+a_O+a_n/\nu)/2$.
Common fixed counting coefficients therefore preserve the terminal
population ranking. Their empirical values are not needed for this
illustrative comparison.

\section{The combined transition}
\label{app:amend_transition}

The proof uses the original equilibrium equations on a region where the
same household constraints remain strictly binding or slack. Young owners
are at their down-payment limits and below their physical caps; young and
old renters are at their rental caps. Young households save positively.
Old owners have slack retention and estate bounds. All consumption,
usable space, fertility, and estate levels are positive. These conditions
hold uniformly over the fixed compact support of entrant types.

\subsection{Equilibrium equations and inherited claims}

Write $u_t=q r_t$, $a=h_R^{\max}$,
$\rho_O=1+\beta(1+\gamma+\omega_B)$, and
$\rho_R=1+\beta(1+\omega_B)$. Housing choices are
$h_i^O=b_i/[(1-\phi)P_t]$ and $h_i^R=a$.
The conditional resource equations and fertility conditions are:
\begin{align}
 \rho_Ox_i^O+\chi n_i^O&=w_i-u_th_i^O,\nonumber\\
 \rho_Rx_i^R+\chi n_i^R&=w_i-u_ta-qu_{t+1}a,\nonumber\\
 \frac{\vartheta}{n_i^m}&=\frac{\chi}{x_i^m}
  +\frac{\alpha\kappa}{h_i^m-\kappa n_i^m},\qquad
 w_i=y_i+b_i+T_t+qT_{t+1}.
 \label{eq:amend_transition_choices}
\end{align}
The future rent affects renter choices and therefore tenure probabilities.
For the cohorts choosing under these prices, old-owner housing is
$\beta\gamma x_i^O/(qu_{t+1})$ and old-renter housing is $a$.

The aggregate state is $Z_t=(P_t,u_t,Y_t,O_t,M_t,R_t)$, where
$M_t/u_t$ and $R_t$ are total old-owner and old-renter housing.
Given this state, rental entry determines $P_{t+1}$ and fiscal balance
determines $T_t$. The two remaining unknowns
$v=(u_{t+1},Y_{t+1})$ solve:
\begin{align}
 F_1(Z_t,v)&=Y_t\int[\pi_i^Oh_i^O+(1-\pi_i^O)a]\,\dd F
              +M_t/u_t+R_t-\bar H=0,\nonumber\\
 F_2(Z_t,v)&=Y_{t+1}-\nu Y_t
              \int[\pi_i^On_i^O+(1-\pi_i^O)n_i^R]\,\dd F=0.
 \label{eq:amend_transition_residuals}
\end{align}
These are housing clearing and the cohort law. The next state is:
\begin{equation}
 G(Z_t,v)=\left(
 P_{t+1},u_{t+1},Y_{t+1},Y_t,
 \frac{\beta\gamma}{q}Y_t\int\pi_i^Ox_i^O\,\dd F,
 aY_t\int(1-\pi_i^O)\,\dd F\right).
 \label{eq:amend_transition_update}
\end{equation}
If $F_v$ is nonsingular, these equations give a smooth local map
$Z_{t+1}=g(Z_t;\vartheta,\phi)$. It is derived from the original equilibrium,
including old-state evolution.

At an unexpected policy date, let $B,A,H$ denote the inherited old-owner
mass, net financial claims, and purchased titles, respectively. Let
$\mathcal I=(Y_{t_p},O_{t_p},B,A,H)$ and
$L=\gamma/(1+\gamma+\omega_B)$. Initial housing demand satisfies:
\begin{equation}
 Y=Y_{t_p},\quad O=O_{t_p},\quad R=a(O_{t_p}-B),\quad
 M=L[A+PH+T(P,Y,O)B].
 \label{eq:amend_actual_boundary}
\end{equation}
The title price and rebate can change, so $M$ is not predetermined wealth.
Financial claims include the original mortgage repayment; the new $\phi$
does not rewrite that debt. The full conditional inherited states are also
retained to check each old household's feasibility.

Previously chosen saving and rent use the forecasts held before the shock.
In particular, rental entry at $t_p-1$ uses the baseline expected price at
$t_p$, and the old cohort's earlier saving condition need not hold at the
surprisingly realized price. Equations~\eqref{eq:amend_transition_choices}--
\eqref{eq:amend_transition_update} apply to the new continuation from $t_p$
onward, with \eqref{eq:amend_actual_boundary} determining the initial old.

\subsection{Local continuation and the two shocks}

Let $Z^*$ be a positive stationary equilibrium satisfying the stated strict
household conditions. Assume $g$ and the four initial restrictions $b$
are smooth, $F_v$ is nonsingular, and $J=g_Z(Z^*)$ has four eigenvalues
strictly inside and two strictly outside the unit circle. Let $B_0=b_Z$
at the initial stationary data. Its restriction to the stable generalized
eigenspace must be invertible. Finally, at this equilibrium the impact
fertility and stationary population derivatives with respect to both
$\vartheta$ and $\phi$ are strictly positive, holding inherited claims fixed
in the impact comparisons.

Here is the local argument. Write $\lambda=(\vartheta,\phi)$ for the two
parameters. Since one is not an eigenvalue of $J$, the stationary state
$Z^*(\lambda)$ is smooth. Center a continuation on it by
$v_k=Z_{t_p+k}-Z^*(\lambda)$. The exact sequence equations are:
\begin{equation}
 b(Z^*(\lambda)+v_0;\mathcal I)=0,\qquad
 v_{k+1}-g(Z^*(\lambda)+v_k;\lambda)+Z^*(\lambda)=0.
 \label{eq:amend_sequence_equations}
\end{equation}
Centering removes the permanent parameter shift from the tail.
Choose $\eta<1$ above the stable spectral radius and use the sequence norm
$\sup_k\eta^{-k}\|v_k\|$. The derivative in $v$ is
$(B_0v_0,\{v_{k+1}-Jv_k\})$. For a forcing $f_k$, its stable and
unstable components are solved by:
\begin{equation}
 v_k^s=J_s^kc+\sum_{j=0}^{k-1}J_s^{k-1-j}f_j^s,
 \qquad
 v_k^u=-\sum_{j=k}^{\infty}J_u^{k-1-j}f_j^u.
 \label{eq:amend_sequence_inverse}
\end{equation}
The boundary determines $c$ uniquely; both sums are bounded in this norm.
Only the unstable block is inverted, so a zero stable root is allowed.
Uniform household margins make the residual continuously differentiable.
The sequence implicit-function theorem therefore gives a smooth,
exponentially converging continuation for nearby parameters and inherited
data. The same inverse on bounded sequences supplies uniqueness among
continuations staying in a sufficiently small common neighborhood.

For the initial preference decline, inherited data are the original
stationary data. Smoothness gives a constant $C$ such that:
\begin{equation}
 \|Z_t^0(\epsilon)-Z^*(\vartheta_0-\epsilon,\phi_0)\|
 \leq C\epsilon\eta^t,\qquad
 \|Z^*(\vartheta_0-\epsilon,\phi_0)-Z^*(\vartheta_0,\phi_0)\|
 \leq C\epsilon.
 \label{eq:amend_uniform_baseline}
\end{equation}
The entire baseline stays uniformly close to the original equilibrium.
The strict preference derivatives give lower initial fertility and a smaller
endpoint for sufficiently small positive $\epsilon$.

Original household budgets reconstruct assets, titles, and tenure
probabilities smoothly on the fixed compact type support. Their deviations
are therefore uniformly of order $\epsilon$, across every date and type.
This controls individual feasibility as well as the aggregate inherited
condition; closeness of aggregate moments alone would not suffice.
At $\delta=0$, local uniqueness identifies the continuation at $t_p$ with
the baseline tail. Continuity preserves the strict policy derivatives for
all these inherited states in one neighborhood. Integrating them from
$\phi_0$ to $\phi_0+\delta$ gives
\eqref{eq:amend_policy_comparison}. The uniform bound permits common positive
upper bounds on $\epsilon$ and $\delta$, independent of the chosen date $t_p$.
This proves Proposition~\ref{prop:amend_combined} once an economy meeting the
strict conditions is supplied.

\subsection{An economy satisfying the conditions}

Use the following primitive values, with homogeneous income and liquid
wealth and a finite ownership-taste scale:
\begin{align}
 q&=\tfrac12,\quad\phi=\tfrac45,\quad b=\tfrac15,\quad
 \alpha=\beta=\omega_B=\tfrac25,\quad\gamma=\tfrac3{10},\nonumber\\
 \chi&=\tfrac3{20},\quad\kappa=\tfrac12,\quad
 \vartheta=\tfrac{141}{400},\quad\nu=2,\quad
 h_R^{\max}=\tfrac14,\quad h_O^{\max}=2,\quad\sigma_\xi=4,\nonumber\\
 \tau^p&=\tfrac{467}{9250},\quad
 y=\tfrac{3521349281}{1677802000},\quad
 \bar H=\tfrac{68104}{68019}.
 \label{eq:amend_example_primitives}
\end{align}
Set $\bar\xi=\sigma_\xi\log(\pi/(1-\pi))-W^O+W^R$ at
$\pi=11/21$ and the allocation below, then hold both taste parameters fixed
during the shocks. This specifies one preference distribution.

There is an exact stationary equilibrium with $P=Y=O=1$,
$u=9717/18500$, $T=3975571/314587875$, and:
\begin{equation}
 (x_O,h_O,n_O)=(1,1,\tfrac34),\qquad
 (x_R,h_R,n_R)=(\tfrac{99}{74},\tfrac14,\tfrac9{40}).
 \label{eq:amend_example_allocation}
\end{equation}
Mean fertility is $1/2$ and $11/21$ of households own. The young owner's
housing value is $16/25$, exceeding the old owner's value $9717/18500$.
All stated household inequalities have strict margins. Thus the allocation
and transition results apply in the same economy.

Differentiating \eqref{eq:amend_transition_choices} gives the two additional
terms needed for the preference shock: $\dd\vartheta/n_i^m$ in the fertility
condition and $\log n_i^m\,\dd\vartheta$ in the conditional value.
Together with \eqref{eq:amend_tenure}, the latter accounts for endogenous
tenure changes. Let $Q_a=g_a$ for $a\in\{\vartheta,\phi\}$; then:
\begin{equation}
 J=G_Z-G_vF_v^{-1}F_Z,\qquad
 Q_a=G_a-G_vF_v^{-1}F_a,\qquad
 Z_a^*=(I-J)^{-1}Q_a.
 \label{eq:amend_example_derivatives}
\end{equation}
These expressions follow directly from the two implicit equations.

Exact root bounds give one zero root, one real stable root in
$(.00319,.00320)$, a nonreal pair with modulus in $(.41,.43)$, and
unstable roots in $(-80.93,-80.92)$ and $(4.148,4.149)$.
The four initial restrictions are independent on the stable space.
Solving the linear recurrence and those restrictions gives:
\begin{center}
\begin{tabular}{@{}lrr@{}}
\toprule
Derivative at the reference equilibrium & $\vartheta$ increases & $\phi$ increases\\
\midrule
Initial mean fertility & $0.84713$ to $0.84714$ & $0.07966$ to $0.07968$\\
Stationary adult households & $4.79728$ to $4.79729$ & $2.33786$ to $2.33788$\\
\bottomrule
\end{tabular}
\end{center}
The bounds are outward enclosures; all four derivatives are positive.
The preference derivative uses $1<\log(10/3)<2$, with a tighter rational
series bound for the displayed intervals. Strict inequalities persist in
an open neighborhood. Uniformly small income and wealth variation on a
fixed compact probability space preserves the argument as well. This does
not impose zero renter mass, zero child goods costs, or zero property tax.

The admissible neighborhood has no stated numerical radius. Larger shocks
can change which constraints bind or eliminate the positive endpoint.
For example, the fertility condition implies:
\begin{equation}
 n<\frac{\vartheta h}{\kappa(\alpha+\vartheta)}
 \leq\frac{\vartheta h_O^{\max}}{\kappa(\alpha+\vartheta)}.
 \label{eq:amend_endpoint_obstruction}
\end{equation}
If the last bound is at most $1/\nu$, a positive closed steady state is
impossible. The local construction does not cover such a decline.

\subsection{The two-panel comparison}
\label{app:amend_figure}

Figure~\ref{fig:amend_credit} uses the analytical first-order response at the
economy above. Take a preference decline of size $\epsilon$, introduce the
credit change at $t_p=1$, and set $\delta=\epsilon/2$. These choices specify
a direction as $\epsilon$ tends to zero; they do not select a finite
reform size or a historical date. If $z^a_k$ is the first-order response to
a permanent unit increase in parameter $a$ from the reference equilibrium,
the two displayed directions are:
\begin{equation}
 z_t^0=-z_t^{\vartheta},\qquad
 z_t^1=\begin{cases}
 z_t^0,&t<1,\\
 z_t^0+\tfrac12 z_{t-1}^{\phi},&t\geq1.
 \end{cases}
 \label{eq:amend_figure_paths}
\end{equation}
Their interaction is of second order. The policy response at date 1
satisfies the homogeneous version of \eqref{eq:amend_actual_boundary},
so it preserves that date's inherited cohorts and claims.

The dashed lines are tangents to housing and fertility schedules evaluated
at constant prices and equal age-cohort masses. Write $N=N_{\rm hh}$. At a
candidate $(P,N)$, set $T=q\tau^pP\bar H/N$, solve the household problems with repeated prices,
and define:
\begin{equation}
 \mathcal H(P,N;\vartheta,\phi)
 =\frac N2(\bar h^Y+\bar h^O)-\bar H.
 \label{eq:amend_figure_schedules}
\end{equation}
The left panel solves $\mathcal H=0$ locally for price as a function of
population; the right panel evaluates mean fertility on that same schedule.
On the underlying exact schedules, an intersection with replacement
fertility is a stationary equilibrium. The drawn tangent intersections give
its first-order location. Other points on the schedules need not satisfy
the demographic law.

The solid arrows connect selected transition equilibria. During adjustment,
cohort proportions, inherited assets, and expected future prices differ
from the stationary schedules, so a transition point need not lie on one
of their dashed lines. The arrows compare initial, intervention-date, and
limiting states; they do not impose monotone motion between them. Displayed coordinates are first-order changes per unit $\epsilon$ from the
old equilibrium; symbolic labels refer to economic levels. Both axes use
affine scales.

\end{document}
```

### File: output/model/simplified_olg_amendments/mixed_transition_proof.md
```md
# A local transition with substantial renting and positive child costs

September 6, 2026. Supporting proof for the illustrative theory; the author’s
manuscript is unchanged. This extends the earlier transition near the all-owner,
zero-child-cost limit. It does not change the household problem or adopt a
public-finance institution. Independent verification is recorded at the end.

## 1. What the result says

There is an open set of economies with endogenous ownership, substantial
renting, positive child goods costs, and a positive rebated property tax in
which a small permanent increase in the financed share has a locally unique
converging competitive-equilibrium transition. Initial mean fertility rises
and the final population is larger. The original initial old keep their
pre-reform financial assets and housing titles.

The construction below has an owner share of $11/21$, so almost half the
households rent. It proves local existence around an interior mixed economy,
not just continuity from a zero-renter limit. The neighborhood and the
permitted reform size are not given numerical lower bounds. It does not prove
global convergence, arbitrary changes of binding constraints, or a general
fertility sign. Its finite small reforms eventually oscillate around the new
steady state: cumulative fertility can rise without fertility increasing at
every date. The short direct-allocation welfare result remains separate.

## 2. Original choices and the maintained branch

Keep the original flow utilities and value functions. Write
$x=c-\chi n$, $s=h-\kappa n$, and $u_t=qr_t$ only to shorten the proof.
The primitives and the fertility weight $\vartheta$ are constant after the
reform. There is no capital-gains tax. Define
\[
 \rho_O=1+\beta(1+\gamma+\omega_B),\qquad
 \rho_R=1+\beta(1+\omega_B),\qquad a=h_R^{\max}.
\]
We maintain the following uniformly strict conditional branches. Young owners
are at their down-payment limit, below their physical owner cap, and save
positively. Young and old renters are at their rental caps and young renters
save positively. Old owners are interior, with slack retention and estate
bounds. Adult goods, adult space, fertility and estates are positive. Every
conditional young owned size exceeds $a$. Conditional choices are evaluated
in both tenures even when one is selected less often.

Let $F$ be the fixed entrant distribution of income and wealth, with compact
support and uniform margins for these inequalities. Additive ownership tastes
have a finite positive logistic scale $\sigma_\xi$. Conditional policies are
independent of that taste draw; ownership probabilities are not. The same
probability must enter housing, births and subsequent old-age demand.

At date $t$, set $w_i=y_i+b_i+T_t+qT_{t+1}$. The original dated budgets and
first-order conditions imply
\[
 h_i^O=\frac{b_i}{(1-\phi)P_t},\quad h_i^R=a,\qquad
 \rho_Ox_i^O+\chi n_i^O=w_i-u_th_i^O,
\]
\[
 \rho_Rx_i^R+\chi n_i^R=w_i-u_ta-qu_{t+1}a,\qquad
 \frac{\vartheta}{n_i^m}
 =\frac{\chi}{x_i^m}+\frac{\alpha\kappa}{h_i^m-\kappa n_i^m}.
\]
Old choices for these young households satisfy
\[
 c_{t+1,i}^{2,m}=\frac{\beta x_{t,i}^m}{q},\quad
 h_{t+1,i}^{2,O}=\frac{\beta\gamma x_{t,i}^O}{q u_{t+1}},\quad
 h_{t+1,i}^{2,R}=a,\quad
 e_{t+1,i}^m=\frac{\beta\omega_Bx_{t,i}^m}{q^2}.
\]
Saving and old net assets are reconstructed from the original young budget
and mortgage repayment. In particular the estate is dated in next-period
goods; the denominator above is $q^2$.

## 3. Six variables suffice for aggregate equilibrium

Use $Z_t=(P_t,u_t,Y_t,O_t,M_t,R_t)$. The variables $Y_t,O_t$ are cohort
masses. The quantity $M_t/u_t$ is total old-owner housing; $M_t$ is not financial
wealth. The variable $R_t$ is total old-renter housing. For $t\geq1$,
\[
 M_t=\frac{\beta\gamma}{q}Y_{t-1}
       \int\pi_{t-1}^O(i)x_{t-1,i}^O\,dF,
 \qquad
 R_t=aY_{t-1}\int(1-\pi_{t-1}^O(i))\,dF.
\]
The user-cost identity and balanced tax rebate give
\[
 P_{t+1}=\frac{(1+q\tau^p)P_t-u_t}{q},\qquad
 T_t=\frac{q\tau^pP_t\bar H}{Y_t+O_t}.
\]
Given $Z_t$, take $v=(u_{t+1},Y_{t+1})$ as two unknowns and set
$T_{t+1}=q\tau^pP_{t+1}\bar H/(Y_{t+1}+Y_t)$. Evaluate the preceding
original conditional choices and the logistic tenure probability. Solve
\[
 F_1(Z_t,v)=Y_t\int[\pi_i h_i^O+(1-\pi_i)a]dF
                    +M_t/u_t+R_t-\bar H=0,
\]
\[
 F_2(Z_t,v)=Y_{t+1}-\nu Y_t
                 \int[\pi_i n_i^O+(1-\pi_i)n_i^R]dF=0.
\]
The exact update is
\[
 Z_{t+1}=G(Z_t,v)=
 \left(P_{t+1},u_{t+1},Y_{t+1},Y_t,
 \frac{\beta\gamma}{q}Y_t\int\pi_i x_i^OdF,
 aY_t\int(1-\pi_i)dF\right).
\]
Where $F_v$ is nonsingular, the implicit function theorem gives a smooth
local map $Z_{t+1}=g(Z_t;\phi)$. This is a representation of the original
equilibrium equations, including the future rental cost in tenure choice.
It is not an assumed adjustment rule.

### A simple sufficient condition for the forward solve

For homogeneous income and wealth, suppose also that each old owner occupies
at least $a$. These conditions imply nonsingularity of $F_v$ even with a
positive property tax; neither a small-tax nor a small-child-cost condition is
needed for this particular step.

To see this, housing clearing fixes the owner probability at the given state:
\[
 \pi=\frac{\bar H-M/u-R-Ya}{Y(h^O-a)}\in(0,1).
\]
Write $z=u_{t+1}$ and let $w$ vary with the future rebate. The envelope
derivatives of $\Delta=W^O-W^R$ are
\[
 \Delta_z=-\frac{\beta\gamma}{z}+\frac{qa}{x_R}<0,
 \qquad \Delta_w=\frac1{x_O}-\frac1{x_R}.
\]
The strict inequality is exactly the old renter's cap condition. Holding
$\pi$ fixed therefore implies
\[
 \frac{dz}{dw}=\frac{1/x_O-1/x_R}
                    {\beta\gamma/z-qa/x_R}.
\]
At fixed housing, fertility is increasing in resources when $\chi>0$:
$n_{m,w}>0$. The renter's dependence on $z$ is
$n_{R,z}=-qa n_{R,w}$; the owner's is zero at fixed $w$. Hence along the
fixed-tenure-probability equation,
\[
 \frac{d\bar n}{dw}=\pi n_{O,w}
 +(1-\pi)n_{R,w}
       \frac{\beta\gamma/z-qa/x_O}{\beta\gamma/z-qa/x_R}\geq0.
\]
The numerator in the fraction is nonnegative precisely when $h^{2,O}\geq a$.
Meanwhile $w$ is nonincreasing in $Y_{t+1}$ because the future tax revenue is
rebated to $Y_{t+1}+Y_t$ households. After solving housing for $z$, the
derivative of the birth-equation residual with respect to $Y_{t+1}$ is
$1-\nu Y_t(d\bar n/dw)(dw/dY_{t+1})\geq1$. This proves local
nonsingularity. At $\chi=0$ the resource derivatives of fertility vanish and
the same conclusion follows. This argument does not assert existence at
distant states or global continuation across a branch change.

## 4. The actual initial old

Keep the complete pre-reform old distribution $G_0$ as exogenous data. Let
$B_0$ be its owner mass, $A_0$ its aggregate owner net financial assets after
the original mortgage repayment, and $H_0$ its aggregate inherited owner
title. The latter is inherited purchased housing, not housing retained after
the old reoptimize. Let $L=\gamma/(1+\gamma+\omega_B)$.

The four initial restrictions are
\[
 Y=Y_0,\quad O=O_0,\quad R=a(O_0-B_0),\qquad
 M=L[A_0+PH_0+T(P,Y_0,O_0)B_0].
\]
The initial owner coefficient varies with the surprising price and rebate;
fixing it would wrongly suppress revaluation of existing title. These four
aggregates suffice for date-zero housing clearing. The full $G_0$ is still
needed to reconstruct and check each old household's consumption and estate.
Uniform slack at the initial equilibrium preserves those inequalities under
a small reform.

## 5. A local nonlinear transition theorem

Suppose the local map and boundary operator are at least $C^2$, $Z^*$ is a
stationary equilibrium on this branch, $F_v$ is
nonsingular, and $J=g_Z(Z^*;\phi)$ has four eigenvalues strictly inside and
two strictly outside the unit circle, counted with multiplicity. Let $B$ be
the derivative of the four initial restrictions with respect to $Z_0$.
Suppose $B$ restricted to the four-dimensional stable generalized eigenspace
is invertible. Then sufficiently small permanent changes of $\phi$, starting
from the original stationary population and claims, have a locally unique
equilibrium sequence converging exponentially to the nearby stationary state.

**Proof.** Since $1$ is not an eigenvalue, $(I-J)$ is invertible and the
stationary state is a smooth function of the reform. Center the sequence on
this new state. Split the linear recurrence into its stable and unstable
generalized eigenspaces. Propagate the stable part forward and obtain its
initial value from $B$; solve the unstable part backward using its inverse
and the requirement that the sequence converge. Choose a weight $\eta<1$
above all stable moduli and above the inverse unstable moduli. The resulting
operator has a bounded inverse on sequences with norm
$\sup_t\eta^{-t}\|v_t\|$. The nonlinear residual and boundary conditions
are smooth on a uniform neighborhood, so the sequence implicit function
theorem applies. A zero stable eigenvalue is harmless: only the unstable
block is inverted. Uniform household margins and differentiation under the
fixed compact-support integral reconstruct the full original equilibrium.
This proves uniqueness among nearby convergent sequences, not global
uniqueness of equilibrium. $\square$

## 6. An exact interior mixed economy

All the following primitive values are rational except the ownership taste
location, which is specified by the exact value functions below:
\[
 q=\tfrac12,\quad\phi=\tfrac45,\quad b=\tfrac15,\quad
 \alpha=\beta=\omega_B=\tfrac25,\quad\gamma=\tfrac3{10},
\]
\[
 \chi=\tfrac3{20},\quad\kappa=\tfrac12,\quad
 \vartheta=\tfrac{141}{400},\quad\nu=2,\quad
 a=\tfrac14,\quad h_O^{\max}=2,\quad\sigma_\xi=1,
\]
\[
 \tau^p=\tfrac{467}{9250},\quad
 y=\tfrac{3521349281}{1677802000},\quad
 \bar H=\tfrac{68104}{68019}.
\]
Set the logistic taste location to
$\bar\xi=\sigma_\xi\log(\pi/(1-\pi))-W^O+W^R$ at the allocations
below, where $\pi=11/21$. This is a primitive specification of a finite taste
distribution, not an equilibrium-dependent policy rule. Hold it fixed during
the credit reform. This construction is an analytical example, not a
calibration.

There is an exact steady state with
\[
 P=Y=O=1,\quad u=\tfrac{9717}{18500},\quad
 T=\tfrac{3975571}{314587875},\quad\pi=\tfrac{11}{21},
\]
\[
 (x_O,h_O,n_O)=(1,1,\tfrac34),\qquad
 (x_R,h_R,n_R)=(\tfrac{99}{74},\tfrac14,\tfrac9{40}).
\]
Mean fertility is exactly $1/2$. Old owner housing is $1480/3239>a$, old
owner consumption is $4/5$, and its estate is $16/25$. The strict original
constraints are checked in the companion receipt. In particular the young
owner's marginal value of space is $.64>u$, whereas the interior old owner's
equals $u$. Thus this same equilibrium also satisfies the short housing
misallocation proposition.

### Exact differentiation

Differentiate the two conditional resource and fertility equations. For
$A_n=\vartheta/n^2+\alpha\kappa^2/s^2$,
\[
 \begin{pmatrix}\rho_m&\chi\\-\chi/x_m^2&A_n\end{pmatrix}
 \binom{dx_m}{dn_m}
 =\binom{dw-h_mdu-u\,dh_m-\mathbf1_{m=R}qa\,du_{t+1}}
              {\alpha\kappa\,dh_m/s_m^2}.
\]
The owner derivative is $dh_O=h_Od\phi/(1-\phi)-h_OdP/P$ and $dh_R=0$.
Envelope derivatives give
\[
 dW_O=\frac{dw-h_Odu}{x_O}
 +(\alpha/s_O-u/x_O)dh_O-\frac{\beta\gamma}{u_{t+1}}du_{t+1},
 \quad dW_R=\frac{dw-a\,du-qa\,du_{t+1}}{x_R},
\]
and $d\pi=\pi(1-\pi)(dW_O-dW_R)/\sigma_\xi$.
All coefficients at the specified equilibrium are rational. Therefore
\[
 J=G_Z-G_vF_v^{-1}F_Z,\qquad
 Q=g_\phi=G_\phi-G_vF_v^{-1}F_\phi
\]
are exact rational matrices. No numerical differentiation is used in their
certificate. Separate complex-step differentiation of the original household
helper checks the derivation.

### Roots, initial conditions and signs

The characteristic polynomial is $z p_5(z)$ up to a positive constant, with
\[
\begin{split}
p_5(z)={}&11624093714366716003479040894500z^5\\
&+258366763596378495144999871594252z^4\\
&-1176628664469803650647365064033229z^3\\
&+408376917290861767134959387855371z^2\\
&-199454780254172995433079109028654z\\
&+694085946210757308455013341344.
\end{split}
\]
An exact Sturm count gives three real roots. Rational isolating intervals
place them respectively near $-26.1500318$, $.003504813$, and $3.5930940$.
Vieta's formulas put the remaining pair's product near $.18132116$ and sum
near $.32660224$. Rational inequalities prove that the pair is nonreal with
modulus in $(.42,.44)$. There are exactly four stable roots including zero,
and two unstable roots.

To check the boundary without rounded complex eigenvectors, let $p_6$ be
the monic characteristic polynomial. A polynomial row $\ell(z)$ constructed
from $J$ satisfies exactly
\[
 \ell(z)(J-zI)=-p_6(z)e_P^\top.
\]
It is nonzero at each unstable root. Those two left eigenvectors annihilate
the stable generalized space. Rational interval Gaussian elimination proves
that the determinant of the matrix whose rows are $B$, $\ell(r_-)$ and
$\ell(r_+)$ excludes zero. After the normalization specified in the driver,
its value is near $.2580033$. This verifies the theorem's boundary condition.

The stationary derivative $Z_\phi^*=(I-J)^{-1}Q$ has
\[
 \frac{dY^*}{d\phi}
 =\frac{38021133978224611827644635285}
 {35866089542871200498536115024}>0.
\]
For the initial derivative $dZ_0$, solve $BdZ_0=0$ and
$\ell(r_\pm)dZ_0=\ell(r_\pm)Z_\phi^*$. Then
$dZ_1=JdZ_0+Q$. A rational enclosure proves
\[
 .0288<\frac{d\bar n_0}{d\phi}<.0289,\qquad
 1.0600<\frac{dY^*}{d\phi}<1.0602.
\]
The full exact enclosures are stored with the driver receipt; the decimal
bounds here are outward bounds for readability. Smoothness gives both signs
for sufficiently small finite positive reforms. Strict margins, root counts,
boundary transversality and strict response signs also persist in an open
neighborhood of the primitives. Small nondegenerate compact-support
perturbations of entrant income and wealth are included, retaining the
uniform conditional branches. This is not a certified rectangular parameter
box or a claim about arbitrary heterogeneity.

## 7. What happens after the initial fertility rise

The constructed transition need not be monotone. In fact every sufficiently
small positive reform eventually places fertility on both sides of
replacement, infinitely often; population similarly crosses its new steady
level infinitely often. The crossings decay exponentially.

Here is a finite-reform argument. At the baseline let
$v=dZ_0-Z_\phi^*$, $w_1=Jv$, and $w_2=J^2v$. If the young-cohort response
had no component from the complex pair, then
$e_Y^\top w_2-r_s e_Y^\top w_1=0$, where $r_s$ is the small real root;
the zero mode has disappeared after one step. The rational interval check
instead places this expression near $-.1346566$, strictly below zero.
Thus the complex modes occur in the young-population response.

For a small reform $\delta$, write the deviations from its new steady state
in smooth local stable coordinates. Let $\lambda(\delta)$ be one of the
complex roots. Uniformly, $|\lambda|\in(.42,.44)$ and the other two stable
roots have modulus below $.01$.

**Uniform tail lemma.** The stable graph and its restricted map may be chosen
$C^2$ jointly in the state and $\delta$. At zero state its nonlinear remainder
and its state derivative are zero. On a sufficiently small common neighborhood
Taylor's theorem therefore bounds the remainder by $C\|v\|^2$, its state
derivative by $C\|v\|$, and its parameter derivative by $C\|v\|^2$.
For the last bound one may use $C^3$ regularity; here the conditional formulas,
implicit map and parameter-dependent stable graph can be taken at least
$C^3$. The fixed compact entrant support allows their first three derivatives
to pass under the integral. Choose a stable-coordinate norm whose linear
contraction is strictly below $.5$. Shrinking the neighborhood preserves a
nonlinear contraction by $\eta=.5$. The initial stable coordinate is a $C^2$
function of $\delta$ and vanishes at zero, so iteration gives
$\|v_t\|\le C|\delta|\eta^t$. Differentiating that iteration, or using the
weighted sequence IFT, gives $\|\partial_\delta v_t\|\le C\eta^t$.
Consequently its nonlinear forcing and derivative are bounded respectively
by $C\delta^2\eta^{2t}$ and $C|\delta|\eta^{2t}$. Constants are uniform
for all sufficiently small $\delta$.

The parameter-dependent graph can also be obtained directly by the one-sided
sequence IFT, prescribing the stable initial coordinate and leaving the
unstable initial coordinate free. That construction requires invertibility
only on the unstable space, so the zero stable root causes no problem here
either. This supplies the regularity used in the lemma. The complex
coordinate therefore satisfies
\[
 a_{t+1}=\lambda a_t+r_t,\qquad
 |r_t|\le C\delta^2\eta^{2t}.
\]
Since $\eta^2<|\lambda|$, the convergent sum
\[
 A(\delta)=a_0+\sum_{s=0}^\infty\lambda^{-s-1}r_s
\]
gives $a_t=A(\delta)\lambda^t+O(\delta^2\eta^{2t})$.
The fast stable coordinates are
$O(|\delta|(.01)^t+\delta^2\eta^{2t})$; unstable coordinates on the local
stable graph are quadratic in the stable coordinates. Consequently
\[
 Y_t-Y^*(\delta)=2\operatorname{Re}
       [C_Y(\delta)\lambda(\delta)^t]
       +O(|\delta|(.01)^t+\delta^2\eta^{2t}).
\]
The preceding nonzero linear signal implies $C_Y'(0)\ne0$, so
$C_Y(\delta)\ne0$ for every sufficiently small nonzero $\delta$.
This differentiability follows from the lemma and differentiation of the
uniform convergent sum. Differentiating $\lambda^{-s-1}$ adds a factor
proportional to $s+1$, which remains summable because
$\eta^2/|\lambda|<1$. The nonlinear summands sum to $O(\delta^2)$.

A nonzero real projection of a rotating nonreal eigenmode has infinitely
many strictly positive and negative values bounded away from zero after
division by its modulus. For a rational rotation this follows over a period;
for an irrational rotation it follows from density on the circle. The
remainder above is smaller along both subsequences. Cohort population
therefore crosses its limiting level infinitely often. Total household
population has the same conclusion, since its leading coefficient is
multiplied by $1+1/\lambda\ne0$. Finally,
\[
 \bar n_t-\frac1\nu=
 \frac{Y_{t+1}-Y_t}{\nu Y_t},
\]
whose leading complex coefficient is multiplied by
$(\lambda-1)/(\nu Y^*)\ne0$. Fertility too approaches replacement from
both sides. This result is compatible with higher initial fertility and a
larger final population. It makes no all-date welfare claim.

## 8. A family of mixed economies

The preceding argument need not rest on one taste scale. Keep the listed
primitives and stationary bundles, let $\sigma=\sigma_\xi>0$, and specify
$\bar\xi(\sigma)=\sigma\log(11/10)-W^O+W^R$. This is a one-parameter
family of primitive taste distributions with the same stationary allocation.
Within each economy, both taste parameters are fixed during the credit
reform. Comparing different $\sigma$ here is not a policy experiment.

**Family result.** Every member with $1\leq\sigma\leq4$ has a locally
unique converging transition after a sufficiently small permanent credit
relaxation. Initial fertility and final population rise. The result persists
under small changes in the other primitives and small entrant heterogeneity
with uniform branch margins. The constructed family's owner share is always
$11/21$; nearby economies need not have that exact share. The oscillation
theorem in section 7 is established near $\sigma=1$, not on the entire
interval. For every $\sigma>0$, the stationary young-cohort derivative is
positive and each tenure's stationary lifetime-value derivative is negative.

The stability part is algebraic for the entire positive half-line. Exact
differentiation gives $J(\sigma)=J(0)+\sigma \mathbf a\mathbf b^\top$ for two constant vectors $\mathbf a,\mathbf b$, and $Q(\sigma)$ affine in $\sigma$. Here $J(0)$ is the algebraic
limit of the matrices; a degenerate taste distribution is not assumed in the
economic theorem. The nonzero-root polynomial is affine in $\sigma$:
\[
 p_{5,\sigma}(z)=p_{5,1}(z)+(\sigma-1)p_1(z),
\]
where $p_{5,1}$ is the integer polynomial in section 6 and
\[
\begin{split}
 p_1(z)={}&209842616525591903123960392500000z^4\\
 &-1024713381846718023989705716419600z^3\\
 &+426827306877372594433334439987600z^2\\
 &-162673088007827714095534626813600z\\
 &+495225810071736959616928665600.
\end{split}
\]
Apply $z=(1+w)/(1-w)$ to this quintic. The first column of the Routh array
for $(1-w)^5p_{5,\sigma}((1+w)/(1-w))$ has signs
\[
 -\,,\quad-\,,\quad-\,,\quad+\,,\quad+\,,\quad-.
\]
Each sign follows because the corresponding rational function of $\sigma$
has numerator coefficients all of the displayed sign and denominator
coefficients all positive. The exact coefficients are saved in the
certificate. No row vanishes. Thus exactly two roots have positive real part
in $w$, and three have negative real part. There are two roots outside the
unit circle in $z$ and three inside, in addition to the zero root of $J$.
The quintic's discriminant is a degree-eight polynomial in $\sigma$ with
all coefficients negative. Together with the three-real-root count at one
point, this proves three distinct real roots and one nonreal pair for every
$\sigma>0$. The signs at $z=-1,1$ and at infinity identify one real root
below $-1$ and one above $1$.

For the boundary and initial fertility sign, the rational interval check
covers the entire closed interval $[1,4]$ with 913 adjacent subintervals.
At each subinterval's endpoints, exact root isolation brackets the two
unstable real roots. A bracket whose polynomial signs agree at both parameter
endpoints brackets the root for every intermediate parameter, since the
polynomial is affine in $\sigma$. This is an interval proof, not a grid of
unchecked intermediate economies.

The rank-one representation provides a polynomial left eigenvector
independent of $\sigma$. For $A=J(0)$, take
$\ell(z)=\mathbf b^\top\operatorname{adj}(zI-A)$. The exact identity
\[
 \ell(z)(J(\sigma)-zI)=-p_{6,\sigma}(z)\mathbf b^\top
\]
follows from the determinant lemma and is also checked coefficient by
coefficient. Interval elimination then establishes a nonzero initial-boundary
determinant and a strictly positive initial fertility derivative throughout
each parameter interval. Adjacent endpoints and full coverage are verified
exactly. The family receipt gives uniform positive lower bounds; the allowed
finite reform size is still not numerically quantified.

There is a simpler direct calculation for the stationary comparison. Write
$D(\sigma)=29694400415045569770000\sigma+7922531612096952274579>0$.
Exact differentiation yields
\[
 \frac{dY^*}{d\phi}=
 \frac{55(625392361132072652215364400\sigma+
                 65900983926556653741810787)}{953456D(\sigma)}>0,
\]
\[
 \frac{dW^{O*}}{d\phi}=
 -\frac{33(568192684411159099406576\sigma+
                 7091548916248680863487)}{224D(\sigma)}<0,
\]
\[
 \frac{dW^{R*}}{d\phi}=
 -\frac{83313(49680061080395802000\sigma+
                 2808938751889938767)}{224D(\sigma)}<0.
\]
These signs hold for every positive $\sigma$. At $\sigma=1$, the value
derivatives are approximately $-2.25302$ for owners and $-.518979$ for
renters. Since both conditional values fall and the additive ownership taste
is fixed, the maximum of the two values falls for every fixed taste draw.
This compares the same entrant in two steady states. It does not rank
different populations socially or describe all transitional cohorts. In
particular, the initial old may benefit from the rise in the value of their
inherited titles. The competitive credit reform does not itself provide the
compensation used in the housing-allocation proposition.

## 9. Verification record

The exact six-variable representation and the actual initial-old boundary
were independently reviewed before the spectral construction. The review
required retaining the full initial distribution for individual feasibility
and using a one-sided sequence theorem because a stable root is zero; both
conditions are included above.

The spectral and sign calculation uses exact rational matrices, exact real
root isolation, and rational interval operations rounded outward to dyadic
multiples of $2^{-180}$. It checks an exact polynomial resolvent identity,
every division, the boundary determinant and strict response signs. All
certificate endpoints are saved. Original-equation finite paths and household
optimizations serve as independent arithmetic checks, not as the convergence
proof. The family certificate adds the exact Routh and discriminant sign
checks, full interval coverage, and original stationary-equation comparisons
at four declared taste scales. Three completed sequential reviews cover the
map, the point certificate/tail argument, and the wider family/welfare result.
The second review's requested explicit uniform tail lemma is included in
section 7 and was checked by the third reviewer. The third review was a
source-and-receipt audit; its model environment lacked `sympy`, so it did not
replay the checker. The lead executed all final modes in the working Python
environment.

Reproduce from the repository root with
`python3 code/model/tools/verify_simplified_olg_mixed_transition.py --smoke`,
then the default command, `--certificate-only --figure`, and `--family-only`.
All four receipts match final source SHA-256
`f22bdc662ca470b763805abe0adec28cb1548c0b921260e57636bd84be08daa4`.
The unchanged old helper's SHA-256 is
`d36f98a38d4f7f39b872830449fffb9aad880f7379bace2ba12265dfbca1c65d`.
The four finite original paths and 24 independent household optimizations
have equilibrium residuals below `2.9e-15` and budget errors below `1.6e-15`.
Initial/final derivative discrepancies are below `3e-9`; original stationary
comparisons across four scales agree with the rational formulas to `6.9e-9`.
The 24/40-date prices agree within `3.2e-15`. All inequalities pass.

The evidence index lists the receipts and review reports. The four-page
reading note and eight-page TeX appendix compile without warnings or
overflows; all twelve pages have been inspected. The author-controlled draft
and all 18 discussion decisions remain unchanged.
```

### File: output/model/simplified_olg_amendments/transition_extensions.md
```md
# Housing allocation and demographic adjustment

Supporting results, September 6–7, 2026. The main theory note remains the
reading copy for discussion. These results retain its preferences, budgets,
tenure constraints and treatment of inherited claims.

There are three additions. The stationary population signs hold over a wider
class of mixed-tenure economies. A simple preference restriction guarantees
local convergence in the all-owner limit. Finally, an explicit finite
fertility decline followed by a later credit reform admits converging paths,
with the intended fertility and population comparisons. The last result is
still restricted to the all-owner limit with zero child goods costs and tax.

## 1. Housing misallocation along the paths

The main allocation argument applies at each date of these transitions.
A young owner with a strictly restrictive down payment values additional
housing more than its current service cost. An old owner with slack retention
and estate bounds values it at that cost. The finite construction below
verifies these inequalities uniformly along both paths.

Consequently, at any chosen date, a small reallocation from old owners to
young owners can improve the young households' welfare while compensating
the old. The original direct-allocation proof preserves fertility, estates
and all future real allocations. This is a separate dated comparison at
each possible intervention date; it does not combine an infinite series of
reallocations or establish the welfare effect of the credit policy.

The same reasoning applies to the existing mixed-tenure local transitions:
the strict household inequalities persist throughout their common
neighborhood. No additional planner power is needed for this corollary.

## 2. Stationary population with both renters and owners

**Conditions.** Entrants have common income and liquid wealth. Set
\(\chi=\tau^p=0\), retaining every other preference and both fixed logistic
ownership-taste parameters. Young owners save positively, have strictly
restrictive down payments and a slack physical cap. Old owners have slack
retention and estate bounds. Young and old renters are strictly at their
rental cap. Suppose also that old owners occupy at least as much housing as
the rental cap.

**Result.** At any stationary equilibrium with these properties, greater
credit availability raises the house price and terminal population. A lower
fertility weight lowers both. This holds for any interior owner share and
any positive taste scale. The stationary root is locally unique and is the
only root on a connected interval retaining the stated household constraints.
The comparisons also hold between finite parameter changes connected by such
an interval.

Here is a short proof. Write \(a=h_R^{\max}\), \(h=h^O>a\),
\(\ell=1-q\), \(d=b/(1-\phi)\), \(w=y+b\), and
\(\rho_O=1+\beta(1+\gamma+\omega_B)\),
\(\rho_R=1+\beta(1+\omega_B)\). The original budgets give:
\[
 P=d/h,\qquad u=\ell P,\qquad
 x_O=\frac{w-\ell d}{\rho_O},\qquad
 x_R=\frac{w-(1+q)ua}{\rho_R}.
\]
Thus owner adult consumption is independent of the stationary price when
the down payment binds. Fertility is proportional to housing in this
zero-child-cost specialization. Replacement therefore fixes mean young
housing, denoted \(B\) only in this proof:
\[
 B=\bar h^Y=\frac{\kappa(\alpha+\vartheta)}{\nu\vartheta}
   =\pi h+(1-\pi)a,\qquad
 C=\frac{h^{2,O}}h=\frac{\beta\gamma x_O}{q\ell d}\in(0,1).
\]
The added housing condition is \(Ch\ge a\). Total housing per young-and-old
pair, and stationary household population, are:
\[
 S=\bar h^Y+\bar h^O=(1+C)B+(1-C)a(1-\pi),\qquad
 N_{\rm hh}^*=\frac{2\bar H}{S}.
\]
An increase in population therefore requires less housing per pair.

The original value difference has a useful exact expression:
\[
 \Delta=W^O-W^R
 =\rho_R\log(x_O/x_R)+(\alpha+\vartheta)\log(h/a)
       +\beta\gamma\log(Ch/a).
\]
The constants common to the two lifetime utilities cancel. Define
\(v=(1+q)ua/x_R\) and \(m=\ell d/x_O=\beta\gamma/(qC)\).
Differentiation with respect to \((h,d,\vartheta)\) gives:
\[
 h\Delta_h=A=\alpha+\vartheta+\beta\gamma-v,\qquad
 -d\Delta_d=F=m+\beta\gamma-v,\qquad
 \Delta_\vartheta=\log(h/a)>0.
\]
The old renter cap and old owner retention imply \(F>0\).
The strictly restrictive down payment implies \(\alpha+\vartheta>m\).
Consequently \(A>F>0\); no new restriction on preference weights is needed.

Replacement and optimal tenure choice form the complete stationary system:
\[
 \pi h+(1-\pi)a=B,\qquad
 \sigma_\xi\log\frac{\pi}{1-\pi}-\bar\xi=\Delta.
\]
Its Jacobian in \((h,\pi)\) has positive determinant:
\[
 \mathcal D=
 \det\begin{pmatrix}
 \pi&h-a\\ -\Delta_h&\sigma_\xi/[\pi(1-\pi)]
 \end{pmatrix}
 =\frac{\sigma_\xi}{1-\pi}+(h-a)\Delta_h>0.
\]
Since \(\pi_h>0\), mean fertility is strictly increasing in \(h\), or
strictly decreasing in \(P=d/h\). This also establishes uniqueness within
a connected admissible branch.

For credit, replacement holds \(B\) fixed. Owner housing rises and the owner
share falls. The elasticity of owner housing is particularly simple:
\[
 E=\frac{d\log h}{d\log d}
   =\frac{F}{A+\sigma_\xi h/(h-B)}\in(0,1).
\]
The price therefore rises, while the ratio of old to young owner housing
falls. Differentiating \(C\) and total housing gives:
\[
 \frac{dC}{d\log d}=-C-\frac{\beta\gamma}{q\rho_O},\qquad
 \frac{dS}{d\log d}
 =\pi h\left[-C-\frac{\beta\gamma}{q\rho_O}
       +\frac{(1-C)a}{h-a}E\right]<0.
\]
The inequality follows from \(Ch\ge a\), which implies
\((1-C)a/(h-a)\le C\). Hence population rises. Larger young-owner homes
coexist with a larger population because fewer households own and old
housing falls; mean young housing is unchanged in this specialization.

For the fertility weight, \(B_\vartheta=-\kappa\alpha/(\nu\vartheta^2)<0\)
and \(C\) is unchanged. Differentiating the same two equations gives
\(h_\vartheta<0\), while the owner-share response can have either sign.
Writing \(L_h=\log(h/a)>0\), old housing satisfies:
\[
 \mathcal D\,\bar h^O_\vartheta
 =B_\vartheta\left[\frac{C\sigma_\xi}{1-\pi}
                   +(Ch-a)\Delta_h\right]
       -(1-C)a\pi L_h<0.
\]
Both ages use less housing as \(\vartheta\) rises. Population and the price
rise. Reversing the change gives the fertility-decline comparison.

**Scope.** Strict signs persist under sufficiently small positive child
goods costs and rebated property tax around each regular economy in this
class. Renters need not be rare. The size of that neighborhood is not
quantified. With general positive child costs, fertility is no longer
proportional to housing and the sign needs additional conditions. The
existing positive-cost counterexample remains valid. These stationary
results do not themselves prove a transition.

An exact half-renter witness and the complete derivative calculations are
preserved in the [review record](transition_extension_reviews.json).
The checker independently differentiates the original lifetime utilities
and verifies the original budgets and inequalities.

## 3. A simple local convergence condition

Consider the all-owner demand limit with \(\chi=\tau^p=0\). Retain the
original strict household constraints. Define the same coefficients as in
the [existing local proof](local_transition_proof.md#10-a-broader-local-transition-result):
\[
 \rho=1+\beta(1+\gamma+\omega_B),\quad D=\frac{\beta\gamma}{\rho},
 \quad L=\frac{\gamma}{1+\gamma+\omega_B},\quad
 r=\frac{y+b}{b/(1-\phi)},\quad
 C=\frac{D[r-(1-q)]}{q(1-q)}.
\]
Here \(r\) is a resource ratio, distinct from dated rent \(r_t\).
The earlier proof gives the exact condition for the required two stable
roots and one unstable root:
\[
 (1-q)+(3+q)C>4D.
\]
A readable sufficient restriction is:
\[
 \boxed{\alpha+\vartheta\le 1+\beta(1+\gamma+\omega_B).}
\]
The combined current housing and fertility weights then do not exceed the
current goods weight plus discounted old-age weights. Indeed, the strict
down-payment constraint implies:
\[
 C>\frac{D\rho}{q(\alpha+\vartheta)}\ge\frac Dq,
 \qquad (3+q)C>4D.
\]
The existing initial-old boundary argument supplies the remaining condition
for a locally unique converging path. This restriction is sufficient, not
necessary. It covers valid economies with \(y+b<b/(1-\phi)\), outside the
earlier convenient resource restriction.

The primitive household restrictions can also be written explicitly:
\[
 \max\left\{\frac{\beta\gamma}{q(\alpha+\vartheta)},
             \frac{L(q-\phi)}{q(1-q)}\right\}<C<1,\quad
 (1-q)\omega_B>q\gamma,\quad
 Kh_O^{\max}>1,
 \qquad K=\frac{\nu\vartheta}{\kappa(\alpha+\vartheta)}.
\]
These enforce restrictive purchase finance, positive saving, slack retention,
slack estate composition, and the owner size limit. Both conditional renter
caps must also remain strictly binding for the extension to positive renter
mass. Exact examples in the receipt verify a convergence margin of
\(847/1160>0\) with \(r=4/5\), and a fully feasible counterexample with
margin \(-1/20\). Household feasibility alone therefore does not imply
convergence.

A small preference decline has the desired initial sign throughout this
convergence region. Normalize the old stationary young cohort and \(K\)
to one, and let \(k\) be the new \(K\). For the actual initial old, the
first housing residual has derivatives:
\[
 F_1=1+\frac{C(1+q)-L}{1-q},\quad
 F_2=-\frac{Cq}{1-q},\quad F_k=\frac L{1-q}-(1+C)<0.
\]
If \(\lambda_1,\lambda_2\) are the two stable roots, the initial cohort
response is:
\[
 \frac{dY_1}{dk}
 =\frac{-F_k-F_2(1-\lambda_1)(1-\lambda_2)}
        {F_1+F_2(\lambda_1+\lambda_2)}>0.
\]
The denominator is positive by the existing boundary argument. The root
product in the numerator is positive for either real stable roots or a
complex conjugate pair; the formula extends continuously to repeated roots.
Since \(K\) rises with \(\vartheta\), a small preference decline lowers
initial fertility and the stationary population. A credit reform's initial
fertility sign still requires its separate condition.

## 4. Finite changes at any later policy date

The following extends the [existing all-owner example](local_transition_proof.md#2-a-limiting-economy-with-strict-household-conditions), with
homogeneous entrants and zero child goods costs and tax. Keep its parameters
and original initial old. Set:
\[
 \vartheta_0=\frac2{15},\qquad
 \vartheta_1=\frac{1998}{15005},\qquad K_1/K_0=\frac{999}{1000}.
\]
The preference decline gives a converging baseline with lower initial
fertility and terminal population. At any later baseline date, every
permanent credit change in
\[
 \frac45<\phi_1\le\frac{81}{101}
 \quad\Longleftrightarrow\quad
 1<d_1/d_0\le\frac{101}{100}
\]
has a converging continuation, higher policy-date fertility than the
continuing baseline, and a larger terminal population. The upper financed
share is about \(0.80198\); the certified change is modest.

**Proof.** Put \(S=K\bar H\) and \(g=(w/d-1)/q>0\). The original fertility
condition, down payment and cohort law imply \(P_t=KdY_t/Y_{t+1}\).
Original housing clearing then becomes, for \(i\ge2\):
\[
 (S-Y_i)\left(Y_{i-1}-qY_i^2/Y_{i+1}\right)
       -D\left(gY_{i-2}Y_i+Y_{i-1}^2\right)=0.
\]
Given its three neighbors, this equation defines \(Y_i\) uniquely in a
population interval \([m,M]\) if:
\[
 0<m<M<S,\quad m^2>qM^2,\quad S-M>2DM,\qquad
 ((1-q)+D)m+DgM\le(1-q)S
       \le((1-q)+D)M+Dgm.
\]
These inequalities give opposite residual signs at the endpoints and a
strictly negative own-coordinate derivative. Its magnitude is bounded below:
\[
 G_{\min}=m-qM^2/m+2qm(S-M)/M+Dgm.
\]
A sufficient upper bound on the sum of the three response magnitudes is:
\[
 \lambda=
 \frac{DgM+S-m-2Dm+q(S-m)M^2/m^2}{G_{\min}}<1.
\]
For an interval of credit changes, use the conservative extrema of \(g\)
in this bound.

The first date uses the actual old, with total financial claims \(A_0\)
and inherited titles \(H_0\):
\[
 (S-Y_1)(Y_0-qY_1^2/Y_2)
      -L[(A_0/d)Y_1+KY_0H_0]=0.
\]
For this equation require opposite signs at \((m,m)\) and \((M,M)\), and:
\[
 G_0=Y_0-qY_1^2/Y_2+2qY_1(S-Y_1)/Y_2+LA_0/d>0,\qquad
 \frac{q(S-Y_1)Y_1^2}{Y_2^2G_0}\le\lambda_0<1.
\]
These ensure a unique first-date root and bound its response to \(Y_2\).
In particular,
the old mortgage is not recalculated at the new financed share. For an
intervention during the baseline, let \(U,Y,V\) be baseline young populations
at the preceding, intervention and following dates. Actual inherited claims
are exactly:
\[
 H_0=Y/K_1,\qquad
 A_0=\frac{(\rho-1)(w-d_0)}{q\rho}U
       -\frac{d_0}{\rho}\frac{Y^2}{V}.
\]
Here \(V\) retains the forecast under which the old chose saving. It is
not replaced by the policy realization.

Updating each population coordinate from its equation is a contraction on
the complete space of infinite sequences in the interval. Hence there is
a unique fixed point there. For its tail distance \(e_t\) from the stationary
population, the three-neighbor bound gives \(e_t\le\lambda e_{t-2}\).
The path therefore converges exponentially. This argument has no terminal
date or imposed terminal allocation.

For \(k=999/1000\), the baseline lower bound is exactly
\(m_B=k-Dg(1-k)/((1-q)+D)\). The certificate checks the baseline interval
\([0.9989614726\ldots,1]\) and a common policy interval \([0.997,1.004]\).
The respective interior derivative bounds are below \(0.660\) and \(0.689\);
the first-date bounds are below \(0.374\) and \(0.392\).
It verifies initial-old feasibility, every original owner inequality, and
the conditional renter caps throughout the intervals. These are fraction
calculations with positive margins, not pointwise sampling.

The initial baseline residual evaluated at \(Y_1=Y_2=1\) is
\((K_1-1)[(1-q)(1+C)-L]<0\), proving the initial fertility decline.
Since \(C(d)\) decreases with credit, stationary population rises with \(d\).
For policy-date fertility, implicit differentiation of the infinite
contraction, keeping \(Y_0,A_0,H_0\) fixed, gives:
\[
 \frac{17}{1000}<\frac{\partial Y_1}{\partial d}<\frac{51}{1000}.
\]
The bound holds uniformly over all inherited baseline states and
\(d\in[1,101/100]\). The checker obtains it by outward interval substitution;
the unretained tail always keeps its unrestricted uniform derivative bound.
It does not substitute a terminal steady state. Integrating the inequality
and using \(n_0=Y_1/(\nu Y_0)\) proves the finite fertility comparison.

**Scope.** This is an explicit finite range in the all-owner limit. A broad
finite transition theorem with material renting and general positive child
costs remains open. The stationary result in section 2 does not close that
gap. Neither result proves an all-date fertility ordering or a welfare gain
from the competitive credit reform.

## Verification

Run from the project root:

~~~sh
PYTHONDONTWRITEBYTECODE=1 python3 code/model/tools/verify_simplified_olg_transition_extensions.py
~~~

The [receipt](transition_extension_checks.json) separates symbolic
original-equation identities, exact finite-interval inequalities, and
floating-point checks against the original household helper.
The [review record](transition_extension_reviews.json) preserves the three
bounded research reports and the separate review of the finite proof's
financial settlement and infinite derivative calculation.

The main TeX/PDF, its two figures, earlier proofs, quantitative code and
author decisions are unchanged.
```

### File: output/model/simplified_olg_amendments/local_transition_proof.md
```md
# A local analytical transition in the simplified OLG model

September 5, 2026. Supporting derivation for the compact assessment. This uses
the original dated household problems with zero gains tax. It is not a
calibration, an implementation of the direct planner, or a global transition
theorem. Original flow utilities and value functions are unchanged.

## 1. What the result can establish

There is a nonempty class of original-model economies with heterogeneous
entrant income and wealth, positive child goods costs, positive renter mass
and ordinary property tax in which a small permanent increase in the financed
share has a locally unique exponentially converging equilibrium. Initial
fertility rises under the stated additional parameter condition, and the final
population is larger. Fertility and tenure remain household choices.

Section 10 gives a general primitive stability inequality for the all-owner,
zero-child-goods-cost and zero-property-tax limit. Section 10.2 gives both an
exact primitive initial-fertility test and a more elementary sufficient one.
The smooth extension admits positive renter mass, costs and taxes near each
regular limiting configuration; section 9 explains heterogeneous entry and
verifies a nondegenerate admissible support. The allowed neighborhood and
reform must be small, but their numerical size has not been certified. This
is a local analytical theorem, not a global continuation result.

The explicit example in sections 2–4 has a closed-form first-order path. Its
population rises at every date and fertility is above replacement throughout
adjustment. Section 6 proves these signs uniformly for sufficiently small
finite reforms in that limiting economy. These stronger all-date signs are
not asserted throughout the broader parameter family or for positive renters.
Indeed section 10 gives an admissible case with lower initial fertility and a
larger final population. Section 8 separately treats the stationary population
sign with positive child goods costs. Welfare is a different comparison:
section 7 records lower stationary household welfare after the uncompensated
credit reform in the plotted limiting example.

## 2. A limiting economy with strict household conditions

Take homogeneous income and wealth, \(q=4/5,\phi=4/5,b=1/5,y=19/20\),
\(\alpha=\beta=2/5,\gamma=3/10,\omega_B=2,\kappa=1/2\),
\(\vartheta=2/15,\nu=2\), owner cap 2 and rental cap \(1/4\).
Initially let \(\chi=\tau^p=0\) and the owner share tend to one. Define
\(d=b/(1-\phi)\), \(w=y+b=23/20\),
\(A=1+\gamma+\omega_B=33/10\), and
\(\rho=1+\beta A=58/25\). The variables \(A,\rho,d\) are temporary
coefficients for this derivation, not replacements for the original value functions.

Choose \(\bar H=1213/928\). The limiting stationary allocation is
\[
P=1,\quad Y=O=1,\quad
(x,h,n,a')=(95/232,1,1/2,627/1160),
\]
\[
a=-301/928,\qquad(c^2,h^2,e)=(95/464,285/928,475/928).
\]
The post-mortgage financial state \(a\) can be negative; the original
nonnegative-saving constraint applies to \(a'\). Saving is positive,
\(h^2<h<2\), \(e>Ph^2\), and
\(MV^Y=19/87>u=1/5\), so purchase rationing is strict. All log arguments
are strictly positive. Concavity and the original first-order conditions
establish global conditional owner optimality.

Conditional renting is also regular at these prices, even though its mass is
zero in the limit. Its exact choices are
\[
(x^R,h^R,n^R,a'_R)=(53/110,1/4,1/8,34/55),
\]
\[
a_R=17/22,\qquad(c^{2R},h^{2R},e^R)=(53/220,1/4,53/88).
\]
The young rental housing value is \(848/825>1/5\); the old value is
\(159/550>1/5\). Thus both rental caps bind strictly and saving remains
positive. In a common neighborhood the conditional owner and renter solutions,
values and old choices are smooth, with all these strict inequalities retained.

## 3. The exact limiting recurrence and the initial old

At \(\chi=0\), the original fertility condition gives \(n_t=h_t/2\).
Purchase rationing and the cohort law therefore give
\[
h_t=d/P_t,\qquad Y_{t+1}=Y_t h_t,\qquad
P_t=dY_t/Y_{t+1}.
\]
The original dated owner budgets imply
\[
\rho x_t=w-d+qd(P_{t+1}/P_t),\qquad
h_{t+1}^2=\frac{\beta\gamma x_t}{q(P_{t+1}-qP_{t+2})}.
\]
Write \(D=\beta\gamma/\rho=3/58\) and \(r=w/d\). Housing clearing
for every \(t\ge1\) is exactly
\[
\bar H=Y_{t+1}+
D\frac{[(r-1)/q]Y_{t-1}+Y_t^2/Y_{t+1}}
 {Y_t/Y_{t+1}-qY_{t+1}/Y_{t+2}}. \tag{T1}
\]
The date-zero old retain their pre-reform asset and title choices. With
\(Y_0=O_0=1\), initial clearing instead requires
\[
F_0=Y_1+\frac1{11}
 \frac{-301/928+d/Y_1}{d/Y_1-qdY_1/Y_2}-\bar H=0. \tag{T2}
\]
Using (T1) at date zero would silently revalue the initial old's predetermined
financial assets. Equation (T2) avoids that error.

The nearby stationary allocation satisfies
\[
P^*=d,\quad x^*=[w-(1-q)d]/\rho,\quad
h^{2*}=\frac{75}{232}[w/d-(1-q)],\quad
Y^*=\frac{\bar H}{1+h^{2*}}. \tag{T3}
\]
In particular \(\partial_dY^*=345/1213>0\) at the displayed baseline.

## 4. Roots, initial conditions, and the analytical response

Linearizing (T1) at the baseline gives the characteristic polynomial
\[
1140z^3-3253z^2+945z-45=0. \tag{T4}
\]
Rational endpoint evaluations bracket its three roots in
\((.059,.060),(.261,.262),(2.53,2.54)\), respectively. Thus
\(0<\lambda_1<\lambda_2<1<\Lambda\). No numerical stability assumption
is required. Their decimal values are .05958565, .26160613 and 2.53231699.

At the baseline, (T2) has derivatives
\[
(F_1,F_2,F_d,F_H)
=(33783/10208,-285/232,1505/10208,-1).
\]
The stable initial-condition matrix is
\[
B=\begin{pmatrix}
1&1\\
F_1\lambda_1+F_2\lambda_1^2&F_1\lambda_2+F_2\lambda_2^2
\end{pmatrix}. \tag{T5}
\]
It is nonsingular: its determinant is
\((\lambda_2-\lambda_1)[F_1+F_2(\lambda_1+\lambda_2)]>0\), using the
rational root brackets. Numerically it is .58886862. The derivative of (T1)
with respect to its last dated quantity is nonzero. The local stable-manifold
argument and (T5) therefore select a unique local converging path, among paths
remaining in the chosen neighborhood, for sufficiently small reforms.

For a unit increase in \(d\), the young-population derivative is
\[
y_t=J+c_1\lambda_1^t+c_2\lambda_2^t,\qquad
J=345/1213,\quad c_1=.87792408\ldots,\quad c_2=-1.16234288\ldots . \tag{T6}
\]
Here \(c_1+c_2=-J\), and the second row of (T5) applied to
\((c_1,c_2)\) is \(-(F_1+F_2)J-F_d\). In particular \(y_0=0\).
The coefficients are defined by those exact equations; their printed decimals
are not the proof. The root brackets give
\[
y_{t+1}-y_t
\ge m\lambda_2^t>0,\quad
m=-c_1(1-\lambda_1)-c_2(1-\lambda_2)=.03265445\ldots . \tag{T7}
\]
Initial fertility has derivative \(m/2>0\), because
\(n_t=Y_{t+1}/(2Y_t)\). Old population follows one period later.
For an explicit rational certificate, write \(S=\lambda_1+\lambda_2\).
The boundary equations imply
\[
m=\frac{-F_d-F_2J(1-\lambda_1)(1-\lambda_2)}{F_1+F_2S}.
\]
The stated rational root brackets bound this below by
\(293918959/9027813150>0\). Thus the initial sign is not inferred
from rounded coefficients. Fertility rises most in the second young cohort
in this example; above-replacement fertility does not mean fertility is
monotone over dates.

For a proportional stock increase, \(\partial_aY^*=1\), and the transient
coefficients are -.01324704 and -.98675296. Both are negative, so the same
limiting first-order population comparison is increasing at every date.

## 5. Extending path existence to the original mixed model

Fix \(\sigma_\xi=1\) and write the original logistic share as
\[
\pi_t^O=\frac1{1+\epsilon\exp[(W_t^R-W_t^O)/\sigma_\xi]},
\qquad \epsilon=\exp(-\bar\xi/\sigma_\xi).
\]
Every \(\epsilon>0\) corresponds to a finite ownership-taste location and
strictly positive renter and owner masses. The zero endpoint is only a demand
limit; no assertion about continuity of unbounded taste-utility levels is needed.

Let \(\mu=(\epsilon,\chi,\tau^p)\). The steady-state equations can be
written in \((P,T)\): replacement fertility and the ordinary per-household
rebate. Their Jacobian at \(\mu=0,P=1,T=0\) is diagonal with entries
\((-1,1)\). Hence pre- and post-reform steady states are smooth near the
limit. The state supplied by the pre-reform old is a finite vector
\[
I_0=(Y^{\rm pre},O^{\rm pre},\pi^{O,\rm pre},
       a^{O,\rm pre},H^{O,\rm pre},a^{R,\rm pre}).
\]
Each component is a smooth function of the stationary prices and the strict
conditional policies in section 2. Initial old housing demand is a smooth
function of this vector and current prices. The normalized distribution of
old states need not be differentiable in total-variation norm; what matters
is its smooth entry into the actual boundary equation. Thus the original
pre-reform renter mass and assets are included rather than held at zero.

To accommodate the added price lead in tenure choice, use sequences instead
of asserting an unchanged-dimensional dynamic map. Center prices and young
population at the post-reform steady state and give the sequences the norm
\(\|v\|_\eta=\sup_{t\ge0}|v_t|/\eta^t\), with \(\eta=2/5\).
Eliminate rebates using
\(T_t=q\tau^pP_t\bar H/(Y_t+O_t)\), where \(O_0\) is predetermined
and \(O_t=Y_{t-1}\) thereafter. Use the original two conditional young
problems, their generated old asset/title states, actual old optimization,
housing clearing, and the demographic law at every date. These define an
operator on two price/population sequences, with the initial population and
old-state boundary included.

**Bounded inverse at the limit.** Linearized demography eliminates the price
sequence by a bounded difference operator. The remaining operator is the
forced version of (T4), plus two initial conditions in (T5). In companion
coordinates its matrix has eigenvalues \(\lambda_1,\lambda_2,\Lambda\).
For any residual sequence \(g\) with finite \(\eta\)-norm, the stable
coordinates satisfy forward sums of the form
\[
s_t=M_s^ts_0+\sum_{j=0}^{t-1}M_s^{t-1-j}b_sg_{j+1},
\]
and the unique weighted-bounded unstable coordinate is
\[
u_t=-\sum_{j=t}^{\infty}\Lambda^{t-1-j}b_ug_{j+1}.
\]
The sums are bounded by geometric constants proportional to
\((\eta-\lambda_2)^{-1}\) and \((\Lambda-\eta)^{-1}\).
The two boundary conditions determine \(s_0\) through (T5), after the
known \(u_0\) and residual terms are included. This constructs a unique
bounded inverse for arbitrary boundary and dated residuals. Recovering prices
uses only another bounded shift/difference.

**Smooth perturbation.** Every finite forward or backward shift is bounded
in this norm. All conditional policies and values are uniformly smooth on
the compact feasible neighborhood from section 2. Changes in \(\mu\),
including the \(P_{t+2}\) terms in tenure choice, therefore perturb the
derivative by a small bounded operator. Pointwise products are continuous in
the weighted space; second-order remainders decay at least as fast as the
product of the sequence deviations. The initial operator is smooth through
the explicitly verified finite vector \(I_0(\mu)\). Ordinary smooth
extensions across the zero parameter boundary may be used for the theorem;
the claimed economies have positive parameters and probabilities.

The bounded inverse, a Neumann-series bound for its perturbation and the
implicit-function theorem give a locally unique infinite exponentially
converging path for all sufficiently small positive \(\mu\) and sufficiently
small reforms. All original inequalities and both tenure masses remain
admissible uniformly in time. This proof does not confuse a terminally closed
finite computation with an infinite path.

Since \(\partial_dY^*=J>0\) and \(\partial_dn_0=m/2>0\) at the limit,
both inequalities persist in a common positive-parameter neighborhood. A
sufficiently small positive credit reform therefore raises initial fertility
and final household population in that original-model neighborhood. Common
stationary adult/child counting conventions give the same final person ranking.

## 6. Uniform finite-reform signs in the limiting economy

This section is deliberately restricted to the all-owner limit. Work at the
post-reform steady state, so its eigenvalues \(\lambda_j(\delta)\) move
with the reform rather than being incorrectly held at baseline values. The
two-dimensional stable manifold has coordinates whose dynamics are
\[
v_{t+1}=\operatorname{diag}(\lambda_1,\lambda_2)v_t+R(v_t),
\qquad R(v)=O(|v|^2).
\]
Its embedding has linear population coordinate \(v_1+v_2\), with a quadratic
remainder. Uniformly for small reforms, the selected initial point is
\(O(\delta)\), \(|v_t|\le C\delta\eta^t\), and
\(\lambda_1<\eta^2<\lambda_2<\eta\), using \(\eta=.4\).

Forward summation gives
\(v_{1t}=v_{10}\lambda_1^t+O(\delta^2\eta^{2t})\).
For the second coordinate the convergent sum
\[
B_2=v_{20}+\sum_{j\ge0}\lambda_2^{-1-j}R_2(v_j)
\]
gives \(v_{2t}=B_2\lambda_2^t+O(\delta^2\eta^{2t})\). The boundary
transversality implies \(v_{10}/\delta\to c_1\) and
\(B_2/\delta\to c_2\). Consequently the population increments equal
the two-mode increments with these continuous coefficients, plus
\(O(\delta^2\eta^{2t})\). By (T7) the two-mode part is at least
\((m/2)\delta\lambda_2^t\) for sufficiently small positive reforms.
Since \(\eta^2<\lambda_2\), the error is smaller uniformly over every
date when \(\delta\) is sufficiently small. Thus
\(Y_{t+1}>Y_t\) at all dates and \(n_t>1/2\), while fertility tends to
replacement. The stock-change argument uses its two negative coefficients
and the same remainder bound.

## 7. Welfare and verification limits

This is a demographic response to an uncompensated credit change. It is not
the direct-allocation Pareto comparison. In the limiting steady state young
housing and replacement fertility are unchanged by \(d\), while adult
consumption and old housing fall. In fact the common owner's lifetime value
has derivative
\[
\partial_dW^O=-\frac{1-q}{x^*}-\frac{\beta\gamma}{d}
=-289/475<0
\]
at the displayed baseline. The ownership taste is constant in this reform
comparison and does not change that derivative. A larger population therefore
does not establish a welfare improvement.

`verify_simplified_olg_local_transition.py` checks the polynomial by symbolic
differentiation and rational sign brackets. Its eight declared cases reconstruct
the dated original problems, initial-old asset states, market/rebate/cohort
equations and strict private conditions. Finite positive/negative perturbations
match (T6); 24- and 40-date terminal closures agree to numerical precision.
These are independent arithmetic checks, not the proof of sections 5–6 or a
certificate of the size of their parameter neighborhood. The new figure uses
(T6) directly, not the finite terminally closed paths. Earlier figures and
the original author manuscript are untouched by the verification driver.

## 8. A stationary population condition with positive child goods costs

The all-owner zero-property-tax comparison also admits an explicit condition
with \(\chi>0\). It is a conditional statement on the same strictly binding
purchase/interior-old branch, not a claim about every tenure configuration.
Let \(\ell=1-q\), \(d=b/(1-\phi)\),
\(\rho=1+\beta(1+\gamma+\omega_B)\),
\(w_\chi=y+b-\chi/\nu\), and \(g=\beta\gamma/(q\ell)\).
Replacement fertility fixes \(n^*=1/\nu\), so the original budget gives
the adult-consumption quantity directly from primitives:
\[
X=\frac{w_\chi-\ell d}{\rho}>0.
\]
The fertility condition and housing budgets then give
\[
h^*=\frac{\kappa}{\nu}
\frac{\nu(\alpha+\vartheta)X-\chi}{\nu\vartheta X-\chi},\quad
P^*=d/h^*,\quad h^{2*}=gXh^*/d,\quad
Y^*=\frac{\bar H}{h^*(1+gX/d)}. \tag{T8}
\]
Require \(\nu\vartheta X>\chi\) and the original private inequalities.
Strict purchase rationing itself reduces to the primitive inequality
\[
(\alpha+\vartheta)X>\ell d+\chi/\nu.
\]
No equilibrium price needs to be solved for in these expressions.

Differentiating (T8), with all other primitives fixed, gives
\[
\frac{\partial\log Y^*}{\partial d}=
\frac1\rho\left[
 \frac{g w_\chi}{d(d+gX)}-
 \frac{\alpha\nu\chi\ell}
 {(\nu\vartheta X-\chi)[\nu(\alpha+\vartheta)X-\chi]}
\right]. \tag{T9}
\]
Hence population rises precisely when the first bracketed term exceeds the
second. Every term is explicit in primitives. The first captures the declining
ratio of old to young housing. The second captures the rise in young housing
needed to maintain replacement fertility as the goods available for children
fall. At \(\chi=0\) the second term vanishes, so the sign is
positive. With \(\chi>0\), it need not be.

For the parameters in section 2 with \(\chi=.015\), (T9) is .2565939892,
and the original purchase/old constraints remain strict. Reducing the old
housing weight can reverse the sign: with \(\gamma=.001\), \(\chi=.015\)
and an adequately small rental cap for the unused renter problem, the same
formula is negative while the owner purchase limit still binds. These are
stationary comparative statics, not a claim that the transition in that second
economy has already been proved. The strict sign and stationary regularity
also permit small positive renter shares by a finite-dimensional implicit
function argument, separately from the dynamic argument in section 5.

## 9. Extending the result to income and wealth heterogeneity

The homogeneous-entrant restriction can be relaxed without changing the
limiting aggregate recurrence. In this section all preferences, including
ownership tastes, retain their original form. Let the fixed entrant
distribution of income and liquid wealth have support in
\[
 y_i\in[.94,.96],\qquad b_i\in[.195,.205],
 \qquad E[y_i]=.95,\quad E[b_i]=.2.
\]
Write this fixed probability distribution as \(F\); the type measure of a
young cohort is \(Y_tF\). It may have atoms or a density and may correlate
income with wealth. Ownership tastes have the original independent logistic distribution.
Only \(\epsilon,\chi,\tau^p\) approach zero, as in section 5; this does not
require the income and wealth dispersion to approach zero.

**Exact aggregation at the limit.** Write \(d_i=b_i/(1-\phi)\) and
\(w_i=y_i+b_i\). As long as all conditional owner branches remain strict,
\[
 h_{it}=d_i/P_t,\qquad n_{it}=h_{it}/2,\qquad
 \rho x_{it}=w_i-d_i+qd_iP_{t+1}/P_t.
\]
Every one of these quantities is linear in \((w_i,d_i)\). Old owner housing
is \(\beta\gamma x_{i,t-1}/(q u_t)\), also linear. Average saving,
post-mortgage assets and inherited title are linear as well. It follows that
housing clearing, mean fertility and the initial-old demand equation are
exactly (T1)–(T2), with \(w_i,d_i\) replaced by their means. The stated mean
restrictions reproduce the same baseline and the same linearized operator.
The initial old cohort is exactly the normalized stationary distribution
generated by pre-reform type-specific choices under \(F\), multiplied by its
predetermined cohort mass. The financial-asset distribution has not been replaced with its mean in a
nonlinear problem: old demand on this branch is linear in those resources.

**Uniform strictness on the support.** At \(P=1\), the owner housing-value
gap is strict when \(w_i/d_i>107/100\), or equivalently
\(y_i>(87/20)b_i\). The worst corner satisfies
\(.94-(87/20).205=.04825>0\). Housing, saving, old consumption,
retention slack and estate slack are affine in \((y_i,b_i)\) on this branch;
positive slack at all four corners establishes it throughout the rectangle.
Conditional renter quantities are likewise affine in resources at \(\chi=0\)
with both rental caps binding, and their marginal housing values exceed user cost at
all four corners. Thus a common compact neighborhood of prices and parameters
preserves every original inequality for all types.

**Smooth aggregation for positive renters and costs.** For each type, the
conditional young solution and value, tenure share, and generated old state
are smooth functions of the finite set of dated prices and rebates used in
the original two-period model. The functions and their first two derivatives
are uniformly bounded on the compact type/price neighborhood above. Their
integrals are therefore differentiable by differentiation under the integral
sign, with the derivative bounded by the supremum of the individual bound.
At date zero, integrate the same smooth old-demand functions over the actual
pre-reform conditional assets, titles and tenure shares. This is a smooth
map of the pre-reform stationary prices and parameters; no total-variation
smoothness of an atomic old-state distribution is needed.

Consequently the equilibrium operator on the weighted sequence space of
section 5 is continuously differentiable. At \(\epsilon=\chi=\tau^p=0\),
it is exactly the same aggregate operator for every distribution in this
class. Its inverse and the two strict comparative-static signs have already
been established. Uniform bounds on individual derivatives make a small
positive \((\epsilon,\chi,\tau^p)\) a small operator perturbation after
integration as well. The same implicit-function argument gives the local
converging path and positive initial-fertility and final-population responses
for these heterogeneous-entrant economies. An arbitrary broad distribution
with households crossing constraint boundaries is outside this result.

This extension broadens the homogeneous-entrant result in section 1; it does not
change the exact plotted first-order path, because the limiting economy
aggregates exactly. It also does not prove all-date fertility signs in the
positive-renter case or provide a numerical size for the admissible parameter
neighborhood.

## 10. A broader local transition result

This section replaces the example-specific stability calculation with a
condition in primitives. It retains \(\chi=\tau^p=0\) and the all-owner
demand limit for the initial derivation. A fixed compact entrant distribution
is allowed when all original conditional owner and renter branches are
uniformly strict, as in section 9. All preferences except these expressly
zero limiting parameters remain original.

Define the following coefficients from primitives and entrant means:
\[
 \ell=1-q,\quad A=1+\gamma+\omega_B,\quad \rho=1+\beta A,
 \quad D=\frac{\beta\gamma}{\rho},\quad L=\frac\gamma A,
\]
\[
 d=\frac{E[b_i]}{1-\phi},\quad w=E[y_i+b_i],\quad r=\frac wd,
 \quad C=\frac{D(r-\ell)}{q\ell},\quad
 K=\frac{\nu\vartheta}{\kappa(\alpha+\vartheta)}.
\]
Here \(C\) equals mean old housing divided by mean young housing at the
limiting steady state, but is calculated directly from primitives. Mean young
housing is \(1/K\), mean old housing is \(C/K\), and
\(Y^*=K\bar H/(1+C)\). Strict individual old retention implies
\(0<C<1\); the strict estate condition is \(\ell\omega_B>q\gamma\),
which implies \(0<D<L<\ell\). All young purchase, saving and physical
housing conditions remain explicit requirements. For homogeneous entrants,
the positive purchase multiplier requires
\((\alpha+\vartheta)(r-\ell)>\rho\ell\).

Assume in addition \(r>1\), meaning \(w>d\). Thus
\(0<D<\ell C\). The following result does **not** require \(q\ge1/2\).

### 10.1 Convergence under the existing strict branches and \(w>d\)

Normalize young population by its initial stationary level. After multiplying
housing clearing by \(K\), the exact recurrence is (T1) with arbitrary
\(q,D,r\). Its linearized polynomial, up to a nonzero factor, is
\[
 \mathcal P(z)=Cqz^3-[\ell-D+C(1+q)]z^2+(C-2D)z+D-\ell C.
 \tag{T10}
\]
It satisfies
\[
 \mathcal P(0)<0,\qquad \mathcal P(1)=-\ell(1+C)<0,
 \qquad \mathcal P'(1)=-\ell(2+C)<0.
\]
The positive leading coefficient and negative derivative at one imply exactly
one simple real root \(\Lambda>1\). On \([1,\infty)\), the cubic first
decreases to its upper critical point and then increases through this root.

The product of the other two roots is \(R/\Lambda\), where
\[
 R=\frac{\ell C-D}{qC}>0,\qquad
 \mathcal P(R)=
 -\frac{\ell(\ell C-D)[qC^2+\ell C-D]}{C^2q^2}<0.
\]
If \(R\le1\), then \(R<\Lambda\). If \(R>1\), the unique root above
one and \(\mathcal P(R)<0\) again imply \(R<\Lambda\). Thus this
product is strictly between zero and one. A complex conjugate pair is
therefore strictly inside the unit circle.

If the two remaining roots are real, they have the same sign. Positive roots
are below one because \(\Lambda\) is the only root above one. For negative
roots, note
\[
 \mathcal P(-1)=4D-\ell-(3+q)C
 <(1-5q)C-\ell<0.
\]
For \(q\ge1/5\), the last inequality is immediate. For \(q<1/5\), it
follows from \(C<1\), giving the upper bound \(-4q<0\). Exactly one
root below minus one would give \(\mathcal P(-1)>0\); two would have
product above one. Both are impossible. Thus the two real roots are also
strictly inside the unit circle. Repeated roots inside the circle are allowed.

The initial-old boundary has derivatives
\[
 F_1=1+\frac{C(1+q)-L}{\ell},\qquad
 F_2=-\frac{Cq}{\ell},\qquad
 F_d=\frac L\ell-C. \tag{T11}
\]
The reform variable here is a proportional increase in \(d\), and initial
young and old cohort masses both normalize to one. If \(S\) is the sum of
the two stable roots, then \(-2<S<2\) and
\[
 F_1+F_2S=1+\frac{C(1+q-qS)-L}{\ell}>1-L/\ell>0.
\]
This is the boundary determinant after removing the distinct-root factor.
For a repeated root, the basis \(\lambda^t,\partial_\lambda\lambda^t\) gives this
same nonzero factor directly. This polynomial-derivative notation includes a
zero stable root when section 10.4 is used. Hence the initial cohort and its predetermined
old state meet the stable solutions transversely. The local stable-manifold
argument gives a unique local converging equilibrium from the pre-reform
steady state for sufficiently small reforms. The weighted inverse can be
constructed with any decay weight between the stable spectral radius and one;
Jordan blocks cause no difficulty for the geometric bounds.

### 10.2 Initial fertility: an exact primitive test and a simple sufficient one

The stationary normalized young-population derivative for a proportional
increase in \(d\) is
\[
 J=\frac{Dr}{q\ell(1+C)}>0.
\]
Write \(m\) for the normalized initial population increment derivative;
initial fertility has derivative \(m/\nu\). The two boundary equations give
\[
 m=\frac{-F_d-F_2J(1-\lambda_1)(1-\lambda_2)}{F_1+F_2S}
 =\frac{C-L/\ell+J(1+C)/(\Lambda-1)}{F_1+F_2S}. \tag{T12}
\]
The second identity uses
\((1-\lambda_1)(1-\lambda_2)=\ell(1+C)/[Cq(\Lambda-1)]\).
The denominator is positive, but initial fertility need not rise.

There is an exact test involving only primitives and a polynomial evaluation.
If \(C\ge L/\ell\), initial fertility rises. If \(C<L/\ell\), define
\[
 Z=1+\frac{Dr}{q(L-\ell C)}>1.
\]
Then \(m>0\) precisely when \(\mathcal P(Z)>0\). Indeed (T12) requires
\(\Lambda<Z\), and the cubic has just one root above one. The test
requires no equilibrium price, borrowing multiplier, or computed dynamic root.

A more elementary sufficient inequality is available when \(q\ge1/2\):
\[
 \ell C[\ell+(1+q)C]>L(\ell-D+C). \tag{T13}
\]
To see this, put \(\Pi=R/\Lambda\). Since \(q\ge\ell\) and
\(C>D/\ell\), the coefficient identity
\(S=[(C-2D)/(Cq)-\Pi]/\Lambda\) implies \(S>0\). Hence
\(\Lambda<[\ell-D+C(1+q)]/(Cq)\). Using this upper bound in (T12)
yields the lower bound
\[
 C-L/\ell+J(1+C)/(\Lambda-1)
 >C-L/\ell+\frac{DrC}{\ell(\ell-D+C)}.
\]
Positivity of the last expression is exactly (T13). At the explicit baseline
its left-minus-right gap is \(229247/47365120>0\). This is a sufficient
condition, not the exact fertility boundary. The existence argument in 10.1
has no restriction \(q\ge1/2\).

### 10.3 Positive renters and costs, and the limit of the conclusion

Under the uniform conditional-policy regularity in section 9, the section 5
argument applies around **every** limiting parameter configuration satisfying
the strict conditions in 10.1. Choose a decay weight above this configuration's
stable spectral radius and below one. Retain the continuously differentiable
actual-pre-reform-old boundary and uniformly bounded operator conditions,
including every dated lead from tenure choice. It gives
local convergence for small positive renter mass, child goods costs and
property tax. The positive terminal-population derivative persists. If the
exact test or strict sufficient condition in 10.2 also holds, positive initial
fertility persists as well. The initial old are always the stationary cohort
generated before the reform. No broad-distribution mean-state approximation
is used once tenure selection matters.

This remains a local result in the added parameters and reform size. Section
6's all-date finite-reform sign uses the particular positive-root example and
is not asserted for the whole family. Stable complex roots can produce
oscillating fertility during convergence.

For a concrete counterexample to an unconditional initial-fertility claim,
use \(q=.4,\phi=.8,b=.2,y=1,\alpha=1.2,\vartheta=.4,\beta=1/12\),
\(\gamma=.3,\omega_B=1.7,\kappa=.5,\nu=2\), rental cap .025,
owner cap 2, \(\bar H=1.05\), and the all-owner \(\chi=\tau^p=0\)
limit. At \(P=Y=O=1\), young choices are \((x,h,n,a')=(.48,1,.5,.52)\)
and old choices are \((c^2,h^2,e)=(.1,.05,.425)\). All original owner
inequalities are strict except the binding purchase limit, with
\(MV^Y=.768>u=.6\); both conditional rental caps can strictly bind at
the stated small cap. Here \((C,D,L)=(.05,.02,.1)\), \(J=2/21>0\),
but \(m=-.11952631\ldots<0\). The local path converges with a larger
final population and lower fertility initially. This is a distinct constructed
economy, not the transition illustrated in Figure 2.

### 10.4 A sharper convergence condition, with a second algebraic proof

The sufficient restriction \(r>1\) can itself be relaxed. Keep the strict
household branch, so \(0<C<1\) and \(0<D<L<\ell\), and require only
\(r>\ell\) for positive adult consumption. The exact condition for (T10)
to have two roots strictly inside the unit circle and one root above one is
\[
 \ell+(3+q)C>4D. \tag{T14}
\]
It is automatic under \(r>1\), by the bound on \(\mathcal P(-1)\)
in 10.1. It is also automatic if \(4D\le\ell\). Neither sufficiency
statement replaces the maintained original household inequalities.

Here is an independent short root proof. Put \(z=(1+v)/(1-v)\), which maps
\(\operatorname{Re}v<0\) onto \(|z|<1\). Direct expansion gives
\[
 (1-v)^3\mathcal P\!\left(\frac{1+v}{1-v}\right)
 =A_3v^3+A_2v^2+\ell(C-1)v-\ell(C+1),
\]
\[
 A_3=\ell+(3+q)C-4D,
 \qquad A_2=\ell+4D+(7q-3)C.
\]
Under (T14), the coefficient signs are positive, unrestricted, negative,
negative. Descartes' rule therefore gives exactly one positive real root
\(v_+\). It lies below one because the polynomial is negative at zero
and equals \(8Cq>0\) at one. The other two roots have positive product:
\(\ell(C+1)/(A_3v_+)>0\). Their sum is negative, because the sum of
pairwise products is \(\ell(C-1)/A_3<0\). Thus both remaining roots
have negative real parts, whether real or complex. Transforming back gives
one \(\Lambda>1\) and exactly two stable roots. This includes repeated
stable roots. The map has no root at its excluded point \(v=1\).

Conversely, if (T14) is an equality, \(\mathcal P(-1)=0\), so there is a
unit root at minus one. If it is reversed, \(\mathcal P(-1)>0\), while
the positive leading cubic tends to minus infinity as \(z\to-\infty\).
There is then a real root below minus one in addition to the root above one.
Thus (T14) is necessary as well as sufficient for precisely two stable roots
under these branch conditions.

The boundary argument (T11) uses only \(S<2\) and \(L<\ell\), so it is
unchanged. Formula (T12), the exact primitive fertility test, positive
stationary population derivative, and the mixed/heterogeneous perturbation
argument likewise need only (T14), not \(r>1\). Choose a decay weight strictly above the stable spectral radius and below
one, retaining the section 5/9 continuously differentiable actual-pre-reform-old
boundary and uniformly bounded-operator assumptions, including every dated
tenure lead. Uniqueness concerns nearby exponentially converging sequences,
not global equilibria. The more convenient sufficient fertility inequality
(T13) is retained with its stated extra restrictions \(q\ge1/2,r>1\). An all-date sign remains restricted to
the particular example in sections 4–6.

## 11. What credit changes in the limiting steady state

The stationary welfare loss in the plotted example is general within the
maintained all-owner, zero-child-goods-cost and zero-property-tax demand
limit. This result is independent of the convergence proof. It compares a
household of the same type entering each steady state, with the world bond
price, preferences and entrant distribution fixed.

Let \(d_i=b_i/(1-\phi)\), \(\bar d=E[d_i]\),
\(\rho=1+\beta(1+\gamma+\omega_B)\), and
\(K=\nu\vartheta/[\kappa(\alpha+\vartheta)]\). Retain positive saving,
strictly restrictive purchase constraints, slack owner physical caps, and
interior old choices with slack retention and estate constraints, uniformly
over the compact type support as in section 9. The original
fertility condition, purchase constraint and replacement fertility give
\[
 P^*=K\bar d,\qquad
 h_i^*=\frac{d_i}{K\bar d},\qquad
 n_i^*=\frac{b_i}{\nu E[b_i]},\qquad
 x_i^*=\frac{y_i+b_i-(1-q)d_i}{\rho}. \tag{T15}
\]
To obtain the price, use \(n_i=\vartheta h_i/[\kappa(\alpha+\vartheta)]\)
and integrate \(h_i=d_i/P\) against the fixed probability distribution
\(F\). The condition \(\nu E[n_i]=1\) then gives \(P=K\bar d\).
A common change in the financed share scales every \(d_i\) and its mean
proportionally. It therefore leaves each young household's housing and
fertility unchanged at the new steady state. It raises the housing price and
reduces adult consumption:
\[
 \frac{\partial\log P^*}{\partial\phi}=\frac1{1-\phi},\qquad
 \frac{\partial x_i^*}{\partial\phi}
 =-\frac{(1-q)d_i}{\rho(1-\phi)}<0.
\]
These signs are stationary comparisons, not statements about prices or housing
at every transitional date.

From the original old budget \(c^2+qe+qr h^2=a+PH\), the interior solution is
\[
 c_i^{2*}=\frac{\beta x_i^*}{q},\qquad
 h_i^{2*}=\frac{\beta\gamma x_i^*}{q(1-q)P^*},\qquad
 e_i^*=\frac{\beta\omega_B x_i^*}{q^2}.
\]
In particular the estate is dated one period after old consumption; its
denominator is \(q^2\), not \(q\). Substitution into the original separate
flow utilities shows that all terms in \(W^{O*}(i)\) that vary with
\(\phi\) are \(\rho\log x_i^*-\beta\gamma\log P^*\). Consequently
\[
 \boxed{\frac{\mathrm d W^{O*}(i)}{\mathrm d\phi}
 =-\frac{(1-q)d_i/x_i^*+\beta\gamma}{1-\phi}<0.} \tag{T16}
\]
The same derivative applies to realized ownership utility when the household's
taste \(\xi_i\) is held fixed. No finite limit of diverging taste levels is
being taken. The comparison concerns conditional owner values and fixed
individual tastes at the all-owner demand limit.

Every old household type also uses less housing:
\[
 \frac{\mathrm d\log h_i^{2*}}{\mathrm d\phi}
 =-\frac{1+(1-q)d_i/(\rho x_i^*)}{1-\phi}<0. \tag{T17}
\]
Mean young housing is unchanged and mean old housing falls. Housing clearing
therefore gives a larger stationary household population. This is how price
capitalization and old downsizing accompany the demographic gain in Figure 2.
It is not an implementation of the compensated Pareto comparison in the
allocation result. Neither a welfare statement about transitional cohorts nor
a ranking of populations under a social objective follows from (T16).

The derivative applies on a sufficiently small interval that preserves the
original household conditions. Positive \(\chi\), positive property tax,
and material tenure switching require separate welfare derivatives. The
existing mixed-economy convergence theorem does not automatically establish
the same uniform type-by-type welfare ordering there.

Verification: `--welfare-only` in the existing transition driver differentiates
the original separate utility terms symbolically, checks all dated budgets
and inequalities for five fixed types at three financed shares, and compares
the derivatives with central differences of the original utility. It performs
ten independent original household optimizations at the central share.
The maximum derivative discrepancy is below \(8.3\times10^{-9}\).
`local_transition_scope_review.md` supplies the independent derivation and
economic assessment. Its discussion of an estate restriction should be read
as a qualification on the interior-demand formula: the original constraint
\(e\geq P_{t+1}h^2\) is not a minimum-retention rule and does not physically
prevent an old owner from selling housing.
```

### File: code/model/tools/verify_simplified_olg_transition_extensions.py
```python
#!/usr/bin/env python3
"""Verify the supporting transition extensions without changing the main note.

These are analytical theory checks, not calibration or a finite-horizon
simulation. The receipt separates exact identities from floating-point checks
against the original household equations.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "output/model/simplified_olg_amendments"
from verify_simplified_olg_local_transition import (
    complex_jacobian,
    parameters,
    young_choices,
)

def symbolic_checks():
    q, beta, gamma, omega, alpha, theta, kappa, w, d, h, a = sp.symbols(
        "q beta gamma omega alpha theta kappa w d h a", positive=True
    )
    ell = 1 - q
    rho_o = 1 + beta * (1 + gamma + omega)
    rho_r = 1 + beta * (1 + omega)
    price = d / h
    x_o = (w - ell * d) / rho_o
    x_r = (w - (1 + q) * ell * price * a) / rho_r

    def value(x, house, old_house):
        n = theta * house / (kappa * (alpha + theta))
        return (
            sp.log(x) + alpha * sp.log(house - kappa * n)
            + theta * sp.log(n)
            + beta * (sp.log(beta * x / q) + gamma * sp.log(old_house)
                      + omega * sp.log(beta * omega * x / q**2))
        )

    delta = value(x_o, h, beta * gamma * x_o / (q * ell * price)) - value(x_r, a, a)
    v = (1 + q) * ell * price * a / x_r
    m = ell * d / x_o
    A = alpha + theta + beta * gamma - v
    F = m + beta * gamma - v
    identities = {
        "owner_value_minus_renter_value_h": sp.diff(delta, h) - A / h,
        "owner_value_minus_renter_value_d": sp.diff(delta, d) + F / d,
        "owner_value_minus_renter_value_theta": sp.expand_log(
            sp.diff(delta, theta) - sp.log(h / a), force=True
        ),
    }
    C = beta * gamma * x_o / (q * ell * d)
    identities["old_housing_ratio_derivative"] = d * sp.diff(C, d) + C + beta * gamma / (q * rho_o)
    x, c, k, C0, L = sp.symbols("x c k C L", positive=True)
    inherited_assets_over_d = ell * C0 / L - 1
    initial = x + L * (inherited_assets_over_d * x + k) / (1 - q * x**2 / c) - k * (1 + C0)
    base = {x: 1, c: 1, k: 1}
    identities["actual_old_initial_x"] = sp.diff(initial, x).subs(base) - (1 + (C0 * (1 + q) - L) / ell)
    identities["actual_old_initial_c"] = sp.diff(initial, c).subs(base) + C0 * q / ell
    identities["actual_old_preference_shift"] = sp.diff(initial, k).subs(base) - (L / ell - 1 - C0)
    # Reconstruct the finite recurrence and surprises from the dated budgets.
    K, stock, lag, prev, nxt, yy, assets, titles, d0 = sp.symbols(
        "K stock lag prev nxt yy assets titles d0", positive=True
    )
    Ppast, Pnow, Pfuture = K*d*lag/prev, K*d*prev/x, K*d*x/nxt
    xp = (w-d+q*d*Pnow/Ppast)/rho_o
    clearing = x/K + lag*beta*gamma*xp/(q*(Pnow-q*Pfuture)) - stock
    residual = (K*stock-x)*(prev-q*x*x/nxt) - beta*gamma/rho_o*((w/d-1)/q*lag*x+prev**2)
    identities["finite_interior_housing_clearing"] = clearing*K*(prev-q*x*x/nxt)+residual
    Pnow, Pfuture = K*d*yy/x, K*d*x/nxt
    clearing0 = x/K+L*(assets+Pnow*titles)/(Pnow-q*Pfuture)-stock
    residual0 = (K*stock-x)*(yy-q*x*x/nxt)-L*(assets*x/d+K*yy*titles)
    identities["finite_actual_old_housing_clearing"] = clearing0*K*(yy-q*x*x/nxt)+residual0
    prior_x = (w-d0+q*d0*yy**2/(lag*nxt))/rho_o
    inherited = lag*(w-d0-prior_x)/q
    identities["baseline_generated_assets"] = inherited-((rho_o-1)*(w-d0)*lag/(q*rho_o)-d0*yy**2/(rho_o*nxt))
    identities["finite_interior_credit_forcing"] = sp.diff(residual,d)-beta*gamma/rho_o*w*lag*x/(q*d**2)
    identities["finite_initial_credit_forcing"] = sp.diff(residual0,d)-L*assets*x/d**2
    for name, expression in identities.items():
        assert sp.simplify(expression) == 0, name
    return {name: "exact zero after differentiation of the original budgets/utilities" for name in identities}

def original_margins(p, price=1.0):
    hh = young_choices([price] * 3, [0.0, 0.0], p)
    q, u = p["q"], (1 - p["q"]) * price
    residuals, margins = {}, {}
    for tenure in ("owner", "renter"):
        own = tenure == "owner"
        x, h, n, saving, c2, h2, estate = hh[tenure]["z"]
        cost = (1 - p["phi"]) * price * h if own else u * h
        residuals[tenure + "_young_budget"] = x + p["chi"] * n + saving + cost - p["y"] - p["b"]
        residuals[tenure + "_old_budget"] = c2 + q * estate + u * h2 - hh[tenure]["assets"] - (price * h if own else 0)
        residuals[tenure + "_fertility"] = p["theta"] / n - p["chi"] / x - p["alpha"] * p["kappa"] / (h - p["kappa"] * n)
        margins[tenure + "_positive_bundle"] = min(x, h - p["kappa"] * n, n, saving, c2, h2, estate)
        margins[tenure + "_young_housing_gap"] = p["alpha"] * x / (h - p["kappa"] * n) - u
        if own:
            residuals["down_payment"] = (1 - p["phi"]) * price * h - p["b"]
            margins.update(owner_physical=p["owner_cap"] - h,
                           owner_retention=h - h2, owner_estate=estate - price * h2)
        else:
            margins["old_renter_housing_gap"] = p["gamma"] * c2 / h2 - u
    assert max(abs(z) for z in residuals.values()) < 1e-10, residuals
    assert min(margins.values()) > 0, margins
    return hh, {"maximum_original_equation_residual": max(abs(z) for z in residuals.values()),
                "strict_margins": margins}

def limiting_witnesses():
    cases = {
        "convergent_with_w_less_than_d": parameters(y=.6, alpha=.9, theta=.3, rental_cap=.05, Hbar=277/232),
        "feasible_but_wrong_stable_dimension": parameters(q=.5, y=.35, beta=.5, gamma=1., omega=2.,
                                                         alpha=24., theta=8., rental_cap=.01, Hbar=31/30),
    }
    result = {}
    for name, p in cases.items():
        hh, checks = original_margins(p)
        q, ell = p["q"], 1 - p["q"]
        rho = 1 + p["beta"] * (1 + p["gamma"] + p["omega"])
        D = p["beta"] * p["gamma"] / rho
        r = (p["y"] + p["b"]) / (p["b"] / (1 - p["phi"]))
        C = D * (r - ell) / (q * ell)
        gap = ell + (3 + q) * C - 4 * D
        exact_gap = sp.Rational(847, 1160) if name.startswith("convergent") else sp.Rational(-1, 20)
        rq, rb, rg, ro, ry, rwealth, rphi = [sp.Rational(str(p[key])) for key in ("q", "beta", "gamma", "omega", "y", "b", "phi")]
        rrho = 1 + rb*(1+rg+ro)
        rD = rb*rg/rrho
        rr = (ry+rwealth)/(rwealth/(1-rphi))
        rC = rD*(rr-(1-rq))/(rq*(1-rq))
        assert sp.simplify((1-rq)+(3+rq)*rC-4*rD-exact_gap) == 0
        assert abs(gap - float(exact_gap)) < 1e-13
        roots = np.roots([C*q, -(ell-D+C*(1+q)), C-2*D, D-ell*C])
        assert sum(abs(z) < 1 for z in roots) == (2 if gap > 0 else 1)
        checks.update(parameters=p, r=r, C=C, exact_convergence_gap=str(exact_gap),
                      stable_root_count=int(sum(abs(z) < 1 for z in roots)))
        result[name] = checks
    return result

def mixed_stationary_checks():
    result = []
    for scale in (.25, 1., 4.):
        p = parameters(q=.5, beta=.4, gamma=.3, omega=.4, alpha=.2, theta=7/15,
                       kappa=7/8, y=99/50, Hbar=99/100, rental_cap=.25, sigma=scale)
        base = young_choices([1.]*3, [0., 0.], p)
        p["taste_weight"] = np.exp((base["owner"]["utility"]-base["renter"]["utility"])/scale)
        hh, checks = original_margins(p)
        assert abs(hh["pi"]-.5) < 1e-13 and abs(hh["fertility"]-.5) < 1e-13

        def original_outputs(v):
            h, d, theta = v
            varied = dict(p, phi=1-p["b"]/d, theta=theta)
            H = young_choices([d/h]*3, [0., 0.], varied)
            old_h = H["pi"]*H["owner"]["z"][5] + (1-H["pi"])*p["rental_cap"]
            return np.array([H["fertility"], H["housing"], old_h, H["pi"], d/h])

        J = complex_jacobian(original_outputs, [1., 1., p["theta"]])
        derivatives = {}
        for shock, col, conversion in (("phi", 1, 5.), ("theta", 2, 1.)):
            h_derivative = -J[0, col] / J[0, 0]
            full = (J[:, col] + J[:, 0]*h_derivative)*conversion
            base_output = original_outputs([1., 1., p["theta"]])
            S = base_output[1] + base_output[2]
            assert abs(S-p["Hbar"]) < 1e-12
            N_derivative = -2*p["Hbar"]*(full[1]+full[2])/S**2
            assert abs(full[0]) < 1e-12 and N_derivative > 0 and full[4] > 0
            if shock == "phi":
                assert 0 < h_derivative < 1 and full[3] < 0 and abs(full[1]) < 1e-12 and full[2] < 0
            else:
                assert h_derivative < 0 and full[1] < 0 and full[2] < 0
            derivatives[shock] = dict(N=float(N_derivative), hOwner=float(h_derivative*conversion),
                                     hYoung=float(full[1]), hOld=float(full[2]), pi=float(full[3]), P=float(full[4]))
        checks.update(taste_scale=scale, taste_location=float(-scale*np.log(p["taste_weight"])),
                      owner_share=float(hh["pi"]), derivatives=derivatives)
        result.append(checks)
    return result

def finite_transition_certificate():
    """Exact infinite-sequence bounds; no finite terminal condition."""
    from fractions import Fraction as F

    q, beta, gamma, omega = F(4,5), F(2,5), F(3,10), F(2)
    alpha, nu, kappa, b, y = F(2,5), F(2), F(1,2), F(1,5), F(19,20)
    ell = 1-q
    A = 1+gamma+omega
    rho = 1+beta*A
    D, L, w, Hbar = beta*gamma/rho, gamma/A, y+b, F(1213,928)
    k = F(999,1000)                         # K after the preference shock; Kpre=1.
    theta = alpha*k/(nu/kappa-k)
    S = k*Hbar
    dlo, dhi = F(1), F(101,100)
    elo, ehi = (w/dhi-1)/q, (w/dlo-1)/q
    mB, MB = k-D*ehi/(ell+D)*(1-k), F(1)
    mP, MP = F(997,1000), F(251,250)
    a_pre, h_pre = -F(301,928), F(1)
    report = {'theta1':str(theta), 'phi1_max':str(1-b/dhi),
              'baseline_box':list(map(float,(mB,MB))),
              'policy_box':list(map(float,(mP,MP)))}

    def full_box(name, dmin, dmax, m, M):
        """Uniform in d, all dates, and every neighbor in the population box."""
        emin, emax = (w/dmax-1)/q, (w/dmin-1)/q
        assert emin>0 and 0<m<M<S
        rmin, rmax = (m/M)**2, (M/m)**2
        # E_b>0 gives exact vertex tests for E(m)>=0 and E(M)<=0.
        slacks = {
          'user_cost': m-q*M*M/m,
          'E_b': S-M-2*D*M,
          'self_lower': ell*S-(ell+D)*m-D*emax*M,
          'self_upper': (ell+D)*M+D*emin*m-ell*S,
        }
        assert all(v>=0 for v in slacks.values())
        assert slacks['user_cost']>0 and slacks['E_b']>0
        # Bounds for G=-E_x and the three absolute row coefficients.
        Gmin = m-q*M*M/m+2*q*m*(S-M)/M+D*emin*m
        Gmax = M-q*m*m/M+2*q*M*(S-m)/m+D*emax*M
        ac = (D*emin*m/Gmax, D*emax*M/Gmin)
        bc = ((S-M-2*D*M)/Gmax, (S-m-2*D*m)/Gmin)
        cc = (q*(S-M)*m*m/(M*M*Gmax), q*(S-m)*M*M/(m*m*Gmin))
        row = ac[1]+bc[1]+cc[1]
        assert Gmin>0 and 0<row<1
        # Original young and generated-old owner conditions, checked over d.
        for d in (dmin,dmax):
            xmin=(w-d+q*d*rmin)/rho
            xmax=(w-d+q*d*rmax)/rho
            vals={
              'adult_consumption':xmin,
              'saving':y-xmax,                       # a'=y-x, since closing cost=b.
              'purchase':(alpha+theta)*xmin-d*(1-q*rmin),
              'physical_owner_cap':F(2)-M/(k*m),
              'old_retention':q*d*rmin*(1-q*rmax)-beta*gamma*xmax,
              'old_estate':omega-q*(omega+gamma)*rmax,
            }
            assert all(v>0 for v in vals.values()), (name,d,vals)
            # Each d-dependent inequality above is affine in d, so endpoints suffice.
            for key,val in vals.items():
                slacks[key]=min(slacks.get(key,val),val)
        # Conditional renter caps also stay strict; not needed for the all-owner map.
        rentcap=F(1,4)
        rhoR=1+beta*(1+omega)
        Pmax=k*dmax*M/m
        U=Pmax*(1-q*rmin)
        cash=w-rentcap*(1+q)*U
        renter_margins=(cash,(alpha+theta)*cash-rhoR*rentcap*U,
                        beta*gamma*cash-rhoR*q*rentcap*U)
        assert all(v>0 for v in renter_margins)
        report[name]={'row_bound':float(row),
                     'slacks':{key:float(val) for key,val in slacks.items()},
                     'renter_cap_margins':list(map(float,renter_margins))}
        return Gmin,Gmax,ac,bc,cc

    base_coeffs=full_box('baseline',dlo,dlo,mB,MB)
    pol_coeffs=full_box('policy',dlo,dhi,mP,MP)

    def E0(x,c,Y,assets,d,titles):
        # assets and titles are TOTAL original old claims, not new-policy choices.
        return (S-x)*(Y-q*x*x/c)-L*((assets/d)*x+k*Y*titles)

    # The first preference shock retains the actual Section2 stationary old.
    base_lower=E0(mB,mB,F(1),a_pre,dlo,h_pre)
    base_upper=E0(MB,MB,F(1),a_pre,dlo,h_pre)
    assert base_lower>0 and base_upper<0
    assert base_upper==(k-1)*(ell*(1+D*(w-ell)/(q*ell))-L)

    # At ANY later baseline date, write U=Y_{t-1}, Y=Y_t, V=Y_{t+1}^{baseline}.
    # Actual claims: Atot=cA U - Y^2/(rho V); Htot=Y/k.
    # V is the OLD baseline forecast, even when policy unexpectedly changes prices.
    cA=(rho-1)*(w-dlo)/(q*rho)
    Atot_min=cA*mB-dlo*MB*MB/(rho*mB)
    Atot_max=cA*MB-dlo*mB*mB/(rho*MB)
    assert Atot_min<Atot_max<0
    # E0 increases in Y once Htot=Y/k is substituted, and decreases in Atot.
    assert S-MP-2*L*MB>0
    policy_lower=E0(mP,mP,mB,Atot_max,dhi,mB/k)
    policy_upper=E0(MP,MP,MB,Atot_min,dlo,MB/k)
    assert policy_lower>0 and policy_upper<0
    report['initial_boundary']={
      'baseline_lower':float(base_lower),'baseline_upper':float(base_upper),
      'policy_lower':float(policy_lower),'policy_upper':float(policy_upper),
      'actual_total_assets_range':list(map(float,(Atot_min,Atot_max)))}

    def boundary_coeffs(name,m,M,Ymin,Ymax,amin,amax,dmin,dmax):
        # Here amin,amax bound TOTAL assets for the first-update derivative.
        assert amin<=amax<0
        Gmin=Ymin-q*M*M/m+2*q*m*(S-M)/M+L*amin/dmin
        Gmax=Ymax-q*m*m/M+2*q*M*(S-m)/m+L*amax/dmax
        cc=(q*(S-M)*m*m/(M*M*Gmax), q*(S-m)*M*M/(m*m*Gmin))
        ff=(L*amin*M/(dmin*dmin*Gmin), L*amax*m/(dmax*dmax*Gmax))
        assert Gmin>0 and 0<cc[0]<=cc[1]<1 and ff[0]<=ff[1]<0
        report[name+'_boundary_row']=float(cc[1])
        return cc,ff

    boundary_coeffs('baseline',mB,MB,F(1),F(1),a_pre,a_pre,dlo,dlo)
    c0,f0=boundary_coeffs('policy',mP,MP,mB,MB,Atot_min,Atot_max,dlo,dhi)

    # Original initial-old choices: c2=R/A, h2=L R/u, e=omega R/(q A).
    # Verify R>0, h2<H, e>Pnext*h2 for every actual inherited household.
    rbmin,rbmax=(mB/MB)**2,(MB/mB)**2
    individual_amin=(w-dlo-(w-dlo+q*dlo*rbmax)/rho)/q
    individual_amax=(w-dlo-(w-dlo+q*dlo*rbmin)/rho)/q
    individual_hmin,individual_hmax=mB/(k*MB),MB/(k*mB)
    for name,m,M,Ymin,Ymax,dmin,dmax,amin,amax,hmin,hmax in [
      ('baseline',mB,MB,F(1),F(1),dlo,dlo,a_pre,a_pre,h_pre,h_pre),
      ('policy',mP,MP,mB,MB,dlo,dhi,individual_amin,individual_amax,
       individual_hmin,individual_hmax)]:
        # Independent bounds intentionally allow incompatible corners: conservative.
        Pmin=k*dmin*Ymin/M; Pmax=k*dmax*Ymax/m
        Pnextmax=k*dmax*M/m
        umin=Pmin-q*Pnextmax
        resource_min=amin+Pmin*hmin
        resource_max=amax+Pmax*hmax
        retention=hmin-L*resource_max/umin
        estate=omega*umin-q*gamma*Pnextmax
        assert umin>0 and resource_min>0 and retention>0 and estate>0
        report[name+'_actual_old']={
          'resource_min':float(resource_min), 'retention_slack':float(retention),
          'estate_slack':float(estate)}

    # Differentiate at FIXED actual inherited Y,A,H. z0=0 and, for i>=2,
    # zi=-ai*z(i-2)+bi*z(i-1)+ci*z(i+1)+fi; z1=c0*z2+f0.
    # Uniform coefficient intervals come directly from E and E0, for all dates/d.
    Gmin,Gmax,ac,bc,cc=pol_coeffs
    fc=(D*w*mP*mP/(q*dhi*dhi*Gmax),D*w*MP*MP/(q*dlo*dlo*Gmin))
    row=max(ac[1]+bc[1]+cc[1],c0[1])
    B=max(fc[1],-f0[0])/(1-row)
    assert row<1 and B>0
    # The inverse (I-DT)^-1 exists on bounded sequences, so every zi is in [-B,B].
    # Repeated interval substitution preserves enclosure of the INFINITE solution.
    # The unretained tail remains [-B,B]; no terminal steady-state value is imposed.
    DYADIC=2**80

    def down(x):return F((x.numerator*DYADIC)//x.denominator,DYADIC)
    def up(x):return -down(-x)
    def add(a,b):return down(a[0]+b[0]),up(a[1]+b[1])
    def mul(a,b):
        products=[x*y for x in a for y in b]
        return down(min(products)),up(max(products))
    def neg(a):return -a[1],-a[0]

    N=100
    z=[(F(0),F(0))]+[(-B,B)]*(N+1)
    for _ in range(100):
        new=[(F(0),F(0)),add(f0,mul(c0,z[2]))]
        for i in range(2,N+1):
            new.append(add(add(fc,neg(mul(ac,z[i-2]))),
                           add(mul(bc,z[i-1]),mul(cc,z[i+1]))))
        new.append((-B,B))
        z=new
    assert z[1][0]>F(17,1000) and z[1][1]<F(51,1000)
    report['policy_impact_derivative']={
      'coefficient_row_bound':float(row), 'derivative_norm_bound':float(B),
      'dY1_dd_interval':list(map(float,z[1])),
      'claimed_exact_enclosure':['17/1000','51/1000'],
      'iterations':100,'retained_dates':100,'outward_dyadic_bits':80}

    C0=D*(w/dlo-ell)/(q*ell)
    C1=D*(w/dhi-ell)/(q*ell)
    assert S/(1+C0)==k and S/(1+C1)>k
    report['stationary_young_population']={
      'pre':1.0,'baseline':float(k),'largest_certified_credit':float(S/(1+C1))}
    report["policy_impact_derivative"]["exact_dY1_dd_interval"] = list(map(str, z[1]))
    report["exact_population_boxes"] = {"baseline": list(map(str, (mB, MB))), "policy": list(map(str, (mP, MP)))}
    return report

def main():
    report = {
        "scope": "Supporting theory only. No calibration, finite-horizon simulation, new planner power, or main-note revision.",
        "symbolic_original_equation_checks": symbolic_checks(),
        "finite_transition_certificate": finite_transition_certificate(),
        "limiting_branch_witnesses": limiting_witnesses(),
        "mixed_stationary_original_equation_checks": mixed_stationary_checks(),
    }
    report["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    evidence = [
        "latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex",
        "output/pdf/simplified_olg_amendment_proposal.pdf",
        "output/model/simplified_olg_amendments/theory_slides_misallocation.pdf",
        "output/model/simplified_olg_amendments/combined_transition_figure.pdf",
        "output/model/simplified_olg_amendments/transition_extensions.md",
        "output/model/simplified_olg_amendments/transition_extension_reviews.json",
        "code/model/tools/verify_simplified_olg_local_transition.py",
    ]
    report["evidence_sha256"] = {
        name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in evidence
    }
    target = OUT / "transition_extension_checks.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"verification": "passed", "receipt": str(target)}, indent=2))

if __name__ == "__main__":
    main()
```
