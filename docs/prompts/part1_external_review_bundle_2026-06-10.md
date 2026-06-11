# External Review Bundle: Part I Analytical Model

Date: 2026-06-10

This is a self-contained bundle for a new ChatGPT Pro or Claude chat. Treat the LaTeX source below as the maintained current draft of standalone Part I. The working audit and figure source are included after the draft.

## Role

You are a senior academic macroeconomist and economic theorist, writing/refereeing at the level of a top JPE/QJE/AER/RESTUD-style theory or quantitative macro paper. Audit this as economics, not as surface prose.

## Project Context

The paper studies intergenerational housing allocation, tenure choice, rental-owner segmentation in family-sized housing, and fertility. The author wants the model section to read like a clear macro theory paper: environment, households, choices, constraints, aggregate equilibrium, planner, CE/planner comparison, then fertility and policy implications.

Important preferences:

- Start from the economy, not from wedges or mysterious technical objects.
- Avoid decorative terminology and avoid “wedge-first” exposition.
- Keep the economics fixed unless you identify a real inconsistency.
- Use Coven et al.-style aggregate housing equilibrium as the reference point: tenure choices, segmented markets, construction/supply, market clearing, then planner comparison.
- The child-input primitive is Stone-Geary: children require both consumption goods and housing services.
- The core fertility mechanism should not depend on old-side property-tax lock-in. Property tax can be an old-side retention force, but young constrained access should stand on its own.

## Current State Of The Draft

Recent changes already incorporated:

1. Young households choose tenure through rental and owner branch problems.
2. Rental and owner housing have different size menus, with large family-sized units more available in the owner segment.
3. Construction and aggregate clearing are explicit by segment.
4. Planner feasibility distinguishes real menus from private financial implementation constraints.
5. Section 9 now contains a local policy sufficient statistic. The central object is the effective family-housing cap \(H_i^m\). A policy affects fertility through this mechanism only if it relaxes the binding component of \(H_i^m\).
6. Two conceptual figure sketches exist: one for the branch mechanism and one for the sufficient-statistic decomposition.

## What I Want From You

Please provide a rigorous audit in six parts:

1. **Executive diagnosis:** Is this version now a coherent Part I for a macro theory/quantitative paper? What are the 5-8 highest-priority remaining issues?
2. **Mathematical audit:** Check the Stone-Geary child input structure, the Section 9 derivative, the equal-scaling condition \(\chi_i=\kappa q^m/\alpha\), and the sufficient-statistic formulas. Correct any mistakes.
3. **Economic audit:** Are tenure choice, segmentation, construction, old retention, and planner feasibility coherently integrated? Are real menu constraints and private financial constraints treated correctly?
4. **Policy statistic:** Is the local sufficient statistic useful and correctly qualified? Can it be made simpler, sharper, or closer to what a referee would accept?
5. **Figures:** Assess whether the two figure concepts are the right two pictures. Suggest better titles, labels, or alternative figures if needed.
6. **Revision patches:** Give compact LaTeX-ready replacement text only where useful. Do not rewrite the whole document.

Be explicit about what is true only locally/fixed-price/fixed-active-set versus what can be stated in full equilibrium.

---

# A. Current Standalone Part I LaTeX Source

```latex
\documentclass[11pt]{article}

\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,amsthm,mathtools}
\usepackage{booktabs,tabularx,array}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{float}
\usepackage{enumitem}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage[round]{natbib}
\usepackage{microtype}

\graphicspath{{figures/}}
\hypersetup{colorlinks=true,linkcolor=blue!50!black,citecolor=blue!50!black,urlcolor=blue!50!black}
\setlength{\parindent}{0em}
\setlength{\parskip}{0.45em}
\setlist[itemize]{leftmargin=1.3em,itemsep=0.12em,topsep=0.12em}
\setlist[enumerate]{leftmargin=1.45em,itemsep=0.12em,topsep=0.12em}
\setlength{\emergencystretch}{3em}
\allowdisplaybreaks

\newtheorem{proposition}{Proposition}
\newtheorem{corollary}{Corollary}
\newtheorem{definition}{Definition}
\newtheorem{remark}{Remark}
\newtheorem{assumption}{Assumption}

\newcommand{\dd}{\mathrm d}
\newcommand{\E}{\mathbb E}
\newcommand{\1}{\mathbf 1}
\newcommand{\KK}{\mathcal K}
\newcommand{\HH}{\mathcal H}
\newcommand{\LL}{\mathcal L}

\title{Intergenerational Housing Misallocation and Fertility\\[-0.15em]
\large Compact Analytical Model}
\author{Draft note}
\date{June 9, 2026}

\begin{document}
\maketitle

\begin{abstract}
This note develops a compact non-spatial theory of intergenerational housing allocation and fertility. Young households decide whether to enter the modeled housing market and, conditional on entry, choose tenure, housing, consumption, and fertility. Rental and owner housing differ in their size menus, and young households face borrowing and payment constraints when they try to buy owner housing. Old households are incumbent owners who choose consumption, retained owner housing, and bequests. The model defines a competitive equilibrium and a constrained planner, then compares marginal allocation rules for family-capable housing.
\end{abstract}

\textbf{Overview.} The compact model first describes the economy, the rental and owner housing segments, and the household problems. It then defines equilibrium and the planner's problem. The main allocation result concerns constrained access to family-capable housing among entered young households. Some young households are unconstrained; others cannot upsize because rental size menus, down-payment constraints, or payment-to-income constraints bind. Property-tax assessment lock-in enters only as a concrete old-side retention mechanism. It can amplify the allocation gain from releasing owner housing, but it is not necessary for the young-access mechanism.

\begin{remark}[Status]
This standalone file contains only the compact analytical model. It is copied from the v4 draft and revised to focus on tenure choice, rental-owner size segmentation, constrained access to family-capable housing, and the associated efficiency results.
\end{remark}

\clearpage
\tableofcontents
\clearpage

% ==========================================================
\part{Compact Analytical Model}
% ==========================================================

\section{Environment}

The economy is populated by two overlapping generations of households, young and old. Young households are potential parents. They decide whether to enter the modeled housing market and, conditional on entry, choose tenure, housing, consumption, and fertility. Old households are incumbent owners. They are no longer fertile and choose how much owner housing to retain, how much to consume, and how much to leave as bequests.

There is one consumption good and two housing segments. A household can rent housing services in the rental segment, denoted $R$, or occupy housing in the owner segment, denoted $O$. The two segments differ in the size menus they provide. Let
\[
\HH^R=\{h:0<h\le \bar h^R\},
\qquad
\HH^O=\{h:0<h\le \bar h^O\},
\]
with
\[
0<\bar h^R<\bar h^O.
\]
Thus there is rental-owner size segmentation: large family-sized units are more available in the owner segment than in the rental segment.
Housing services are divisible within each segment; $\bar h^R$ and $\bar h^O$ are individual occupancy caps, not aggregate stocks.

A renter pays flow rental cost $q^R h$. An owner pays flow user cost $q^O h$ and faces owner asset price $P^O$. The owner asset price is
\begin{equation}
P^O=\frac{q^O}{\rho_p},
\qquad
\rho_p=r+\delta^O+\tau^p,
\label{eq:compact_owner_price}
\end{equation}
where $r$ is the interest rate, $\delta^O$ is owner-housing depreciation, and $\tau^p$ is a recurring property tax. For a fixed owner user cost $q^O$, a higher $\tau^p$ lowers the asset price and hence lowers the down payment required to buy a given amount of owner housing.

Let $\bar H^R$ and $\bar H^O$ denote inherited rental and owner-segment stocks. New net construction in segment $m\in\{R,O\}$ is $I^m\ge0$, so total segment stocks are
\begin{equation}
H^m=\bar H^m+I^m,
\qquad m\in\{R,O\}.
\label{eq:compact_segment_stocks}
\end{equation}
Construction uses the consumption good. Costs are
\[
\Phi^m(I^m;Z^m),
\qquad
\Phi_I^m(I^m;Z^m)>0,
\qquad
\Phi_{II}^m(I^m;Z^m)>0,
\qquad m\in\{R,O\},
\]
where $Z^m$ shifts the ease of construction in segment $m$. Competitive construction implies, for interior construction,
\begin{equation}
q^R=\Phi_I^R(I^R;Z^R),
\qquad
q^O=\Phi_I^O(I^O;Z^O).
\label{eq:compact_construction_focs}
\end{equation}
Housing payments are transfers across households, old owners, passive landlords, and fiscal accounts. Lump-sum rebates from these accounts are included in disposable resources. The aggregate resource cost of housing expansion is the cost of net construction, not the value of the whole existing stock.

Let
\[
x=1-\tau^b
\]
denote the after-tax bequest share. The rule assigning estate-tax revenue or rebates is part of fiscal policy.

\section{Young households}

There is a mass $\bar M>0$ of potential young households. A potential young household has type
\[
i=(y_i,a_i,B_i),
\]
drawn from $G_Y$. Here $y_i>0$ is disposable nonhousing income, $a_i\ge0$ is liquid wealth, and $B_i\ge0$ is gross inheritance exposure. Effective liquid wealth for purchase constraints is
\[
A_i(x)=a_i+xB_i+T_i^b(x),
\]
where $T_i^b(x)$ is the estate-tax rebate or transfer assigned to type $i$.

Children require consumption and housing services. Let $C_i$ denote total nonhousing consumption expenditure. If a young household has fertility $n_i>0$, the adult residual consumption and housing services are
\[
c_i=C_i-\chi_i n_i>0,
\qquad
s_i=h_i-\kappa n_i>0,
\]
where $\chi_i>0$ is the consumption requirement per child and $\kappa>0$ is the housing requirement per child. Preferences over total nonhousing consumption have the Stone--Geary form
\[
U_i^Y(C_i,h_i,n_i)
=
\log(C_i-\chi_i n_i)+\alpha\log(h_i-\kappa n_i)+\beta\log n_i,
\qquad
\alpha,\beta>0.
\]
Equivalently, using adult residual consumption $c_i=C_i-\chi_i n_i$, the branch problems below use the residual utility
\[
u_i^Y(c_i,h_i,n_i)=
\log c_i+\alpha\log(h_i-\kappa n_i)+\beta\log n_i
\]
and put the child consumption requirement in the budget. In a branch with user cost $q^m$, the equal-scaling benchmark is $\chi_i=\kappa q^m/\alpha$: one child then uses consumption and housing in the same proportions as the adult residual Cobb--Douglas bundle in that branch.

Conditional on entry, household $i$ chooses tenure. The renter branch is
\begin{equation}
W_i^R(q^R)
=
\max_{c_i,h_i,n_i}
\left\{
\log c_i+\alpha\log(h_i-\kappa n_i)+\beta\log n_i
\right\}
\label{eq:compact_renter_problem}
\end{equation}
subject to
\begin{align}
c_i+q^R h_i+\chi_i n_i&=y_i,
\label{eq:compact_renter_budget}\\
h_i&\in\HH^R,
\label{eq:compact_renter_menu}\\
h_i-\kappa n_i&>0.
\notag
\end{align}

The owner branch is
\begin{equation}
W_i^O(q^O,P^O,x)
=
\max_{c_i,h_i,n_i}
\left\{
\log c_i+\alpha\log(h_i-\kappa n_i)+\beta\log n_i
\right\}
\label{eq:compact_owner_problem_young}
\end{equation}
subject to
\begin{align}
c_i+q^O h_i+\chi_i n_i&=y_i,
\label{eq:compact_owner_budget_young}\\
(1-\phi)P^O h_i&\le A_i(x),
\qquad
\phi\in(0,1),
\label{eq:compact_owner_dp}\\
q^O h_i&\le \psi y_i,
\qquad
\psi\in(0,1),
\label{eq:compact_owner_pti}\\
h_i&\in\HH^O,
\notag\\
h_i-\kappa n_i&>0.
\notag
\end{align}
Using $P^O=q^O/\rho_p$, owner housing is bounded above by
\[
h_i\le
\bar h_i^O(q^O,\tau^p,x)
=
\min\left\{
\bar h^O,
\frac{A_i(x)\rho_p}{(1-\phi)q^O},
\frac{\psi y_i}{q^O}
\right\}.
\]

The inside value of entry is
\begin{equation}
W_i^Y(q^R,q^O,P^O,x)
=
\max\{W_i^R(q^R),W_i^O(q^O,P^O,x)\}.
\label{eq:compact_tenure_value}
\end{equation}
Let $o_i\in\{R,O\}$ denote an optimal tenure.

\begin{remark}[Rental premium formulation]
The cap in \eqref{eq:compact_renter_menu} is a compact representation of size segmentation. A smooth formulation lets renters pay
\[
R(h)=q^R h+\varrho(h),
\]
where $\varrho_h(h)\ge0$ for large units. Then the renter marginal conditions below replace $q^R$ by $R_h(h)$.
\end{remark}

\section{Old households}

Old households are indexed by $j$ and distributed according to $F_O$. An old household has income $m_j>0$, initial owner housing $H_j^0$, and a warm-glow bequest utility function $\mathcal B_j$, possibly equal to zero.

Old household $j$ chooses consumption $c_j^O$, retained housing $h_j^O$, and bequest expenditure $e_j\ge0$. Preferences are
\[
u_j^O(c_j^O,h_j^O,e_j)
=
\log c_j^O+\gamma\log h_j^O+\mathcal B_j(e_j),
\qquad
\gamma>0.
\]

The baseline young-access mechanism below does not require any old-side retention advantage. To make old retention concrete when it is present, the compact model uses property-tax assessment lock-in. The market property-tax liability per unit of owner housing is $\tau^pP^O$. Old owner $j$ has protected assessed liability $\tilde\tau_j\tilde P_j$ per unit. Define the per-unit assessment advantage
\[
L_j^\tau(q^O,\tau^p)
=
\tau^pP^O-\tilde\tau_j\tilde P_j
\ge0.
\]
Assume
\[
0\le L_j^\tau(q^O,\tau^p)<q^O.
\]
The old budget is
\begin{equation}
c_j^O+e_j+
\left[q^O-L_j^\tau(q^O,\tau^p)\right]h_j^O
=
m_j+q^O H_j^0,
\qquad
0\le h_j^O\le H_j^0.
\label{eq:compact_old_budget}
\end{equation}
The term $L_j^\tau$ is not a preference. It is a private carrying-cost reduction created by the tax system. The planner counts true old housing utility and true warm-glow bequest utility, but it does not count assessment advantages as utility or as real resource savings.

\section{Entry}

A potential young household can enter the modeled housing market or take an outside option. The outside value is $\bar W^E$. Entry is subject to Type-I extreme-value taste shocks with scale $\kappa_E>0$. The probability that type $i$ enters is
\begin{equation}
\pi_i^E(q^R,q^O,P^O,x)
=
\frac{
\exp\{W_i^Y(q^R,q^O,P^O,x)/\kappa_E\}
}{
\exp\{W_i^Y(q^R,q^O,P^O,x)/\kappa_E\}
+
\exp\{\bar W^E/\kappa_E\}
}.
\label{eq:compact_entry_prob}
\end{equation}
The endogenous measure of young entrants is
\[
\dd F_Y^E(i;q^R,q^O,P^O,x)
=
\bar M\,\pi_i^E(q^R,q^O,P^O,x)\,\dd G_Y(i).
\]
The limit $\kappa_E\downarrow0$ gives deterministic entry. Conditional on entry, the household chooses tenure according to \eqref{eq:compact_tenure_value}.

\section{Competitive equilibrium}

The equilibrium object is aggregate market clearing. Households take prices as given, choose tenure and housing, and their choices aggregate into segment demands. Construction and inherited stocks determine segment supplies. There is no pairwise assignment of old units to young households in the competitive equilibrium.

Given policies $(\tau^p,x)$, construction shifters $(Z^R,Z^O)$, inherited stocks $(\bar H^R,\bar H^O)$, and transfer rules, aggregate rental demand at candidate prices is
\[
H^{R,D}
=
\bar M\int \pi_i^E\1\{o_i=R\}h_i^R\,\dd G_Y(i),
\]
and aggregate owner-segment demand is
\[
H^{O,D}
=
\bar M\int \pi_i^E\1\{o_i=O\}h_i^O\,\dd G_Y(i)
+
\int h_j^O\,\dd F_O(j).
\]
The corresponding segment supplies are
\[
H^{R,S}=\bar H^R+I^R,
\qquad
H^{O,S}=\bar H^O+I^O,
\]
with construction governed by the marginal cost conditions in \eqref{eq:compact_construction_focs}.

\begin{definition}[Competitive equilibrium]
Given $(\tau^p,x,Z^R,Z^O,\bar H^R,\bar H^O)$, a competitive equilibrium consists of user costs $(q^R,q^O)$, owner asset price $P^O$, construction $(I^R,I^O)$, segment stocks $(H^R,H^O)$, young branch allocations, tenure choices, entry probabilities, old allocations, and fiscal transfers such that:
\begin{enumerate}
\item conditional on entry, young households solve \eqref{eq:compact_renter_problem} and \eqref{eq:compact_owner_problem_young}, then choose tenure by \eqref{eq:compact_tenure_value};
\item old households solve \eqref{eq:compact_old_budget};
\item entry is given by \eqref{eq:compact_entry_prob};
\item the owner asset price satisfies \eqref{eq:compact_owner_price};
\item segment stocks satisfy \eqref{eq:compact_segment_stocks} and competitive construction satisfies \eqref{eq:compact_construction_focs};
\item aggregate rental and owner-segment markets clear,
\begin{align}
H^{R,D}&=H^{R,S}=\bar H^R+I^R,
\label{eq:compact_rental_clearing}\\
H^{O,D}&=H^{O,S}=\bar H^O+I^O;
\label{eq:compact_owner_clearing}
\end{align}
\item passive housing-income, property-tax, estate-tax, and policy-instrument accounts are closed by lump-sum transfers.
\end{enumerate}
\end{definition}

Let $\Xi^\Pi\ge0$ denote real resource costs of policy instruments. The transfer rule is chosen so that, after summing household budgets and using \eqref{eq:compact_rental_clearing}--\eqref{eq:compact_owner_clearing}, the aggregate goods constraint is
\begin{align*}
&\bar M\int
\pi_i^E
\left[
\1\{o_i=R\}(c_i^R+\chi_i n_i^R-y_i)
+
\1\{o_i=O\}(c_i^O+\chi_i n_i^O-y_i)
\right]\dd G_Y(i)
\notag\\
&\qquad
+
\int(c_j^O+e_j-m_j)\,\dd F_O(j)
+
\Phi^R(I^R;Z^R)+\Phi^O(I^O;Z^O)+\Xi^\Pi
\le0.
\end{align*}
Housing quantities are disciplined separately by segment feasibility. The goods constraint charges only net construction costs, and the old endowment term $q^O H_j^0$ is a claim on the inherited owner stock rather than new output.

Aggregate fertility is
\[
N
=
\bar M\int \pi_i^E
\left[
\1\{o_i=R\}n_i^R+\1\{o_i=O\}n_i^O
\right]\dd G_Y(i).
\]

\section{Constrained planner}

The planner has the same preferences, type distributions, entry technology, inherited segment stocks, construction technologies, and rental-owner size menus as the competitive economy. The planner chooses entry-relevant inside allocations, tenure assignments, young consumption, housing and fertility, old consumption, housing and bequests, net construction, and policy instruments. The planner respects segment feasibility and the rental and owner menus.

The planner counts true old housing utility and true warm-glow bequest utility. It does not count $L_j^\tau$ as utility or as a resource saving. The planner also cannot relax a down-payment or payment-to-income constraint merely with a notional transfer unless that transfer, payment subsidy, credit guarantee, or contract intervention is in the feasible instrument set and its real cost is included in $\Xi^\Pi$.

Let $\omega_i^Y>0$ and $\omega_j^O>0$ be Pareto weights. Let $\mathcal G_i(n_i)$ denote an optional social payoff from births beyond parental utility. The fertility-neutral planner has $\mathcal G_i\equiv0$.

If the planner assigns young type $i$ tenure $o_i^P\in\{R,O\}$ and allocation $(c_i,h_i,n_i)$, the inside value relevant for entry is
\[
V_i^P=u_i^Y(c_i,h_i,n_i).
\]
The induced entry probability is
\[
\pi_i^P
=
\frac{
\exp\{V_i^P/\kappa_E\}
}{
\exp\{V_i^P/\kappa_E\}
+
\exp\{\bar W^E/\kappa_E\}
}.
\]
The corresponding expected entry value, up to an additive constant, is
\[
\mathcal W_i^E(V_i^P)
=
\kappa_E
\log
\left[
\exp\{V_i^P/\kappa_E\}
+
\exp\{\bar W^E/\kappa_E\}
\right].
\]

The planner solves
\begin{align}
\max\; &
\bar M\int
\omega_i^Y\mathcal W_i^E(V_i^P)\,\dd G_Y(i)
+
\bar M\int
\pi_i^P\mathcal G_i(n_i)\,\dd G_Y(i)
+
\int
\omega_j^O u_j^O(c_j^O,h_j^O,e_j)\,\dd F_O(j)
\label{eq:compact_planner_problem}
\end{align}
subject to segment feasibility,
\begin{align*}
\bar M\int \pi_i^P\1\{o_i^P=R\}h_i\,\dd G_Y(i)&\le \bar H^R+I^R,
\\
\bar M\int \pi_i^P\1\{o_i^P=O\}h_i\,\dd G_Y(i)+\int h_j^O\,\dd F_O(j)&\le \bar H^O+I^O,
\end{align*}
menu feasibility,
\[
h_i\in\HH^R\text{ if }o_i^P=R,
\qquad
h_i\in\HH^O\text{ if }o_i^P=O,
\]
financial implementation, and goods feasibility. If the planner assigns type $i$ to the owner segment and does not use an instrument that relaxes owner financial access for that type, then the assignment must also satisfy
\[
(1-\phi)P^O h_i\le A_i(x),
\qquad
q^O h_i\le \psi y_i.
\]
If a down-payment subsidy, credit guarantee, payment subsidy, or contract intervention relaxes these constraints, the instrument must be feasible and its real resource cost is included in $\Xi^\Pi$. Goods feasibility is
\begin{align*}
&\bar M\int
\pi_i^P
(c_i+\chi_i n_i-y_i)\dd G_Y(i)
\notag\\
&\qquad
+
\int(c_j^O+e_j-m_j)\,\dd F_O(j)
+
\Phi^R(I^R;Z^R)+\Phi^O(I^O;Z^O)+\Xi^\Pi
\le0.
\end{align*}

Let $\Lambda$ be the multiplier on goods feasibility. Let $\Lambda Q^R$ and $\Lambda Q^O$ be the multipliers on rental and owner segment feasibility. Conditional on entry probabilities and tenure assignments, and away from menu corners, an interior planner allocation satisfies
\begin{align*}
\frac{u^Y_{h,i}}{u^Y_{c,i}}&=Q^R
&&\text{for renter assignments},
\\
\frac{u^Y_{h,i}}{u^Y_{c,i}}&=Q^O
&&\text{for young owner assignments},
\\
\frac{u^O_{h,j}}{u^O_{c,j}}&=Q^O
&&\text{for old owner housing},
\\
Q^R&=\Phi_I^R(I^R;Z^R),
\\
Q^O&=\Phi_I^O(I^O;Z^O),
\\
\frac{\beta c_i}{n_i}+\nu_i^g&=\chi_i+\kappa Q^{o_i^P},
\\
\frac{\mathcal B_j'(e_j)}{u^O_{c,j}}&=1
\qquad \text{if }e_j>0.
\end{align*}
The term $\nu_i^g$ is the goods-equivalent social value of a marginal birth beyond parental utility. If $\mathcal G_i\equiv0$, then $\nu_i^g=0$.

\section{Marginal conditions in equilibrium}

The household optimality conditions imply the marginal value of relaxing the relevant housing limits and financial access constraints.

For a renter, let $\lambda_i^R$ be the multiplier on \eqref{eq:compact_renter_budget} and $\theta_i^R\ge0$ the multiplier on the rental size cap $h_i\le\bar h^R$. At an interior allocation with respect to consumption and fertility,
\begin{align}
\frac{1}{c_i^R}&=\lambda_i^R,
\notag\\
\frac{\alpha}{h_i^R-\kappa n_i^R}&=\lambda_i^R q^R+\theta_i^R,
\notag\\
\frac{\beta}{n_i^R}&=\lambda_i^R\chi_i+\frac{\alpha\kappa}{h_i^R-\kappa n_i^R}.
\label{eq:compact_renter_foc_n}
\end{align}
Define the goods-equivalent value of relaxing the rental size limit by
\[
\zeta_i^R=\frac{\theta_i^R}{\lambda_i^R}\ge0.
\]
Then
\begin{equation}
MV_i^R
\equiv
\frac{u^Y_{h,i}}{u^Y_{c,i}}
=
\frac{\alpha c_i^R}{h_i^R-\kappa n_i^R}
=
q^R+\zeta_i^R.
\label{eq:compact_renter_mv}
\end{equation}

For a young owner, let $\lambda_i^O$ be the multiplier on \eqref{eq:compact_owner_budget_young}, $\mu_i\ge0$ the multiplier on the down-payment constraint \eqref{eq:compact_owner_dp}, $\eta_i\ge0$ the multiplier on the payment-to-income constraint \eqref{eq:compact_owner_pti}, and $\theta_i^O\ge0$ the multiplier on the owner menu cap $h_i\le\bar h^O$. At an interior allocation with respect to consumption and fertility,
\begin{align}
\frac{1}{c_i^O}&=\lambda_i^O,
\notag\\
\frac{\alpha}{h_i^O-\kappa n_i^O}
&=
\lambda_i^O q^O+
\mu_i(1-\phi)P^O+
\eta_i q^O+
\theta_i^O,
\notag\\
\frac{\beta}{n_i^O}
&=
\lambda_i^O\chi_i+\frac{\alpha\kappa}{h_i^O-\kappa n_i^O}.
\label{eq:compact_owner_foc_n}
\end{align}
Define the goods-equivalent value of relaxing owner access by
\[
\zeta_i^O
=
\underbrace{\frac{\mu_i}{\lambda_i^O}(1-\phi)P^O+
\frac{\eta_i}{\lambda_i^O}q^O}_{\zeta_i^{O,F}\text{, financial component}}
+
\underbrace{\frac{\theta_i^O}{\lambda_i^O}}_{\text{owner-menu component}}
\ge0.
\]
If the owner menu cap does not bind, then $\theta_i^O=0$ and $\zeta_i^O=\zeta_i^{O,F}$. The young owner's marginal value of owner housing is
\begin{equation}
MV_i^O
\equiv
\frac{u^Y_{h,i}}{u^Y_{c,i}}
=
\frac{\alpha c_i^O}{h_i^O-\kappa n_i^O}
=
q^O+\zeta_i^O.
\label{eq:compact_owner_mv}
\end{equation}

For an old owner at an interior retained-housing choice,
\begin{equation}
MV_j^O
\equiv
\frac{u^O_{h,j}}{u^O_{c,j}}
=
\frac{\gamma c_j^O}{h_j^O}
=
q^O-L_j^\tau(q^O,\tau^p).
\label{eq:compact_old_mv}
\end{equation}
If bequests are interior, then
\[
\mathcal B_j'(e_j)=\frac{1}{c_j^O}.
\]
This is a true preference margin, not an allocation distortion.

\begin{proposition}[Constrained family-housing access and fertility]
Fix an entered young household and its optimal tenure branch. If the household rents, then
\[
\frac{\beta c_i^R}{n_i^R}
=
\chi_i+\kappa(q^R+\zeta_i^R).
\]
If the household owns, then
\[
\frac{\beta c_i^O}{n_i^O}
=
\chi_i+\kappa(q^O+\zeta_i^O).
\]
Thus rental size limits and owner access constraints enter fertility through the housing services required per child.
\end{proposition}

\begin{proof}
For renters, multiply \eqref{eq:compact_renter_foc_n} by $1/\lambda_i^R=c_i^R$ and use \eqref{eq:compact_renter_mv}. For owners, multiply \eqref{eq:compact_owner_foc_n} by $1/\lambda_i^O=c_i^O$ and use \eqref{eq:compact_owner_mv}.
\end{proof}

If no access or menu constraint binds in the relevant tenure branch, the corresponding access value is zero. For example, an unconstrained young owner satisfies
\begin{align*}
c_i^u&=\frac{y_i}{1+\alpha+\beta},
\\
n_i^u&=\frac{\beta y_i}{(1+\alpha+\beta)(\chi_i+\kappa q^O)},
\\
h_i^u&=\frac{\alpha y_i}{(1+\alpha+\beta)q^O}
+
\kappa\frac{\beta y_i}{(1+\alpha+\beta)(\chi_i+\kappa q^O)}.
\end{align*}
The renter analog replaces $q^O$ by $q^R$ and is feasible only if $h_i^u\le\bar h^R$.

\section{Efficiency}

The efficiency question is whether the competitive allocation that clears the aggregate housing markets also solves the planner's allocation problem. Market clearing alone is not enough. In equilibrium, prices make aggregate rental and owner-segment demand equal aggregate supply. Efficiency requires, in addition, that the households who receive the marginal housing units are the households whose marginal values match the social marginal cost of those units.

\begin{proposition}[Planner-equilibrium comparison]
Consider an interior competitive equilibrium and a planner with the same preferences, construction technologies, inherited segment stocks, and rental-owner size menus. Suppose policy accounts are resource-neutral and the planner has no social birth payoff, so $\nu_i^g=0$. Conditional on the equilibrium entry measure and tenure assignments, and away from real menu corners, the competitive allocation satisfies the planner's conditional housing, fertility, and construction first-order conditions only if the financial access gap is zero for entered young owners whose financial constraints are private implementation constraints,
\[
\zeta_i^{O,F}=0,
\]
and the old assessment-lock-in gap is zero,
\[
L_j^\tau(q^O,\tau^p)=0
\]
for old owners at an interior retained-housing choice.
When these conditions fail, the competitive allocation differs from the planner's allocation in the following directions:
\[
MV_i^O=q^O+\zeta_i^{O,F}>q^O=Q^O
\]
for a young owner constrained by down-payment or payment-to-income requirements that the planner can relax with feasible instruments, and
\[
MV_j^O=q^O-L_j^\tau(q^O,\tau^p)<q^O=Q^O
\]
for an old owner with a property-tax assessment advantage. Conversely, if these implementation-relevant gaps are zero for all relevant households, then, conditional on entry, tenure assignments, and real menu constraints, the competitive allocation satisfies the planner's conditional housing, fertility, and construction first-order conditions for suitable Pareto weights.
\end{proposition}

\begin{proof}
Aggregate market clearing and competitive construction imply $Q^R=q^R=\Phi_I^R(I^R;Z^R)$ and $Q^O=q^O=\Phi_I^O(I^O;Z^O)$ at an interior allocation. The planner's conditional housing rules therefore require young renters away from rental-menu corners to satisfy $MV_i^R=q^R$, young owners away from owner-menu corners to satisfy $MV_i^O=q^O$, and old owners to satisfy $MV_j^O=q^O$. The competitive marginal conditions are \eqref{eq:compact_renter_mv}, \eqref{eq:compact_owner_mv}, and \eqref{eq:compact_old_mv}. Comparing these conditions gives the gaps above. The fertility conditions follow from the housing requirement per child: when the relevant access gap is zero and $\nu_i^g=0$, the competitive fertility rules coincide with the planner's rules. Pareto weights can then be chosen to support the competitive consumption allocation. If an implementation-relevant gap is nonzero for a household with positive weight at an interior choice, the competitive allocation violates a necessary planner first-order condition.
\end{proof}

The broader young-access mechanism combines real family-size scarcity in the segment menus with private financial limits on owner access. The implementation-relevant gap in the proposition is the financial component $\zeta_i^{O,F}$: some young households cannot upsize into owner-segment housing services required for children because purchase constraints bind. For them, the marginal value of owner housing exceeds the segment user cost even though the owner segment clears in the aggregate. The old-side mechanism is assessment lock-in: it lowers the old owner's private carrying cost below the social opportunity cost of owner housing. The old-side mechanism can amplify the young-access problem, but it is not required for the young-access mechanism.

Old occupancy is not inefficient by itself. If old housing demand reflects true old housing utility or true warm-glow bequest utility, those motives are part of $u_j^O$ and are counted by the planner. Inefficiency arises when young households face constrained access to family-capable housing, when old retention is supported by a property-tax assessment advantage, or both. The young-access result is present even when the property-tax assessment advantage is zero.

\begin{remark}[Real menus and instruments]
The rental-owner size menus are real constraints of the compact environment. A positive multiplier on a rental or owner menu cap is a shadow value of relaxing that real constraint, not by itself a violation of constrained efficiency when the planner faces the same menus. If the planner can use construction, conversion, assignment, or tenure-access instruments to change family-capable housing access, those instruments and their resource costs must appear in $\Xi^\Pi$. By contrast, down-payment and payment-to-income constraints are private financial constraints. If the planner has feasible instruments that relax them, the financial component $\zeta_i^{O,F}$ is an efficiency-relevant access gap. If the planner must instead respect the same financial constraints, those constraints must appear in the planner problem and their multipliers enter the planner's first-order conditions.
\end{remark}

\section{Fertility and entry}

The fertility implication does not require a planner who directly values births. It follows from the private household problem. Children require consumption and housing services, so a constraint on family-capable housing can change fertility even when fertility enters only parental utility.

It is useful to separate a branch-level comparative static from the equilibrium comparison. Fix a tenure branch and suppress the type index. Let $q$ be the branch user cost and let $H$ be the effective upper bound on occupied housing in that branch. For renters, $H=\bar h^R$. For young owners,
\[
H=
\min\left\{
\bar h^O,
\frac{A_i(x)\rho_p}{(1-\phi)q^O},
\frac{\psi y_i}{q^O}
\right\}.
\]
Holding $(y,\chi,q,\alpha,\beta,\kappa)$ fixed, the branch problem is
\[
\max_{c,h,n}
\left\{
\log c+\alpha\log(h-\kappa n)+\beta\log n
\right\}
\]
subject to
\[
c+qh+\chi n=y,
\qquad
h\le H,
\qquad
h-\kappa n>0.
\]
The unconstrained branch allocation is
\[
c^u=\frac{y}{1+\alpha+\beta},
\qquad
n^u=
\frac{\beta y}{(1+\alpha+\beta)(\chi+\kappa q)},
\]
\[
h^u=
\frac{\alpha y}{(1+\alpha+\beta)q}
+
\kappa
\frac{\beta y}{(1+\alpha+\beta)(\chi+\kappa q)}.
\]
Suppose $0<H<h^u$, so the housing cap binds. Then $h(H)=H$. Let
\[
c(H)=y-qH-\chi n(H),
\qquad
s(H)=H-\kappa n(H).
\]
The constrained fertility choice is the unique solution to
\[
\frac{\beta}{n(H)}
=
\frac{\chi}{c(H)}
+
\frac{\alpha\kappa}{s(H)}.
\]
Implicit differentiation gives
\[
\frac{\dd n}{\dd H}
=
\frac{
\frac{\alpha\kappa}{s(H)^2}
-
\frac{\chi q}{c(H)^2}
}{
\frac{\beta}{n(H)^2}
+
\frac{\chi^2}{c(H)^2}
+
\frac{\alpha\kappa^2}{s(H)^2}
}.
\]
Therefore a relaxation of the binding housing cap raises fertility at $H$ if and only if
\[
\alpha\kappa c(H)^2>\chi q s(H)^2.
\]

This condition has a direct Stone--Geary interpretation. At an unconstrained interior branch allocation,
\[
\frac{q s^u}{c^u}=\alpha.
\]
The child input bundle has the same housing-to-consumption expenditure ratio as the adult residual bundle when
\[
\frac{q\kappa}{\chi}=\alpha,
\qquad\text{equivalently}\qquad
\chi=\frac{\kappa q}{\alpha}.
\]
More generally, fertility is strictly increasing in every binding cap relaxation if and only if
\[
\chi\le \frac{\kappa q}{\alpha},
\]
so the child consumption requirement is no larger than the requirement implied by a scaled copy of the adult residual consumption-housing bundle. A binding cap implies $u_h/u_c=\alpha c(H)/s(H)>q$, so the condition above holds strictly under equal scaling whenever the cap binds. Hence, under equal scaling, every relaxation of a strictly binding family-housing cap raises fertility until the cap ceases to bind. If $\chi>\kappa q/\alpha$, fertility need not rise from every cap relaxation; in particular, $\dd n/\dd H<0$ for binding caps sufficiently close to the unconstrained optimum. The reason is that a larger cap raises housing slack, which lowers the crowding cost of children, but it also induces the household to buy more housing at price $q$, reducing consumption. Fertility rises when the slack effect dominates the consumption-cost effect.

This branch result maps into the model's access margins. A relaxation of the rental size cap raises $H=\bar h^R$. If the owner down-payment constraint is the unique binding component of the effective owner cap, then
\[
\frac{\partial H}{\partial A_i}
=
\frac{\rho_p}{(1-\phi)q^O}>0.
\]
If the payment-to-income constraint is the unique binding component, then
\[
\frac{\partial H}{\partial \psi}
=
\frac{y_i}{q^O}>0.
\]
The derivative above signs the fertility response to these primitive cap relaxations locally, holding tenure, prices, income, and child costs fixed.

The competitive-equilibrium/planner comparison is a different statement because prices, construction, consumption, tenure assignments, and entry can all change. In the competitive equilibrium, an entered young household in branch $m$ satisfies
\[
\frac{\beta c_i^{CE,m}}{n_i^{CE,m}}
=
\chi_i+\kappa(q^{CE,m}+\zeta_i^{CE,m}).
\]
A fertility-neutral planner has no direct social birth payoff. If the planner allocation for the same type and branch is interior with respect to any remaining access constraint, then
\[
\frac{\beta c_i^{P,m}}{n_i^{P,m}}
=
\chi_i+\kappa Q^{P,m}.
\]
Hence
\[
n_i^{P,m}>n_i^{CE,m}
\quad\Longleftrightarrow\quad
\frac{c_i^{P,m}}{c_i^{CE,m}}
>
\frac{\chi_i+\kappa Q^{P,m}}
{\chi_i+\kappa(q^{CE,m}+\zeta_i^{CE,m})}.
\]
Thus a fertility-neutral planner does not mechanically deliver higher fertility merely because the competitive allocation has a positive access value. The branch-level result gives a primitive sufficient condition for a signed fertility response when a policy or planner locally raises a type's effective housing cap holding the branch primitives fixed.

If the planner directly values births through $\mathcal G_i(n_i)$, the interior planner condition becomes
\[
\frac{\beta c_i^P}{n_i^P}+\nu_i^g
=
\chi_i+\kappa Q^{P,m}.
\]
A positive $\nu_i^g$ is a separate fertility motive. Holding consumption, segment scarcity costs, and any remaining access constraints fixed, it raises desired fertility. In the full equilibrium, aggregate fertility also depends on prices, consumption, tenure assignments, and entry.

Aggregate fertility is
\[
N
=
\bar M\int \pi_i^E
\left[
\1\{o_i=R\}n_i^R+\1\{o_i=O\}n_i^O
\right]\dd G_Y(i).
\]
Across two allocations,
\[
N^P-N^{CE}
=
\bar M\int
\left[
n_i^P(\pi_i^P-\pi_i^E)
+
\pi_i^E(n_i^P-n_i^{CE})
\right]\dd G_Y(i),
\]
where \(n_i^{CE}=\1\{o_i=R\}n_i^R+\1\{o_i=O\}n_i^O\). The second term is the intensive fertility margin among the households who would enter in the competitive equilibrium. The first term is the entry margin induced by the change in inside value. If entry is held fixed, the aggregate fertility effect is just
\[
N^P-N^{CE}
=
\bar M\int
\pi_i^E(n_i^P-n_i^{CE})\,\dd G_Y(i).
\]
This is the aggregate analog of the planner-equilibrium comparison: the planned allocation delivers higher fertility exactly to the extent that it raises conditional fertility for constrained young types and, possibly, increases entry into the modeled housing market.

\subsection{A local policy statistic}

The branch calculation also gives a simple local statistic for the fertility effect of a policy. Let $z$ be a policy parameter and let
\[
H_i^m(z)
\]
denote the effective family-housing cap faced by type $i$ in branch $m$. For example, $H_i^R=\bar h^R$ for renters, while for owners
\[
H_i^O
=
\min\left\{
\bar h^O,
\frac{A_i(x)\rho_p}{(1-\phi)q^O},
\frac{\psi y_i}{q^O}
\right\}.
\]
Let
\[
P_i^m(z)=\frac{\partial H_i^m}{\partial z}
\]
be the policy pass-through to the effective family-housing cap. If a policy does not relax the component of $H_i^m$ that binds for type $i$, then $P_i^m(z)=0$ for this mechanism.

For a type whose cap binds in branch $m$, define the fixed-branch fertility response
\[
\Psi_i^m
=
\frac{
\frac{\alpha\kappa}{(s_i^m)^2}
-
\frac{\chi_i q^m}{(c_i^m)^2}
}{
\frac{\beta}{(n_i^m)^2}
+
\frac{\chi_i^2}{(c_i^m)^2}
+
\frac{\alpha\kappa^2}{(s_i^m)^2}
}.
\]
Under equal scaling, $\chi_i=\kappa q^m/\alpha$, and a strictly binding cap, $\Psi_i^m>0$.

Holding prices, tenure branches, and active sets fixed, the local aggregate fertility effect among entered households is
\[
\frac{\dd N^{cond}}{\dd z}
=
\bar M
\int_{\mathcal B(z)}
\pi_i^E
\Psi_i^m
P_i^m(z)
\dd G_Y(i),
\]
where $\mathcal B(z)$ is the set of entered types for whom the effective family-housing cap binds in the relevant branch. This expression says that conditional fertility responds to the mass of constrained households, the pass-through of the policy to their effective cap, and the fertility response to that cap.

If entry is allowed to respond, the logit entry block adds an entry term. Since
\[
\frac{\partial \pi_i^E}{\partial W_i^Y}
=
\frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
\]
and the utility value of relaxing the effective cap is $\lambda_i^m\zeta_i^m$, the local fixed-price, fixed-branch statistic becomes
\[
\frac{\dd N}{\dd z}
=
\bar M
\int_{\mathcal B(z)}
\left[
\pi_i^E\Psi_i^m
+
n_i^m
\frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
\lambda_i^m\zeta_i^m
\right]
P_i^m(z)
\dd G_Y(i).
\]
Thus policy effectiveness is exposure to the binding family-space constraint, times pass-through to the effective cap, times the fertility response, plus the induced entry response. The statistic is local: a full general-equilibrium policy experiment also changes prices, construction, tenure assignments, old retention, fiscal transfers, and possibly the resource cost $\Xi^\Pi$.

\clearpage

\begin{thebibliography}{99}
\small

\bibitem[Barro and Becker(1988)]{barro1988}
Barro, Robert J., and Gary S. Becker. 1988. ``A Reformulation of the Economic Theory of Fertility.'' \emph{Quarterly Journal of Economics} 103(1): 1--25.

\bibitem[Barro and Becker(1989)]{barro1989}
Barro, Robert J., and Gary S. Becker. 1989. ``Fertility Choice in a Model of Economic Growth.'' \emph{Econometrica} 57(2): 481--501.

\bibitem[Coven et al.(2024)]{coven2024}
Coven, Joshua, Sebastian Golder, Arpit Gupta, and Abdoulaye Ndiaye. 2024. ``Property Taxes and Housing Allocation Under Financial Constraints.'' Working paper / CESifo Working Paper No. 11203.

\bibitem[Diamond(1965)]{diamond1965}
Diamond, Peter A. 1965. ``National Debt in a Neoclassical Growth Model.'' \emph{American Economic Review} 55(5): 1126--1150.

\bibitem[Fonseca et al.(2024)]{fonseca2024}
Fonseca, Julia, Lu Liu, and Pierre Mabille. 2024. ``Unlocking Mortgage Lock-In: Evidence from a Spatial Housing Ladder Model.'' Working paper.

\bibitem[Galor and Weil(1996)]{galorweil1996}
Galor, Oded, and David N. Weil. 1996. ``The Gender Gap, Fertility, and Growth.'' \emph{American Economic Review} 86(3): 374--387.

\bibitem[Gervais(2002)]{gervais2002}
Gervais, Martin. 2002. ``Housing Taxation and Capital Accumulation.'' \emph{Journal of Monetary Economics} 49(7): 1461--1489.

\bibitem[Sommer and Sullivan(2018)]{sommersullivan2018}
Sommer, Kamila, and Paul Sullivan. 2018. ``Implications of US Tax Policy for House Prices, Rents, and Homeownership.'' \emph{American Economic Review} 108(2): 241--274.

\end{thebibliography}

\end{document}

```

---

# B. Working Audit And Sufficient-Statistic Memo From Codex

```markdown
# Part I Deep Audit And Policy Sufficient Statistic

Date: 2026-06-10

Scope: standalone Part I analytical model in `latex/intergenerational_housing_fertility_part1.tex`. This note is a working audit, not a replacement draft.

## 1. Deep Audit

The current Part I is on the right path. It now starts from an economy: young potential parents, old incumbent owners, rental and owner housing segments, construction, aggregate market clearing, a planner, and then marginal comparisons. That is the right macro-theory order. It no longer reads primarily as a list of wedges.

The strongest part is the young household block. Tenure is now primitive, and the child-input structure is economically coherent: children require adult residual consumption and adult residual housing services through
\[
c_i=C_i-\chi_i n_i,\qquad s_i=h_i-\kappa n_i.
\]
The branch representation with residual consumption \(c_i\) is the right way to write the FOCs. The text must keep reminding the reader that \(c_i\) is adult residual consumption, not total nonhousing consumption.

The Section 9 derivative is now doing useful work. Under equal scaling,
\[
\chi_i=\frac{\kappa q^m}{\alpha},
\]
one child uses consumption and housing in the same expenditure proportions as the adult residual Cobb-Douglas bundle. Then any relaxation of a strictly binding effective family-housing cap raises fertility until the cap stops binding. This is the cleanest fertility statement in the current draft.

The main remaining economic issue is the distinction between three objects that are still too easy to blur:

1. Real size menus: \(h\le \bar h^R\) and \(h\le \bar h^O\). These are real constraints unless the planner has construction, conversion, or assignment instruments that change the menu.
2. Financial implementation constraints: down-payment and payment-to-income constraints. These can create a CE/planner gap if the planner has feasible instruments that relax them at resource cost.
3. Aggregate stock scarcity: \(H^m=\bar H^m+I^m\), cleared by prices and construction.

The current text now handles this better, but the paper should keep these objects separate in every efficiency and policy statement. A binding rental size cap is not by itself an inefficiency if the planner faces the same cap. It is a shadow value of relaxing a real constraint. A binding down-payment constraint is a private implementation gap only if the planner can relax it through an explicit instrument.

There is also a subtle construction/menu issue. Current construction \(I^m\) expands aggregate segment services and affects \(q^m\). It does not, by itself, raise the individual size cap \(\bar h^m\). If the policy is "more family-sized rental units" or "conversion of owner units into family-capable rentals," the model needs an instrument that changes the effective cap or menu, not merely an increase in aggregate stock. Otherwise the policy changes prices but not the family-size access constraint that drives the clean Section 9 result.

The planner section is improved but still should be sold cautiously. It is a conditional marginal comparison, not a full welfare theorem. It is conditional on entry, tenure assignments, interior construction, and being away from real menu corners. That is acceptable, but the prose should keep saying "conditional marginal comparison" rather than "the CE is inefficient" without qualifications.

The old-side property-tax mechanism is conceptually clean. \(L_j^\tau\) is not old utility and not a real resource saving. It lowers the old owner's private carrying cost below the social opportunity cost. The old-side channel can amplify young access problems, but the young-access/fertility mechanism does not depend on property tax.

Notation cleanup still matters before this becomes paper-facing. Superscript \(O\) means young owner tenure and old household variables. That is readable locally but fragile. A later pass should use \(B\), \(\mathrm{own}\), or \(b\) for young buyers and reserve "old" notation for old households. Also, \(A_i(x)\) should be described as pledgeable liquid resources for purchase constraints, not full wealth, unless it is put into the lifetime budget.

## 2. Sufficient Statistic For Policy Effectiveness

The simplest sufficient-statistic object is local and branch-specific. It should not pretend to solve the full equilibrium. Fix a policy parameter \(z\). Fix a type \(i\), tenure branch \(m\), and active set. Let the policy change the household's effective family-housing cap
\[
H_i^m(z).
\]
For renters, \(H_i^R=\bar h^R\). For owners,
\[
H_i^O=
\min\left\{
\bar h^O,
\frac{A_i(x)\rho_p}{(1-\phi)q^O},
\frac{\psi y_i}{q^O}
\right\}.
\]

For a type whose cap binds, define the branch fertility response to a cap relaxation:
\[
\Psi_i^m
=
\frac{
\frac{\alpha\kappa}{s_i^2}
-
\frac{\chi_i q^m}{c_i^2}
}{
\frac{\beta}{n_i^2}
+
\frac{\chi_i^2}{c_i^2}
+
\frac{\alpha\kappa^2}{s_i^2}
}.
\]
This is \(\partial n_i^m/\partial H_i^m\) holding the branch, prices, income, child costs, and preferences fixed. Under equal scaling and a strictly binding cap, \(\Psi_i^m>0\).

Define policy pass-through to the effective cap:
\[
P_i^m(z)=\frac{\partial H_i^m}{\partial z}.
\]
Examples:

- Rental cap relaxation: if \(z=\bar h^R\), then \(P_i^R=1\) for capped renters.
- Down-payment liquidity grant: if the down-payment component uniquely binds and \(z=A_i\), then
\[
P_i^O=\frac{\rho_p}{(1-\phi)q^O}.
\]
- Payment-to-income relaxation: if \(z=\psi\), then
\[
P_i^O=\frac{y_i}{q^O}.
\]
- A policy that makes the wrong units cheaper has \(P_i^m=0\) for fertility if it does not relax the binding family-space cap.

With fixed entry, tenure, and prices, the local aggregate fertility effect is
\[
\frac{\dd N^{cond}}{\dd z}
=
\bar M
\int_{\mathcal B(z)}
\pi_i^E\,
\Psi_i^m\,
P_i^m(z)
\,
\dd G_Y(i),
\]
where \(\mathcal B(z)\) is the set of entered types in the active branch whose effective family-housing cap binds.

This is the clean fertility sufficient statistic:
\[
\text{fertility effect}
=
\text{mass exposed}
\times
\text{policy pass-through to effective family space}
\times
\text{fertility response to family space}.
\]

If entry is allowed to respond but tenure and prices are still held fixed, the logit entry block adds an entry term. Since
\[
\frac{\partial \pi_i^E}{\partial W_i^Y}
=
\frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
\]
and the utility value of relaxing the cap is \(\lambda_i^m\zeta_i^m\), the local effect becomes
\[
\frac{\dd N}{\dd z}
=
\bar M
\int_{\mathcal B(z)}
\left[
\pi_i^E\Psi_i^m
+
n_i^m
\frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
\lambda_i^m\zeta_i^m
\right]
P_i^m(z)
\,
\dd G_Y(i).
\]
The first term is the intensive fertility effect among entrants. The second is the entry effect: policies that relax binding family-space constraints raise inside value and can draw in additional households.

This statistic is sufficient for a local, fixed-price, fixed-branch policy experiment. It is not sufficient for the full general equilibrium response. A full GE policy effect also needs price pass-through, construction responses, tenure switching, changes in old retention, and the resource cost of instruments.

For welfare, the analogous local statistic is even simpler:
\[
\frac{1}{\Lambda}\frac{\dd W}{\dd z}
=
\bar M
\int
\pi_i^E
\zeta_i^{relaxed}
P_i^m(z)
\dd G_Y(i)
-
\frac{\dd \Xi^\Pi}{\dd z}
\]
plus any old-side term if the policy changes old retained housing:
\[
\int L_j^\tau
\left(-\frac{\partial h_j^O}{\partial z}\right)
\dd F_O(j).
\]
Thus the birth statistic uses \(\Psi_i^m\), while the welfare statistic uses \(\zeta_i\). They are related but not the same object.

## 3. Simple Paper Statement

The clean statement is:

> A policy affects fertility only to the extent that it relaxes the effective family-housing cap of households for whom that cap binds. Under equal scaling of child consumption and housing requirements, the fertility effect of a local cap relaxation is positive for every strictly constrained household. Aggregate effectiveness is the mass of constrained households times the policy pass-through to their effective cap times their branch fertility response, plus any entry response.

That sentence avoids two common overclaims. It does not say every housing policy raises fertility. It says the policy must relax the binding family-space constraint. It also does not say the planner always raises fertility in GE. It says the local fertility response is signed when the policy relaxes the effective cap under the stated primitive condition.

## 4. Two Figures To Use

Figure 1 should show the branch mechanism. Put the effective cap \(H\) on the horizontal axis and fertility \(n(H)\) on the vertical axis. Under equal scaling, the curve rises over the binding region \(H<h^u\) and flattens once the cap no longer binds. Label the shaded binding region "family-space constrained" and annotate the policy arrow \(\dd H>0\), \(\dd n>0\). This picture should be the visual anchor for Section 9.

Figure 2 should show the sufficient-statistic decomposition. The picture should not be a supply-demand graph. It should be a pipeline:
\[
\text{policy}
\to
\text{pass-through to } H_i^m
\to
\text{fertility response } \Psi_i^m
\to
\text{aggregate births}
\]
with a separate entry branch. The bottom of the figure should display the formula
\[
\dd N
\approx
\bar M\int
\left[
\pi_i^E\Psi_i^m
+
n_i^m
\frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
\lambda_i^m\zeta_i^m
\right]
\dd H_i^m
\dd G_Y(i).
\]
The visual point is that the policy works only through \(\dd H_i^m\), the relaxation of the effective family-space cap. If \(\dd H_i^m=0\), the policy may affect ownership labels or prices but not this fertility mechanism.


```

---

# C. Figure Sheet Source

The compiled PDF is `latex/part1_policy_sufficient_stat_figures_20260610.pdf`. The source below defines two conceptual figures: the branch cap/fertility mechanism and the local policy-statistic pipeline.

```latex
\documentclass[11pt]{article}

\usepackage[landscape,margin=0.6in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{tikz}
\usepackage{xcolor}
\usepackage{microtype}

\usetikzlibrary{arrows.meta,positioning,calc,patterns}

\newcommand{\dd}{\mathrm d}
\newcommand{\1}{\mathbf 1}

\definecolor{ink}{RGB}{35,35,35}
\definecolor{accent}{RGB}{38,98,135}
\definecolor{muted}{RGB}{105,105,105}
\definecolor{softblue}{RGB}{231,242,248}
\definecolor{softgreen}{RGB}{232,244,236}
\definecolor{softred}{RGB}{249,235,232}

\pagestyle{empty}

\tikzset{
  axis/.style={->,thick,draw=ink},
  curve/.style={very thick,draw=accent},
  arrow/.style={-{Latex[length=3mm]},thick,draw=accent},
  redarrow/.style={-{Latex[length=3mm]},thick,draw=red!65!black},
  note/.style={font=\small,align=center,text=ink},
  smallnote/.style={font=\scriptsize,align=center,text=muted},
  box/.style={draw=ink,rounded corners=2pt,thick,align=center,inner sep=7pt,minimum height=1.0cm,fill=white},
  bluebox/.style={box,fill=softblue},
  greenbox/.style={box,fill=softgreen},
  redbox/.style={box,fill=softred}
}

\begin{document}

\begin{center}
{\Large Figure 1. Binding Family-Space Caps And Fertility Under Equal Scaling}
\end{center}

\vspace{0.3cm}

\begin{tikzpicture}[x=1cm,y=1cm]
  \fill[softblue] (0.7,0.45) rectangle (7.8,5.25);
  \node[note,text=accent] at (4.2,4.9) {family-space constrained region};

  \draw[axis] (0.7,0.55) -- (10.8,0.55);
  \node[note] at (5.7,-0.2) {effective housing cap $H$};
  \draw[axis] (0.7,0.55) -- (0.7,5.55) node[above left,note] {fertility $n(H)$};

  \draw[curve]
    plot[smooth,tension=0.72] coordinates {
      (0.95,0.82) (1.55,1.25) (2.2,1.75) (3.0,2.35)
      (3.8,2.92) (4.75,3.42) (5.75,3.82) (6.8,4.05)
      (7.8,4.16) (9.5,4.16)
    };

  \draw[dashed,thick,draw=muted] (7.8,0.55) -- (7.8,4.16);
  \node[smallnote,below] at (7.8,0.48) {$h^u$};
  \node[smallnote,above right] at (7.8,4.2) {cap no longer binds};

  \draw[dashed,draw=muted] (2.0,0.55) -- (2.0,1.62);
  \node[smallnote,below] at (2.0,0.48) {$\bar h^R$};
  \draw[dashed,draw=muted] (3.45,0.55) -- (3.45,2.7);
  \node[smallnote,below] at (3.45,0.48) {$\widehat h_i^O$};

  \draw[redarrow] (2.0,1.72) -- (3.45,2.7)
    node[midway,above,sloped,note,text=red!65!black] {$\dd H>0,\ \dd n>0$};

  \node[bluebox,text width=5.0cm] at (14.0,2.05) {
    Equal scaling:
    \[
      \chi=\frac{\kappa q}{\alpha}
    \]
    child inputs have the same expenditure mix as the adult residual bundle.
  };

  \node[greenbox,text width=5.0cm] at (14.0,4.35) {
    While the cap binds,
    \[
      \frac{u_h}{u_c}>q
      \quad\Rightarrow\quad
      \frac{\dd n}{\dd H}>0.
    \]
  };

  \node[smallnote,text width=12.0cm] at (6.1,-1.35) {
    The curve is conceptual. Section 9 gives the exact local slope:
    \[
    \frac{\partial n}{\partial H}
    =
    \frac{\alpha\kappa/s^2-\chi q/c^2}
    {\beta/n^2+\chi^2/c^2+\alpha\kappa^2/s^2}.
    \]
  };
\end{tikzpicture}

\clearpage

\begin{center}
{\Large Figure 2. A Local Sufficient Statistic For Policy Effectiveness}
\end{center}

\vspace{0.2cm}

\begin{tikzpicture}[x=1cm,y=1cm]
  \node[bluebox,text width=2.7cm] (policy) at (0,1.6) {
    Policy $z$\\
    \scriptsize grant, PTI rule, cap, conversion
  };
  \node[bluebox,text width=3.0cm] (pass) at (3.9,1.6) {
    Pass-through\\
    \[
      P_i^m=\frac{\partial H_i^m}{\partial z}
    \]
  };
  \node[greenbox,text width=3.0cm] (cap) at (7.8,1.6) {
    Effective cap\\
    \[
      \dd H_i^m=P_i^m\,\dd z
    \]
  };
  \node[greenbox,text width=3.1cm] (intensive) at (11.9,3.1) {
    Intensive fertility\\
    \[
      \Psi_i^m=\frac{\partial n_i^m}{\partial H_i^m}
    \]
  };
  \node[greenbox,text width=3.1cm] (entry) at (11.9,0.1) {
    Entry margin\\
    \[
      \frac{\partial\pi_i}{\partial W_i}
      =
      \frac{\pi_i(1-\pi_i)}{\kappa_E}
    \]
  };
  \node[redbox,text width=2.8cm] (agg) at (15.8,1.6) {
    Aggregate births\\
    \[
      \dd N
    \]
  };

  \draw[arrow] (policy) -- (pass);
  \draw[arrow] (pass) -- (cap);
  \draw[arrow] (cap.east) -- (intensive.west);
  \draw[arrow] (cap.east) -- (entry.west);
  \draw[arrow] (intensive.east) -- (agg.west);
  \draw[arrow] (entry.east) -- (agg.west);

  \node[box,text width=16.0cm,fill=white] at (7.9,-3.65) {
    Local fixed-price, fixed-branch statistic:
    \[
    \dd N
    \approx
    \bar M
    \int
    \left[
    \pi_i^E\Psi_i^m
    +
    n_i^m
    \frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
    \lambda_i^m\zeta_i^m
    \right]
    \dd H_i^m
    \,\dd G_Y(i).
    \]
    Effectiveness is exposure $\times$ pass-through to effective family space $\times$ fertility response, plus entry.
  };
\end{tikzpicture}

\end{document}

```
