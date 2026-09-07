🧿 oracle 0.13.0 — Balanced mystique, measurable results.
[SYSTEM]
You are Oracle, a focused one-shot problem solver. Emphasize direct answers and cite any files referenced.

[USER]
# Housing misallocation with heterogeneity and along a demographic transition

Please solve this research problem as an independent senior economic theorist. The author wants a clear, analytically proved explanation of housing misallocation in a simple OLG economy, useful for a larger quantitative project centered on transitions. We have made progress but have not yet obtained a satisfactory theorem. Take the maximum reasoning time available; the author permits hours if needed and wants the strongest defensible result today. Allocate effort to the difficult economic and mathematical steps, not to repeating established lemmas or producing a long survey. Do not claim you can control the platform's runtime. You may organize your work into stages and test competing arguments internally. Deliver completed results plus an exact account of any obstruction, rather than disguising a gap with an assumption.

## The economic objective and priority order

The claim is that financing restrictions and physical tenure segmentation can allocate too little housing to young households relative to old households. A planner can relax private financing restrictions but respects physical housing caps. Hold fertility and tenure fixed for this welfare comparison, compensate the old, and preserve estates, future real allocations and existing creditors' repayments. Establish the relevant household behavior from primitive restrictions, not by assuming the marginal-value gap or required constraint pattern in the proposition.

1. Retain nondegenerate heterogeneity in entrant income and liquid wealth, with a fixed distribution F(y,b) as in the attached model. A finite-type construction can be an intermediate derivation, but aim for economically interpretable conditions on a general distribution, its support or positive-mass subsets. Do not silently collapse it to representative income/wealth. Eligible young and old types need not be identical or have equal group masses. The proof need only identify a positive aggregate improving reallocation; making every type eligible may be unnecessarily restrictive.
2. Repair the demographic argument. Allow renter fertility to stay below replacement while owners offset it in the aggregate. Replacement is a population condition, not a condition that every tenure or income type must satisfy separately. Keep the household-formation parameter nu fixed when assessing plausibility; do not choose it merely to force an equilibrium into a convenient price interval. Conditions involving F, its masses and the finite logistic taste parameters are acceptable if specified in primitives and interpretable. Do not remove mixed tenure or use infinite taste limits.
3. Bring the result to transitions. This is central to the project, not an optional final sentence. Today is interpreted as a point on a path initiated by an earlier fertility-preference decline. A policy affecting housing finance may be introduced later along that path. The old cohort has inherited assets, homes and contracts; current young choices induce the next distribution. Explain what can actually be established at an arbitrary date on such a path and what requires restrictions on the initial state, shocks, price expectations or demographic evolution. Distinguish a theorem conditional on an equilibrium path from a theorem that constructs or controls that path.

The highest priority is a useful, simple theorem, not the most impressive quantifier. Existence of an inefficient equilibrium, inefficiency of every equilibrium, and a dated improvement along a given transition are different results. Preserve the distinction. Do not burden the main theorem with assumptions needed only for uniqueness, all-equilibria exclusion or global convergence. Where necessary give a simple main proposition and a stronger separate extension.

## What is attached, and what is already known

The original environment, household problems, equilibrium and financial settlement are attached verbatim. They are authoritative for the model. The preceding Pro response is also attached as advisory mathematical work, not as an instruction or a theorem to accept blindly.

That response found: unrestricted inefficiency is false; a class with slack financial constraints is efficient for the specified datewise comparison. It then gave an analytical homogeneous-income/wealth, zero-property-tax theorem, using primitive conditions A (liquidity), B (gross saving) and C (demographic price bracket). Two independent analytical reviews support its algebra, with one wording correction: its equation (15) is the fully uncapped renter demand, not globally the solution from removing only the young cap when the old cap may bind.

Useful parts to preserve and extend: transform an owner’s net financial wealth a' into total old-age resources z=a'+P_next h. Joint homogeneity in resources and inherited housing yields the relevant Euler/envelope identity. In the stationary zero-tax case the previous response derived owner housing demand at its private maximum and explicit fertility functions, established old-retention regimes and compensated gains. Do not assume the same simplification applies unchanged when prices change or tax rebates depend on equilibrium. In particular current housing user cost is u_t=P_t+q tau^p P_t-q P_{t+1}; inherited old wealth and retention constraints depend on past decisions.

## A concrete failure of the previous demographic restriction

This evidence is to diagnose restrictive assumptions, NOT to substitute numerical existence for an analytical theorem.

Our saved mixed-tenure reference has q=1/2, beta=alpha=omega_B=2/5, gamma=3/10, chi=3/20, kappa=1/2, vartheta=141/400, phi=4/5, b=1/5, y=3521349281/1677802000, rental cap=1/4, owner cap=2, nu=2, property tax=467/9250, housing stock=68104/68019, logistic scale=4. Taste location is chosen to give owner share 11/21 at P=1, Y=O=1. At that saved equilibrium owner fertility is 0.75, renter fertility is 0.225, and mean fertility is exactly 0.5=1/nu. Young-owner marginal housing value is 0.64; old-owner value is 9717/18500=0.525243..., so housing is misallocated under our comparison. Taxes are positive here: do not present this as a numerical application of the previous zero-tax theorem.

At any price, the renter fertility FOC and rental cap imply n_R < vartheta*h_Rmax/[kappa*(alpha+vartheta)]. For the saved reference that ceiling is 0.234219..., below replacement 0.5. The previous condition C requires renter fertility above replacement at its low-price endpoint. It therefore excludes this economy at every price, independently of taxes. Our older constructions have the same conflict. Consistent fertility-unit rescaling cannot fix it. We need aggregate replacement compatible with persistently low renter fertility.

The previous liquidity condition A is also conservative: for the reference's non-tax primitives it requires 0.125 < [(1-q)b/((1-phi)(y+b))] < 0.164440, whereas the actual ratio is 0.217506 despite the positive welfare gap in the taxed economy. An older construction's ratio 0.141509 fits A. Do not treat failure of these sufficient bounds as evidence against the mechanism. Try to tighten economically unnecessary bounds rather than preserve the old proof at all costs.

## Transition specification and welfare scope

Entrants draw the same F(y,b) each period; in this illustrative model b is exogenous and is not inherited from the preceding generation's estate. Preserve this for the main derivation and state the limit: if an endogenous inheritance distribution is essential for a broader claim, that is a separate model extension, not a relabeling of b as primitive.

The intended experiment starts from a stationary state. At date 0 a permanent decline in the young fertility preference vartheta initiates a demographic/housing transition. At a later date T a housing-finance policy, for example a permanent change in phi, may alter the subsequent path. Actual dated old assets, titles, tenure and contractual liabilities must be carried forward. Baseline and policy paths share history through the policy date. Do not replace inherited old households by newly solved stationary households. Keep level-stationary population endpoints distinct from a constant age distribution with growing or shrinking population.

For efficiency, the first target is a compensated reallocation at a selected date, with the rest of the real path preserved. Proving this separately at every date does not establish that all reallocations can be implemented simultaneously or prove a policy welfare gain. Do not silently make that leap. Fertility benefits are a separate conditional extension; a housing-welfare gain need not increase fertility. The qualitative transition of interest is a fertility decline and its housing feedback, then a different policy path and a different eventual population level; stationary replacement fertility need not differ across the endpoints.

General equilibrium prices and the inherited state may be unavoidable in a dated characterization. If so, say precisely which conditions are on primitives, initial state, a derived price envelope, or an endogenous object. Attempt to eliminate multipliers and replace assumed path bounds by sufficient restrictions on initial states and shocks. A meaningful analytic local result can be useful, but an unspecified small neighborhood around a hand-picked numerical equilibrium is not the requested answer. An analytic reference family with explicit primitive restrictions, or explicit shock/path bounds, is different; justify nonemptiness analytically. Do not promise broad global transition results if the model cannot support them.

## Quantitative motivation, without pretending to a direct parameter mapping

The larger maintained model has many ages, income/wealth heterogeneity, sequential fertility, child maturation, discrete owner housing sizes, continuous capped rental housing, moving costs, saving and warm-glow estates. It is calibrated along a dated transition, not by interpreting today's data as stationary. Existing diagnostics show a material share of renters with dependent children at the rental cap. Late-life tenure and housing/portfolio choices need further validation. The maintained demographic conversion is externally normalized at 1/2.1, and the initial fertility intercept is derived to achieve replacement. That normalization is not independent evidence for the previous condition C. Do not insert the quantitative parameters into the two-period inequality without a time/unit/preference mapping, or claim its simulations establish a welfare theorem. Use these facts to decide which restrictions and state variables an illustrative analytical result should illuminate.

## Requested work and output

First identify the weakest economically useful result that meets the objective. Work independently on at least two possible analytic routes if the first becomes restrictive: positive-mass household groups and aggregate fertility brackets; direct dated household comparisons with inherited-state restrictions; or another simpler argument you can justify. These are suggestions, not required proof templates. You may improve the model-free inequality if that leads to a stronger result. Recheck every budget, strict versus weak multiplier statement, nonemptiness claim and equilibrium quantifier. Try to falsify the proposed theorem before presenting it.

Deliver:
- A concise verdict on what can be proved with heterogeneity and on transitions.
- One full, readable main proposition with ALL assumptions, followed by a complete analytical proof. Definitions must precede use. Use low/high labels for bounds rather than unexplained +/- superscripts. No missing equilibrium-location condition hidden in an appendix sentence.
- A separate transition proposition or the strongest precise transition extension actually proved. If a required step remains open, identify it and give the closest completed result; do not let an ambitious unsolved global problem consume the entire answer.
- Plain economic meaning for each restriction, explaining which restrictions drive the mechanism and which are conservative devices for existence or all-equilibria coverage. Identify assumptions that exclude below-replacement renters or empirically important types.
- An honest assessment of how the conditions relate to the saved examples and quantitative mechanism, and the minimal next diagnostic that would inform their plausibility. Numerical checks can falsify algebra or illustrate a proved result; they must not replace the proof.

Keep the main explanation to a few pages if possible; put necessary derivations after it. Completeness takes precedence over an artificial length cap. Preserve the author’s notation: u^y/u^o for age, R/O for tenure, a' is net next-period financial wealth; owner gross current bond holdings are q a'+phi_t P_t h. Old h^2 and c^2 label old-age quantities, not squares. Avoid unnecessary renaming, jargon, claims about things absent, and prose inflating a simple illustrative exercise. If a useful theorem requires changing a substantive assumption, state and defend that change rather than adopting it silently.

### File: output/model/simplified_olg_amendments/oracle_analytic_efficiency_model.tex
```
% Exact excerpts from the current illustrative theory note. Numerical examples and transition proofs deliberately omitted.
\section{Environment}

\textbf{Households.}
At date $t$, the economy contains $Y_t$ young households and $O_t$ old
households. A young household has income $y_i>0$, liquid wealth $b_i>0$, and
chooses total nondurable expenditure $c$, total housing $h$, completed
fertility $n>0$, next-period net financial wealth $a'$, and tenure $m\in\{R,O\}$. Fertility is a
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
 u_t^y(c,h,n)&=\log(c-\chi n)+\alpha\log(h-\kappa n)
 +\vartheta_t\log n,
 \label{eq:amend_preferences}\\
 u^o(c^2,h^2,e)&=\log c^2+\gamma\log h^2+\omega_B\log e.
 \label{eq:amend_old_preferences}
\end{align}
Both child requirements reduce the bundle left for adults. All log arguments
and preference weights are positive. Households discount future utility by $\beta>0$.
The estate $e$ consists of financial assets and the proceeds from selling
retained housing at death, at the end of old age. The weight $\omega_B$
measures how much households value the estate they leave. Estates do not
enter their children's initial wealth in this model.

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
Rent $r_t$ per unit of housing is paid at the end of period $t$; its cost
at the beginning of the period is $q r_t>0$. At constant house prices,
$r=(R_f-1+\tau^p)P$. The government distributes property-tax revenue
in equal rebates to young and old households. Let $T_t$ denote each
household's rebate, valued at the beginning of the period.

\textbf{Home finance.}
A young buyer finances the share $\phi_t\in(0,1)$ and posts
$(1-\phi_t)P_th$ at purchase. Closing occurs before current income and rebates
are available, and unsecured borrowing is unavailable. Thus only liquid
wealth $b_i$ can satisfy the down payment. The owner cannot buy more housing
when old or rent its retained housing to another household. At death the
estate sells it back into the fixed stock, and the proceeds enter warm-glow
estate expenditure.

\section{Household problems}

Households take prices and transfers as given. In both tenures, $a'$ denotes
financial wealth on entering old age, net of any mortgage repayment.

\textbf{Young renters.}
Conditional on renting, the young household solves
\begin{align}
 W_t^R(i)=\max_{c,h,n,a'}\;&u_t^y(c,h,n)+\beta V_{t+1}^R(a') \nonumber\\
 \text{s.t.}\quad&
 c+qa'+q\,r_th=y_i+b_i+T_t,\nonumber\\[-2pt]
 &h\leq h_R^{\max},\qquad a'\geq0.
 \label{eq:amend_young_renter}
\end{align}
The renter sets aside $qa'$ to enter old age with financial wealth $a'$.

\textbf{Young owners.}
Conditional on owning, the young household solves
\begin{align}
 W_t^O(i)=\max_{c,h,n,a'}\;&u_t^y(c,h,n)+\beta V_{t+1}^O(a',h) \nonumber\\
 \text{s.t.}\quad&
 c+qa'+(1+q\tau^p)P_th=y_i+b_i+T_t,\nonumber\\[-2pt]
 &(1-\phi_t)P_th\leq b_i,\qquad h\leq h_O^{\max},\nonumber\\[-2pt]
 &qa'+\phi_tP_th\geq0.
 \label{eq:amend_young_owner}
\end{align}
The owner sets aside $qa'+\phi_tP_th$ in bonds. After subtracting mortgage
debt, its old-age financial wealth is $a'$, alongside title
$H=h$. The last constraint requires nonnegative bond holdings; $a'$ itself
can be negative. Below, positive saving means positive bond holdings,
before subtracting mortgage debt. The property tax is also set aside in a bond.

\textbf{Old households.}
An old household has financial wealth $a$. An owner also carries housing $H$
from young age and retains $h^2\leq H$ for use while old. Let $a^e\geq0$
denote financial saving during old age. For a household old at date $t$,
the estate left at death is
\begin{equation}
 e=\begin{cases}
 R_f a^e, & \text{renter},\\
 R_f a^e+P_{t+1}h^2, & \text{owner}.
 \end{cases}
 \label{eq:amend_estate_components}
\end{equation}
Financial saving earns the bond return, and the owner's retained housing
is sold at the end of old age. Housing sold at the beginning of old age,
$H-h^2$, instead provides current funds. The owner's budget before
substituting for financial saving is
\begin{equation}
 c^2+a^e+q\tau^pP_t h^2=a+P_t(H-h^2)+T_t.
 \label{eq:amend_old_cash_budget}
\end{equation}
Substituting $a^e=q(e-P_{t+1}h^2)$ and using
\eqref{eq:amend_usercost} gives the owner's budget below.
\par\needspace{12\baselineskip}
Renters and owners solve, respectively,
\begin{align}
 V_t^R(a)&=\max_{c^2,h^2,e}u^o(c^2,h^2,e),\nonumber\\
 &c^2+qe+q\,r_th^2=a+T_t,\quad 0<h^2\leq h_R^{\max},
 \label{eq:amend_old_renter}\\
 V_t^O(a,H)&=\max_{c^2,h^2,e}u^o(c^2,h^2,e),\nonumber\\
 &c^2+qe+q\,r_th^2=a+P_tH+T_t,\quad
 0<h^2\leq H,\quad e\geq P_{t+1}h^2.
 \label{eq:amend_old_budget}
\end{align}
The owner's estate bound $e\geq P_{t+1}h^2$ is equivalent to
nonnegative financial saving $a^e\geq0$.

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

% Prior conditional settlement argument; not the desired equilibrium characterization.
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
```

### File: output/model/simplified_olg_amendments/oracle_analytic_efficiency_response.md
```md
# Pro response: analytical equilibrium inefficiency

Source: https://chatgpt.com/c/6a9f2beb-c8b0-83e9-b649-3304f9801ad3

Retrieved 7 September 2026. External model response, preserved for independent review; this is not an adopted manuscript result. Equations extracted from the rendered math labels. Attachment citation buttons omitted.

## Verdict

**The unrestricted claim is false.** The model admits stationary equilibria that are efficient for the packet’s fixed-fertility, fixed-tenure, datewise welfare comparison.

However, a stronger result is available than existence of one inefficient equilibrium: **under explicit primitive inequalities, every stationary equilibrium is inefficient.** The result below retains endogenous fertility, positive child goods and space costs, finite housing caps, and finite-scale logistic tenure choice. It uses homogeneous income and wealth and zero property tax. It also covers old owners whose estate constraint binds—not only those leaving positive financial assets.

---

## 1. Setup and an analytical efficient class

Set $\tau^p=0$ , let all entrants have the same $(y,b)$ , and keep all preference weights and child costs strictly positive. Write

$$
\begin{gathered} w=y+b,\qquad d=1-q,\qquad u=qr=dP,\\ A=1+\gamma+\omega_B,\qquad M=1+\beta A,\qquad D=M+\alpha+\vartheta,\\ \rho=h_R^{\max},\qquad \bar h=h_O^{\max},\qquad B=\frac{b}{1-\phi}. \end{gathered}
$$

Thus $Ph\le B$ is the down-payment constraint. Importantly, the owner’s gross bond holdings are $qa'+\phi Ph$ , not $qa'$ .

I use **stationarity in levels** : $Y=O=N>0$ and

$$
\nu\bar n=1.
$$

Individual fertility remains a household choice; this is an equilibrium consistency condition. The stationary population $N$ is part of the stationary state, with its scale determined by housing clearing.

### A class of efficient stationary equilibria

Define, entirely from primitives,

$$
u_*=\frac{\nu\vartheta w/D-\chi}{\kappa}, \qquad P_*=\frac{u_*}{d}, \qquad x_*=\frac wD, \qquad n_*=\frac1\nu,
$$

and

$$
\begin{aligned} h_*^Y&=\frac{\alpha x_*}{u_*}+\frac{\kappa}{\nu}, & z_*&=\frac{\beta A x_*}{q},\\ c_*^2&=\frac{\beta x_*}{q}, & h_*^2&=\frac{\beta\gamma x_*}{q u_*}, & e_*&=\frac{\beta\omega_B x_*}{q^2}. \end{aligned}
$$

Suppose the following explicit inequalities hold:

$$
\boxed{ \begin{gathered} u_*>0,\qquad d\omega_B>q\gamma,\\ h_*^2<h_*^Y<\rho<\bar h,\\ P_*h_*^Y<B,\qquad qz_*+(\phi-q)P_*h_*^Y>0. \end{gathered}} \tag{E}
$$

These formulas solve the relaxed household problems, and (E) verifies that the omitted constraints are slack. Both tenures therefore choose the displayed real allocations, with

$$
a_R'=z_*, \qquad a_O'=z_*-P_*h_*^Y.
$$

Their conditional values coincide, so finite logistic tastes produce strictly positive shares of both tenures. Setting

$$
N=\frac{\bar H}{h_*^Y+h_*^2}
$$

completes the stationary equilibrium.

This allocation is **globally efficient within the specified datewise feasible set** , not merely immune to one old-to-young transfer. At fixed fertility and estates,

$$
\frac{\alpha x_*}{h_*^Y-\kappa n_*} = \frac{\gamma c_*^2}{h_*^2} =u_*.
$$

Concavity gives, for each young household,

$$
x_*\Delta U_i^Y\le \Delta c_i+u_*\Delta h_i,
$$

and, for each old household,

$$
c_*^2\Delta U_j^O\le \Delta c_j^2+u_*\Delta h_j^2.
$$

Summing over households makes the right-hand side zero for any consumption- and housing-preserving reallocation. A Pareto improvement on a positive-measure set is therefore impossible.

Thus the **presence** of credit and rental-size limits does not establish inefficiency. The next result makes both young-household limits bind with strictly positive multipliers.

---

## 2. Primitive sufficient conditions for inefficiency of every stationary equilibrium

Define

$$
v=\max\left\{d,\frac{\gamma}{\gamma+\omega_B}\right\}, \qquad \ell=dB,\qquad E=w-\ell,
$$

and

$$
\underline\ell=\frac{\beta\gamma w}{qM+\beta\gamma}.
$$

The central restrictions are

$$
\boxed{ \underline\ell<\ell< \frac{\alpha d\,w}{\alpha d+v(M+\vartheta)} } \tag{A}
$$

and

$$
\boxed{ \frac{\beta A E}{M+\vartheta}>(q-\phi)_+B. } \tag{B}
$$

Condition (B) is automatic when $\phi\ge q$ , but that restriction on $\phi$ is not required.

We next give a demographic restriction using only explicit functions of primitives.

### An explicit fertility function

For $m,E,H>0$ , define

$$
\mathcal N_m(E,H)= \frac{2\vartheta EH} {K_m(E,H)+ \sqrt{K_m(E,H)^2- 4\chi\kappa(m+\alpha+\vartheta)\vartheta EH}},
$$

where

$$
K_m(E,H) =\kappa E(\alpha+\vartheta)+\chi H(m+\vartheta).
$$

It is the unique feasible solution of

$$
mx+\chi n=E, \qquad \frac{\vartheta}{n} =\frac{\chi}{x}+\frac{\alpha\kappa}{H-\kappa n}. \tag{1}
$$

It is strictly increasing in both $E$ and $H$ .

### Primitive price endpoints

Set

$$
t_-=\max\left\{\underline\ell,\frac{\ell\rho}{\bar h}\right\}, \qquad g=\frac{\chi\rho}{\kappa},
$$

and

$$
t_+= \frac{(\alpha+\vartheta)w-Dg+ \sqrt{\big[(\alpha+\vartheta)w-Dg\big]^2+4D\alpha wg}} {2D}.
$$

Define

$$
P_-=\frac{t_-}{d\rho}, \qquad P_+=\frac{t_+}{d\rho},
$$

and the four endpoint fertilities

$$
n_O^\pm=\mathcal N_M\!\left(E,\frac{\ell\rho}{t_\pm}\right), \qquad n_R^\pm=\mathcal N_M(w-t_\pm,\rho).
$$

The final restriction is

$$
\boxed{ \max\{n_O^+,n_R^+\} < \frac1\nu < \min\{n_O^-,n_R^-\}. } \tag{C}
$$

This is an **analytical root-bracketing condition** , rather than a closed-form equilibrium price. Its endpoints and inequalities contain no unknown equilibrium quantities.

### Proposition

Under (A)–(C), with $0<\rho<\bar h<\infty$ , $\bar H>0$ , finite $\bar\xi$ , and $0<\sigma_\xi<\infty$ :

**A stationary state exists. Every stationary equilibrium has $P\in(P_-,P_+)$ , and every such equilibrium is Pareto inefficient for the packet’s datewise comparison.**

In every such equilibrium:

- Young owners have a strictly binding down-payment constraint, a slack physical cap, and positive gross bond holdings.
- Young renters have a strictly binding rental-size cap.
- Old owners sell part of their inherited housing; old renters’ housing caps are slack.
- Both tenures have positive measure.

The old owner’s financial estate is positive when $d\omega_B>q\gamma$ , zero with a strictly binding estate constraint when $d\omega_B<q\gamma$ , and zero with a zero estate-constraint multiplier at equality.

---

## 3. Proof: establish the constraints rather than assume them

### Old owners and the relevant housing value

Let

$$
z=a'+Ph
$$

be the young owner’s total wealth entering old age, including the inherited house. Its lifetime budget becomes

$$
c+qz+uh=w. \tag{2}
$$

The old problem has resources $z$ , housing-retention limit $h^2\le h$ , and estate bound $e\ge Ph^2$ , exactly as obtained from the packet’s budgets.

When retention is slack, solving the old problem gives

$$
c^2=\frac zA, \qquad h^2=\frac{\gamma z}{AvP}, \tag{3}
$$

and

$$
e= \begin{cases} \displaystyle\frac{\omega_B z}{Aq}, &d\omega_B\ge q\gamma,\\[5pt] \displaystyle Ph^2, &d\omega_B<q\gamma. \end{cases} \tag{4}
$$

Consequently,

$$
MV^O=\frac{\gamma c^2}{h^2}=vP. \tag{5}
$$

Thus the donor’s fixed-estate housing value equals user cost only in the financially slack branch. When the estate bound binds strictly, $v>d$ : compensation must cover a higher housing value.

### The owner chooses the largest privately admissible house

Temporarily drop the nonnegative-gross-bond constraint, retaining all other constraints. Write $K=Ph\le B$ .

The fertility first-order condition implies

$$
\chi n<\vartheta x. \tag{6}
$$

Homogeneity of the old problem, together with the saving first-order condition, gives

$$
qz\le\beta A x,
$$

with equality when old housing retention is slack. Combining this with (2),

$$
w-dK=x+\chi n+qz<(M+\vartheta)x.
$$

Therefore

$$
x>\frac{w-dK}{M+\vartheta}. \tag{7}
$$

Since $s<h$ , condition (A) now implies

$$
\begin{aligned} MV^Y =\frac{\alpha x}{s} &> \frac{\alpha(w-dK)P}{(M+\vartheta)K}\\ &\ge \frac{\alpha E P}{(M+\vartheta)B} >vP\ge u. \end{aligned} \tag{8}
$$

Increasing current housing also weakly relaxes the future retention limit. Hence the household’s optimized value is strictly increasing in $h$ up to its current limits:

$$
\boxed{h_O^Y(P)=\min\{\bar h,B/P\}.} \tag{9}
$$

The omitted gross-bond constraint is indeed slack. If retention is slack,

$$
qa'+\phi Ph =qz+(\phi-q)K =\beta A x+(\phi-q)K>0
$$

by (7) and (B). If retention binds, $h^2=h$ , the old budget and $e\ge Ph$ imply

$$
z=c^2+qe+uh>Ph=K,
$$

so $a'=z-K>0$ , which also gives positive gross bonds. Thus dropping that constraint changed nothing.

### Policies inside the bracket

For $P\ge B/\bar h$ , consider

$$
h_O^Y=B/P,\qquad n_O=\mathcal N_M(E,B/P),\qquad x_O=\frac{E-\chi n_O}{M}. \tag{10}
$$

The lower inequality in (A) verifies old retention slack:

$$
\frac{h_O^2}{h_O^Y} =\frac{\beta\gamma x_O}{qvB} \le\frac{\beta\gamma x_O}{q\ell} < \frac{\beta\gamma E}{qM\ell}<1. \tag{11}
$$

Hence (10) is the actual owner solution, with

$$
z_O=\frac{\beta A x_O}{q}, \qquad c_O^2=\frac{\beta x_O}{q}, \qquad h_O^2=\frac{\beta\gamma x_O}{qvP}. \tag{12}
$$

For renters inside the proposed bracket,

$$
h_R^Y=\rho,\qquad n_R=\mathcal N_M(w-dP\rho,\rho),\qquad x_R=\frac{w-dP\rho-\chi n_R}{M}. \tag{13}
$$

Their old housing is

$$
h_R^2=\frac{\beta\gamma x_R}{qdP}<\rho, \tag{14}
$$

because $dP\rho\ge t_-\ge\underline\ell$ .

Finally, without the young rental cap, renter housing demand is

$$
h_{R,\mathrm{free}}^Y(P) =\frac wD\left[ \frac{\alpha}{dP} +\frac{\kappa\vartheta}{\chi+\kappa dP} \right]. \tag{15}
$$

The equation $h_{R,\mathrm{free}}^Y=\rho$ has the explicit solution $dP\rho=t_+$ . Thus the young rental cap binds with positive multiplier whenever $P<P_+$ .

At an interior-bracket price, the owner’s physical cap is slack, its gross bonds are positive, and old retention is slack. Its down-payment multiplier $\mu$ therefore satisfies

$$
MV_O^Y-u=(1-\phi)P\mu x_O.
$$

Equation (8) proves $\mu>0$ . **No assumption that the future estate bound is slack enters this envelope calculation.**

---

## 4. Why the bracket captures every stationary equilibrium

The crucial global fact is:

$$
\boxed{\text{Conditional owner and renter fertility each decrease strictly with }P.}
$$

Their logistic mixture need not be monotone.

For owners above $B/\bar h$ , this follows from (10): $E$ is constant and $B/P$ decreases. Below that price, (9) gives $h=\bar h$ . Depending on the old constraints, fertility is one of

$$
\begin{aligned} &\mathcal N_M(w-dP\bar h,\bar h),\\ &\mathcal N_{M_0}(w-(1+q)dP\bar h,\bar h), \qquad M_0=1+\beta(1+\omega_B),\\ &\mathcal N_{1+\beta}(w-P\bar h,\bar h). \end{aligned} \tag{16}
$$

These correspond respectively to slack retention; binding retention with slack estate bound; and both old constraints binding. Each decreases strictly on its applicable range. The branches join continuously. The relevant wealth-to-house-value ratios decrease along each branch, so the old regimes cannot switch back and forth.

For renters, (A) implies $\alpha>\beta\gamma/q$ . An old rental cap therefore cannot bind while the young rental cap is slack. The successive fertility schedules are

$$
\mathcal N_{M_0}(w-(1+q)dP\rho,\rho),\qquad \mathcal N_M(w-dP\rho,\rho),\qquad \frac{\vartheta w}{D(\chi+\kappa dP)}. \tag{17}
$$

Again, each decreases strictly and the branches join continuously.

Now observe

$$
t_-<\ell<\frac{\alpha w}{D}<t_+<w.
$$

At the intermediate price

$$
P_0=\frac B\rho,
$$

both tenures choose the same current housing, consumption, and fertility:

$$
n_O(P_0)=n_R(P_0)=n_0:=\mathcal N_M(E,\rho).
$$

Consequently,

$$
\max\{n_O^+,n_R^+\} <n_0< \min\{n_O^-,n_R^-\}. \tag{18}
$$

This also proves that the demographic interval in (C) is nonempty.

Let $\pi(P)$ be the logistic owner share and set

$$
\bar n(P)=\pi(P)n_O(P)+(1-\pi(P))n_R(P).
$$

It is continuous. Condition (C) gives

$$
\bar n(P_-)>\nu^{-1}>\bar n(P_+),
$$

so a stationary fertility root exists. Global conditional monotonicity also gives

$$
P\le P_-\implies\bar n(P)>\nu^{-1}, \qquad P\ge P_+\implies\bar n(P)<\nu^{-1}.
$$

**Every stationary root is therefore inside the bracket.**

For any root, let $\pi=\pi(P)$ and set

$$
N= \frac{\bar H} {\pi(h_O^Y+h_O^2)+(1-\pi)(\rho+h_R^2)}. \tag{19}
$$

Together with the induced asset distribution, this satisfies housing clearing and stationary demography. Finite logistic tastes give $0<\pi<1$ .

---

## 5. The compensated Pareto improvement

At every covered equilibrium,

$$
MV_O^Y>vP=MV_O^O. \tag{20}
$$

Select a common positive submass of young and old owners. This is a proof device, not an assumption about matching or cohort shares.

Give each selected young owner $\epsilon$ additional housing and take it from a selected old owner. Keep fertility, tenure, and the old estate fixed. The donor’s exact consumption compensation is

$$
D(\epsilon) =c^2\left[ \left(\frac{h^2}{h^2-\epsilon}\right)^\gamma-1 \right].
$$

It preserves the donor’s utility exactly, and

$$
D'(0)=\frac{\gamma c^2}{h^2}=vP.
$$

The young owner’s gain is

$$
\Delta U^Y = \log\!\left(1-\frac{D(\epsilon)}x\right) +\alpha\log\!\left(1+\frac{\epsilon}{s}\right),
$$

so

$$
\left.\frac{d\Delta U^Y}{d\epsilon}\right|_{\epsilon=0} =\frac{MV_O^Y-vP}{x}>0.
$$

A sufficiently small positive transfer therefore benefits the young strictly while respecting their physical cap and keeping all consumption positive.

The packet’s settlement applies unchanged. Set

$$
L(\epsilon)=D(\epsilon)-u\epsilon.
$$

The old owner receives $P\epsilon+L(\epsilon)$ , consumes $D(\epsilon)$ , and buys a bond costing $qP\epsilon$ , replacing the title removed from its estate. The young owner’s new borrowing proceeds are also $qP\epsilon$ ; next period it sells the additional title for $P\epsilon$ and repays that new debt. Existing estates, future real allocations, and promised repayments on existing loans remain unchanged.

This also works when the old owner initially has **zero** financial saving: the settlement creates positive financial saving to replace the transferred estate housing.

---

## 6. Economic meaning and scope

The lower bound on $\ell$ prevents the inherited owner house from being so small that the old household wants to retain all of it. The upper bound makes liquid wealth sufficiently scarce relative to income and young housing preferences that young owners value additional housing more than eligible old donors do. Condition (B) ensures positive **gross** saving; it does not require positive net $a'$ .

The factor $v$ is important. When the estate constraint binds, old housing is worth more than ordinary user cost at fixed estate. The upper liquidity bound automatically becomes tighter rather than incorrectly treating such donors as having marginal value $u$ .

The region is demonstrably nonempty. For arbitrary positive old-age preference weights, choose

$$
\alpha> \frac{\beta\gamma v(M+\vartheta)}{dqM}.
$$

Then the interval in (A) is nonempty. Choosing $\phi\ge q$ makes (B) automatic; for any admissible $\ell$ and $w>0$ , set

$$
b=\frac{(1-\phi)\ell}{d}, \qquad y=w-b>0.
$$

Positive $\chi,\kappa$ and any finite ordered caps are admissible, and (18) supplies a nonempty interval for $\nu$ .

For **existence alone** , (C) can be weakened to the taste-dependent endpoint condition

$$
\bar n(P_+)<\nu^{-1}<\bar n(P_-),
$$

using the explicit endpoint policies in the logistic rule. That weaker condition guarantees an inefficient stationary state but does not exclude other stationary equilibria.

The proved result is for homogeneous entrants and zero property tax. It does not establish price uniqueness, convergence from arbitrary initial states, or inefficiency outside the displayed region. If the stationary cohort mass $N$ is prescribed independently, housing-clearing compatibility in (19) must additionally hold. No unresolved step remains inside the proposition; heterogeneous entrants and positive, price-dependent tax rebates require a separate equilibrium argument.
```
Copied markdown to clipboard (~14.7k tokens; 2 files).
