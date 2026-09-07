🧿 oracle 0.13.0 — From shruggy agents to shippable PRs.
[SYSTEM]
You are Oracle, a focused one-shot problem solver. Emphasize direct answers and cite any files referenced.

[USER]
We need a simple analytical economic theorem, not another numerical existence result. Please work independently from the attached original OLG household and equilibrium equations.

The objective is to establish that competitive equilibrium can allocate too little housing to young households, relative to a planner who respects the same physical tenure limits but can relax private financing constraints. Fertility is held fixed for this welfare comparison; valuing a larger population is not part of the argument. The author wants a clear problem-set-style result, ideally broad and with interpretable restrictions on primitive parameters.

Proceed in this order:
1. First test whether equilibrium must be inefficient generally. Distinguish the mere presence of borrowing and rental-size limits from cases where both actually bind. If a general theorem is false, establish an analytical counterexample or a class of efficient allocations, with the precise welfare scope. Failure of one old-to-young transfer test is not proof of efficiency.
2. Then derive explicit, economically interpretable sufficient conditions on primitives under which the stationary equilibrium is inefficient because housing can be moved from old to young. Establish the household constraint pattern as an equilibrium outcome; do not put the desired pattern or an MRS inequality directly into the theorem's assumptions and call the problem solved. Work from primitives such as income, liquid wealth, credit limits, housing caps, preferences, and household formation. An implicit equilibrium price/multiplier condition is not the requested final condition. Analytical root bracketing with endpoints specified entirely in primitives is acceptable if transparent; label that distinction. Try to establish the largest useful class, and identify when you prove existence of an inefficient equilibrium versus inefficiency of every equilibrium.
3. A homogeneous income-and-wealth specialization is acceptable as a first theory result, but keep endogenous fertility in the competitive household problem, positive child goods/space costs, finite housing caps, and mixed tenure with finite logistic taste scale where possible. If a simplification such as zero property tax or a different tenure assumption is necessary, state it explicitly and explain what it buys. Do not silently change the model or the planner.

Planner scope: at a chosen date, reallocate consumption and housing, keeping fertility and tenure fixed and respecting physical housing limits. Can relax private financing restrictions. Seek a compensated improvement preserving each existing estate, all future real allocations, and creditors' promised repayments, as in the attached settlement. Positive aggregate reallocation needs eligible positive-measure households, but equal group masses, individual matching, or uniform gaps should be proof devices if needed rather than arbitrary economic hypotheses. Do not recast this as a market policy, tax implementation, constrained-efficiency-with-markets result, or a new intergenerational welfare criterion. If broader reallocations are required to assess universal efficiency, state clearly the feasible set and distinguish it from the datewise result.

What is already known, but does not answer our question: given a strictly positive young-vs-old marginal housing value gap, exact compensation proves a Pareto improvement. We also constructed one numerical reference equilibrium and a small open neighborhood where the needed constraints hold, and proved local transitions with computer-verified bounds. The author specifically rejects that approach as the main illustrative theorem. Do NOT return a calibrated example, a numerical sweep, a computer-assisted interval certificate, or existence by continuity around chosen numbers. We need analytical inequalities in primitives and their economic meaning.

Useful algebra to verify, not assume blindly: for an owner with positive gross bond holdings, a slack current physical cap, and slack future housing-retention bound, the young housing value exceeds current user cost if the down-payment limit has positive shadow value. Future estate slack is not needed for that envelope step. Old housing value equals user cost when the old owner sells part of the inherited house and leaves financial assets as well as housing. The difficulty is deriving those statuses from primitives. A constraint holding with equality need not have positive shadow value.

Notation in the attachment: young a' is NET next-period financial wealth, after mortgage repayment. Owner gross current bond purchases equal q a'+phi_t P_t h >= 0, so positive saving does not mean a'>0. In the old-age estate definition, a^e is current financial saving and e=R_f a^e+P_{t+1}h^2 for owners. h^2 is old-age housing, not a square. The attached settlement uses equation references from the full note: u_t=q r_t=P_t+q tau^p P_t-qP_{t+1}; x_i=c_i-chi n_i; the donor's exact compensation is D_j(epsilon)=c_j^2[(h_j^2/(h_j^2-epsilon))^gamma-1]. No real moving costs in the core question.

Deliver a concise verdict first, then the strongest analytical proposition you can actually prove, with a complete readable derivation, an explicit list of primitive restrictions, their economic interpretation, and any unresolved obstruction. Aim for a few pages, not a long survey. If you cannot establish the desired region, say exactly which step remains unproved. Preserve the author's notation and plain prose. This is pure theory: there is no repository implementation task or need to run software.

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
Copied markdown to clipboard (~8.5k tokens; 1 files).
