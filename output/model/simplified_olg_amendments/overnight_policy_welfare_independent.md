# Independent check: living-household welfare from a rebated property tax

September 9, 2026 UTC. This note changes no household equation or welfare weight. It uses the uniformly slack finance/cap/estate branch and local transition construction in the independent transition memo. No live Pro output, model solve, or numerical search was used.

**Result.** An exact analytical family has a positive intervention-date utilitarian welfare derivative for every interior ownership share and every bounded proportional endowment distribution specified below. Young aggregate welfare rises; initial-old aggregate welfare can fall at high ownership shares. The capital price falls on impact even though its stationary endpoint rises. Outside rental financiers bear an explicit capital loss and are excluded from the maintained utility sum. This is a market-feasible example, not a mortgage-mechanism or Pareto-inefficiency theorem.

## 1. The intervention-date envelope

Normalize the stationary type law \(F\) to one; each living age group has mass \(N\). Write \(w_i=y_i^y+b_i\), \(v_i=y_i^o\), \(K=1+\gamma+\omega_B\), and \(D=1+\alpha+\vartheta+\beta K\). At the zero-tax stationary reference,
\[
x_i=(w_i+qv_i)/D,\qquad
c_i^o=\beta x_i/q,\qquad
A_y=\frac{\alpha}{p}+\frac{\kappa\vartheta}{\chi+\kappa p}.
\]
Here \(A_yx_i\) is the house bought by the preceding young cohort, which becomes the initial old owner's inherited title. It is not its current old-age housing choice. Put
\[
I_x=\int x_i^{-1}\,dF(i).
\]
Owner probability \(\pi\) is constant and independent of endowments in this branch because the two unconstrained conditional value functions coincide.

Young continuation utility and old warm-glow utility remain in the objective and are allowed to change. The original reduced budgets and the envelope theorem give
\[
\frac{d\mathcal W_y}{N}
=I_x(dT_0+q\,dT_1)-A_y\,dp_0-\beta\gamma\,dp_1/p. \tag{1}
\]
Initial-old resources change by \(dT_0+H_i^{own}dP_0\), holding their inherited financial obligations fixed. Since their marginal utility of resources is \(q/(\beta x_i)\),
\[
\frac{d\mathcal W_o}{N}
=\frac q\beta I_x\,dT_0-\gamma\,dp_0/p
+\pi\frac q\beta A_y\,dP_0. \tag{2}
\]
The estate term is included: optimizing \(\log c^o+\gamma\log h^o+\omega_B\log e\) gives marginal value \(1/c^o\). There is no omitted fixed-estate assumption.

Therefore the proposed envelope is correct:
\[
\boxed{\begin{aligned}
\frac{d\mathcal W}{N}
={}&I_x[(1+q/\beta)dT_0+q\,dT_1]\\
&-\left[\frac{\alpha+\gamma}{p}
+\frac{\kappa\vartheta}{\chi+\kappa p}\right]dp_0
-\beta\gamma\,dp_1/p
+\pi(q/\beta)A_y\,dP_0 .
\end{aligned}} \tag{3}
\]
Initial cohort masses are fixed, so their welfare weights acquire no population derivative. Later entrants are not added to this intervention-date sum.

All terms can be evaluated from primitives in the specified branch. The rebate benefit is weighted by reciprocal adult consumption; higher current and next-period housing-service costs reduce utility; the initial owner title has a capital gain or loss. A positive terminal population derivative alone does not sign (3).

## 2. Exact analytical family and admissibility

Set
\[
q=\beta=\tfrac12,\quad
\alpha=\gamma=\vartheta=\kappa=\nu=1,\quad
\omega_B=2,\quad\chi=\tfrac32.
\]
Let \(w\) have any bounded nondegenerate positive distribution with mean eight and let \(v=w/2\). Then
\[
D=5,\quad x_i=w_i/4,\quad X=2,\quad
p=\tfrac12,\quad P=1,\quad N=\bar H/9.
\]
The reference choices are
\[
n_i=w_i/8,\quad h_i^y=5w_i/8,\quad
c_i^o=w_i/4,\quad h_i^o=w_i/2,\quad e_i=w_i.
\]
Take \(\phi\ge q\), any interior logistic owner share, and finite \(H_R<H_O\) with \(H_R>5\sup w/8\). The renter saves \(w_i/2>0\). A young owner has net financial wealth \(-w_i/8\), but its mortgage inequality has strict slack; at \(\phi=q\), current cash spending is \(3w_i/4<w_i\). Old owner financial saving is \(w_i/4>0\). Thus all required margins are uniformly strict on the bounded support. Heterogeneity and both tenures remain.

## 3. The exact first-order equilibrium path

A prime denotes the derivative with respect to a permanent tax introduced unexpectedly at date zero. Set \(\ell_t=p_t'/p\), \(y_t=Y_t'/N\), and \(g_X=X_t'/X\). At the zero-tax reference, the fiscal base derivative is multiplied by zero, so
\[
T_t'=\frac{qP\bar H}{2N}=\frac94,\qquad
g_X=\frac{(1+q)T_t'}{DX}=\frac{27}{80}
\]
at every date. The transition coefficients are
\[
d=33/4,\quad \eta=1/4,\quad
\mathcal A=28/33,\quad k=4/33,\quad f=8/11,\quad G=20/33.
\]
For \(t\ge1\),
\[
y_{t+1}=\mathcal A y_t-ky_{t-1}+f g_X,\qquad
y_0=0,\quad y_1=g_X-\ell_0/4.
\]
To check the price jump independently, let \(S_y=\sum_{t\ge0}q^t y_t\) and \(S_\ell=\sum_{t\ge0}q^t\ell_t\). Summing the stable recurrence and demographic equation gives
\[
S_y=\frac{q y_1+f g_Xq^2/(1-q)}G
=\frac{57}{40}g_X-\frac{33}{160}\ell_0,
\]
\[
S_\ell=\frac{23}{10}g_X+\frac{33}{40}\ell_0,\qquad
P_0'=\tfrac12S_\ell-1
=\frac{33}{80}\ell_0-\frac{979}{1600}. \tag{4}
\]
The inherited owner exposure is \(5\pi\) per old household. Date-zero housing clearing therefore requires
\[
\frac{33}{4}\ell_0=\frac{45}{16}+\frac{5\pi}{2}P_0'.
\]
Its coefficient is \(\mathcal J=33(8-\pi)/32>0\). Solving yields
\[
\boxed{
\ell_0=\frac{1800-979\pi}{660(8-\pi)},\qquad
P_0'=-\frac{377}{100(8-\pi)},\qquad
\ell_1=\frac{63}{110}-\frac5{33}\ell_0.
} \tag{5}
\]
Thus impact capital prices fall for every \(0<\pi<1\). Separately, the stationary formulas give \(P_\infty'=7/20>0\) and \(N_\infty'/N=9/10>0\). No date-by-date fertility ordering is inferred.

## 4. Welfare signs, including the initial old

Define the arithmetic-to-harmonic mean ratio
\[
B_w=8\int w_i^{-1}\,dF(i)\ge1.
\]
Then \(I_x=B_w/2\). Substituting (5) and the funded rebates into (3) gives
\[
\frac{\mathcal W'}N
=\frac{45}{16}B_w-\frac94\ell_0-\frac12\ell_1
+\frac{5\pi}{2}P_0',
\]
or, exactly,
\[
\boxed{
\frac{\mathcal W'}N
=\frac{45}{16}(B_w-1)
+\frac{622008-380105\pi}{43560(8-\pi)}
>0.} \tag{6}
\]
The second term is decreasing in \(\pi\) and bounded below by
\(241903/304920>0\). Equal rebates weighted by marginal utility therefore outweigh service-cost increases and incumbent-owner capital losses even in the equal-resource benchmark. Dispersion increases the first-order aggregate gain within this family; no income-risk mechanism is being asserted.

The separate cohort derivatives are
\[
\boxed{
\frac{\mathcal W_y'}N
=\frac{27}{16}(B_w-1)
+\frac{348768+14839\pi}{43560(8-\pi)}>0,
} \tag{7}
\]
\[
\boxed{
\frac{\mathcal W_o'}N
=\frac98(B_w-1)
+\frac{1035-1496\pi}{165(8-\pi)}.
} \tag{8}
\]
Young aggregate welfare is uniformly positive; its equal-resource lower bound is \(1211/1210\). Initial-old aggregate welfare is positive below ownership share \(1035/1496\), approximately 0.692, without requiring dispersion. At higher shares it can be negative when dispersion is small. Young gains still dominate in (6). At the threshold equality, dispersion makes the old gain strict; otherwise their derivative is zero. These are cohort aggregates, not individual Pareto comparisons.

## 5. Outside asset account and finite local interpretation

Under the explicit outside-residual-owner convention, let \(H^{rent}_{-1}\) be the entire inherited rental stock. Its surprise revaluation is
\[
(A^{rent}_{outside})'=H^{rent}_{-1}P_0'<0.
\]
If outside intermediaries own all reference rental units, then
\[
H^{rent}_{-1}=(1-\pi)\bar H=9(1-\pi)N,\qquad
\boxed{\frac{(A^{rent}_{outside})'}N
=-\frac{3393(1-\pi)}{100(8-\pi)}.} \tag{9}
\]
Using only the rentals occupied by the preceding young would omit the rentals occupied by the exiting old. This monetary capital loss belongs to the outside residual owner's account. It is not erased, added to domestic resources, or included in the maintained living-household utility sum. A welfare comparison including outside claimants would need their preferences and weights; monetary loss cannot simply be added to utility.

The fiscal ledger holds exactly along the underlying paths; at the reference its derivative is \(2NT'=qP\bar H\). Existing creditor payments are not written down. The result inherits the memo's explicit residual-financier convention and does not assume that a zero-equity intermediary can absorb losses without an account.

Combined with the memo's local implicit-function theorem, the strict derivative implies a welfare increase for sufficiently small finite permanent tax increases. The same conclusion holds at intervention dates on sufficiently nearby small preference-decline transitions, using their actual inherited state and the same remaining preference path on both tax and baseline paths. This extension follows from differentiability of the actual equilibrium and welfare; equation (3)'s stationary pairing must not be reused literally at a nonstationary date.

Uniform slackness, positive estate saving, bounded-price selection, the outside ownership account, and local policy size remain qualifications. The example establishes neither a binding-mortgage mechanism, an optimal tax, nor a Pareto improvement across all agents and generations.

**Verification:** exact symbolic identities confirm the date-zero housing residual, discounted asset-price equation, next service price, cohort decomposition, and total welfare expression. No model simulation was used.

**One Pro question:** does this positive living-household welfare derivative, together with the locally convergent permanent-tax path, survive a positive mass of binding young financing constraints when initial capital revaluation, endogenous tenure, future rebates, and outside asset losses are all retained?
