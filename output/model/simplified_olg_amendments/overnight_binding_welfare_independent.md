# Independent welfare theorem with binding young finance

September 9, 2026. Analytical diagnostic; no model runs or numerical-reference-point argument. The model remains frozen. This note uses the exact local equilibrium path and explicit inverse bound in `overnight_binding_transition_independent.md`. The criterion is the sum of remaining utility of households alive at intervention, with each household weighted one. Young continuation utility retains its private discount factor \(\beta\); future estates are endogenous. Outside residual rental financiers are excluded from this utility sum and their capital loss is recorded below.

**Result.** A primitive mean-income inequality, together with an explicit positive upper bound on the ownership share, gives a welfare-improving small permanent rebated property tax in the branch where all young financial constraints bind strictly. An exact heterogeneous family satisfies the inequality and has positive young and initial-old cohort welfare effects, higher initial fertility, and a larger limiting population. Neither the ownership bound nor the endowment restrictions are claimed to be empirically mild.

## 1. Branch and household envelopes

Let \(w_i=y_i^y+b_i\), \(v_i=y_i^o\), \(W=\int w_i\,dF\), \(V=\int v_i\,dF\), \(E=1+\alpha+\vartheta\), and \(K=1+\gamma+\omega_B\). The probability distribution \(F\) has compact positive endowment support. Maintain
\[
\phi=q,\quad \tau=0,\quad
\frac{v_i}{w_i}>\frac{\beta K}{qE}\ \text{uniformly},\quad
\omega_B(1-q)>q\gamma,\quad
\chi<\frac{\nu\vartheta W}{E},
\tag{1}
\]
with finite tenure caps uniformly above both conditional housing choices. The transition memo constructs the stationary equilibrium explicitly and proves a nearby convergent equilibrium path for the ownership interval used below. No standalone restriction on the sign of \(\beta-q\) is needed; (1) still imposes an income-dependent upper bound on \(\beta\).

At this reference, define adult consumption \(x_i=c_i^y-\chi n_i\), adult space \(s_i=h_i^y-\kappa n_i\), and \(A=\alpha+\vartheta\eta\), where \(\eta=\kappa p/(\chi+\kappa p)\). Both tenures choose
\[
x_i=w_i/E,\quad s_i=\alpha x_i/p,\quad
n_i=\vartheta x_i/(\chi+\kappa p),\quad
h_i^y=Ax_i/p,\quad z_i=v_i.
\]
Old marginal utility of resources is \(m_i=K/v_i\). The young cash multiplier is
\(\mu_i=1/x_i-\beta K/(qv_i)>0\).

For a renter, the envelope of remaining lifetime utility is
\[
dW^y_{Ri}=\frac{dT_0-h_i^y\,dp_0}{x_i}
 +\beta m_i\,dT_1-\beta\gamma\frac{dp_1}{p}.
\tag{2}
\]
For an owner, its difference from (2) is
\[
dW^y_{Oi}-dW^y_{Ri}
=q\mu_i h_i^y(dP_0-dP_1).
\tag{3}
\]
Indeed, the owner's current cash price is \(L_t=(1-q+q\tau_t)P_t\), while \(L_t-p_t=q(P_{t+1}-P_t)\). Its full housing envelope is
\(-h_i^y[\beta m_i\,dp_0/q+\mu_i\,dL_0]\), which gives (3). This term retains the binding-mortgage channel.

An initial old household instead has
\[
dW^o_i=m_i(dT_0+\mathbf1_O h_i^{y,\mathrm{inherited}}dP_0)
 -\gamma\frac{dp_0}{p}.
\tag{4}
\]
The inherited title is its own preceding young house, and its inherited face debt remains fixed. The estate utility term is already incorporated in \(m_i=K/v_i\); it has not been frozen or omitted.

Conditional tenure values and real choices coincide at the reference, so the ownership probability \(\pi\) is independent of endowments. For young welfare the expected-maximum envelope, including the ownership taste, weights (2)–(3) by their reference choice probabilities. At an indifferent switching household the two values including taste coincide, so the switching boundary contributes zero at first order. This justifies the aggregate formula without freezing tenure probabilities off the reference.

## 2. Primitive sufficient condition and positive ownership

Let
\[
p=\frac{\nu\vartheta W/E-\chi}{\kappa},\quad
a=\frac{\alpha W}{Ep},\quad b=\frac{\gamma V}{Kp},\quad
c=\frac\kappa\nu,\quad h_y=a+c,\quad d=a+b+\eta c,
\]
\[
t=T_\tau=\frac{qp(a+b+c)}{2(1-q)},\qquad
L_0=\frac{h_y/W+b/V}{d},\qquad
L_1=L_0+\frac{h_y}{d}\left(\frac1W-\eta L_0\right).
\tag{5}
\]
Here \(L_0t\) and \(L_1t\) are the respective impact and next-period proportional rent derivatives at the zero-ownership limit. These quantities follow from housing clearing and demographic evolution, with both initial cohorts fixed. At zero initial tax, the permanent policy has \(T'_0=T'_1=t\).

Summing (2) and (4) at zero ownership gives the exact expression
\[
\frac{W'_0}{Nt}
=\int\frac E{w_i}\,dF+(1+\beta)\int\frac K{v_i}\,dF
 -(A+\gamma)L_0-\beta\gamma L_1.
\tag{6}
\]
Consequently the explicit mean-income test
\[
\boxed{J:=\frac EW+(1+\beta)\frac KV
 -(A+\gamma)L_0-\beta\gamma L_1>0}
\tag{7}
\]
implies \(W'_0/N\ge tJ>0\) by Jensen's inequality. Every object in (7) is given by primitive means and preferences through (5); no equilibrium multiplier or assumed marginal-utility ordering appears in this test.

A simpler conservative income-ratio certificate is available. Define \(\rho=EV/(KW)\), \(B=\gamma E/K\), and \(Q=\alpha+\vartheta\eta^2\). Then
\[
WL_0=\frac{A+B}{Q+\gamma\rho},\qquad
WL_1\le\frac{2A+B}{Q+\gamma\rho}.
\]
The denominator is at least \(\gamma\rho\), so
\[
WJ\ge E+\frac{E(1+\beta)
 -[(A+\gamma)(A+B)+\beta\gamma(2A+B)]/\gamma}{\rho}.
\]
It therefore suffices that
\[
\boxed{\rho>\max\left\{0,
\frac{(A+\gamma)(A+B)+\beta\gamma(2A+B)}{E\gamma}
 -(1+\beta)\right\}.}
\tag{8}
\]
Replacing \(A\) on the right by \(\alpha+\vartheta\) gives a more conservative restriction using only preference weights and the old/current income ratio. This says that future income is sufficiently large relative to current purchasing resources, with a threshold that accounts for the housing and continuation-utility costs of the policy. It is sufficient, not necessary.

For completeness the extension to positive ownership is quantified. Take the explicit \(C_0,C_K,R_\tau,\pi_{\rm inv}\) from equations (12), (13), and (16) of the transition memo, and define
\[
M=2C_0^2C_KR_\tau,\quad U=2C_0R_\tau,\quad
m=\frac KE\int\frac{w_i}{v_i}\,dF<q/\beta.
\]
That memo proves, for \(0<\pi<\pi_{\rm inv}\),
\[
\sup_t|P'_{t,\pi}/p-P'_{t,0}/p|\le M\pi,
\quad\sup_t|P'_{t,\pi}/p|\le U,
\quad\sup_t|p'_{t,\pi}/p-p'_{t,0}/p|\le(1+q)M\pi.
\]
The owner coefficients in (3)–(4), multiplied by \(p\), are \(A(q-\beta m)>0\) and \(Am>0\). Hence
\[
\left|\frac{W'_\pi-W'_0}{N}\right|\le\pi L_W,\quad
L_W=[A+(1+\beta)\gamma](1+q)M
 +A[2(q-\beta m)+m]U.
\tag{9}
\]
Thus (1), (7), the finite cap margins, and
\[
\boxed{0<\pi<\min\{\pi_{\rm inv},\ tJ/(2L_W)\}}
\tag{10}
\]
give \(W'_\pi/N>tJ/2>0\). A finite logistic taste location generates any such strictly positive share. The nonlinear local equilibrium theorem and differentiability then imply a welfare gain for sufficiently small positive permanent taxes. This is a general analytical sufficient theorem, not a positive numerical point followed by an unquantified continuity claim.

## 3. Exact heterogeneous compatibility family

Set
\[
q=\beta=1/2,\quad \alpha=\gamma=\vartheta=\kappa=\nu=\chi=1,
\quad\omega_B=2,\quad W=6,\quad v_i=6w_i.
\tag{11}
\]
Any nondegenerate distribution of \(w_i\in[3,9]\) with mean six is allowed. For example take finite caps \(h_R^{\max}=14<h_O^{\max}=15\). Then \(E=3,K=4,p=1,P=2,N=\bar H/12\), and
\[
x_i=w_i/3,\quad n_i=w_i/6,\quad
c_i^y=h_i^y=w_i/2,\quad c_i^o=h_i^o=3w_i/2.
\]
All cap margins are strict. Old estates are \(6w_i\), with positive financial estates \(3w_i\). The young finance multiplier is \(7/(3w_i)>0\). Thus this is an exact analytical equilibrium family with heterogeneity, positive estates, and genuinely binding young finance.

Its transition constants are
\[
a=2,\ b=9,\ c=1,\ \eta=1/2,\ d=23/2,\ t=6,
\quad \ell_0=9/23,\quad y_1=37/46,\quad \ell_1=318/529.
\tag{12}
\]
The zero-ownership demographic recurrence is
\(y_{t+1}=(20/23)y_t-(9/23)y_{t-1}+37/46\).
Discounted summation of this recurrence and the bounded asset-pricing equation gives
\[
P'_0=-26/61,\qquad P'_1=512/1403.
\tag{13}
\]
These are exact fractions. If \(B_w=6\int w_i^{-1}\,dF\ge1\), the separate cohort derivatives are
\[
\frac{W'_{o,0}}N=\frac23B_w-\frac9{23}\ge\frac{19}{69}>0,
\]
\[
\frac{W'_{y,0}}N=\frac{10}{3}B_w-\frac{27}{46}-\frac{159}{529}
\ge\frac{7763}{3174}>0,
\]
\[
\frac{W'_0}N=4B_w-\frac{1353}{1058}\ge\frac{2879}{1058}>0.
\tag{14}
\]
The conservative price-free certificate (8), with \(A\) replaced by \(2\), requires \(\rho>49/24\); here \(\rho=9/2\). Thus compatibility does not depend on evaluating the sign only at the fractions in (14).

An explicit common ownership interval preserves both cohort signs. For this family, set
\[
C_y=\frac{3/2}{(1-3/\sqrt{23})^2},\quad
C_0=2+\frac{48}{23}C_y,\quad C_K=1/2,\quad R_\tau=14/23,
\]
\[
M=14C_0^2/23,\quad U=28C_0/23,\quad
L_y=3M+7U/6,\quad L_o=3M/2+U/3.
\]
Choose
\[
0<\pi<\min\left\{\frac12,\frac1{C_0},
\frac{7763}{6348L_y},\frac{19}{138L_o},
\frac{13}{61M},\frac{37}{92M}\right\}.
\tag{15}
\]
Equations (9) and the separate envelopes preserve at least half each positive cohort margin in (14). The final two bounds preserve \(P'_0<-13/61<0\) and \(y_1>37/92>0\). All constants are elementary functions of specified primitives. The resulting bound is conservative and small; it is not a claim about ordinary empirical ownership shares.

This family also satisfies the separate dated planner comparisons with caps slack. At the market allocation, dated consumption and housing per young-old pair are both 12. Holding each fertility choice fixed, the full consumption-and-housing planner gives each age adult consumption and adult space \(11/2\), so aggregate young housing rises from 3 to \(13/2\). Allowing fertility jointly gives \(x=c_o=s=h_o=24/5\), \(n=12/5\), and young housing \(36/5\). These choices fit the stated finite caps. The tax transition has a larger limiting population because the transition memo's sign expression is strictly positive when \(V>W\); it also has higher initial fertility under (15). Limiting fertility remains replacement fertility, so this is not a permanent positive fertility-level gap.

## 4. Sharper family result from the retained-owner boundary solution

The transition agent separately establishes existence and initial-boundary invertibility for every \(\pi\in[0,1]\) in family (11), using the characteristic cubic rather than the Neumann bound. The following welfare algebra uses its two inputs: a unique bounded derivative and an inverse unstable root \(s\in(0,1/2]\). Verification of those inputs belongs to the transition memo's all-share extension; the welfare inequalities below are independently derived here. They remove (15) for this exact family once that extension is combined with this note.

Write \(F=P'_0\), \(D=(207-41\pi)/36\), \(A_0=(414-68\pi)/36\), and \(n=(1-\pi)/4\). The boundary formula is
\[
F=\frac{-7s+\tfrac12s^2(3+9s)/(1-s)}
{(1-s)D+ns(3+9s)}<0,\qquad
P'_1=\frac{A_0F+7}{D}.
\tag{16}
\]
It also implies the uniform bound \(F>-3/4\). To see this, the denominator is at least \((1-s)83/18\). Since the numerator is negative, it suffices to show
\[
83-334s+287s^2+108s^3>0\quad(0\le s\le1/2).
\]
Putting \(u=1/2-s\) turns the expression into
\(5/4-34u+449u^2-108u^3\ge5/4-34u+395u^2>0\); the final quadratic has negative discriminant and positive leading coefficient.

The initial housing and demographic rows give
\[
\ell_0=F-P'_1/2+1,\quad
y_1=1-\ell_0/2-\pi(P'_1-F)/4,
\]
\[
\ell_1=\frac{9/2+3y_1-(41/18)\pi-(7/18)\pi P'_1-(3/4)\pi F}{2D}.
\tag{17}
\]
Substitution into the separate envelopes gives, at \(B_w=1\),
\[
\frac{W'_o}{N}=
\frac{171+41\pi+\pi(228-41\pi)F}{3(207-41\pi)}
\ge\frac{287}{2484}>0.
\tag{18}
\]
For the bound, use \(F>-3/4\): the numerator is at least
\(171-130\pi+(123/4)\pi^2\ge287/4\).
Young welfare is affine in \(F\), with coefficient
\[
-\frac{\pi(2583\pi^2-28334\pi+74691)}{4(207-41\pi)^2}\le0.
\]
Since \(F<0\), dropping this nonnegative contribution yields
\[
\frac{W'_y}{N}\ge
\frac{49610\pi^2-361539\pi+628803}{6(207-41\pi)^2}
\ge\frac{158437}{128547}>0.
\tag{19}
\]
The numerator decreases on \([0,1]\) and is at least \(316874\); the denominator is at most \(6\cdot207^2\). Heterogeneity adds exactly \((2/3)(B_w-1)\) to (18) and \((10/3)(B_w-1)\) to (19). Thus both cohort averages rise for every interior ownership share in the analytically solved family, including shares bounded away from zero. This is a uniform symbolic inequality, not a numerical-root calculation.

## 5. Incidence and scope

Under the explicit external residual-owner convention, let \(H^{rent}_{-1}\) denote the **entire inherited rental stock**, including rentals occupied by both ages. Its surprise revaluation is \(H^{rent}_{-1}P'_0\). With the stationary convention \(H^{rent}_{-1}=(1-\pi)\bar H=12N(1-\pi)\), (15), or the all-share result (16), implies a strict outside capital loss at every covered interior share. At the zero-ownership limit it is exactly \(-312N/61\). Using only the preceding young renters' homes would undercount this account. Existing creditor face claims are not written down, and this monetary loss cannot be added to household utility without specifying the outside claimants' preferences and welfare weights.

The fiscal ledger remains balanced: at the reference \(2Nt=qP\bar H\). The result is a market-feasible living-household utilitarian improvement, not a Pareto improvement over every household or outside financier. Positive cohort aggregates do not imply every incumbent owner gains. The all-renter limit already gains, so the theorem includes binding mortgages without establishing that the mortgage channel is necessary or dominant. Outside the analytically parameterized family, no general welfare sign at large ownership shares, with binding size caps, or with binding old estate floors is proved.

The strict result extends to sufficiently nearby ongoing preference-decline transitions by the local equilibrium theorem, comparing the same inherited state and the same subsequent preference path with and without tax. The stationary envelope pairing must not be reused literally at a nonstationary intervention date.

**One precise follow-up for independent review:** can the primitive mean-income condition (7), together with (1), deliver a useful non-small ownership interval by solving the retained-owner derivative system analytically and bounding (2)–(4), while keeping the entire outside rental-stock loss explicit?
