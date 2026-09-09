# Independent binding-finance transition theorem

September 9, 2026. Analytical diagnostic in the exact two-age model; no model
runs, numerical roots, Pro output, or household-equation changes. This report
is separate from the earlier slack-finance transition memo. All claims below
are local. The ownership interval is explicit and may be conservative.

## 1. Regime and stationary reference

Write \(w_i=y_i^y+b_i\), \(v_i=y_i^o\), \(W=\mathbb E w_i\),
\(V=\mathbb E v_i\), and
\[
E=1+\alpha+\vartheta,\qquad K=1+\gamma+\omega_B.
\]
Use compact, genuinely heterogeneous positive endowment support, \(\phi=q\),
and zero reference tax. Impose uniform strict inequalities
\[
\frac{v_i}{w_i}>\frac{\beta K}{qE},\qquad
\omega_B(1-q)>q\gamma,\qquad
\chi<\frac{\nu\vartheta W}{E}.                                      \tag{1}
\]
Take finite tenure caps strictly above both conditional young and old housing
demands. The renter cap can remain strictly smaller than the owner cap.
The first inequality makes every young household strictly cash constrained
in both tenures. The second makes old owners' financial-estate floor slack.

At the stationary reference, both tenures have identical real choices:
\[
x_i=\frac{w_i}{E},\quad
n_i=\frac{\vartheta w_i}{E(\chi+\kappa p)},\quad
h_i^y=\frac{\alpha w_i}{Ep}+\kappa n_i,
\quad c_i^o=\frac{v_i}{K},\quad h_i^o=\frac{\gamma v_i}{Kp}.
\tag{2}
\]
Young renters save zero. Young owners borrow the maximum principal
\(qPh_i^y\) and owe \(Ph_i^y\) on entering old age. The old house sale
exactly offsets this face repayment at stationary prices; old resources are
\(v_i\). This offset does not hold along the transition.

Conditional tenure values coincide, so the stationary ownership probability
is \(\pi=\operatorname{logistic}(\bar\xi/\sigma_\xi)\), independent of type.
Any strictly positive share below the explicit bound below is attainable
with a finite taste location and positive logistic scale.

Replacement fertility and market clearing determine
\[
p=\frac{\nu\vartheta W/E-\chi}{\kappa}>0,\quad
P=\frac{p}{1-q},\quad
N=\frac{\bar H}{\kappa/\nu+(\alpha W/E+\gamma V/K)/p}.
\tag{3}
\]
Neither aggregate external wealth nor population is independently fixed.

## 2. Exact transition accounting

Let \(p_t=(1+q\tau_t)P_t-qP_{t+1}\),
\(L_t=(1-q+q\tau_t)P_t\), and \(\Delta P_t=P_{t+1}-P_t\).
Thus \(L_t-p_t=q\Delta P_t\), even when the tax is permanent.
In the maintained binding regime, the owner solves
\[
c+L_th=w_i+T_t,\quad
z=v_i+T_{t+1}+\Delta P_t h,
\]
\[
\frac{\alpha}{s}=\frac{L_t}{x}
                  -\frac{\beta K\Delta P_t}{z},\qquad
\frac{\vartheta_t}{n}=\frac\chi x+\frac{\alpha\kappa}s.
\tag{4}
\]
The renter has \(c+p_th=w_i+T_t\) and \(z=v_i+T_{t+1}\).
These exact conditional maps include future rebates and prices. Strict
concavity and the uniform inequalities in (1) make them smooth locally.
The endogenous logistic tenure probabilities must also be recomputed.

Put \(A^O_t=\mathbb E[\pi^O_t(i)h^O_t(i)]\). For \(t\ge1\), mean old
resources and housing are
\[
Z_t=V+T_t+(P_t-P_{t-1})A^O_{t-1},\qquad
\bar h^o_t=\frac{\gamma Z_t}{Kp_t}.                                \tag{5}
\]
At date zero, instead use the actual inherited financial claims and titles:
\[
Z_0=\bar a_0+P_0\bar H_0+V+T_0.
\tag{6}
\]
At the reference \(\bar H_0=\pi\bar h^y\) and
\(\bar a_0=-P\bar H_0\); these inherited face claims stay fixed when policy
is announced. Complete the equations with demographic evolution, housing
clearing, and \(T_t=q\tau_tP_t\bar H/(Y_t+O_t)\).

## 3. Linearization retaining positive ownership

Define positive housing quantities and a cost share
\[
X=W/E,\quad a=\alpha X/p,\quad b=\gamma V/(Kp),\quad c=\kappa/\nu,
\quad h_y=a+c,\quad H_*=a+b+c,
\]
\[
\eta=\frac{\kappa p}{\chi+\kappa p}\in(0,1),\quad
d=a+b+\eta c,\quad j=\vartheta/E,\quad A_h=\alpha+\vartheta\eta.
\]
Heterogeneity enters the first-order owner response through the moment
\[
B_f=\frac{\beta K}{E}\frac{\mathbb E(w_i^2/v_i)}W\in(0,q).
\]
Set
\[
d_n=\frac{A_h B_f}{E}+\eta(q-B_f),\quad
d_h=\frac{h_y A_h B_f}{E}+(a+\eta c)(q-B_f),\quad
d_o=\frac{\gamma h_y}{K}.                                         \tag{7}
\]
All three are strictly positive. For example, the owner's effective adult
space price \(\rho_i=\alpha x_i/s_i\) satisfies
\(d\rho_i=dp+(q-\beta K w_i/(Ev_i))d\Delta P\). This gives (7) after
integrating. Derivatives of tenure shares cancel from real aggregates at
the reference because conditional quantities coincide. The cancellation
does not freeze tenure probabilities away from the reference.

Use normalized variations \(y_t=dY_t/N\), \(F_t=dP_t/p\),
\(\ell_t=dp_t/p\), and \(\zeta_t=d\vartheta_t/\vartheta\). Let
\(\xi_t=dT_t/W-j\zeta_t\), \(\xi_t^o=dT_t/V\).
With initial populations fixed, set \(y_{-1}=y_0=0\) and \(F_{-1}=0\).
The last convention represents fixed initial titles and face debt through
(6), not an endogenous pre-policy price. The exact derivative system is
\[
y_{t+1}-y_t+\eta\ell_t+\pi d_n(F_{t+1}-F_t)=\xi_t+\zeta_t,
\]
\[
d\ell_t-h_yy_t-by_{t-1}
 +\pi d_h(F_{t+1}-F_t)-\pi d_o(F_t-F_{t-1})
 =h_y\xi_t+c\zeta_t+b\xi_t^o,
\tag{8}
\]
\[
\ell_t=F_t-qF_{t+1}+\frac q{1-q}d\tau_t,\qquad
dT_t=T' d\tau_t,\quad T'=\frac{qpH_*}{2(1-q)}>0.
\tag{9}
\]
In particular, the date-zero housing row contains the initial-owner capital
gain or loss \(\pi d_oF_0\). It is not discarded.

For reference, homogeneous modes obey the explicit cubic
\[
(1-qr)\{rd(r-1)+\eta(h_yr+b)\}
 +\pi(r-1)\{d_n(h_yr+b)+(r-1)(d_hr-d_o)\}=0.                     \tag{10}
\]
The following operator argument avoids numerical roots or an unquantified
appeal to continuity in the owner share.

## 4. Explicit inverse and nonlinear local path

Eliminate the fiscal and price rows, divide the housing row by \(d\), and
let \(S\) be the forward shift. The reduced operator on
\(u=(y_{t+1},F_t)_{t\ge0}\) is \(\mathcal L_\pi=\mathcal L_0+\pi\mathcal K\):
\[
\begin{split}
y_{t+1}-y_t+\eta(I-qS)F_t+\pi d_n\Delta^+F_t&=f_t,\\
(I-qS)F_t-(h_y/d)y_t-(b/d)y_{t-1}
 +\pi(d_h/d)\Delta^+F_t-\pi(d_o/d)\Delta^-F_t&=g_t.
\end{split}                                                       \tag{11}
\]
Use the maximum of the two supremum norms. At \(\pi=0\), elimination gives
\[
y_{t+1}=\mathcal A y_t-ky_{t-1}+f_t-\eta g_t,\quad
\mathcal A=1-\eta h_y/d,\quad k=\eta b/d.
\]
Here \(0<\mathcal A<1\), \(0<k<1\), and the three strict quadratic Jury
inequalities hold. If
\[
r_*=
\begin{cases}(\mathcal A+\sqrt{\mathcal A^2-4k})/2,&\mathcal A^2\ge4k,\\
\sqrt{k},&\mathcal A^2<4k,
\end{cases}
\]
then \(r_*<1\). The recurrence impulse coefficients satisfy
\(|\psi_m|\le(m+1)r_*^m\), including repeated roots. Consequently, with
\[
C_y=\frac{1+\eta}{(1-r_*)^2},\quad
C_F=\frac{1+(H_*/d)C_y}{1-q},\quad C_0=\max\{C_y,C_F\},
\quad C_K=\max\{2d_n,2(d_h+d_o)/d\},                              \tag{12}
\]
we have \(\|\mathcal L_0^{-1}\|\le C_0\) and
\(\|\mathcal K\|\le C_K\). To recover the asset price use the unique
bounded inverse \((I-qS)^{-1}=\sum_{m\ge0}q^mS^m\).
Every constant in (12) is an explicit finite function of primitives and
the stationary formulas (3).

If
\[
0<\pi<\pi_{\rm inv}:=\min\{1/2,(2C_0C_K)^{-1}\},               \tag{13}
\]
the Neumann series gives
\(\|\mathcal L_\pi^{-1}\|\le2C_0\). This argument works on both the
Banach space \(c\) of convergent sequences and \(\ell^\infty\).
Stable convolution and discounted forward summation preserve convergence;
the exceptional initial row is a bounded finite-rank operation.

For the nonlinear system, take four unknown sequences
\((Y_{t+1},p_t,P_t,T_t)\) and four residual sequences: demography, housing,
asset pricing and fiscal balance. Equations (4)–(6) give its conditional
choices and initial row. Compact support and uniform regime margins make
this residual map continuously differentiable on open neighborhoods in
both \(c^4\) and \((\ell^\infty)^4\). At zero tax, the fiscal derivative
first pins \(dT\); the price row and its arbitrary residual can then be
eliminated. The remaining derivative is (11), including arbitrary bounded
residual sequences. Hence (13) proves full derivative invertibility, not
only a homogeneous stability test.

The Banach implicit-function theorem on \(c\) gives an exact convergent
path for sufficiently small convergent tax/taste perturbations and nearby
admissible inherited states. Applied on \(\ell^\infty\), it gives local
uniqueness among nearby bounded paths. These are the same solution. All
young constraints remain strictly binding, and old floors and tenure caps
remain slack. A surprise intervention followed by perfect foresight is
therefore covered. No global uniqueness, distant-path exclusion, or large
policy claim follows.

## 5. Signed population and initial-fertility effects

At a stationary tax, \(L=p\) again. Define
\[
A_0=\frac{(\alpha+\vartheta)W}{E}+\frac{\gamma V}{K}-\frac\chi\nu,
\quad B_0=\frac{\alpha+\vartheta}{E}+\frac\gamma K,
\quad g_\tau=\frac{q\tau}{2(1-q+q\tau)}.
\]
The exact endpoint is independent of the stationary owner share:
\[
T(\tau)=\frac{g_\tau A_0}{1-g_\tau B_0},\quad
p(T)=\frac{\nu\vartheta(W+T)/E-\chi}{\kappa},\quad
N(\tau)=\frac{\bar H p(T)}{A_0+B_0T}.                             \tag{14}
\]
The local sign is
\[
S_{\rm tax}:=\frac{\gamma\nu\vartheta(V-W)}{EK}
       +\chi\left(\frac\alpha E+\frac\gamma K\right),\qquad
N_\tau(0)=\frac{\bar H S_{\rm tax}}{\kappa A_0^2}T'.              \tag{15}
\]
Thus \(V\ge W\) is an easily read sufficient condition for a larger
limiting population; \(S_{\rm tax}>0\) is the exact sign condition in
this branch. The limiting fertility of both paths is still \(1/\nu\).

There is also an explicit positive-ownership interval preserving the
initial signs. For a permanent unit tax perturbation, define
\[
m_\tau=\frac{T'}d\left[\frac{a(1-\eta)+b}{W}-\frac{\eta b}{V}\right],
\]
\[
R_\tau=\max\left\{\left|\frac{T'}W-\frac{\eta q}{1-q}\right|,
\left|\frac{T'}d\left(\frac{h_y}W+\frac bV\right)-\frac q{1-q}\right|\right\}.
\tag{16}
\]
At zero ownership, \(y_1/ d\tau=m_\tau\), which is positive exactly when
\(S_{\rm tax}>0\). For a permanent normalized taste change
\(\zeta=d\vartheta/\vartheta\), define
\[
m_\vartheta=\frac{(a+b)(1-j)+\eta aj}{d}>0,\quad
R_\vartheta=\max\{1-j,|c(1-j)-aj|/d\}.                           \tag{17}
\]
The resolvent identity gives, for either normalized forcing \(r\),
\[
\|u_\pi-u_0\|\le2\pi C_0^2 C_K\|r\|.
\]
Therefore choose the entirely explicit bound
\[
\boxed{\quad
0<\pi<\pi_*:=\min\left\{\pi_{\rm inv},
\frac{m_\tau}{4C_0^2C_KR_\tau},
\frac{m_\vartheta}{4C_0^2C_KR_\vartheta}\right\}.\quad}          \tag{18}
\]
Under \(S_{\rm tax}>0\), this interval has positive length. Within it,
a small permanent rebated tax raises initial average fertility, and a small
permanent decline in the fertility taste lowers initial average fertility.
The nonlinear signs follow from strict first-order signs. Since \(Y_0\)
is fixed, the sign of \(y_1\) is the sign of initial fertility.

At zero tax, a larger taste parameter raises the endpoint population:
\(p_\vartheta=\nu W(1+\alpha)/(\kappa E^2)>0\), while
\(\alpha W/E+\gamma V/K\) falls, so housing per stationary household
pair falls. Thus a small taste decline lowers the endpoint. Nearby ongoing
convergent taste-decline paths are also parameters of the theorem. At a
policy date, compare two solutions from exactly the same inherited state
and same subsequent taste sequence. A small permanent tax yields higher
initial fertility locally and a larger limiting population. It need not
restore the pre-decline population or increase fertility at every date.

## 6. Compatibility and economic limits

The simpler shared region for the uncapped planner comparisons is
\[
\frac{v_i}{w_i}>\max\left\{1,\frac KE,\frac{\beta K}{qE}\right\}.
\tag{19a}
\]
It gives strict young finance, \(V>W\), \(c_i^o>x_i\), and
\(\alpha h_i^o>\gamma s_i\). These are the consumption and adult-space
inequalities used in the separate fixed-fertility and joint-fertility
planner argument. The transition proof does not require old total housing
to exceed young total housing. If that stronger matched-type comparison
is desired, a sufficient pointwise restriction is
\[
\frac{v_i}{w_i}>K\max\left\{\frac{\alpha+\vartheta}{E\gamma},
                             \frac\beta q,\frac1K\right\}
\tag{19b}
\]
which implies old housing exceeds young housing by matched type, strict
young finance, and \(V>W\). A stronger mean-income restriction used in a separate
joint-fertility proposition can also be imposed alongside the pointwise
finance bound. The set is nonempty with heterogeneous compact support.
An explicit finite-cap requirement is
\[
h_R^{\max}>\max\left\{\frac{w_{\max}}W h_y,
                              \frac{\gamma v_{\max}}{Kp}\right\},
\qquad h_O^{\max}>h_R^{\max}.
\tag{20}
\]
Strict gaps preserve cap slackness for the local path. This establishes joint
logical compatibility, not empirical plausibility. For fixed incomes,
the allowed \(\beta\) is bounded by (1); required future-income thresholds
grow with \(\beta\). No separate assumption \(\beta R_f<1\) or
\(\beta R_f>1\) is made.

Unlike the slack-finance memo, this result has positively binding mortgages
for every owner and a positive owner mass. The \(\Delta P\) terms in (8)
retain that finance channel. The theorem does not identify this channel as
the sole or necessary cause of the tax effect: the limiting renter economy
already has the same strict sign. The explicit owner interval can be small.

The physical housing stock and domestic household budget equations are
respected. Outside rental financiers bear the intermediary capital losses
under the maintained external-financing convention; this report neither
adds domestic landlord wealth nor treats tax rebates as free resources.
Initial household creditor claims are fixed; estate proceeds and subsequent
financial positions are endogenous. A fixed-debt, zero-equity intermediary
interpretation would need an explicit loss-bearing account. No current-
household welfare theorem follows from the population result.

For the two common-state paths,
\[
\frac{Y_T^{\rm pol}}{Y_T^{\rm base}}
=\prod_{t=0}^{T-1}\frac{\bar n_t^{\rm pol}}{\bar n_t^{\rm base}}.
\]
The theorem's finite positive limits make the limiting cumulative log gap
equal to \(\log(N^{\rm pol}/N^{\rm base})>0\). Total adult-household
population converges to \(2N\); resident-person population requires a
separate counting definition.

One precise further review target is independent verification of (7)–(13)
and the nonlinear sequence-space construction. The general model with
large owner shares or binding size/old-finance restrictions remains open.
