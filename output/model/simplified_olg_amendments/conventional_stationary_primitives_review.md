# Primitive stationary existence and eligible-owner mass

Independent analytic review, September 8, 2026. This concerns exactly the
conventional-finance discussion proposal in
`docs/model/simplified_olg_conventional_finance_proposal.md`; it changes no model,
manuscript, target, or calibration. The proof is analytic and uses no numerical
equilibrium construction. The final section gives an explicit nonempty family,
including any fixed positive discount factor.

## 1. Primitives and the statement

Write \(\theta=\vartheta\), \(\omega=\omega_B\), and \(\tau=\tau^p\). Assume
\(0<q<1\), \(0\le\phi\le1\), \(0\le\tau<2\), and positive
\(\alpha,\theta,\beta,\gamma,\omega,\chi,\kappa,\nu,\bar H\).
The stationary entering-type distribution \(F\) has bounded support and
\[
w_0=y^y+b\ge\underline w>0,\qquad v_0=y^o\ge0.
\]
The distribution is exogenous across successive entering cohorts; no estate
feedback changes it. The ownership taste has a logistic distribution with
strictly positive, finite scale and finite location, independently of type
(the same proof permits a continuous conditional distribution with uniformly
positive finite scale on the compact type support). Tenure is chosen once and
retained when old. Owner and renter physical caps are finite and satisfy
\(h_O^{\max}>h_R^{\max}>0\). All log arguments are strictly positive.

Define the stationary housing service cost, owner cash cost, and preference sums:
\[
d_p=1-q+q\tau>0,\quad d_L=1-\phi+q\tau\ge0,\quad
p=d_pP,\quad L=d_LP,
\]
\[
K=1+\gamma+\omega,\quad B=\beta K,\quad
D=1+\alpha+\theta+B,\quad a=\max\{d_p,d_L\}.
\]
Let
\[
\bar M_0=\int(w_0+qv_0)\,dF,\qquad
\bar T=\frac{\tau\bar M_0}{(1-q)(2-\tau)},\qquad
\bar M=\bar M_0+(1+q)\bar T.
\]
The following entirely primitive inequalities are sufficient for a positive
stationary equilibrium:
\[
h_R^{\max}>\frac\kappa\nu,\qquad
\theta\nu>
\frac{\chi D}{\underline w}
+\frac{\alpha\kappa}{h_R^{\max}-\kappa/\nu}.
\tag{E}
\]
They are sufficient restrictions, not necessary conditions or empirical claims.
Define two positive price bounds:
\[
P_- =\frac{\alpha\underline w}{D a h_R^{\max}},\qquad
P_+ =\frac{\nu\bar M}{d_p\kappa}.
\tag{P}
\]
Every positive stationary equilibrium satisfies
\(P_-<P<P_+\) and \(0\le T\le\bar T\).
Existence does not require a fiscal contraction, uniqueness, an assumed smooth
active set, or a numerically checked price rectangle.

## 2. Full household problem and uniform adult-consumption bound

At fixed \((P,T)\), write \(w=w_0+T\), \(v=v_0+T\),
\(M=w+qv\), \(x=c^y-\chi n\), and \(s=h^y-\kappa n\).
Combining both dated household budgets gives, for either tenure,
\[
x+ps+(\chi+p\kappa)n+qc^o+qp h^o+q^2e=M.
\tag{B}
\]
The young owner's mortgage restriction is
\[
x+Ls+(\chi+L\kappa)n\le w.
\]
For a renter replace \(L\) by \(p\); this is exactly the no-borrowing
restriction \(a'\ge0\). Young housing satisfies
\(s+\kappa n\le h_m^{\max}\), and old housing satisfies
\(h^o\le h_m^{\max}\). The owner's old-age restriction is the homogeneous
estate floor \(Ph^o-e\le0\). These constraints retain both physical caps and
the original estate floor.

Let \(\lambda>0\) multiply (B), \(\mu\ge0\) multiply current cash, and
\(\eta_y,\eta_o\ge0\) multiply the two housing caps. Scaling all six
positive choice coordinates and using complementary slackness gives
\[
D=\lambda M+\mu w+\eta_y h_m^{\max}+\eta_o h_m^{\max}.
\]
The estate multiplier contributes zero because its restriction is homogeneous.
The adult-consumption first-order condition is \(1/x=\lambda+\mu\).
Since \(M\ge w\),
\[
x\ge\frac wD\ge\frac{\underline w}{D}.
\tag{X}
\]
This bound remains valid when either physical cap or the estate floor binds.

Conditional household policies exist, are unique in real allocations, and are
continuous in \((P,T,\text{type})\) on compact positive-price rectangles.
For completeness: the closed resource set is compact; a strictly positive
feasible allocation is obtained by choosing all current real uses small and
retaining positive old-age resources. The log objective rules out zero log
arguments. Its strict concavity gives uniqueness. Feasible-set continuity,
including at changing active sets, follows from the linear restrictions and a
strict feasible allocation; the maximum theorem then gives continuity. The
owner's gross mortgage/bond portfolio need not be unique, which does not affect
these real policies. Logistic integration makes aggregate tenure-weighted
policies continuous even if \(F\) has atoms.

## 3. Low and high price fertility signs

If the young housing cap is slack, its first-order condition gives
\[
\frac{\alpha x}{s}
=\frac{\lambda p+\mu L_m}{\lambda+\mu}
\le aP,
\qquad L_O=L,\quad L_R=p.
\]
If a young cap binds, housing is at least \(h_R^{\max}\). Therefore, in
either tenure,
\[
h^y\ge\min\left\{h_R^{\max},
\frac{\alpha\underline w}{DaP}\right\}.
\]
For \(P\le P_-\), young housing is at least \(h_R^{\max}\). Were
\(n\le1/\nu\), the fertility first-order condition and (X) would imply
\[
\theta\nu\le\frac\theta n
=\frac\chi x+\frac{\alpha\kappa}{s}
\le\frac{\chi D}{\underline w}
+\frac{\alpha\kappa}{h_R^{\max}-\kappa/\nu},
\]
contrary to (E). Thus every conditional type/tenure policy has
\(n>1/\nu\), uniformly over \(T\in[0,\bar T]\) and \(P\le P_-\).
Consequently aggregate fertility exceeds replacement regardless of the
logistic tenure mix.

At high prices, (B) and positive adult/old uses imply
\[
n<\frac{M}{\chi+d_pP\kappa},\qquad
\bar n(P,T)<\frac{\bar M_0+(1+q)T}{\chi+d_pP\kappa}.
\]
Hence \(\bar n<1/\nu\) whenever \(P\ge P_+\) and
\(T\le\bar T\). The price interval is automatically nonempty:
\[
\frac{P_-}{P_+}
=\frac\alpha D\frac{\underline w}{\bar M}
\frac{d_p}{a}\frac\kappa{\nu h_R^{\max}}<1.
\]
No separate primitive condition \(P_-<P_+\) is needed.

## 4. Fiscal self-map and stationary equilibrium

Let \(\bar h=E[h^y+h^o]\), including the conditional tenure weights.
Since \(q<1\), the household lifetime budget implies
\[
qd_pP\bar h<\bar M_0+(1+q)T.
\]
Define the balanced-rebate map
\[
\mathcal T(P,T)=\frac{q\tau P}{2}\bar h(P,T).
\]
It is nonnegative and satisfies
\[
\mathcal T(P,T)\le
\frac\tau{2d_p}\{\bar M_0+(1+q)T\}\le\bar T
\quad\text{for }T\in[0,\bar T].
\]
The final inequality is precisely
\(2d_p-\tau(1+q)=(1-q)(2-\tau)>0\). At any stationary equilibrium the
same inequality, rearranged, proves \(T\le\bar T\) without first assuming
that the equilibrium lies in the rectangle. For positive tax the resource
inequality is strict, so its upper bound is also strict; for zero tax \(T=0\).

For any fixed \(\eta>0\), the continuous map
\[
(P,T)\longmapsto
\left(\operatorname{clip}_{[P_-,P_+]}
\{P+\eta(\nu\bar n(P,T)-1)\},\ \mathcal T(P,T)\right)
\]
maps \([P_-,P_+]\times[0,\bar T]\) into itself. Brouwer yields a fixed
point. The strict fertility signs exclude either price boundary, so its
interior price satisfies \(\nu\bar n=1\). This also works when \(\tau=0\):
the rectangle then reduces to a compact interval, or one can directly use the
intermediate value theorem.

Set the common young and old cohort mass to
\[
N=\frac{\bar H}{\bar h(P,T)}>0.
\]
Housing clears, the replacement law \(N'=\nu\bar nN\) preserves cohort
mass, and \(2NT=q\tau P\bar H\). The stationary old distribution is the
pushforward of the stationary entering distribution and tenure draw through
the unique young policies. Thus the household, population, housing, and fiscal
conditions all hold. This is an existence result, not a uniqueness or dynamic
stability result.

## 5. Fully primitive conditions for positive eligible-owner mass

The following additions make a positive mass of strictly constrained owners
have slack young and old housing caps and a strictly slack old estate floor at
every stationary equilibrium. They also apply to their comparison allocations
with the young mortgage constraint removed.

Assume \(d_L>0\), set \(\ell=d_L/d_p\), and define
\[
C_{\min}=1+\alpha\ell+\theta\min\{1,\ell\},\qquad
k=\frac{D}{C_{\min}}-1.
\]
Let \(S\) be a specified positive-\(F\)-mass subset of primitive types such
that
\[
qv_0>kw_0+(k-q)_+\bar T\quad\text{on }S,
\tag{I}
\]
where \((z)_+=\max\{z,0\}\). Define a finite bound
\[
M_S\ge\sup_{i\in S}\{w_{0i}+qv_{0i}+(1+q)\bar T\}.
\]
Impose
\[
\omega d_p>q\gamma,
\qquad
h_O^{\max}>
\frac{M_S}{d_pP_-}\max\left\{1,\frac\gamma{qK}\right\}.
\tag{C}
\]
The cap condition is not circular: \(P_-\) depends on the renter cap, not
the owner cap. The bound may be conservative but is explicit and finite.

**Cap and estate verification.** In both actual and mortgage-relaxed
allocations, the lifetime budget gives \(h^y<M/p\) and old resources
\(z<M/q\). The unrestricted old choice under the strict estate inequality in
(C) is
\[
c^o=z/K,\quad h^o=\gamma z/(Kp),\quad e=\omega z/(Kq).
\]
The owner cap in (C) makes this old housing strictly feasible for every
\(z\le M_S/q\). It also makes the young cap strictly slack by the raw
lifetime-budget bound. These arguments apply to the actual and relaxed
problems, so no use of an uncapped formula assumes the conclusion being proved.
A separate young bound \(h^y<w/L\) is optional and unnecessary here.

**Strict finance verification.** With these caps slack, the mortgage-relaxed
current cash demand divided by current resources is
\[
\frac{M}{Dw}
\left[1+\alpha\ell+
\theta\frac{\chi+\ell p\kappa}{\chi+p\kappa}\right].
\]
The square bracket is at least \(C_{\min}\), because its last ratio is a
convex combination of \(1\) and \(\ell\). Condition (I) implies
\(qv>kw\) for every \(T\in[0,\bar T]\), and hence
\(M/w>D/C_{\min}\). Relaxed current cash demand therefore strictly exceeds
\(w\), so the actual mortgage multiplier is strictly positive by strict
concavity. This argument remains valid if \(k<0\).

Both conditional values are finite. The logistic owner probability is
strictly positive for every type; therefore
\(\int_S\pi_O(i)\,dF(i)>0\). The physical eligible-owner mass is this
integral times \(N>0\). Pointwise strict conditions can be restricted to a
positive-mass subset with uniform slack if the separate finite-mass local
Pareto argument needs such margins.

For the independently established fixed-fertility housing direction, add
\[
B d_L\ge\phi-q.
\tag{H}
\]
This is exactly \(BL\ge p-L\). Existence and positive constrained-owner mass
do not themselves prove a welfare comparison; (H) is the separate sufficient
restriction used in the fixed-fertility complete-finance comparison.

## 6. Explicit analytic compatibility, with no upper bound on beta

Here is an explicit family rather than an unspecified neighborhood. It is a
logical nonemptiness construction, not a proposed calibration. Fix arbitrary
\(0<q<1\), \(0\le\tau<2\), \(0\le\phi<1\), and any \(\beta>0\).
Set \(\alpha=\gamma=\kappa=\nu=\underline w=1\),
\(\chi=1/4\), and \(h_R^{\max}=2\). Choose
\[
\omega=1+\frac q{d_p}+\frac{(\phi-q)_+}{\beta d_L},
\qquad K=2+\omega,\quad B=\beta K,\quad\theta=6+B.
\]
Then (H) holds, \(\omega d_p>q\), and \(D=8+2B\). Condition (E) is
\(6+B>3+B/2\), which is strict for every \(\beta>0\).

Compute \(k\) from these primitives and write
\[
c_\tau=\frac\tau{(1-q)(2-\tau)},\qquad
V=1+2k_++4c_\tau(k-q)_+,\qquad
\epsilon=\frac1{2(V+1)}>0.
\]
Give mass \(1-\epsilon\) to a uniform rectangle
\(w_0\in[1,2],\ qv_0\in[0,1]\), and mass \(\epsilon\) to a uniform
rectangle \(S\) with \(w_0\in[1,2],\ qv_0\in[V,V+1]\).
Both rectangles have positive width. If separate \((y^y,b)\) heterogeneity
is desired, draw \(r\) uniformly on \([1/4,1/2]\), independently, and set
\(b=rw_0\), \(y^y=(1-r)w_0\); this produces positive, heterogeneous liquid
wealth and young income without changing \(w_0\).

This bounded, nondegenerate distribution has
\[
\bar M_0\le3+\epsilon V<4,
\qquad \bar T<4c_\tau\quad\text{if }\tau>0.
\]
For zero tax, \(\bar T=0\). Every type in \(S\) satisfies
\[
qv_0\ge V>
2k_++4c_\tau(k-q)_+
\ge kw_0+(k-q)_+\bar T.
\]
Thus (I) holds on an explicitly specified positive-mass rectangle. Take
\[
M_S=V+3+4(1+q)c_\tau,\qquad
h_O^{\max}=3+\frac{M_S}{d_pP_-}
\max\left\{1,\frac1{qK}\right\}.
\]
This finite owner cap exceeds the renter cap and satisfies (C). Any
\(\bar H>0\) and logistic scale \(\sigma_\xi>0\) complete the construction.
All inequalities are primitive, strict where required, and valid for any
fixed positive \(\beta\). The construction establishes that the sufficient
conditions are jointly satisfiable; it does not establish their empirical
plausibility at a particular parameterization.

## Review conclusion

The proposed transfer bound, KKT consumption bound, low/high fertility signs,
clipped-price fixed point, and income-timing restriction are valid. Two useful
simplifications are proved above: the price interval is automatically
nonempty, and the raw lifetime budget eliminates the need for a separate
actual-owner cash-based cap bound. Physical cap assumptions must cover both
actual and finance-relaxed allocations. The strict estate-floor condition is
an explicit primitive restriction. No transition existence, convergence,
unconditional policy sign, or endogenous-fertility welfare theorem follows
from this stationary result.
