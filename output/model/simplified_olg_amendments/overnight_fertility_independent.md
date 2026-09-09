# Independent check: fertility after the dated allocation

September 8, 2026. **Result:** average fertility after private rechoice can
fall even when the fixed-fertility planner raises aggregate young housing
and the joint planner raises fertility. The analytical family below is a
positive stationary competitive equilibrium of the maintained household model.
It retains heterogeneous endowments, positive shares of both tenures, strictly
binding young finance, and slack old estate floors. All caps are finite but
slack in this counterexample. This is a dated parental-welfare comparison,
not a market reform or demographic transition.

## 1. The two experiments

Normalize each age's type law to probability \(Q\), with cohort mass \(N\).
Let \(m=\int n_i\,dQ\) denote reference fertility. For shifted-log utility,
the uncapped fixed-fertility planner assigns
\[
c_i^A=\chi n_i+X,\quad h_i^A=\kappa n_i+S,\qquad
X=\frac{C/N-\chi m}{2},\quad
S=\frac{\alpha}{\alpha+\gamma}(\bar H/N-\kappa m).
\]
Experiment A then fixes these gross bundles and lets each young privately
rechoose fertility, denoted \(n_i^A\). Its subsequent adult resources differ
from \(X,S\).

Experiment B jointly optimizes current goods, housing, and fertility.
With slack joint-planner caps, it chooses common \(n^B\) solving
\[
\frac{\vartheta}{n^B}
=\frac{2\chi}{C/N-\chi n^B}
+\frac{\kappa(\alpha+\gamma)}{\bar H/N-\kappa n^B}.
\]
Both experiments preserve incumbent continuation opportunities and old
net estates. Changing births does not preserve every future population.

## 2. An explicit stationary competitive family

Set \(\alpha=\gamma=\vartheta=1\), \(\beta=\phi=q\in(0,1)\), and
\(\tau^p=0\). Choose \(\chi,\kappa,\nu>0\), and \(u>0\), \(u\ne1\). Define
\[
p=\chi u/\kappa,\quad P=p/(1-q),\quad m=1/\nu,\quad
\bar x=(\chi+p\kappa)m,\quad K=2+\omega_B,
\]
where \(\omega_B(1-q)>q\) ensures positive old financial estates.
There are two equally likely endowment types, indexed by signs:
\[
x_\pm=\bar x(1\pm d),\quad 0<d<1,\qquad
w_\pm=3x_\pm,\qquad y_\pm^o=kKx_\pm.
\]
Split \(w_i\) into positive current income and entrant wealth, for example
\(y_i^y=b_i=w_i/2\). The scalar \(k>1\) will lie in the explicit interval below.

In either tenure, conditional competitive choices are
\[
n_i=m(1\pm d),\quad
c_i^y=x_i+\chi n_i,\quad h_i^y=x_i/p+\kappa n_i,\quad
z_i=kKx_i,
\]
\[
c_i^o=kx_i,\qquad h_i^o=kx_i/p,\qquad e_i=\omega_Bkx_i/q.
\]
These are optimal, not merely feasible: \(c_i^y+ph_i^y=3x_i=w_i\);
\(\Lambda_i=1/(kx_i)\) and
\(\mu_i=(k-1)/(kx_i)>0\) satisfy the young conditions, including
\(\vartheta/n_i=\chi/x_i+\kappa/(x_i/p)\). Concavity suffices for optimality.

A renter saves zero when young. An owner has \(a_i'=-Ph_i^y\), borrows
principal \(qPh_i^y\), and repays it on entering old age; hence
\(qa_i'+\phi Ph_i^y=0\) and \(z_i=y_i^o\). Old owner financial saving is
\[
a_i^e=q(e_i-Ph_i^o)
=kx_i\left[\omega_B-\frac q{1-q}\right]>0.
\]
Both old tenure problems therefore have the same unconstrained optimum.
Conditional lifetime values coincide, so any finite logistic taste location
and positive scale produce a constant ownership share strictly between zero
and one.

Choose \(H_R<H_O\) large enough for the displayed choices and both dated
optima, retaining finite caps. For example,
\[
H_R>2(1+d)\left[\kappa m+(1+k)\bar x/p\right],\qquad H_O=2H_R
\]
is sufficient: the bracket is total housing per paired reference household,
and bounds each joint-planner home; the remaining displayed demands also
satisfy this bound. Housing clearing sets
\[
N=\frac{\bar H}{\kappa m+(1+k)\bar x/p}>0.
\]
Mean fertility equals \(1/\nu\), rents satisfy \(qr=p\), and the zero-tax
rebate is zero. This completes stationarity and the original finance,
estate, tenure, and housing checks.

## 3. Positive housing and joint-fertility effects

Put \(a=(1+k)/2>1\). The fixed-fertility planner has
\[
X=a\bar x,\qquad S=a\bar x/p.
\]
Therefore aggregate young housing increases by
\(N(a-1)\bar x/p>0\). Young aggregate gross consumption also increases.

At \(n=m\), the right side of the joint-planner fertility equation is
\(1/(am)<1/m\). Its left side decreases and its right side increases in
\(n\); consequently \(n^B>m\). Thus both desired aggregate signs hold
throughout the family \(k>1\).

## 4. Nevertheless private rechoice can lower average fertility

Measure goods and space in child units:
\[
R=1+u,\qquad T=1+u^{-1},\qquad
W=R+T=RT>4,\qquad M=\max\{R,T\}.
\]
Let \(\delta_i(a)=n_i^A-n_i\). Its exact condition is
\[
\frac1{n_i+\delta_i}
=\frac1{aRm-\delta_i}+\frac1{aTm-\delta_i}.
\]
At the auxiliary value \(a=1\), writing \(z=n_i/m\), the relevant quadratic
root gives
\[
\frac{\delta_i(1)}m
=\frac{W-z-D(z)}3,\qquad
D(z)=\sqrt{z^2+Wz+W^2-3W}.
\]
Since
\[
D''(z)=\frac{3W(W-4)}{4D(z)^3}>0,\qquad D(1)=W-1,
\]
heterogeneity implies the strictly positive, explicit number
\[
J=\frac{D(1-d)+D(1+d)-2D(1)}6>0,
\qquad \int\delta_i(1)dQ=-mJ.
\]

This strict sign is extended by a bound, not by a numerical point or an
unspecified continuity neighborhood. Implicit differentiation gives
\[
0<\partial_a\delta_i
=\frac{Rm/(aRm-\delta_i)^2+Tm/(aTm-\delta_i)^2}
{1/(n_i+\delta_i)^2+1/(aRm-\delta_i)^2+1/(aTm-\delta_i)^2}
<mM.
\]
Hence the **analytically nonempty primitive interval**
\[
\boxed{1<k<1+\frac JM}
\]
implies
\[
\int(n_i^A-n_i)dQ
<-mJ+mM\frac{k-1}{2}
<-\frac{mJ}{2}<0.
\]
Every equilibrium in this interval has strictly constrained young households,
higher planner young housing, and \(n^B>m>\int n_i^A\,dQ\).
The sufficient interval is narrow and conservative; no quantitative relevance
is claimed.

## 5. Sufficient restrictions and the general-preference lesson

For the log model, a clean restriction eliminates this obstruction. If the
fixed-fertility planner's common residuals satisfy
\[
X/\chi=S/\kappa=:r,
\]
then private rechoice is affine in reference fertility:
\[
n_i^A=\frac{\vartheta}{1+\alpha+\vartheta}(n_i+r).
\]
Its average rises exactly when \(\vartheta r>(1+\alpha)m\), which is also
the sign condition for \(n^B>m\) when the joint optimum is uncapped. This
balances adult goods and space measured in child units; it is an additional
restriction, not a general property of the planner.

For unspecified gross utility \(a(c,n)+b(h,n)+v(n)\), local complementarity
\(a_{cn},b_{hn}>0\) signs individual resource responses but does not aggregate
them through a redistribution. Even a concave, degree-one homogeneous
fertility choice function is insufficient. In the counterfamily, the exact
private rule is
\[
\mathfrak n(c,h)=\frac{c/\chi+h/\kappa
-\sqrt{(c/\chi)^2-(c/\chi)(h/\kappa)+(h/\kappa)^2}}3.
\]
It is increasing, homogeneous of degree one, and concave because the square
root is a norm. Both young mean resources rise, yet mean fertility falls.
Concavity provides an upper Jensen bound, not the needed lower bound.

A transparent general sufficient restriction is an affine conditional rule
\(\mathfrak n(c,h)=r_0+r_c c+r_hh\), \(r_c,r_h\ge0\): weak increases in both
young mean resources then cannot lower mean fertility. It follows, for example,
from \(U_n=\lambda(r_0+r_c c+r_hh-n)\), \(\lambda>0\), with suitable
concavity and monotonicity on the specified domain. Alternatively, if all
reference and new bundles lie on the same ray, a homogeneous rule makes
aggregate fertility depend only on aggregate scale. Neither restriction
follows from joint concavity.

Checked analytically against the supplied budgets, fertility first-order
conditions, planner equations, and demographic law. No model run, policy
implementation, or claim for binding-cap regimes is included.

## 6. Joint fertility with binding physical caps

**A proved extension.** Let \(\bar x_d\) be reference mean young adult goods
conditional on retained tenure \(d\), \(\bar x\) their overall mean, and
\(\bar c^o\) old mean consumption. Suppose \(\alpha\ge\gamma\), reference
old aggregate housing is at least young aggregate housing, and
\(\bar H<2N\int H_d\,dQ\). If
\[
\boxed{\frac{\bar x+\bar c^o}{2}\ge
\max\{\bar x_R,\bar x_O\},} \tag{C}
\]
the joint planner strictly raises average fertility, including when its
physical caps bind. The retained tenure weights \(\pi_R,\pi_O>0\) stay fixed,
sum to one, and agree across the two current cohorts. All fertility choices
are interior. Reference endowments and tenure remain heterogeneous.
Condition (C) concerns equilibrium averages, not primitive parameters.

Paired \(c_i^o\ge x_i\) guarantees \(\bar c^o\ge\bar x\). It therefore suffices
for (C) when adult-goods means are equal across retained tenures.
With tenure sorting, (C) requires the old mean to cover that difference:
\(\bar c^o\ge2\max_d\bar x_d-\bar x\). This is sufficient, not necessary.

**Proof.** For \(X,h>0\), define \(\mathcal N(X,h)\) as the positive solution
of the fertility first-order condition
\[
\frac{\vartheta}{n}=\frac{\chi}{X}
+\frac{\alpha\kappa}{h-\kappa n}.
\]
This is a first-order-condition solution map, not Experiment A's private
choice holding a gross bundle fixed. Writing \(t=h/\kappa\),
\[
\mathcal N(X,h)=
\frac{\chi t+(\alpha+\vartheta)X
-\sqrt{[\chi t+(\alpha-\vartheta)X]^2+4\alpha\vartheta X^2}}
{2\chi}.
\]
It is increasing in both arguments, homogeneous of degree one, jointly
concave, and strictly concave in \(h\) at fixed \(X\). Concavity follows
because the square root is a norm of a nonsingular linear transformation.
Reference fertility satisfies \(n_i=\mathcal N(x_i,h_i^y)\).

At the joint optimum, adult consumption equals a common \(X^J\) for both
ages. Denote mean fertility by \(m^J\). Goods clearing gives
\[
2X^J+\chi m^J=\bar x+\bar c^o+\chi m.
\]
Suppose, for contradiction, \(m^J\le m\). Then
\(X^J\ge(\bar x+\bar c^o)/2\ge\bar x_d\).
Conditional Jensen inequalities imply
\[
\bar n_d^{eq}
\le\mathcal N(\bar x_d,\bar h_d^{y,eq})
\le\mathcal N(X^J,\bar h_d^{y,eq}). \tag{7}
\]

For the joint housing multiplier \(\lambda_H>0\), uncapped young
fertility equals \(\vartheta/(\chi/X^J+\kappa\lambda_H)\), and adult space
equals \(\alpha/\lambda_H\). Thus the desired gross house is common:
\[
h_*=\frac{\kappa\vartheta}{\chi/X^J+\kappa\lambda_H}
+\frac{\alpha}{\lambda_H},\qquad
h_d^{y,J}=\min\{H_d,h_*\}.
\]
Clipping follows from the strictly concave joint problem in housing and
fertility at fixed resource multipliers.

This allocation also maximizes \(\int\mathcal N(X^J,h_i)dQ\) at its own
young housing total. To verify, set
\(\zeta=\mathcal N_h(X^J,h_*)\) when some young household is uncapped.
Every uncapped marginal equals \(\zeta\); every capped marginal is at least
\(\zeta\), because \(\mathcal N_h\) decreases in \(h\).
These are exactly the sufficient conditions for that concave maximization
problem. If all young households are capped, their allocation is forced.

Moreover, the joint housing conditions give
\(h_d^{y,J}\ge h_d^{o,J}\), strictly wherever old housing is uncapped:
uncapped young housing equals \(\kappa n+\alpha/\lambda_H\), whereas old
housing equals \(\gamma/\lambda_H\). Unused aggregate capacity therefore gives
\[
H_Y^J>\bar H/2\ge H_Y^{eq}.
\]
Starting from the tenure-mean reference housing allocation in (7), this
strictly larger young total permits a positive feasible housing addition
somewhere below the retained caps. Since \(\mathcal N\) is increasing,
maximization at the joint allocation yields
\[
m^J=\int\mathcal N(X^J,h_d^{y,J})dQ
>\sum_d\pi_d\mathcal N(X^J,\bar h_d^{y,eq})
\ge m,
\]
a contradiction. This proof uses the actual optimum, not merely a feasible
paired resource average.

**A conservative primitive certificate.** In the zero-tax subcase of the
independent housing theorem, let \(w=y^y+b\), \(v=y^o\), and
\(K=1+\gamma+\omega_B\). Define its resource-loss bound
\[
\ell=\min\{1,(1-\phi)/(1-q)\},\qquad
\delta_0=\frac{\alpha+\vartheta}{q(1+\alpha+\vartheta)}(\ell^{-1}-1).
\]
That lemma gives \(z_i\ge v_i-\delta_0w_i\),
\(c_i^o\ge z_i/K\), and \(x_i<w_i\), allowing both estate regimes and
binding caps. Hence
\[
\boxed{\bar v-\delta_0\bar w\ge2K w_{\max}} \tag{8}
\]
implies (C). Add (8) to the housing theorem's primitive conditions.
It is a strong future-income restriction relative to current-cash dispersion,
not a proposed calibration. It imposes neither side of \(\beta/q=1\).
Raising bounded old incomes can satisfy it without changing the maintained
planner or eliminating physical caps.

**What remains open.** Paired \(c_i^o\ge x_i\) alone gives only
\(X^J\ge\bar x\) under the contradiction; it does not give the tenure-specific
inequalities needed in (7). Sufficiency of that weaker ordering in a positive
stationary equilibrium with unused capacity remains unproved here.
