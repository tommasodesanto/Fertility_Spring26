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
