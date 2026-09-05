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
