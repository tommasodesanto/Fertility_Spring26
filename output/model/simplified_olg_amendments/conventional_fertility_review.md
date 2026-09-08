# Fertility under conventional mortgage finance

Independent mathematical review, September 8, 2026. This is supporting work for
`docs/model/simplified_olg_conventional_finance_proposal.md`; the proposal and
the protected manuscript are unchanged. No model simulation or parameter search
is used. All statements below hold at a fixed price and rebate path.

## Main findings

There are two useful results, with different economic scopes.

1. **Housing and fertility together.** If the cash cost of owner housing is at
   least its lifetime service cost, \(L\ge p\), a strictly binding mortgage
   constraint implies that a marginal relaxation increases both housing and
   fertility. This holds for every strictly positive goods cost and space cost
   of children, every positive discount factor, and every nonnegative future
   income. Young and old housing caps must be slack for the strict housing
   derivative; a young cap allows weak housing growth and still strict fertility
   growth while finance binds.
2. **Fertility when \(L<p\), including tenure switching.** Write
   \(B=\beta(1+\gamma+\omega_B)\). If
   \[
   \boxed{\qquad \frac{p}{\alpha+B}\le\frac{\chi}{\kappa}
                    \le\frac{p}{\alpha},\qquad L\le p,\qquad}
   \]
   conditional-owner fertility rises with credit access, and owner fertility
   is at least renter fertility even when the rental size cap is smaller.
   Consequently a credit relaxation weakly increases fertility averaged over
   the initial heterogeneous cohort, including its change in tenure shares.
   Old-age housing caps must be slack at the relevant conditional optima.
   **Housing itself can decrease in this second result.**

These are primitive restrictions, not assumptions about unknown multipliers.
Neither requires a small \(\beta\), a numerical reference equilibrium, or a
planner preference for population. They are fixed-price household and
fixed-cohort statements, not general-equilibrium or transition theorems.

## 1. Reduced problem and exact constraint test

Let \(x=c-\chi n\) be adult consumption, \(s=h-\kappa n\) adult space,
\(S=qz\) discounted old-age resources, and \(V=qv\) discounted future income.
Set \(\eta=\chi/\kappa>0\), the goods cost of a child per unit of required
space, and \(E=1+\alpha+\vartheta\). With uncapped old-age housing, conditional
old utility is \(K\log z+C_m\), where \(K=1+\gamma+\omega_B\) and
\(B=\beta K>0\). The constants \(C_O,C_R\) may differ because of the owner's
estate floor; they affect tenure choice but not conditional allocations.

The young owner solves
\[
\max\;\log x+\alpha\log s+\vartheta\log n+B\log S
\]
subject to
\[
x+ps+(\chi+p\kappa)n+S=w+V,\qquad
x+Ls+(\chi+L\kappa)n\le w,\qquad s+\kappa n\le H_O.
\]
Assume \(p>0,L\ge0,w>0,V\ge0\) and positive preference weights. Strict
concavity gives a unique real allocation. Let \(\Delta=L-p\). At a strictly
binding financial constraint,
\[
c=w-Lh,\qquad S=V+\Delta h.
\]
In the original mortgage convention, \(dL/d\phi=-P\). Every sign result below
is more generally a result about lowering \(L\) at fixed \(p,w,V\); a different
loan covenant can therefore use the same theorem if its reduction is valid.

With slack young and old caps, define
\[
H(p)=\frac\alpha p+\frac\vartheta{\eta+p},\qquad D=E+B.
\]
Finance strictly restricts the household exactly when
\[
\boxed{\quad \frac{V}{w}>\frac{D}{E+(L-p)H(p)}-1.\quad}
\]
This is the proposal's income-timing test in shorter notation. The denominator
is positive because it equals
\(1+\alpha L/p+\vartheta(\eta+L)/(\eta+p)\).
If a cap binds, use the proposal's exact test \(c^*+Lh^*>w\), retaining that cap
when computing the allocation with unrestricted finance.

## 2. A sharp primitive fertility theorem

**Proposition 1.** Suppose young and old housing caps are slack and finance is
strictly binding.

- If \(L\ge p\), then \(\partial n/\partial\phi>0\) for every
  \(\chi,\kappa,\alpha,\vartheta,B>0\) and \(V\ge0\). The proposal's housing
  proof also gives \(\partial h/\partial\phi>0\).
- If \(0\le L<p\), the condition
  \[
  \boxed{\quad Bp(L+\eta)\ge(p-L)(p-\alpha\eta)\quad}                 \tag{1}
  \]
  guarantees \(\partial n/\partial\phi>0\) for every strictly constrained
  income profile. A nonpositive right-hand side makes the condition automatic.
  Condition (1) is sharp for a guarantee uniform over current and future
  income at this fixed \(L,p\): if it fails, there are strictly constrained
  positive-income households whose fertility falls when credit relaxes.

Consequently,
\[
\boxed{\quad (\alpha+B)\chi\ge p\kappa\quad}                         \tag{2}
\]
is sufficient for increasing fertility for every nonnegative \(L\). It is also
necessary for a guarantee uniform over all \(L\ge0\) and income profiles in
the uncapped model. For a finite permitted range \(L\ge L_{\min}\), use the
less restrictive condition (1) evaluated at \(L_{\min}<p\), together with the
automatic \(L\ge p\) result. The fertility weight \(\vartheta\) affects levels
and who is constrained, but cancels out of this uniform sign restriction.

**Proof.** Rescale children by their space requirement, \(k=\kappa n\); this
changes utility only by a constant. Work with \(\eta k\) goods expenditure.
At the uncapped constrained optimum, define the effective marginal price
\(\rho=\alpha x/s\) and the auxiliary number \(t=Bx/S\). The first-order
conditions imply
\[
0<t<1,\qquad \rho=(1-t)L+tp,\qquad
\frac hx=H(\rho):=\frac\alpha\rho+\frac\vartheta{\eta+\rho}.
\]
The strict inequality \(t<1\) is the positive marginal value of relaxing
finance; these auxiliary objects are used only in the proof, not in (1)-(2).

Put \(A=\rho^2/\alpha\). Implicit differentiation of the strictly concave
two-variable problem in \((h,k)\) shows that the fertility derivative has the
sign of
\[
\mathcal N=(1-t)(A-L\eta)
 +H(\rho)\left[A(L+\eta)
       +\frac{t^2(L-p)}B(A-p\eta)\right].                            \tag{3}
\]
For completeness, the Hessian entries before this simplification are
\[
\begin{aligned}
a&=L^2/x^2+\alpha/s^2+B\Delta^2/S^2,\\
b&=-L\eta/x^2+\alpha/s^2,\\
d&=\eta^2/x^2+\alpha/s^2+\vartheta/k^2,\\
F&=(x+Lh)/x^2-BV/S^2.
\end{aligned}
\]
Their determinant \(ad-b^2\) is positive, and
\[
\frac{\partial k}{\partial\phi}
 =\frac{P}{ad-b^2}\left(bF+a\eta h/x^2\right),\qquad
bF+a\eta h/x^2=\mathcal N/x^3.
\]

First remove the last term in (3), calling the remainder \(\mathcal N_0\).
Since \(H(\rho)\ge\alpha/\rho\),
\[
\mathcal N_0\ge(1-t)A+\rho L+\eta tp>0.                            \tag{4}
\]
If \(L\ge p\) and \(A\ge p\eta\), the last term is nonnegative. If
\(A<p\eta\), the accounting identity \(S=V+\Delta h\) and \(V\ge0\) give
\(t\Delta H(\rho)\le B\). Therefore the negative last term is bounded below
by \(t(A-p\eta)\), and
\[
\mathcal N\ge A+H(\rho)A(L+\eta)-\eta\rho
            \ge A+\rho L>0.
\]
This proves the first claim without restricting either child cost or \(B\).

If \(L<p\), the last term is nonnegative when \(A\le p\eta\). Otherwise
\(A>p\eta>L\eta\), and the first term of (3) is strictly positive. Since
\(\rho<p\) and \(t<1\),
\[
t^2\left(1-\frac{p\eta}{A}\right)
 <1-\frac{\alpha\eta}{p}.
\]
Condition (1) then makes the expression in square brackets in (3) positive.
To see sharpness, let \(t\) approach one from below. Then \(\rho\) approaches
\(p\), and the limiting sign expression divided by \(H(p)\) is
\[
\frac{p^2}{\alpha}(L+\eta)
 -\frac{p-L}{B}\left(\frac{p^2}{\alpha}-p\eta\right).
\]
It is negative exactly when (1) fails. Such \(t\)'s correspond to admissible
positive-income households, rather than to arbitrary shadow prices: choose
any \(x>0\), set
\[
w=x\{E+t(L-p)H(\rho)\},\qquad
V=x\{B/t-(L-p)H(\rho)\}>0.
\]
The resulting allocation satisfies both budgets and all first-order conditions.
Finally, when \(\eta<p/\alpha\), the right-hand threshold for \(B\) in (1)
is decreasing in \(L\) and is largest at \(L=0\), where it is
\(p/\eta-\alpha\). This yields (2). \(\square\)

## 3. Finite reforms and physical caps

**Young cap.** If the young owner cap binds, housing is fixed at \(H_O\).
While finance also binds, the fertility first-order condition gives
\[
\frac{\partial n}{\partial\phi}
 =\frac{P\chi H_O/x^2}
 {\chi^2/x^2+\alpha\kappa^2/s^2+\vartheta/n^2}>0.                   \tag{5}
\]
If finance is slack, changing \(\phi\) leaves the conditional allocation
unchanged. Strict concavity gives continuity at changes of constraint status.
Thus Proposition 1 extends to finite reforms satisfying its restriction at
every intervening \(L\), with old caps remaining slack and an arbitrary young
cap. Fertility weakly
increases; it strictly increases if a nonzero part of the reform relaxes a
strictly binding constraint. The joint housing-and-fertility result over a
range \(L\ge p\) allows housing to stop rising when its cap is reached.

**Old caps.** A binding old housing cap changes the continuation function, so
the homogeneous \(B\) must not be used without checking it. Each active
old-owner regime instead has
\[
\mathcal V^O(z)=A_o\log(z-z_0)+C,
\]
where

| Active old regime | \(A_o\) | \(z_0\) |
|---|---:|---:|
| Housing cap slack; estate floor either status | \(K\) | \(0\) |
| Housing cap \(H_O\) binding; estate floor slack | \(1+\omega_B\) | \(u_oH_O\) |
| Housing cap and estate floor binding | \(1\) | \((u_o+qP_{o+1})H_O\) |

Within a regime replace \(B\) by \(\beta A_o\), \(S\) by
\(q(z-z_0)\), and \(V\) by \(q(v-z_0)\). Condition (1) remains sufficient
when \(L<p\); the shifted \(V\) is automatically positive since
\(S=V+(L-p)h>0\). Thus
\[
(\alpha+\beta)\chi\ge p\kappa
\]
is a convenient sufficient condition for owner fertility to increase through
all young and old cap transitions over a range \(L\le p\).

For \(L\ge p\), the shifted future income can be negative, so the first part
of Proposition 1 cannot simply be quoted. A sufficient alternative is
\(\eta\le p/\alpha\): then \(\rho\ge p\) implies \(A\ge p\eta\), and
(3)-(4) give the sign without an income assumption. Consequently the narrower
interval
\[
\frac{p}{\alpha+\beta}\le\frac\chi\kappa\le\frac p\alpha
\]
guarantees conditional-owner fertility monotonicity over all nonnegative
\(L\), including every old and young cap regime. This is sufficient, not a
claimed necessary restriction. It does not by itself order owner and renter
fertility, since their binding old caps change continuation incentives.

## 4. Tenure selection and heterogeneous households

Conditional fertility growth alone does not establish aggregate fertility
growth. If \(\pi\) is a type's ownership probability, its mean fertility is
\[
\bar n=\pi n_O+(1-\pi)n_R,
\qquad
\frac{d\bar n}{d\phi}
 =\pi\frac{dn_O}{d\phi}
  +\frac{d\pi}{d\phi}(n_O-n_R).                                    \tag{6}
\]
The owner's feasible set expands, so \(d\pi/d\phi\ge0\); this uses only a
fixed distribution of additive tenure tastes, not specifically logistic
tastes. The second term nevertheless needs a fertility ordering.

**Proposition 2.** Suppose old-age housing caps are slack at both relevant
conditional optima, \(H_O\ge H_R\), \(L\le p\), and
\[
\frac p{\alpha+B}\le\eta\le\frac p\alpha.                          \tag{7}
\]
Then \(n_O\ge n_R\), conditional-owner fertility is nondecreasing in a finite
credit relaxation, and the type's mean fertility is nondecreasing. The same
holds after integrating over any fixed heterogeneous distribution whose types
satisfy these inequalities. There is no common-income or common-taste
requirement.

**Proof.** A renter's no-borrowing restriction is exactly
\(c+ph\le w\). Hence, ignoring the allocation-irrelevant constants
\(C_O,C_R\), renting is the problem with \(L=p\) and cap \(H_R\).
First increase this young cap from \(H_R\) to \(H_O\), holding \(L=p\).
When a cap binds, its housing first-order inequality gives
\(\rho=\alpha x/s\ge p\). If finance also binds, the sign of the fertility
response to raising the cap is the sign of
\[
\alpha\kappa/s^2-p\chi/x^2
 =\frac\kappa{x^2}(\rho^2/\alpha-p\eta)\ge0.
\]
If finance is slack, optimizing saving first replaces this expression by
\[
\frac\kappa{x^2}\left(\rho^2/\alpha-
                            \frac{p\eta}{1+B}\right)>0.
\]
The upper bound in (7) therefore makes fertility nondecreasing as the cap
expands. Next lower \(L\) from \(p\) to its owner value. The lower bound in
(7), Proposition 1, and (5) make fertility nondecreasing on this second leg.
The endpoint is exactly the conditional owner allocation. Thus
\(n_O\ge n_R\). Equation (6), or its finite-change counterpart, completes
the argument. \(\square\)

A strictly positive mass of households whose finance strictly relaxes gives a
strict aggregate increase when these households have positive ownership
probability. A logistic tenure taste supplies a strictly positive probability
for every finite conditional-value difference. Binding rental space limits
can make the tenure fertility gap strict as well. No fertility conclusion
follows from the tenure-share response alone outside these restrictions.

## 5. Exact examples separating housing from fertility

The proposal's existing rational counterexample already provides the desired
separation. Its strictly constrained owner has \(h=2,n=1\), and the exact
derivatives are
\[
\frac{\partial h}{\partial\phi}=-\frac{8280}{25271}<0,
\qquad
\frac{\partial n}{\partial\phi}=\frac{26900}{25271}>0.
\]
Thus a fertility benefit does not require that total housing occupied increase.
For that example \(p=1,\alpha=31/40,B=1,\eta=1\), so (7) holds. Old caps
are slack. This is also an exact example consistent with the tenure-ordering
conditions above, not merely an unexplained numerical sign.

Fertility can fall when the primitive condition fails. Keep the same financial
prices as that example: \(q=1/2,P=P_{next}=2,\tau^p=0,\phi=19/20\), hence
\(p=1,L=1/10\). Set \(\alpha=\vartheta=B=\kappa=1\),
\(\chi=1/10\), and choose
\[
x=1,\quad s=\frac{100}{91},\quad n=\frac{100}{101},\quad
S=\frac{10}{9},\quad
w=\frac{12021}{9191},\quad V=\frac{247430}{82719}.
\]
The first-order conditions hold with lifetime multiplier \(9/10\) and cash
multiplier \(1/10\); both budgets hold exactly. An owner cap of three is
slack. Taking \(\beta=1/2,\gamma=1/4,\omega_B=3/4\) yields \(B=1\)
and slack old estate and housing constraints. Exact differentiation gives
\[
\frac{\partial h}{\partial\phi}
 =-\frac{182499868000}{92226644393}<0,\qquad
\frac{\partial n}{\partial\phi}
 =-\frac{417283358000}{645586510751}<0.
\]
This uses ordinary positive goods and space costs and no vanishing discount
factor. Additive ownership tastes can select this conditional owner without
changing any real choice. It disproves an unconditional fertility theorem.

Both examples' budgets and derivatives were checked using exact rational
arithmetic. Symbolic simplification independently verifies the identity in
(3): its residual is exactly zero. These checks did not solve or search a
numerical model.

## 6. Policy, compensation, and general equilibrium are distinct

At given \(c,h\), the fertility first-order condition gives the differential
\[
\left(\chi^2/x^2+\alpha\kappa^2/s^2+\vartheta/n^2\right)dn
 =\frac\chi{x^2}\,dc+\frac{\alpha\kappa}{s^2}\,dh.                 \tag{8}
\]
Providing extra housing while holding current total consumption fixed raises
fertility. If consumption pays for the extra house, fertility rises only when
\(-dc/dh<\rho^2/(\alpha\eta)\). Along a reallocation compensated to hold
young utility fixed, \(dc/dh=-\rho\), so fertility rises exactly when
\(\rho>\alpha\eta\). Compensation in old age instead produces a different
current-consumption path. None of these is automatically the credit-policy
derivative: that derivative must obey both financial budgets.

An allocation or policy can improve household welfare even if fertility falls;
fertility is chosen inside the retained utility function. No additional social
value of population has been introduced. Conversely, higher fertility does not
establish a welfare gain or an equilibrium housing misallocation.

In equilibrium, a credit policy can change \(p,L,w,V\), tenure values and the
distribution of future cohorts. The propositions hold these objects fixed
apart from the direct \(L\) change. They establish neither the sign of total
GE fertility nor a population transition, and they should not be represented
as doing so. A separate primitive equilibrium argument would have to control
those price, transfer, income and composition responses.
