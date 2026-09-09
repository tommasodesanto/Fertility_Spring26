# General preferences: housing allocation and fertility

Working memo, September 8, 2026. This separate branch starts from unspecified
preferences over **gross** goods, housing, and fertility. It changes neither
the household model nor the dated planner agreed in
`simplified_olg_utilitarian_work.md`. Linear child needs enter only
the final special case.

## 1. Unspecified utility and the maintained benchmark

Let young utility be \(U^y(c,h,n)\). Consumption \(c\) and housing \(h\) are
the gross resources in the existing budgets. Specify its domain as a primitive;
do not initially impose \(c>\chi n\) or \(h>\kappa n\). Assume twice continuous
differentiability, joint concavity, and \(U^y_c,U^y_h>0\). Old utility
\(U^o(c_o,h_o,e)\) is increasing and jointly concave. All comparisons retain
the specified cardinal utility scales.

**Do not assume \(U^y_n>0\) everywhere at fixed \(c,h\).** Without a separate
fertility resource cost or an upper domain bound, that assumption rules out
a finite interior fertility optimum. Net marginal utility of children may
become negative because of costs or crowding represented inside \(U^y\).
Concavity alone does not guarantee existence on an unbounded fertility domain.
The derivative results below concern an existing regular interior optimum;
active domain restrictions would add their own multipliers.

Each current age group has mass \(N\) and the same stationary probability law
\(Q\) of paired types and retained tenures. The physical cap is
\(H_i\in\{H_R,H_O\}\), with \(H_R<H_O\). First fix individual fertility and
tenure. The planner chooses **all current consumption and housing**, fixes
future real continuation opportunities and old net estates, and relaxes
individual financing restrictions while honoring obligations. It maximizes
\[
N\int[U^y(c_i^y,h_i^y,n_i)+U^o(c_i^o,h_i^o,e_i)]\,dQ
\]
subject to
\[
N\int(c_i^y+c_i^o)dQ=C,\qquad
N\int(h_i^y+h_i^o)dQ=\bar H,\qquad h_i^a\le H_i.
\]
Here \(C,\bar H\) are the reference resource totals. Write
\(H_Y=N\int h_i^y\,dQ\), \(H_O=N\int h_i^o\,dQ\); superscripts \(eq,*\)
denote the competitive reference and planner allocation.
Continuation and fixed-tenure taste terms are constant. The owner's financial
estate floor is relaxed; the promised net estate is unchanged.

With resource multipliers \(\lambda_C,\lambda_H\) and planner cap multipliers
\(\eta_i^{a*}\ge0\), interior-domain optimality requires
\[
U^y_c=U^o_c=\lambda_C,\qquad
U^y_h=\lambda_H+\eta_i^{y*},\qquad
U^o_h=\lambda_H+\eta_i^{o*}.
\]
Together with feasibility and complementary slackness these characterize an
optimum under concavity, if one exists. Strict concavity ensures uniqueness.
Without separability, consumption changes housing marginal utilities.

## 2. The competitive wedge is fully general

At a positive stationary equilibrium, let \(q=1/R_f\),
\(p=(1-q+q\tau^p)P\), \(L_R=p\), and
\(L_O=(1-\phi+q\tau^p)P\). Here \(p\) is housing's service cost and \(L_d\)
its coefficient in the current financing constraint. Current cash is
\(w_i=y_i^y+b_i+T\), old income including rebate is \(v_i^o=y_i^o+T\), and
\(z_i\) is total old-age resources after repayment of the young mortgage.
The reduced budgets remain
\[
c+ph+qz=w_i+qv_i^o,\qquad c+L_dh\le w_i.
\]
The young objective is \(U^y(c,h,n)+\beta V_d(z)\). The maintained continuation
value has **no direct dependence on \(n\)**; \(\beta\) discounts the parent's
own old age. For budget and financing multipliers \(\Lambda_i,\mu_i\), and
competitive young cap multiplier \(\eta_i^y\),
\[
q\Lambda_i=\beta m_i,\qquad
U^y_c=\Lambda_i+\mu_i,\qquad
U^y_h=p\Lambda_i+L_d\mu_i+\eta_i^y,
\quad m_i=V_d'(z_i)>0.
\]
For the old budget \(c_o+ph_o+qe=z\), put \(\rho_i\ge0\) on the owner's
floor \(e-Ph_o\ge0\), and set \(\rho_i=0\) for renters. Let \(\eta_i^o\)
be the competitive old housing-cap multiplier. The envelope and old
first-order conditions give
\[
U^o_c=m_i,\qquad U^o_e=qm_i-\rho_i,\qquad
U^o_h=pm_i+P\rho_i+\eta_i^o.
\]
Stationarity matches the young household's future old allocation with its
current old counterpart. Consequently
\[
\boxed{U^y_h-U^o_h
=\left(\frac{\beta}{q}-1\right)pm_i
+L_d\mu_i+\eta_i^y-P\rho_i-\eta_i^o.}
\]
This uses neither logarithms nor linear child costs. The old estate floor
raises old direct housing marginal utility. No restriction on \(\beta R_f\)
has been imposed. A bound used to sign this identity is an additional
sufficient restriction, not a consequence of logarithmic utility.

A positive gap permits a small improving transfer to an **uncapped** young
recipient, holding consumption fixed. It does not establish the full
optimum's aggregate direction. Young finance need not dominate old finance
or the age-weight term.

Cardinal comparison matters independently of curvature: replacing \(U^o\)
by \(A U^o\) and \(\beta\) by \(\beta/A\) leaves all competitive choices
unchanged but multiplies the planner's weight on old utility by \(A>0\).
Thus concavity alone cannot establish a universal age direction. This is a
family of different social comparisons, not a normalization of one fixed
criterion.

## 3. Fertility and a useful intermediate class

Given a gross bundle and fixed continuation opportunities, an interior
fertility choice satisfies \(U^y_n=0\). If \(U^y_{nn}<0\), then
\[
dn=-\frac{U^y_{nc}\,dc+U^y_{nh}\,dh}{U^y_{nn}}.
\]
Positive consumption and housing responses therefore require positive
cross-partials with fertility. Joint concavity does not sign them.

For a counterexample with no linear resource offsets, take
\[
U^y(c,h,n)=\sqrt c+\sqrt{h+n}+A\sqrt n-kn,\qquad A,k>0.
\]
This is strictly concave and increasing in \(c,h\). At each positive bundle
there is a unique interior fertility optimum: \(U_n\) decreases from
\(+\infty\) to \(-k\). Yet
\(U_{nh}=-[4(h+n)^{3/2}]^{-1}<0\), so more housing lowers fertility.

An intermediate specification separates the two interactions:
\[
\boxed{U^y(c,h,n)=a(c,n)+b(h,n)+v(n).}
\]
Assume \(a_c,b_h,v'>0\); \(a_n,b_n\) may be negative due to goods costs and
crowding. Require joint concavity and an existing interior optimum. Define
\(D=a_{nn}+b_{nn}+v''<0\). Its fertility condition and response are
\[
a_n+b_n+v'=0,\qquad
dn=-\frac{a_{cn}\,dc+b_{hn}\,dh}{D}.
\]
Thus \(a_{cn}>0\) and \(b_{hn}>0\) suffice for both resource effects to be
positive, without a shift specification. They say that additional resources
raise the marginal utility of children. If consumption falls while housing
rises, this numerator gives the local tradeoff.

At fixed fertility this class also separates the planner's consumption and
housing problems. A useful **additional age-comparison restriction** is
\[
U^o=f_o(c_o)+b(h_o,0)+B(e),\qquad b_{hh}<0,\quad b_{hn}>0.
\]
The old housing component is cardinally matched to the childless young
component. Suppose the cross-partial restriction holds between \(0\) and
each \(n_i>0\), on a common feasible housing domain, and positive solutions
exist. Then \(b_h(h,n_i)>b_h(h,0)\). A common planner housing multiplier
therefore yields
\[
h_i^{y*}\ge h_i^{o*},
\]
strictly wherever the old counterpart is uncapped. If
\(\bar H<2N\int H_i\,dQ\), not all old households can be capped, so
\(H_Y^*>\bar H/2\). Hence \(H_Y^{eq}\le H_O^{eq}\) suffices for a strict
aggregate young housing gain. This is a theorem under explicit preference
and reference-allocation restrictions, not a generic implication of finance.
Neither this result nor \(b_{hn}>0\) proves that every young household gains.

## 4. Joint fertility choice by the dated planner

If the planner also chooses \(n_i\), the gross-variable conditions are simply
\[
U^y_c=\lambda_C,\qquad
U^y_h=\lambda_H+\eta_i^{y*},\qquad U^y_n=0
\]
at an interior choice. There is no additional fertility resource term in
the maintained gross budgets. The planner values children through currently
living parents; it attaches no new welfare weight to future people.
Conditional fertility responses do not establish greater average fertility
at this joint optimum. Changing births also changes future entrant masses,
so this is not a completed dynamic allocation holding every future
population fixed.

Separate explicit child-resource constraints would constitute another model
architecture. They require their own definitions and are not adopted here.
The same child goods or space already included in gross \(c,h\) must not be
charged again as an additional aggregate resource requirement.

## 5. Restricted mapping: the existing shifted specification

Only now impose the existing form
\[
a(c,n)=\log(c-\chi n),\qquad
b(h,n)=\alpha\log(h-\kappa n),\qquad
v(n)=\vartheta\log n.
\]
Its domain \(c>\chi n,\ h>\kappa n,\ n>0\) belongs to **this special case**.
Old utility is \(\log c_o+\gamma\log h_o+\omega_B\log e\).
Here \(a_{cn}=\chi/(c-\chi n)^2>0\) and
\(b_{hn}=\alpha\kappa/(h-\kappa n)^2>0\), so the intermediate fertility result
reproduces the existing positive conditional responses. The gross condition
\(U^y_n=0\) becomes exactly
\[
\frac{\vartheta}{n}
=\frac{\chi}{c-\chi n}+\frac{\alpha\kappa}{h-\kappa n}.
\]

A still-restricted nonlogarithmic extension replaces the two logarithms by
\(f(c-\chi n)\) and \(\alpha g(h-\kappa n)\), with increasing, strictly
concave functions. With additively separable old housing utility
\(\gamma g(h_o)\), fixed-fertility planner housing is
\[
h_i^{y*}=\min\{H_i,\kappa n_i+r_y\},\quad
h_i^{o*}=\min\{H_i,r_o\},\quad
r_y=(g')^{-1}(\lambda_H/\alpha),\quad
r_o=(g')^{-1}(\lambda_H/\gamma).
\]
Assume positive solutions and the appropriate derivative range, extending
the inverse to infinity when desired housing exceeds every finite cap.
Let \(\bar n=\int n_i\,dQ\). For \(\alpha=\gamma\) and slack planner caps,
\[
r_y=r_o=\frac{\bar H/N-\kappa\bar n}{2},\qquad
H_Y^*=\frac{\bar H+\kappa N\bar n}{2}.
\]
This equalizes **adult space**, a concept specific to the shifted class.
Young household \(i\) gains housing exactly when
\(h_i^{y,eq}-\kappa n_i<r_y\). Aggregate young housing increases exactly when
old mean reference housing exceeds young mean reference adult space.
For \(\alpha\ge\gamma\), the paired planner housing ordering survives binding
caps. These formulas are not conclusions for unspecified \(U^y(c,h,n)\).

The fully general competitive age ordering, equilibrium existence, aggregate
fertility comparison, and funded policy transition remain unproved.
