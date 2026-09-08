# Ownership, financing, and family housing

Discussion proposal, September 7, 2026. This document develops the household
problem and the conditions under which financing restricts housing. It does
not amend the existing theory note or establish equilibrium inefficiency.

## 1. Proposed environment

**Households.** Households live for two periods. A young household has liquid
wealth \(b>0\), income \(y^y>0\) when young, and income \(y^o\ge0\) when old. The joint
distribution of these endowments can be heterogeneous. Future income is known
when the young household chooses housing, fertility, saving, and tenure.

**Preferences.** Retain the existing preferences:
\[
u^y(c,h,n)=\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta\log n,
\qquad
u^o(c^2,h^2,e)=\log c^2+\gamma\log h^2+\omega_B\log e.
\]
Here \(c\) includes spending on children, each child requires \(\chi\) goods
and \(\kappa\) units of space, and \(e\) is the estate left at death.
All log arguments are positive. Fertility is one continuous choice, and
households discount future utility by \(\beta>0\). The superscript 2 labels
old-age quantities, as in the existing note.

**Markets.** Retain the fixed housing stock, external bond market, rental
intermediary, and property-tax rebate. The bond price is \(q\), its gross return
is \(1/q\), and the house price is \(P_t\). The beginning-of-period cost of
renting is \(u_t=qr_t\), where:
\[
u_t=(1+q\tau^p)P_t-qP_{t+1}>0.
\]
Thus a house's purchase price includes a resale claim as well as its housing
services. Renters can occupy at most \(h_R^{\max}\), and owners at most
\(h_O^{\max}>h_R^{\max}\).

**Finance.** Current income and liquid wealth are available when the household
buys. A mortgage can finance at most \(\phi_t\) of the purchase price.
Unsecured borrowing is unavailable. The mortgage is repaid when the household
enters old age. The tax payment is reserved at purchase, preserving the
existing timing convention.

**Old-age housing.** The proposed baseline retains tenure chosen when young,
but lets old owners buy or sell housing within the owner size limit. Their
choice is no longer bounded by the size of their inherited home. Old
households do not borrow. Retained owner housing is sold at death and enters
the estate; estates do not enter the next cohort's initial wealth.

This last paragraph is a substantive proposal, not an accounting correction.
Keeping lifetime tenure avoids adding another discrete choice to the present
exercise. Allowing old owners to resize removes the inherited-home restriction
that previously entered the housing-value comparison. Old-age tenure switching
remains a separate possible extension.

## 2. Household problems

**Old households.** Let \(z\) be total resources available at the beginning of
old age, including old-age income. Conditional on tenure, define:
\[
\begin{aligned}
\mathcal V_t^R(z)&=\max_{c^2,h^2,e}u^o(c^2,h^2,e),\\
&c^2+qe+u_th^2=z,\qquad 0<h^2\le h_R^{\max},\\[3pt]
\mathcal V_t^O(z)&=\max_{c^2,h^2,e}u^o(c^2,h^2,e),\\
&c^2+qe+u_th^2=z,\qquad 0<h^2\le h_O^{\max},\quad e\ge P_{t+1}h^2.
\end{aligned}
\]
The owner's resources are \(z=a+P_tH+y^o+T_t\); the renter's are
\(z=a+y^o+T_t\). The notation \(\mathcal V\) expresses the existing value
function in total resources: allowing owners to resize makes inherited
financial wealth and housing wealth interchangeable at the prevailing price.

The owner's estate still consists of financial saving and the house sold at
death:
\[
e=q^{-1}a^e+P_{t+1}h^2,\qquad a^e\ge0.
\]
Substituting this expression gives the equivalent cash budget
\(c^2+a^e+(1+q\tau^p)P_th^2=a+P_tH+y^o+T_t\).
The estate floor is therefore the restriction against borrowing in old age.

**Young households.** Write \(w_t=y^y+b+T_t\) for current resources and
\(v_{t+1}=y^o+T_{t+1}\) for future income and the future rebate. The variable
\(a'\) remains financial wealth on entering old age, net of mortgage repayment.
Young owners solve:
\[
\begin{aligned}
W_t^O=\max_{c,h,n,a'}\;&u^y(c,h,n)
 +\beta\mathcal V_{t+1}^O(a'+P_{t+1}h+v_{t+1}),\\
&c+qa'+(1+q\tau^p)P_th=w_t,\\
&qa'+\phi_tP_th\ge0,\qquad h\le h_O^{\max}.
\end{aligned}
\]
Current income can pay for either consumption or the down payment. There is
no separate requirement that the down payment come from \(b\) alone.

To see the mortgage explicitly, let \(d\) be the amount borrowed at purchase
and \(k\) the amount invested in bonds. Then the same feasible set is:
\[
c+k+(1+q\tau^p)P_th=w_t+d,\qquad
k\ge0,\quad 0\le d\le\phi_tP_th,\quad a'=(k-d)/q.
\]
Eliminating \(k,d\) gives the preceding net-asset constraint. With equal
borrowing and saving rates the gross portfolio need not be unique, but the
net financial position and real choices are well defined.

Young renters solve:
\[
\begin{aligned}
W_t^R=\max_{c,h,n,a'}\;&u^y(c,h,n)
 +\beta\mathcal V_{t+1}^R(a'+v_{t+1}),\\
&c+qa'+u_th=w_t,\qquad a'\ge0,\quad h\le h_R^{\max}.
\end{aligned}
\]
Retain the existing ownership taste: a household owns if
\(W_t^O+\xi\ge W_t^R\). The existing logistic distribution gives the same
ownership probability formula with these conditional values.

## 3. What the mortgage limit does

Hold the price and rebate path fixed. Suppress time subscripts and define
\(p=u_t\), the cost of current housing services, and
\(L=(1-\phi_t+q\tau^p)P_t\), the cash required per unit of owner housing,
including the tax reserve. The owner's budget and borrowing restriction imply:
\[
c+Lh\le w.
\]
The household must pay current consumption and the cash portion of housing
from current resources. In particular, \(h<w/L\). Future income supports
repayment but does not remove this initial cash requirement.

Let \(z=a'+P_{t+1}h+v\) be old-age resources, and let \(M=w+qv\) be the
present value of current resources and future income. The complete young
owner problem can equivalently be written:
\[
\max_{c,h,n,z}\ u^y(c,h,n)+\beta\mathcal V_{t+1}^O(z)
\quad\text{subject to}\quad
c+ph+qz=M,\quad c+Lh\le w,\quad h\le h_O^{\max}.
\]
This separates the lifetime resource cost of housing from the cash needed
to acquire it. Both restrictions come directly from the ordinary mortgage
budget.

**An exact test, including housing caps.** First solve this conditional owner
problem without \(c+Lh\le w\), keeping all other restrictions. Denote its
unique real allocation by \((c^*,h^*,n^*,z^*)\). The mortgage limit strictly
restricts the household's choice if and only if:
\[
c^*+Lh^*>w.
\]
If the inequality is reversed the desired allocation is affordable. At
equality the limit can bind with zero marginal value. This test compares the
actual problem with one precisely stated financing relaxation; it is not
an efficiency theorem.

### An explicit income condition

When the relevant old-age housing caps are slack, old-age indirect utility
has the form \(\mathcal V^m(z)=K\log z+C_m\), with
\(K=1+\gamma+\omega_B\). The owner's estate floor may bind; its homogeneity
preserves this expression. The constants \(C_m\) depend on tenure, preferences,
and old-age prices and matter for tenure choice.

Define \(B=\beta K\) and \(D=1+\alpha+\vartheta+B\).
The expression \(\chi+p\kappa\) is the goods-and-space cost of a child.
If the young owner's housing cap is also slack at the relaxed optimum,
the solution without the mortgage limit is:
\[
\begin{aligned}
c^*&=\frac{M}{D}\left(1+\frac{\vartheta\chi}{\chi+p\kappa}\right),\\
h^*&=\frac{M}{D}\left(\frac{\alpha}{p}
                 +\frac{\vartheta\kappa}{\chi+p\kappa}\right),\\
n^*&=\frac{\vartheta M}{D(\chi+p\kappa)},\qquad
qz^*=\frac{BM}{D}.
\end{aligned}
\]
These expressions follow by allocating lifetime resources across adult
consumption, adult space, children, and old age.

Consequently the mortgage limit strictly restricts the household if and only
if:
\[
\boxed{\quad
\frac{w+qv}{w}>
\frac{1+\alpha+\vartheta+\beta(1+\gamma+\omega_B)}
{1+\alpha L/p+\vartheta(\chi+L\kappa)/(\chi+p\kappa)}.
\quad}
\]
The household has enough lifetime resources to want an allocation whose
current consumption and down payment exceed its current resources. This
is an explicit income-timing condition at given prices. It accommodates
heterogeneous \((y^y,b,y^o)\) and does not impose an upper bound on \(\beta\).
Prices and rebates still have to be determined in equilibrium.

**Housing rather than borrowing alone.** A binding mortgage limit need not,
by itself, imply \(h<h^*\). A simple sufficient condition that does is
\(h^*>w/L\): the desired house alone exceeds the maximum affordable size.
Another useful sufficient case is:
\[
L\ge p\quad\Longleftrightarrow\quad qP_{t+1}\ge\phi_tP_t.
\]
In the uncapped regime just described, strict financial restriction then
implies \(h<h^*\). The cash required for a unit of housing is at least its
lifetime service cost, so the restriction reduces housing relative to the
allocation with unrestricted finance. At constant prices this condition is
\(q\ge\phi\). It is sufficient, not a condition for every possible housing
restriction. Its empirical relevance depends on the period length and loan
repayment convention; it must not be assumed solely to obtain the desired sign.

Under the same condition, at a strictly constrained owner allocation with
slack young and old housing caps and \(v\ge0\), a small increase in \(\phi\)
strictly increases housing. This statement allows fertility to adjust and
holds prices, incomes, and rebates fixed.

## 4. Renting and ownership

For renters the lifetime budget has the same form, but the financial
restriction is \(z\ge v\), and housing is bounded by \(h_R^{\max}\):
\[
c+ph+qz=M,\qquad z\ge v,\qquad h\le h_R^{\max}.
\]
Renting avoids the owner's down payment, but cannot provide space beyond
the rental size limit. The rental problem and its value must therefore be
compared with the constrained owner problem, not with an assumed owner
allocation.

If both renter housing caps are slack, the unconstrained formulas above also
apply to renting, with its own old-age constant \(C_R\). The saving restriction
is slack when \(BM/D\ge qv\). Otherwise \(a'=0\), old-age resources equal
\(v\), and the young household allocates its current resources as follows:
\[
c^R=\frac{w}{E}\left(1+\frac{\vartheta\chi}{\chi+p\kappa}\right),\qquad
h^R=\frac{w}{E}\left(\frac{\alpha}{p}
                  +\frac{\vartheta\kappa}{\chi+p\kappa}\right),\qquad
n^R=\frac{\vartheta w}{E(\chi+p\kappa)},
\quad E=1+\alpha+\vartheta.
\]
If this housing choice exceeds the rental cap, solve the same problem with
\(h=h_R^{\max}\). If an old-age cap binds, use the old-age value function
in Section 2; the log-homogeneous formulas no longer apply.

**Ownership response.** At fixed prices and rebates, a higher \(\phi\) expands
the owner's feasible set and leaves the renter's problem unchanged. Therefore
\(W^O\) and the ownership probability weakly increase. With a strict financial
restriction and a positive marginal value of relaxing it, the conditional
owner value increases strictly. The household continues to compare this
value with renting, including its ownership taste.

**Role of rental segmentation.** If a household seeks housing above
\(h_R^{\max}\), renting cannot reproduce that housing choice. This identifies
the joint role of mortgage finance and the rental menu. Whether the household
buys a smaller home or remains a renter follows from the value comparison.
Without a minimum owner size, this model restricts the size of a purchase;
it does not make every owner home infeasible for a household with positive
current resources.

## 5. Changes requiring an author decision

| Feature | Existing note | Proposal |
|---|---|---|
| Current income at purchase | Excluded from the separate down-payment test | Available for consumption and down payment |
| Income when old | No separate endowment | Add heterogeneous, known \(y^o\) |
| Young mortgage restriction | Net-asset limit plus a separate \(b\)-only test | Ordinary mortgage limit, expressed in net assets |
| Old owner's housing | Can retain or downsize inherited home | Can resize within the owner cap, without borrowing |
| Tenure | Chosen when young | Retained |
| Fertility and estate preferences | Goods and space costs; warm-glow estate | Retained |
| Physical rental and owner limits | \(h_R^{\max}<h_O^{\max}\) | Retained |

The extra old-age income is a lifetime-income endowment, not a transfer policy.
The loan's maturity still spans one model period. These choices need to be
acceptable economically before proceeding to an equilibrium welfare claim.

The next work is to define the planner's intertemporal resource constraints
and check whether a housing reallocation improves welfare, then derive the
fertility response including its financing. None of the household results
above alone establishes that equilibrium housing is misallocated.

## Appendix: analytical checks

Write \(x=c-\chi n\) and \(s=h-\kappa n\). In the uncapped old-age regime,
the young objective differs by a constant from
\(\log x+\alpha\log s+\vartheta\log n+B\log(qz)\). Its constraints are:
\[
x+ps+(\chi+p\kappa)n+qz=M,\qquad
x+Ls+(\chi+L\kappa)n\le w.
\]
The relaxed problem allocates fractions \(1/D,\alpha/D,\vartheta/D,B/D\)
of \(M\) to these four expenditures, giving the closed-form allocation and
the binding test in Section 3.

For a binding mortgage limit, let \(\lambda>0\) and \(\mu>0\) be multipliers
on lifetime resources and current cash, respectively. The first-order
conditions are:
\[
x=\frac1{\lambda+\mu},\quad
s=\frac{\alpha}{\lambda p+\mu L},\quad
n=\frac{\vartheta}{\lambda(\chi+p\kappa)+\mu(\chi+L\kappa)},\quad
qz=\frac B\lambda.
\]
The two binding resource equations determine \(\lambda,\mu\). Strict
concavity makes a feasible solution sufficient and unique. With housing caps,
retain the corresponding constraint and multiplier; the two-budget reduction
and the relaxed-allocation binding test continue to hold.

For the housing comparison, set \(r=\mu/(\lambda+\mu)\) and
\(\rho=p+r(L-p)\). Define
\(A(\rho)=1+\vartheta\chi/(\chi+\kappa\rho)\) and
\(H(\rho)=\alpha/\rho+\vartheta\kappa/(\chi+\kappa\rho)\).
The first-order conditions and lifetime budget imply:
\[
\frac hM=\left[p+\frac{A(\rho)}{H(\rho)}
                         +\frac{B}{(1-r)H(\rho)}\right]^{-1}.
\]
The function \(H\) is decreasing, and \(A/H\) is increasing on positive
prices. If \(L\ge p\), then \(\rho\ge p\), while \(r>0\); the denominator
is strictly larger than at \(r=0\). Hence \(h<h^*\).

Finally, the homogeneity of the old-age problem can be checked directly.
Without an old-age size cap, an owner with a slack estate floor chooses
\(c^2=z/K\), \(h^2=\gamma z/(Ku)\), and \(e=\omega_Bz/(Kq)\).
The floor is nonrestricting when \(\omega_Bu\ge q\gamma P_{t+1}\), and strictly
slack when the inequality is strict. With a strictly binding
floor, the choices are
\(c^2=z/K\), \(h^2=(\gamma+\omega_B)z/[K(u+qP_{t+1})]\), and
\(e=P_{t+1}h^2\). Both cases give \(K\log z+C_O\).
The displayed young-household formulas require their implied old-age housing
to satisfy the physical cap, in the constrained and comparison allocations.

### Local housing response

With binding finance, substitute \(c=w-Lh\) and
\(S=qz=qv+(L-p)h\). Apart from a constant, the objective becomes:
\[
\log(w-Lh-\chi n)+\alpha\log(h-\kappa n)+\vartheta\log n
 +B\log\{qv+(L-p)h\}.
\]
Let \(a=-U_{hh}>0\), \(b_0=U_{hn}\), and \(d_0=-U_{nn}>0\) for this
reduced objective. With \(x=w-Lh-\chi n\), \(s=h-\kappa n\), and
\(\Delta=L-p\), differentiation gives:
\[
\begin{aligned}
a&=L^2/x^2+\alpha/s^2+B\Delta^2/S^2,\\
b_0&=-L\chi/x^2+\alpha\kappa/s^2,\\
d_0&=\chi^2/x^2+\alpha\kappa^2/s^2+\vartheta/n^2,\\
\frac{\partial h}{\partial\phi}
&=\frac{P}{ad_0-b_0^2}
\left[d_0\left(\frac{w-\chi n}{x^2}-\frac{Bqv}{S^2}\right)
                   +\frac{b_0\chi h}{x^2}\right].
\end{aligned}
\]
Strict concavity gives \(ad_0-b_0^2>0\). Define
\(J=L+b_0\chi/d_0\); substitution shows
\(J=[L\vartheta/n^2+\alpha\kappa(\chi+L\kappa)/s^2]/d_0>0\).
The sign of the derivative is therefore the sign of
\((x+hJ)/x^2-Bqv/S^2\). Strict financial restriction gives
\(B/S<1/x\), so for \(qv\ge0\):
\[
\frac{x+hJ}{x^2}-\frac{Bqv}{S^2}
\ge \frac{hJ}{x^2}+\frac{h\Delta}{xS}>0
\qquad\text{when }\Delta\ge0.
\]
This proves the stated sufficient sign with endogenous fertility. Caps must
remain slack in a neighborhood for this derivative formula to apply.

### Why a binding mortgage limit is not enough for a housing sign

The following exact counterexample concerns a household at given prices;
it does not assert a market-clearing equilibrium. Set
\(q=1/2\), constant house prices \(P=2\), \(\tau^p=0\),
\(\phi=19/20\), \(\chi=\kappa=1\), \(\alpha=31/40\),
\(\vartheta=71/40\), \(\beta=1/2\), \(\gamma=1/4\),
\(\omega_B=3/4\), \(w=11/5\), and \(v=94/15\).
Let \(h_O^{\max}=3\) and \(h_R^{\max}=1\). Then \(p=1\), \(L=1/10\),
and \(B=1\). The strictly constrained owner chooses:
\[
x=s=n=1,\qquad c=h=2,\qquad qz=4/3,\qquad
\lambda=3/4,\quad\mu=1/4.
\]
These choices satisfy both budgets and all four first-order conditions.
Old choices are \(c^2=4/3\), \(h^2=1/3\), \(e=2\); the estate floor and
housing cap are slack. Strict concavity establishes the conditional optimum.
The rental housing cap prevents the renter from attaining this allocation;
with zero ownership taste, ownership is strictly preferred here.

Removing the mortgage limit gives \(h^*=76/39<2\). Implicit differentiation
of the two budgets and four first-order conditions also gives:
\[
\frac{\partial h}{\partial\phi}=-\frac{8280}{25271}<0.
\]
The borrowing limit restricts how the household spends lifetime resources,
but this household already devotes more of them to housing than it would
with unrestricted finance. A loan against housing is available while a loan
against income alone is not. The example explains why the direction of the
housing response needs its own condition.
