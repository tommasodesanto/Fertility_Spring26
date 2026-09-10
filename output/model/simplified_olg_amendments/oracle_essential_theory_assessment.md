# Theory decision brief

**September 10.** We have a conditional result for the proposed two-period
model: at fixed fertility, the full planner gives the young more housing and
can take consumption from them. When it also chooses fertility, fertility rises
if children's needs are sufficiently housing intensive. Its allocation is
explicit. Pro and independent analytical reviews support the result. The income restriction remains substantive, and
this does not yet establish a fertility effect of policy along a transition.

## 1. The economy

Young households differ in resources \(w\) and receive common retirement income
\(y>0\). Their preferences and old households' preferences are

\[
u^y=\log(c^y-\chi n)+\alpha\log(h^y-\kappa n)+\vartheta\log n,
\qquad
u^o=\log c^o+\alpha\log h^o+\omega_B\log e.
\]

Children require \(\chi\) goods and \(\kappa\) housing; \(e\) is the
estate. Households discount old-age utility at \(\beta\). Rentals have maximum
size \(r\); owners face a mortgage limit. The proposed financing rule gives
rent \(p=(1+\tau-q)P\) and an equity requirement
\(q(z-y)\geq\Gamma ph\), where
\(\Gamma=q(1-\phi)/(1+\tau-q)>0\). Here \(P\) is the house price,
\(q\) prices next-period goods, \(\tau>0\) is property tax, \(\phi<1\)
is the financed share, and \(z\) is resources upon retirement. The budgets
and mortgage payments are specified below.

## 2. What the planner does

Start with a stationary competitive allocation, with equal young and old cohort
masses. Let \(\mathcal C\) and \(\mathcal H\) denote its total consumption
goods and housing per young household, including both ages. The planner gives
equal weight to each living household's current utility and chooses all current
consumption and housing. It holds individual fertility, continuation resources
and estates fixed and can relax private finance. It preserves the rental-size
limit by assigning ownership to larger homes. Because continuation utility and
estates are fixed, maximizing current utility also maximizes the equally
weighted sum of remaining lifetime utilities in this comparison.

Its problem and complete solution, with \(\bar n\) denoting competitive mean
fertility, are

\[
\max_{\{c^y,h^y,c^o,h^o\}}
\{\overline{u^y}+\overline{u^o}\}
\quad\text{subject to}\quad
\bar c^y+\bar c^o=\mathcal C,\qquad
\bar h^y+\bar h^o=\mathcal H;
\]
\[
\boxed{
\begin{gathered}
X=\frac{\mathcal C-\chi\bar n}{2},\qquad
S=\frac{\mathcal H-\kappa\bar n}{2},\\
(c_i^{y,F},h_i^{y,F})=(X+\chi n_i,S+\kappa n_i),\\
(c_j^{o,F},h_j^{o,F})=(X,S).
\end{gathered}}
\]

The planner divides adult goods and space equally, then adds each family's
children's needs. Each parent therefore receives \(\kappa n_i\) more
housing than each old household in the planner allocation. To establish a
reallocation **relative to equilibrium**, we must also establish that the young
initially receive too little adult space as a group.

## 3. When housing moves toward the young

Put \(K=1+\alpha+\omega_B\), \(J=1+\alpha+\vartheta\). The following
conditions are sufficient:

- Old households can finance their unrestricted housing and estate choices:
  \(\omega_B\geq\alpha\Gamma\).
- A positive share of young households wants to borrow against retirement
  income and wants a freely sized rental larger than \(r\):

\[
\boxed{
w<\frac{Jqy}{\beta K},\qquad
\frac wJ\left(\frac\alpha p+
\frac{\vartheta\kappa}{\chi+\kappa p}\right)>r.
}
\]

- The retirement-income floor satisfies

\[
\boxed{\frac yK\geq\frac{\bar w}{J}.}
\]

**Proposition. Under these conditions, the full planner allocates more aggregate
housing to the young.** Young consumption need not rise. The comparison covers
optimal young renters and owners; it does not assume their tenure choices or
the desired young–old housing ranking. The second condition has a wider exact
income interval in the supporting proof.

The argument is short. The housing restrictions imply that affected young
households value extra space more than the old, in consumption units:

\[
\frac{\alpha(c_i^y-\chi n_i)}{h_i^y-\kappa n_i}
>p=\frac{\alpha c_j^o}{h_j^o}.
\]

This is the June argument's core: the young could compensate the old for a
small housing transfer and still gain. The extra income restriction makes the
**full maximizing allocation** move housing toward the young overall. It
compares a floor on old consumption with young current resources divided by
their utility weights. It is conservative because it omits accumulated old
wealth and young saving. We have not shown that it is empirically mild.

The exact mean changes under the full planner are

\[
\bar h^{y,F}-\bar h^{y,E}
=\frac{\bar h^o-\overline{(h^y-\kappa n)}}2>0,
\qquad
\bar c^{y,F}-\bar c^{y,E}
=\frac{\bar c^o-\overline{(c^y-\chi n)}}2.
\]

An individual parent gains housing precisely when its initial adult space is
below \(S\). The aggregate theorem does not say that every financially
constrained parent receives more under equal-weight redistribution.

## 4. What follows for fertility

Now let the same planner choose parental fertility as well, valuing existing
parents' utility. Its complete allocation has common fertility \(n^F\) and

\[
(c^{y,F},h^{y,F})=
\left(\frac{\mathcal C+\chi n^F}{2},
      \frac{\mathcal H+\kappa n^F}{2}\right),\qquad
(c^{o,F},h^{o,F})=
\left(\frac{\mathcal C-\chi n^F}{2},
      \frac{\mathcal H-\kappa n^F}{2}\right).
\]

It chooses \(n^F\) to maximize
\(2\log(\mathcal C-\chi n)+2\alpha\log(\mathcal H-\kappa n)
+\vartheta\log n\). The solution is the feasible root of a quadratic,
written explicitly in the supporting proof. Under the allocation conditions,

\[
\boxed{\kappa p\geq\alpha\chi\quad\Longrightarrow\quad n^F>\bar n.}
\]

The condition says that children require housing at least as intensively as
unrestricted adults, in expenditure terms. **Each parent voluntarily chooses
\(n^F\) at the final joint-planner bundle.** Its private fertility condition
is exactly the planner's condition there. The planner can therefore be described
as assigning fixed consumption and housing while anticipating parents' choices.
The assigned bundles are fixed before those choices and are not conditional on
realized births.

A compensated local housing transfer also raises an affected parent's private
fertility under this condition. Mean private fertility can fall if households
instead receive bundles calculated by first holding fertility at its competitive
level. We have an exact counterexample for that different
experiment. None of these direct allocations establishes implementation through
a cash-transfer policy with freely reoptimizing housing markets.

## 5. The morning decision

**My recommendation is to present this as a conditional static allocation
exercise, with the full planner solution at its center.** The economics and
proof are complete for that claim in the proposed model. The income restriction
and planner's ability to relax finance should be visible. Neither children nor
a binding mortgage limit alone guarantees the desired aggregate reallocation.

The remaining decision is whether to adopt this proposed financing and closure:
two periods, a common divisible housing stock, mortgage interest paid during the
age, common positive retirement income, external finance, and entrant endowments
independent of estates. Property taxes fund fixed public services here. These
are changes from parts of the earlier model; the calculation does not validate
those older versions.

**We do not yet have a theorem that a specified housing policy raises fertility
along an equilibrium transition.** The static result can introduce that question;
it cannot stand in for its answer. No empirical mildness or global equilibrium
uniqueness is established. The proof, complete planner formula, proposed six
slides and verification are below; the main deck and manuscript are unchanged.

<details>
<summary>Detailed model, verified results, and work record</summary>

# Housing allocation and fertility: the shorter argument

**September 10, overnight review.** The main result now covers young renters and
owners without prescribing which type chooses each tenure. It allows the planner
to give the young more housing while taking some consumption from them. The
household, allocation and joint-fertility arguments have passed separate
analytical cross-checks. The completed [Pro review](oracle_all_tenure_simplification_response.md)
confirms the all-tenure theorem and the private-fertility counterexample. The
full joint planner allocation is also explicit, as shown below.
This is a proposed theory statement; the main deck and manuscript are unchanged.

The result remains conditional. The borrowing and rental restrictions explain
why a housing transfer can improve the allocation. A separate restriction on
resources determines whether the **full equally weighted planner** gives more
housing to the young as a group. Children alone do not establish that direction.

## The local argument to lead with

**Author clarification, September 10.** Recover the simple June–July argument
about marginal housing values and build the exposition from it. This is a
change in emphasis, not an instruction to abandon the full planner or fertility
work. The June 9 source is
[the local allocation note](../../../latex/archive/theory_note_cleanup_20260622/local_allocation_result_note_20260609.tex),
especially its local reallocation experiment and fertility condition. The later
[short note](../../../latex/archive/theory_note_cleanup_20260622/intergen_housing_fertility_short_note_corina_plus_planner.tex)
contains the same argument with richer tenure and tax accounting. These sources
are historical evidence, not certifications of all their equations.

The useful core measures housing values in units of consumption. In the
corrected model, the income-and-access conditions below imply that a young
household's valuation exceeds the old household's unrestricted valuation:

\[
M_i^y=\frac{\alpha(c_i^y-\chi n_i)}{h_i^y-\kappa n_i}
>p=
\frac{\alpha c_j^o}{h_j^o}=M_j^o.
\]

Move a little housing from old to young and transfer enough goods back to
compensate the old. The young value the extra space more than the compensation
costs. This is a feasible improvement when the authority can relax finance,
with individual continuation resources and estates preserved. It therefore
proves that the competitive allocation does not maximize the full equally
weighted dated objective. The planner is still allowed to choose consumption
and housing; exhibiting one improving direction suffices to show the original
allocation is not its optimum.

The average-income bound in the proposition is needed only for the additional
claim that the **maximizing allocation gives the young more housing in
aggregate**. It is not needed for this compensated local gain. The old
unrestricted-finance condition and the income/access conditions establishing a
strict young distortion remain essential and must be shown rather than assumed.
Historical tax/retention terms must not be imported into the repaired budgets
without rederiving their fiscal and estate counterparts.

The same marginal valuation enters the parental fertility condition:

\[
\frac{\vartheta(c_i^y-\chi n_i)}{n_i}=\chi+\kappa M_i^y.
\]

Children require space, so limited housing access increases its private
marginal cost. A signed fertility response still requires the specified
resource change; the compensated-transfer result below supplies one explicit
condition. Lead the next Pro follow-up and the final exposition with this local
argument, then make the full-allocation and joint-fertility statements central,
as the author requested in his subsequent clarification. Pro has completed its
review; no further broad prompt is pending.

## Model

Young households differ in resources \(w>0\), and all receive income \(y>0\)
in old age. There are two periods. Preferences are

\[
u^y=\log(c^y-\chi n)+\alpha\log(h^y-\kappa n)+\vartheta\log n,
\qquad
u^o=\log c^o+\alpha\log h^o+\omega_B\log e.
\]

Here \(n\) is fertility, \(\chi\) and \(\kappa\) are children's goods and
space needs, and \(e\) is the net estate. Households discount old-age utility
at \(\beta\). All utility and child-cost coefficients are positive.

Housing comes from a common divisible stock. Rentals cannot exceed \(r\).
Owners may choose any size. Let \(P\) be the house price, \(q\) the price of
next-period goods, \(\phi\) the financed share, and \(\tau>0\) the property
tax. Under the proposed serviced-interest convention, owners satisfy

\[
c+(1+\tau)Ph+a=w,\qquad
z=y+a/q+Ph,\qquad a\geq-q\phi Ph.
\]

If saving is \(s\geq0\) and mortgage principal is \(L\leq\phi Ph\), then
\(a=s-qL\). The present cost of mortgage interest payments is \((1-q)L\),
paid from young resources; principal \(L\) is repaid at sale. This gives the
budget above. Retirement resources are \(z\), and owners must retain equity
\(z-y\geq(1-\phi)Ph\). Renters satisfy \(c+ph+a=w\),
\(z=y+a/q\), \(a\geq0\), and \(h\leq r\).

An old owner has the same financial budget, with its net estate as the payoff:

\[
c^o+(1+\tau)Ph^o+a^o=z^o,\qquad
\boxed{e=a^o/q+Ph^o},\qquad a^o\geq-q\phi Ph^o.
\]

The estate consists of net financial assets and the house sale proceeds at
death. An old renter instead has \(c^o+ph^o+a^o=z^o\),
\(e=a^o/q\), \(a^o\geq0\), and \(h^o\leq r\). The bequest weight
\(\omega_B\) determines how much the old value the estate they leave.

Competitive finance determines stationary rent and the equity requirement:

\[
p=(1+\tau-q)P,\qquad
\Gamma=\frac{q(1-\phi)}{1+\tau-q},\qquad
q(z-y)\geq\Gamma ph.
\]

Property taxes finance public services for this calculation. External finance
and entrants' endowments independent of parental estates remain the proposed
closure. A stationary competitive equilibrium satisfies household optimality,
housing clearing and replacement, \(\nu\bar n=1\), where \(\nu n\) is the
number of entering households per parent. Old resources are generated by the
preceding young cohort's saving and equity, so \(z^o\geq y\).

Write \(K=1+\alpha+\omega_B\) and \(J=1+\alpha+\vartheta\). We maintain
\(\omega_B\geq\alpha\Gamma\), which makes the unrestricted old allocation
feasible:

\[
c^o=\frac{z^o}{K},\qquad
h^o=\frac{\alpha z^o}{Kp},\qquad
qe=\frac{\omega_Bz^o}{K}.
\]

Old households can sell and resize; this condition concerns their desired
housing–estate mix, rather than an inability to sell.

## Allocation statement

The planner chooses all current consumption and housing, giving equal weight to
each living household. It preserves each young household's continuation
resources, each old household's estate, continuation prices, and the current
goods and housing totals. Initially it holds each parent's competitive fertility
fixed. It can relax private finance and assign ownership to homes exceeding
the rental ceiling.

For an elementary sufficient test, define the home a young household would rent
with no saving and no rental ceiling:

\[
h_0(w,p)=\frac{w}{J}
\left[\frac\alpha p+
\frac{\vartheta\kappa}{\chi+\kappa p}\right].
\]

The following restrictions identify households that want to bring retirement
income forward and want more space than the largest rental:

\[
\boxed{\quad w<\frac{Jqy}{\beta K},\qquad h_0(w,p)>r.\quad}
\tag{Housing access}
\]

They imply a strict housing distortion whether the household ultimately rents
or owns. They are sufficient income-and-price tests, not assumptions that the
equilibrium already has a particular housing ordering. The exact, wider income
interval is given below.

**Proposition.** Consider any stationary competitive equilibrium of this economy.
Suppose a positive share of young households satisfies the housing-access test.
A small transfer of housing from an old household to any of these young
households can compensate the old in goods and make the young strictly better
off, after relaxing finance.

If average young resources also satisfy

\[
\boxed{\qquad \frac yK\geq\frac{\bar w}{J},\qquad}
\tag{Resources}
\]

the full equally weighted planner assigns **more aggregate housing to the
young**. If, in addition,

\[
\boxed{\qquad \kappa p\geq\alpha\chi,\qquad}
\tag{Children's needs}
\]

the joint planner choosing consumption, housing and fertility chooses **higher
mean fertility**. The objective values existing parents' utility.

The resource condition compares a floor on old consumption with young current
resources per unit of utility weight. It is conservative because it omits old
accumulated wealth and young saving. The child-cost condition says that
children's housing-to-goods expenditure ratio is at least the unrestricted
adult ratio \(\alpha\). Neither condition orders \(\beta\) and \(q\).
The housing-access and child-cost tests still involve equilibrium rent;
the resource bound is solely a restriction on primitives.

## Short proof

Let \(x_i=c_i^y-\chi n_i\) be adult consumption, and express adult space in
goods as \(v_i=p(h_i^y-\kappa n_i)/\alpha\). Private optimality gives
\(x_i\geq v_i\), strictly under the housing-access test. Unrestricted old
households have \(ph_i^o/\alpha=c_i^o\). Thus a constrained young household
values extra space more highly in consumption units than an old household:

\[
\frac{\alpha x_i}{h_i^y-\kappa n_i}>p
=\frac{\alpha c_j^o}{h_j^o}.
\]

A small transfer costs the old \(p\) goods per unit to first order and leaves
a strictly positive surplus for the young. The exact compensation appears
below.

For the full planner, the useful accounting identity concerns current spending.
Define \(\widehat x_i=(c_i^y+ph_i^y)/J\). The budget and household first-order
conditions imply

\[
v_i\leq\widehat x_i\leq x_i,\qquad
\widehat x_i=\frac{w_i-q(z_i-y)}J\leq\frac{w_i}{J}.
\]

The first two inequalities are strict for a household with a strict housing
distortion. If \(B=\bar c^o\), the resource condition gives
\(B\geq y/K\geq\bar w/J\geq\overline{\widehat x}>\bar v\).

The planner equalizes adult consumption and adult space. For equal cohort
masses its allocation is

\[
X=\frac{\bar x+B}{2},\qquad
S=\frac{\alpha}{p}\frac{\bar v+B}{2},\qquad
(c_i^{y,F},h_i^{y,F})=(X+\chi n_i,S+\kappa n_i),\qquad
(c_j^{o,F},h_j^{o,F})=(X,S).
\]

Consequently the mean young housing gain is
\(\alpha(B-\bar v)/(2p)>0\). These formulas also show why an individual
young household need not gain, and why young consumption need not rise.
In fact, the weaker condition \(B\geq\overline{\widehat x}\) suffices.
At equality, \(\bar x>\overline{\widehat x}=B\), so young consumption falls.

This makes the economic roles of redistribution and housing constraints
explicit. The mean housing change decomposes as

\[
\bar h^{y,F}-\bar h^{y,E}
=\frac{\alpha}{2p}
\left[\underbrace{B-\overline{\widehat x}}_{\text{resource difference}}
+\underbrace{\overline{\widehat x}-\bar v}_{\text{housing distortion}}\right].
\]

The second term is strictly positive when some young households are distorted.
The resource assumption makes the first term nonnegative. At
\(B=\overline{\widehat x}\), the planner gives the young more housing and
takes consumption from them. If old resources are sufficiently low, the first
term can overturn the second. That is why the full allocation's direction needs
more than the local valuation gap.

An individual parent gains housing exactly when \(s_i<S\), and gains
consumption exactly when \(x_i<X\). The constrained group and the group
receiving housing under equal-weight redistribution need not coincide.

The primitive resource floor has a direct interpretation in income ratios:

\[
\frac{y}{\bar w}\geq
\frac{1+\alpha+\omega_B}{1+\alpha+\vartheta}
=1+\frac{\omega_B-\vartheta}{J}.
\]

If the bequest weight does not exceed the fertility weight, retirement income
at least as large as mean young resources is sufficient. If bequests receive
more weight, the required floor is higher because the old devote more resources
to estates. The two old-resource restrictions can hold together precisely when

\[
\alpha\Gamma\leq\omega_B
\leq\frac{Jy}{\bar w}-1-\alpha.
\]

These are statements about the endowments and utility weights of this two-age
exercise; no direct mapping to annual retirement income has been established.
High patience does not automatically rule out affected young households. Their
simple borrowing test is \(\beta<qJy/(Kw_i)\). Under the resource floor,
this upper bound is at least \(q\bar w/w_i\), which exceeds \(q\) for a
below-average-resource household. The rental-size test must still hold: the
very poorest may want a home small enough to rent without restriction.

For fertility, put \(C=\chi+\kappa p\). The fertility first-order condition
and the child-cost restriction give the stronger comparison

\[
v_i\leq\frac{Cn_i}{\vartheta}\leq\widehat x_i\leq x_i,
\qquad
n_i=\vartheta\left(\frac\chi{x_i}+\frac{\kappa p}{v_i}\right)^{-1}.
\]

The middle inequality is strict for a strictly distorted household, so
\(B/C>\bar n/\vartheta\). The function
\(g(x,v)=(\chi/x+\kappa p/v)^{-1}\) is concave. With
\(V=pS/\alpha=(\bar v+B)/2\), Jensen's inequality yields

\[
g(X,V)\geq\frac{g(\bar x,\bar v)+g(B,B)}2
\geq\frac{\bar n/\vartheta+B/C}{2}
>\frac{\bar n}{\vartheta}.
\]

This says that, after allocating consumption and housing optimally, the
marginal value of increasing fertility is positive at the competitive mean.
Strict concavity of the joint planner's objective gives \(n^F>\bar n\).
The formulas generalize directly to any predetermined positive ratio of old
to young households.

## The complete joint allocation

Let \(\mu>0\) be the predetermined ratio of old to young households, and let
\(\mathcal C,\mathcal H\) include both age groups' current resources per
young household. At the full joint optimum every parent chooses the same
fertility. All households receive common adult goods \(X\) and space \(S\):

\[
X(n)=\frac{\mathcal C-\chi n}{1+\mu},\qquad
S(n)=\frac{\mathcal H-\kappa n}{1+\mu},\qquad
\frac{\vartheta}{n}=\frac{\chi}{X(n)}+\frac{\alpha\kappa}{S(n)}.
\]

The last equation equates the marginal utility of children to the goods and
space they require, evaluated at the planner's marginal resource values.
The objective is strictly concave and tends to minus infinity at each end of
\(0<n<\min\{\mathcal C/\chi,\mathcal H/\kappa\}\). It therefore has
one interior solution. The resulting complete bundles are

\[
\boxed{
\begin{aligned}
c^{y,F}&=\frac{\mathcal C+\mu\chi n^F}{1+\mu},&
h^{y,F}&=\frac{\mathcal H+\mu\kappa n^F}{1+\mu},\\
c^{o,F}&=\frac{\mathcal C-\chi n^F}{1+\mu},&
h^{o,F}&=\frac{\mathcal H-\kappa n^F}{1+\mu}.
\end{aligned}}
\]

Young continuation resources and old estates remain individually fixed. The
planner equalizes current adult resources; it does not equalize future wealth.
The comparison with competitive consumption is

\[
\bar c^{y,F}-\bar c^{y,E}
=\frac{\mu}{1+\mu}
\left[B-\bar x+\chi(n^F-\bar n)\right].
\]

Consequently the consumption loss established for a fixed-fertility allocation
does not by itself establish a loss after the joint fertility choice. For
housing the corresponding expression adds
\(\mu\kappa(n^F-\bar n)/(1+\mu)>0\) to the fixed-fertility gain.

<details>
<summary>Explicit fertility formula and feasibility proof</summary>

Write \(M=1+\mu\), \(D=\vartheta+M(1+\alpha)\), and define the positive
coefficient

\[
\mathcal B=\kappa\mathcal C(\vartheta+M\alpha)
+\chi\mathcal H(\vartheta+M).
\]

Multiplying the first-order condition by its positive denominator gives
\(\chi\kappa Dn^2-\mathcal Bn+\vartheta\mathcal C\mathcal H=0\).
The unique feasible root is

\[
\boxed{
n^F=\frac{2\vartheta\mathcal C\mathcal H}
{\mathcal B+\sqrt{\mathcal B^2-
4\chi\kappa D\vartheta\mathcal C\mathcal H}}.
}
\]

To verify the root, put \(N_c=\mathcal C/\chi\) and
\(N_h=\mathcal H/\kappa\). The normalized polynomial is

\[
Q(n)=Dn^2-[(\vartheta+M\alpha)N_c+(\vartheta+M)N_h]n
+\vartheta N_cN_h.
\]

It satisfies \(Q(0)>0\), \(Q(N_c)=MN_c(N_c-N_h)\), and
\(Q(N_h)=M\alpha N_h(N_h-N_c)\). When the resource bounds differ,
the smaller root is below both; the larger root lies between them and gives a
negative adult resource. When both bounds equal \(N\), the roots are
\(\vartheta N/D\) and \(N\), and the second leaves nothing for adults.
The discriminant is positive because it equals

\[
[\kappa\mathcal C(\vartheta+M\alpha)
-\chi\mathcal H(\vartheta+M)]^2
+4M^2\alpha\chi\kappa\mathcal C\mathcal H.
\]

In the proportional-resource case \(\mathcal C/\chi=\mathcal H/\kappa=N\),
the complete allocation simplifies to

\[
n^F=\frac{\vartheta N}{D},\qquad
(c^{y,F},h^{y,F})=\frac{J}{D}(\mathcal C,\mathcal H),\qquad
(c^{o,F},h^{o,F})=\frac{1+\alpha}{D}(\mathcal C,\mathcal H).
\]

This special case retains positive goods and housing costs of children. The
condition \(\kappa p=\alpha\chi\) alone does not imply proportional
aggregate resources. The closed form solves the full conditional planner;
it does not imply a closed form for every heterogeneous competitive equilibrium.

</details>

## Voluntary fertility at the optimal allocation

The joint optimum can be implemented by giving each parent its final fixed
consumption and housing totals and then letting it choose fertility. At the
assigned bundle \((X^F+\chi n^F,S^F+\kappa n^F)\), its private condition is

\[
\frac{\vartheta}{n^F}-\frac{\chi}{X^F}
-\frac{\alpha\kappa}{S^F}=0.
\]

This is the planner's condition. The private objective is strictly concave in
fertility, so the parent's unique choice is \(n^F\). The transfers do not
respond to actual births: \(n^F\) is computed before assigning fixed totals.
Continuation resources are unchanged.

This also establishes equivalence of the two planning problems. Allocations
induced by fixed bundles and voluntary fertility form a subset of the joint
planner's feasible allocations. The joint optimum itself belongs to that subset.
Consequently choosing bundles while anticipating parental fertility attains
exactly the joint maximum. The planner and parent agree because a fertility
change at a fixed bundle reallocates that parent's resources between adults and
children, and the welfare criterion values the existing parent's utility.

The result concerns direct allocations of goods and housing. A cash transfer
followed by free housing, tenure and saving choices is a different experiment.

## Which fertility comparison is established?

The joint planner chooses fertility together with parents' resources. That
result does **not** imply that mean fertility rises if the planner first
allocates bundles based on existing fertility and then lets every parent
rechoose fertility at its assigned total consumption and housing.

For that second experiment, an individual parent's fertility rises exactly when

\[
n_i<n_c,\qquad
n_c=\frac{\vartheta}{\chi/X+\alpha\kappa/S}.
\]

The proposition implies \(n_c>\bar n\), so every parent with fertility at or
below the competitive mean increases fertility. Under the child-cost condition,
a conservative income-only sufficient test is \(w_i\leq Jy/(2K)\).
Mean private fertility can nevertheless fall: an explicit stationary
counterexample satisfies all the proposition's conditions. This negative
result has also passed independent review.

There is a separate, direct private-fertility implication for the compensated
local housing transfer. Its sign is positive exactly when

\[
\kappa p\left(\frac{x_i}{v_i}\right)^2>\alpha\chi.
\]

Thus the weak child-cost condition suffices for any strictly distorted
recipient, despite the consumption needed to compensate the old. This is a
local allocation comparison, not a tax or mortgage-policy implementation.

<details>
<summary>Six-slide exposition for author review</summary>

These are proposed frames in the current model's notation. They have not been
inserted into the main deck. The graph frames should follow once their precise
comparisons are agreed; an old transition diagram is not evidence for a theorem
under the revised financing convention.

### 1. Economy

Two-period households differ in young resources \(w\) and receive common
retirement income \(y\). Children require \(\chi\) goods and \(\kappa\) space.
Preferences are

\[
u^y=\log(c^y-\chi n)+\alpha\log(h^y-\kappa n)+\vartheta\log n,
\qquad
u^o=\log c^o+\alpha\log h^o+\omega_B\log e.
\]

Households choose consumption, housing, saving, tenure and fertility, and
discount old-age utility at \(\beta\). The estate is \(e\).

### 2. Housing finance

Housing can be rented up to size \(r\), or owned with financed share
\(\phi<1\). At house price \(P\) and rent \(p\), the young satisfy

\[
\begin{array}{ll}
\text{Rent:}&c+ph+a=w,\quad z=y+a/q,\quad a\geq0,\quad h\leq r;\\[1ex]
\text{Own:}&c+(1+\tau)Ph+a=w,\quad z=y+a/q+Ph,\quad a\geq-q\phi Ph.
\end{array}
\]

Here \(q\) prices next-period goods, \(\tau\) is the property tax, and
\(a=s-qL\) is saving less discounted mortgage principal. Interest is paid
from current resources. Ownership requires net housing equity \(z-y\geq(1-\phi)Ph\); retirement income is unavailable for that equity.

### 3. Competitive equilibrium

A stationary equilibrium consists of household choices, prices and cohort mass
\(N\). Optimal choices clear the common stock \(\bar H\) and reproduce the
cohort:

\[
p=(1+\tau-q)P,\qquad
N(\bar h^y+\bar h^o)=\bar H,\qquad \nu\bar n=1.
\]

Here \(\nu n\) is the number of entering households per parent. Put
\(K=1+\alpha+\omega_B\), \(J=1+\alpha+\vartheta\), and
\(\Gamma=q(1-\phi)/(1+\tau-q)\). When
\(\omega_B\geq\alpha\Gamma\), old choices satisfy
\(c^o=z^o/K\), \(h^o=\alpha z^o/(Kp)\), and \(qe=\omega_Bz^o/K\).

### 4. Planner allocation

The planner gives equal weight to living households and reallocates all current
consumption and housing, holding fertility and individual future resources fixed:

\[
\begin{gathered}
\max_{\{c^y,h^y,c^o,h^o\}}
\left\{\overline{u^y(c^y,h^y,n)}+\overline{u^o(c^o,h^o,e)}\right\},
\qquad
\bar c^y+\bar c^o=\mathcal C,\quad
\bar h^y+\bar h^o=\mathcal H;\\[1ex]
X=\frac{\mathcal C-\chi\bar n}{2},\qquad
S=\frac{\mathcal H-\kappa\bar n}{2},\qquad
(c_i^{y,F},h_i^{y,F})=(X+\chi n_i,S+\kappa n_i),\quad
(c_j^{o,F},h_j^{o,F})=(X,S).
\end{gathered}
\]

Private finance can be relaxed. The planner divides adult resources equally,
then adds each family's children's needs; \(\mathcal C,\mathcal H\) include
both cohorts' resources per young household.

### 5. Housing reallocation

Suppose old households attain the unrestricted choice on frame 3 and a positive
share of young households satisfies the borrowing and rental-size tests below:

\[
\begin{gathered}
w<\frac{Jqy}{\beta K},\qquad
\frac wJ\left[\frac\alpha p+
\frac{\vartheta\kappa}{\chi+\kappa p}\right]>r,\qquad
\frac yK\geq\frac{\bar w}{J}\\[1ex]
\Longrightarrow\qquad
\bar h^{y,F}-\bar h^{y,E}
=\frac{\bar h^o-\overline{(h^y-\kappa n)}}2>0.
\end{gathered}
\]

The constraints raise young households' valuation of space. The income condition
ensures the full planner gives the young more housing overall; it can take
consumption from them in return.

### 6. Parental fertility

The planner chooses new fixed consumption and housing totals, anticipating
parents' fertility choices. With \(\mathcal C,\mathcal H\) unchanged, its
joint solution satisfies

\[
\begin{gathered}
\frac{\vartheta}{n^F}
=\frac{2\chi}{\mathcal C-\chi n^F}
+\frac{2\alpha\kappa}{\mathcal H-\kappa n^F},\\[1ex]
(c^{y,F},h^{y,F})=
\left(\frac{\mathcal C+\chi n^F}{2},
      \frac{\mathcal H+\kappa n^F}{2}\right),\\[1ex]
\kappa p\geq\alpha\chi\quad\Longrightarrow\quad n^F>\bar n^E.
\end{gathered}
\]

The conclusion uses the preceding allocation conditions. Children require
housing at least as intensively as unrestricted adults. The new totals are fixed
before parents choose fertility; parents voluntarily choose \(n^F\) there.

</details>

<details>
<summary>Exact household interval, proof details, and overnight work record</summary>

The elementary housing-access test in the proposition deliberately avoids
piecewise cutoffs. The exact affected-income interval is larger. Write
\(\delta=\beta K\), \(Y=qy\),
\(A(p)=\alpha+\vartheta\kappa p/(\chi+\kappa p)\), and
\(x_r=pr/A(p)\). Define

\[
w_r(p)=Jx_r+\max\{\delta x_r-Y,0\},\qquad
w_O(p)=
\begin{cases}
Y[J+\Gamma A(p)]/[\delta-\Gamma A(p)],&\delta>\Gamma A(p),\\
\infty,&\delta\leq\Gamma A(p).
\end{cases}
\]

Every optimal tenure has a strict housing distortion exactly for
\(w_r(p)<w<w_O(p)\). Both choices are distorted at a tenure tie inside
that interval. The interval is nonempty iff
\(\delta pr<A(p)[Y+\Gamma pr]\). Outside it, the hypothetical uncapped
rental optimum or the unrestricted lifetime optimum is implementable through
one of the actual tenure choices. The result compares all feasible fertility,
saving and housing deviations. It does not assume a tenure ranking.

For the spending identities let \(t_i=x_i/v_i\geq1\) and
\(\lambda=\chi/(\kappa p)\). Substituting the fertility condition gives

\[
\widehat x_i-v_i=\frac{v_i(t_i-1)}J
\left(1+\frac{\vartheta\lambda}{\lambda+t_i}\right),\qquad
x_i-\widehat x_i=\frac{v_i(t_i-1)}J
\left(\alpha+\frac{\vartheta t_i}{\lambda+t_i}\right),
\]
\[
\frac{\widehat x_i}{C}-g(x_i,v_i)=
\frac{v_i(t_i-1)(t_i-\alpha\lambda)}{JC(\lambda+t_i)},\qquad
Cg(x_i,v_i)-v_i=\frac{\lambda v_i(t_i-1)}{\lambda+t_i}.
\]

The child-cost condition is \(\alpha\lambda\leq1\); the required strictness
survives equality when \(t_i>1\). Concavity follows from

\[
d^2g=-\frac{2\chi\kappa p(v\,dx-x\,dv)^2}
{(\chi v+\kappa p x)^3}\leq0.
\]

For the joint planner, total goods and housing per initial young household are
\(\mathcal C=\bar x+B+\chi\bar n\) and
\(\mathcal H=\alpha(\bar v+B)/p+\kappa\bar n\). Optimizing consumption
and housing first gives

\[
X(n)=\frac{\mathcal C-\chi n}{2},\quad
S(n)=\frac{\mathcal H-\kappa n}{2},\quad
\max_n\{2\log X(n)+2\alpha\log S(n)+\vartheta\log n\}.
\]

Its derivative \(\vartheta/n-\chi/X(n)-\alpha\kappa/S(n)\) is strictly
decreasing from positive infinity to negative infinity on the feasible domain.
Identical parental preferences give common chosen fertility. Young housing
rises by an additional \(\kappa(n^F-\bar n)/2\) after this joint choice.

For a local transfer of \(\epsilon\) housing from old household \(j\), its
exact goods compensation is

\[
\Delta c_j^o=c_j^o\left[
\left(\frac{h_j^o}{h_j^o-\epsilon}\right)^\alpha-1\right]
=p\epsilon+O(\epsilon^2).
\]

Taking those goods from the recipient gives utility gain
\([\alpha/(h_i^y-\kappa n_i)-p/x_i]\epsilon+O(\epsilon^2)>0\).
The derivative of private fertility has the sign of
\(\alpha\kappa/(h_i^y-\kappa n_i)^2-p\chi/x_i^2\), yielding the
condition in the text. Financial settlement uses
\(T_i=\Delta c_i+p\Delta h_i\), \(\sum_iT_i=0\), and
\(\Delta a_i=-qP\Delta h_i^{\mathrm{owned}}\), including intermediaries.
Young continuation resources, old estates and total public spending are fixed.
Since \(z_i\geq y\) and estates are positive, a 100%-balance mortgage can
implement the required portfolio positions.

Analytical nonvacuity was checked against the earlier mixed-tenure stationary
construction, including old wealth generated by prior young decisions and the
replacement condition. Pro adds a separate open all-renter family satisfying
the strict primitive resource floor and a loss of young consumption under the
fixed-fertility full planner. Its global tenure and saving comparisons,
replacement, and housing clearing passed an independent review. The simple
construction chooses \(B_0=\bar w/J\), obtains \(\bar x>B_0\) from a capped
richer renter, and then sets \(KB_0<y<K\bar x\). Consequently
\(\bar w/J<B=y/K<\bar x\): the housing theorem applies and young consumption
falls at fixed fertility. Increasing \(y\) leaves young choices unchanged in
this regime because they save zero; endogenous cohort mass absorbs the change
in old housing demand. This is a comparison across stationary equilibria,
not a fixed-population price-invariance result.

The construction's extra conditions, including its \(\beta/q\) restriction,
are existence-witness restrictions and are not assumptions of the general
allocation theorem. The construction does not sign consumption at the joint
fertility optimum. No numerical reference point supports these existence
claims. Generic global uniqueness across tenure regimes is not established.

**Private-fertility counterexample.** The complete symbolic calculation is in
the all-tenure proof note, Section 9. It fixes \(\nu=1/2\) and competitive
mean fertility at two. For every \(0<\epsilon\leq1/10000\), all assumptions
of the broader exact-interval theorem hold, and mean private fertility after
assigned bundles is strictly below two, while joint-planner fertility exceeds
two. Global ownership deviations, saving, estate finance, housing clearing and
replacement were checked analytically. The simpler housing-access sufficient
test also covers the cap-distorted group in this example. The counterexample
disproves a universal private-fertility claim; it is not the basis for claiming
an empirically relevant equilibrium region.

**Work record.** Original bounded agents began about03:50UTC on September10.
The all-tenure interval and stationary witness were independently checked by
the fertility agent; the full allocation/joint-fertility proof was independently
checked by the necessity agent. Both passed. The subsequent private-fertility
counterexample and income cutoff were derived by the all-tenure agent; the
necessity agent independently passed the counterexample at about04:39UTC.
The lead checked the spending identities, Jensen argument, compensation and
limiting private-fertility roots. No model runs were used.

The three complete derivations are backed up in the [verification record](oracle_all_tenure_verification.md). Their original working files are
`tmp/theory_designs/overnight_all_tenures.md`,
`tmp/theory_designs/overnight_fertility_sharp.md`, and
`tmp/theory_designs/overnight_necessity.md`. Estate notation in those scratch
notes sometimes collides with a spending index; this reader uses
\(\widehat x\) for spending and \(e\) only for estates.

The author explicitly authorized Pro overnight. The
[focused prompt](../../../docs/prompts/oracle_all_tenure_simplification.md) and
[exact packet](oracle_all_tenure_simplification_packet.md) were submitted to
GPT-6 Pro in the existing conversation. Generation was verified at approximately
00:30EDT. Packet SHA-256:
`01237412c91a7db282d9a01eb6bf16407a8b9cc48966f9df0718de6f058a0ef7`.

The [completed Pro answer](oracle_all_tenure_simplification_response.md) was
captured with all 121 mathematical source objects. Its normalized browser and
local text agree: 15,058 characters, FNV-1a 1433091282; file SHA-256
`0a61e39cbf3bf1d367f0c689ab38f2c867444c8e27c379b546ecc1373358b3f1`.

The lead and an independent reviewer checked Pro's new strict-floor equilibrium
construction. The full joint allocation, its quadratic formula, and its implementation with
voluntary fertility at fixed final bundles were derived and checked separately.
All new conclusions are consolidated above; original
working derivations remain in the verification record. No main deck, protected
author draft, model code, transition proof or numerical run was changed. The
quantitative overnight work and sleep-prevention process remain separate.
A final independent reading pass found the recommendation honest and required
no mathematical correction. Its two clarifications were incorporated: fixing
future utility makes the dated objective equivalent to remaining lifetime
utility in this comparison, and the voluntary-fertility result uses newly
optimized fixed bundles. It does not simply release fertility at the earlier
fixed-fertility allocation. The morning decision is whether to adopt this
conditional static exercise; no further mathematical review is pending.

</details>

<details>
<summary>Earlier checked Pro result and discussion history</summary>

## Two periods with positive retirement income — current result

**The revised endowments work.** Both types can receive the same positive retirement income. Under explicit income and financing conditions, poorer young households rent at the size ceiling, richer young households own, and a planner choosing consumption and housing allocates more housing to the young. Every poorer young household receives more goods and housing, and chooses higher fertility at that bundle. A planner choosing fertility too raises mean fertility.

The essential distinction is between the housing inefficiency and the direction of the full allocation. The constraints generate a gap in households' willingness to give up consumption for housing. An additional restriction on the distribution of resources makes the equal-weight planner move both goods and space toward the young. That restriction remains economically substantive; it is not evidence that the result is universal or empirically mild.

The lead checked the planner, financial settlement, fertility and analytical parameter ranges. A separate Astra/max review passed the revised household choices, all tenure and fertility deviations, price cutoffs and existence proof. No consequential gap was found, so no further correction was sent to Pro. The proposal remains outside the main deck and manuscript.

[Exact new Pro response](oracle_two_period_positive_retirement_response.md) · [Question](../../../docs/prompts/oracle_two_period_positive_retirement.md) · [Conversation](https://chatgpt.com/c/6aa1e33d-d328-83ea-971e-0541c9e3e359)

### Model and statement

There are two periods and two types. A fraction \(f\) receives \(w_L\) when young and the remainder receives \(w_H>w_L>0\); both receive \(y>0\) when old. Preferences are

\[
u^y=\log(c^y-\chi n)+\alpha\log(h^y-\kappa n)+\vartheta\log n,
\qquad
u^o=\log c^o+\alpha\log h^o+\omega_B\log e,
\]

with lifetime utility \(u^y+\beta u^o\). Here \(n\) is children per parent household, \(\chi\) and \(\kappa\) are their goods and space needs, and \(e\) is the net estate. All coefficients are positive. Rentals cannot exceed \(r\); the common divisible housing stock can be owned in any quantity.

Let \(P\) be the house price, \(q\) the price of a unit of goods next period, \(\phi\) the financed share, and \(\tau>0\) the property-tax rate. Take \(0<q,\phi,f<1\) and a positive rental ceiling and housing stock. Write \(p=(1+\tau-q)P\) for rent. Let \(z\) be total resources upon retirement and \(a\) the present-value net financial position after interest service. A young owner chooses

\[
c+(1+\tau)Ph+a=w_i,\qquad
z=y+a/q+Ph,\qquad a\geq-q\phi Ph.
\]

Thus housing equity must satisfy \(z-y\geq(1-\phi)Ph\). Known retirement income is not available to finance that equity. Renters satisfy \(c+ph+a=w_i\), \(a\geq0\), and \(h\leq r\). Old households use the same financial rule, with their current wealth replacing \(w_i\), their estate replacing \(z\), and no subsequent income. They can sell and resize.

The following quantities are functions of primitives. Define

\[
K=1+\alpha+\omega_B,\quad A=\alpha+\vartheta,\quad J=1+A,\quad
\Gamma=\frac{q(1-\phi)}{1+\tau-q},
\]
\[
x=\frac{w_H+qy}{J+\beta K},\qquad
\ell=\frac{w_L}{J},\qquad m=\frac{qy}{\beta K},\qquad
b=\frac{(1+\Gamma)\ell}{1+\Gamma\ell/m},\qquad
B=f\frac yK+(1-f)\frac{\beta x}{q}.
\]

Here \(x\) is richer young consumption after children's goods needs, \(\ell\) is the corresponding unrestricted zero-saving rental demand of the poorer type, and \(m\) is its marginal-saving threshold. The quantity \(b\) bounds poorer young adult consumption in the regime below; \(B\) is average old consumption.

**Proposition.** Suppose

\[
\underbrace{\omega_B\geq\alpha\Gamma,\qquad
(\beta K-\Gamma A)x\geq qy}_{\text{old and richer young choices can be financed}},
\qquad
\underbrace{\ell<m}_{\text{poorer young want to bring retirement resources forward}}.
\]

Suppose also that replacement fertility occurs in the price interval derived below, where the rental ceiling binds and poorer households strictly prefer renting to every feasible purchase. Then the constructed competitive equilibrium is inefficient relative to an authority that can relax finance: a small housing transfer from an old household to a poorer young household can compensate the old in goods and leave the young better off.

If, in addition,

\[
\boxed{\quad B\geq fb+(1-f)x,\quad}
\]

the equally weighted planner choosing all current consumption and housing gives the young more housing in aggregate and every poorer young household more consumption and housing. At that bundle each poorer parent chooses higher fertility. The joint planner, valuing existing parents and choosing fertility too, chooses mean fertility above the competitive mean.

This theorem constructs and characterizes one equilibrium regime. The price and real allocations are unique in that regime; it does not exclude equilibria in other regimes. Old tenure and gross mortgage/saving positions can tie.

### What the conditions mean

The two young types lie on opposite sides of the saving decision. The poorer type would like to use future income now. The richer type wants to carry additional wealth into retirement and can finance a larger owned home. The income restrictions can be written directly as

\[
\frac{w_L}{y}<\frac{qJ}{\beta K},
\qquad
\frac{w_H}{y}\geq
\frac{q(J+\Gamma A)}{\beta K-\Gamma A},\qquad
\beta K>\Gamma A.
\]

The old estate condition ensures their desired housing–estate combination is feasible. These are sufficient restrictions supporting the stated tenure pattern, not necessary restrictions for every housing inefficiency.

The additional boxed condition ensures that the old have enough consumption resources for the full planner to increase young adult consumption as well as space. It is automatically satisfied under the preceding bounds when \(\beta\geq q\). When \(\beta<q\), it is

\[
f\left(\frac yK-b\right)
\geq(1-f)\left(1-\frac\beta q\right)x.
\]

The poorer type's rise in consumption resources at retirement must offset the richer type's consumption decline. A very large rich-young endowment or a small poorer share can make this fail. There is no imposed ordering between \(\beta\) and \(q\); the financing restrictions still depend on patience.

The price restriction is also a primitive test, not an assumed equilibrium housing ordering. At a candidate rent \(p\), solve the poorer renter's strictly concave fertility choice

\[
\frac{\vartheta}{n_L(p)}
=\frac{\chi}{w_L-pr-\chi n_L(p)}
+\frac{\alpha\kappa}{r-\kappa n_L(p)}.
\]

Set \(x_L=w_L-pr-\chi n_L\), \(s_L=r-\kappa n_L\), and \(M(p)=\alpha x_L/s_L\). Let \(p_-\) and \(p_+\) solve \(M(p_-)=kp_-\) and \(M(p_+)=p_+\), where \(k=b/\ell>1\). Define

\[
F(p)=fn_L(p)+(1-f)\frac{\vartheta x}{\chi+\kappa p}.
\]

If one parent generates \(\nu n\) entering households, the required test is
\(F(p_+)<1/\nu<F(p_-)\). Both cutoffs are quadratic roots, and \(F\) is strictly decreasing. Hence exactly one rent in this interval satisfies replacement; housing clearing determines the common cohort mass. Household choices and the cutoffs are elementary functions; the equilibrium rent is a unique scalar root.

### The planner and the fertility argument

The planner gives equal weight to each currently living household. It chooses current consumption and housing, preserves each young household's continuation wealth and each old household's estate, and holds continuation prices fixed. Initially it also fixes individual fertility. Let \(\bar x\) and \(\bar s\) be young average consumption and housing after subtracting children's needs, and let \(\bar h_o\) be average old housing. At the competitive allocation,

\[
\bar x=fx_L+(1-f)x,\qquad
\bar s=fs_L+(1-f)\alpha x/p,\qquad
\bar h_o=\alpha B/p.
\]

The resource bound and the binding rental ceiling imply
\(B>\bar x\) and \(\bar h_o>\bar s\). The planner equalizes adult goods and space across all households:

\[
X=\frac{\bar x+B}{2},\qquad S=\frac{\bar s+\bar h_o}{2},
\qquad
(c_i^{y,F},h_i^{y,F})=(X+\chi n_i,S+\kappa n_i),\qquad
(c_i^{o,F},h_i^{o,F})=(X,S).
\]

Therefore

\[
\boxed{\quad \bar h_y^F-\bar h_y^E
=\frac{\bar h_o-\bar s}{2}>0.\quad}
\]

Also \(X>x_L\) and \(S>s_L\), so each poorer young household receives more total goods and housing. Its fertility derivative at the old choice becomes

\[
\frac{\vartheta}{n_L}-\frac{\chi}{X}-\frac{\alpha\kappa}{S}>0,
\]

which implies higher chosen fertility at the assigned bundle. If the planner chooses fertility too, it chooses a common \(n^F\). Concavity of
\((\chi/x+\alpha\kappa/s)^{-1}\), together with \(X>\bar x\) and \(S>\bar s\), makes its fertility derivative positive at competitive mean fertility. Thus \(n^F>\bar n^E\). No independent welfare weight is assigned to unborn people, and no small goods-cost assumption is needed for these two conclusions.

The compensated efficiency improvement is a separate comparison: there the young give up goods to obtain housing, so fertility need not increase. Outside the stronger resource condition, the exact fertility test for the full assigned bundle is

\[
\chi\left(\frac1{x_L}-\frac1X\right)
+\alpha\kappa\left(\frac1{s_L}-\frac1S\right)>0.
\]

### Scope and decisions still belonging to the author

This is an analytical conditional result for two periods with unequal young endowments and common positive retirement income. Its parameter region is established through explicit inequalities, without a numerical reference point. It does not require every old household to occupy a larger home than every young household. The individual housing gain concerns poorer young households; the age-wide gain is aggregate.

The calculation retains equal housing tastes across ages, old secured borrowing, a common stock that can change tenure, no owner size bound, external finance, and estates that do not determine entrant endowments. Positive property taxes finance an additively valued public service, held fixed in the comparison. That fiscal closure and these other proposed simplifications have not been adopted into the paper. The dated fertility results do not establish a policy implementation, transfers-only constrained inefficiency, or a transition to a new population level.

<details>
<summary>Lead verification: planner resources, fertility, and nonempty primitive ranges</summary>

Normalize by the common cohort mass. Available current goods and housing are
\(\mathcal C=\bar x+B+\chi\bar n\) and
\(\mathcal H=\bar s+\bar h_o+\kappa\bar n\). With fixed fertility and continuation values, strict concavity and the two resource constraints give the displayed unique planner allocation. Young continuation utility and old estate utility are constant in this dated problem.

For financial settlement use \(T_i=\Delta c_i+p\Delta h_i\) and
\(\Delta a_i=-qP\Delta h_i^{\mathrm{owned}}\), with intermediary counterpart positions. Transfers sum to zero because aggregate consumption and housing are unchanged. Every young continuation target satisfies \(z_i\geq y\), and every old estate is positive. At a new owned home of value \(Ph\), set outstanding principal to \(\max\{Ph-(z-y),0\}\) and terminal liquid assets to \(\max\{z-y-Ph,0\}\); for the old replace \(z-y\) with \(e\). Principal is at most \(Ph\), so allowing a fully financed balance suffices. The consolidated and dated budgets both hold. Public expenditure is fixed because \(P\) and the physical stock are fixed.

For the separate compensated transfer, preserve the old household's utility and estate by giving it
\(\Delta c_o=c_o[(h_o/(h_o-\epsilon))^\alpha-1]\).
Its leading term is \(p\epsilon\), since the unrestricted old MRS is \(p\). Taking these goods from a poorer young household and giving it \(\epsilon\) housing yields a first-order gain
\((\alpha/s_L-p/x_L)\epsilon>0\). Positivity and strict preference hold for sufficiently small transfers. It is specifically the additional bound on \(B\), not the condition \(\ell<m\), that is unnecessary for this efficiency comparison.

For joint fertility, symmetry and strict concavity imply common fertility. The reduced objective, apart from constants, is
\(2\log X(n)+2\alpha\log S(n)+\vartheta\log n\), where
\(X(n)=(\mathcal C-\chi n)/2\) and
\(S(n)=(\mathcal H-\kappa n)/2\). Its derivative is
\(\vartheta/n-\chi/X(n)-\alpha\kappa/S(n)\), strictly decreasing from positive infinity to negative infinity over the feasible interval. The second differential of
\(g(x,s)=xs/(\chi s+\alpha\kappa x)\) is
\(-2\chi\alpha\kappa(s\,dx-x\,ds)^2/(\chi s+\alpha\kappa x)^3\leq0\).
Since private fertility equals \(\vartheta g(x_i,s_i)\), Jensen gives
\(\vartheta/\bar n\geq\chi/\bar x+\alpha\kappa/\bar s\).
The strict resource inequalities make this exceed the planner's marginal child cost at \(\bar n\), proving the claimed sign.

For analytical nonvacuity, fix \(\omega_B>\alpha\Gamma\) and
\(\beta K>\Gamma A\). If \(\beta\geq q\), any
\(0<w_L<Jm\), any
\(w_H>qy(J+\Gamma A)/(\beta K-\Gamma A)\), and any \(0<f<1\)
satisfy the strict financial and resource bounds. If \(\beta<q\), write
\(\rho=\beta/q\) and instead choose

\[
0<w_L<\frac{J\rho m}{1+\Gamma(1-\rho)},\qquad
w_H>\frac{qy(J+\Gamma A)}{\beta K-\Gamma A},\qquad
\frac{(1-\rho)x}{y/K-b+(1-\rho)x}<f<1.
\]

The first bound implies \(b<\rho m=y/K\), so the interval for \(f\) is nonempty. Both cases imply \(w_H>w_L\). Strictly decreasing \(F\) and distinct price cutoffs provide a nonempty open interval for \(1/\nu\). This proves an open primitive region, not that an externally fixed replacement conversion necessarily passes the test. For a specified \(\nu\), the displayed market inequalities must be checked.

Capture verification: all 139 math objects retained. Browser and local whitespace-normalized text match exactly: 18,582 characters, FNV-1a 1053618808. No numerical model run or PDF build was used.

</details>

<details>
<summary>Independent verification of the revised household and equilibrium argument</summary>

# Independent audit: common positive retirement income

**Verdict: PASS within the stated regime. No consequential algebraic error found.**

Reviewed the new household, tenure, cutoff, and analytical nonvacuity arguments in `output/model/simplified_olg_amendments/oracle_two_period_positive_retirement_response.md`, especially equations (3)–(16), (24)–(27), and the two concluding quadratics. This review takes the stipulated two-period serviced-interest contract as given. The lead owns the dated planner and parental-fertility review; those claims are not independently certified here. No simulations, parameter search, source edits, or lifecycle redesign were used.

## 1. The finance conditions and consumption ordering are correct

Write \(d=1+\tau-q>0\), \(p=dP\), \(\eta=q(1-\phi)>0\), and \(\Gamma=\eta/d\). An owner satisfies

\[
c+ph+qz=w_i+qy,
\qquad q(z-y)\geq\eta Ph=\Gamma ph.
\]

For an old household, the unrestricted logarithmic optimum is

\[
c^o=z/K,\qquad ph^o=\alpha z/K,\qquad qe=\omega_Bz/K.
\]

Its finance condition is exactly \(\omega_B\geq\alpha\Gamma\), independently of \(z>0\). This makes \(V^o(z)=K\log z+C(p)\) an attainable upper bound against every old-tenure choice.

For the richer young, the unrestricted allocation has adult consumption \(x=(w_H+qy)/(J+\beta K)\), adult housing \(\alpha x/p\), fertility \(\vartheta x/(\chi+\kappa p)\), and \(qz_H=\beta Kx\). Its **exact** finance condition is

\[
\beta Kx-qy\geq
\Gamma x\left(\alpha+
\vartheta\frac{\kappa p}{\chi+\kappa p}\right).
\]

Thus equation (11) correctly subtracts \(qy\). The proposed sufficient condition

\[
(\beta K-\Gamma A)x\geq qy
\]

implies \(\beta K>\Gamma A\), \(x>m=qy/(\beta K)\), and feasibility at every relevant price. Indeed, because \(\chi>0\), the exact housing-finance restriction is strictly slack under that sufficient condition even if its displayed weak bound holds at equality. No ordering between \(\beta\) and \(q\) is used.

## 2. The sharper \(b\) bound and the global tenure argument hold

Let \(u=\ell/m\in(0,1)\). The construction gives

\[
k=\frac{1+\Gamma}{1+\Gamma u}>1,
\qquad
\frac bm=ku=\frac{(1+\Gamma)u}{1+\Gamma u}<1,
\]

and, for \(\sigma=1-b/m>0\), exactly \(k-1=\Gamma\sigma\).

For a candidate capped renter at \(p\in(p_-,p_+)\), put \(t=M(p)/p\in(1,k)\). The fertility first-order condition and \(r=s_L+\kappa n_L\) imply

\[
\vartheta x_L=(\chi+t\kappa p)n_L,
\qquad
tpr=Ax_L-\chi n_L.
\]

Combining these with \(w_L=x_L+pr+\chi n_L\) yields the response's identity

\[
tw_L-Jx_L=(t-1)(x_L+\chi n_L)>0,
\]

so \(x_L<t\ell<b<m\). In particular, the no-saving conclusion is derived, not imposed.

Define \(\xi=\alpha/s_L-p/x_L=(t-1)p/x_L>0\) and \(\zeta=1/x_L-1/m>0\). Then

\[
\zeta>\frac{\sigma}{x_L},
\qquad
\xi<\frac{\Gamma\sigma p}{x_L}<\eta P\zeta.
\]

The response uses a weak inequality for the first relation; that weaker statement is valid. At the candidate, the gradient of optimized lifetime utility in \((c,h,n,z)\) is

\[
\left(\frac1{x_L},\frac p{x_L}+\xi,0,\frac qm\right).
\]

Concavity and the common consolidated budget therefore imply, for every feasible alternative including different fertility,

\[
U-U_L\leq\xi(h-r)-q\zeta(z-y).
\]

Every rental alternative has \(h\leq r\), \(z\geq y\), so it cannot improve. Every ownership alternative has \(q(z-y)\geq\eta Ph\), giving

\[
U_O-U_L\leq(\xi-\eta P\zeta)h-\xi r<0.
\]

All feasible plans have \(h>\kappa n>0\). This proves strict preference against **all** owned-home sizes and fertility choices, not merely a comparison at the original \(n_L\). Strict concavity also gives uniqueness of the real renter allocation.

Finally,

\[
n_L=\frac{\vartheta x_L}{\chi+t\kappa p}
<\frac{\vartheta x}{\chi+\kappa p}=n_H,
\qquad
s_L=\frac{\alpha x_L}{tp}<\frac{\alpha x}{p}.
\]

Consequently \(h_H^y>r\). The richer type's feasible unrestricted optimum is unattainable by renting, so ownership is strict for this young type as well.

## 3. The monotonicity and both explicit quadratics are correct

For each \(p\in(0,w_L/r)\), the poorer fertility first-order condition is strictly decreasing in \(n\), tending to opposite infinities at the endpoints of

\[
0<n<\min\{(w_L-pr)/\chi,r/\kappa\}.
\]

There is exactly one feasible solution. Set

\[
Q=\frac{\vartheta}{n_L^2}+
\frac{\chi^2}{x_L^2}+
\frac{\alpha\kappa^2}{s_L^2}>0.
\]

Implicit differentiation gives

\[
n_L'=-\frac{\chi r}{x_L^2Q}<0,
\quad
x_L'=-r\frac{\vartheta/n_L^2+\alpha\kappa^2/s_L^2}{Q}<0,
\quad
s_L'=-\kappa n_L'>0.
\]

Thus \(M=\alpha x_L/s_L\) strictly decreases, has a positive finite limit at \(p=0\), and tends to zero at \(p=w_L/r\). This proves the two cutoff roots and their strict ordering. Both components of \(F\) strictly decrease, which proves the unique interior replacement root under Market.

At \(M=tp\), eliminating \(x_L,n_L\) gives exactly

\[
\kappa tr(t+A)p^2+
\left[\chi r\{\alpha+t(1+\vartheta)\}
-\kappa tAw_L\right]p-\alpha\chi w_L=0.
\]

The leading coefficient is positive and the constant negative, so precisely one root is positive. The polynomial at \(p=w_L/r\) equals

\[
\frac{\kappa t^2w_L^2}{r}+\chi tw_L(1+\vartheta)>0,
\]

which places that root inside the valid price domain. It corresponds to the unique feasible candidate and introduces no extraneous positive cutoff.

Multiplying the fertility first-order condition by its positive denominators gives exactly

\[
\chi\kappa Jn^2-
\big[(\vartheta+\alpha)\kappa v+(\vartheta+1)\chi r\big]n
+\vartheta rv=0,\qquad v=w_L-pr.
\]

The smaller root is the unique feasible one. If the two subsistence boundaries coincide, the larger root lies exactly at the excluded boundary; otherwise it lies beyond the smaller boundary.

## 4. Equations (24)–(26) genuinely prove analytical nonvacuity

Put \(\rho=\beta/q\), so \(y/K=\rho m\).

- If \(\rho\geq1\), (24) is \(0<\ell/m<1\), which gives \(b<m\leq y/K\).
- If \(0<\rho<1\), solving \(b/m<\rho\) gives exactly
  \[
  \frac{\ell}{m}<\frac{\rho}{1+\Gamma(1-\rho)}.
  \]
  This is (24), and also implies \(\ell<m\).

Hence \(\delta=y/K-b>0\) in both cases. Equation (25) is algebraically equivalent to the strict richer-finance condition:

\[
(\beta K-\Gamma A)\frac{w_H+qy}{J+\beta K}>qy.
\]

Its lower threshold exceeds \(Jm\), since

\[
\frac{qy(J+\Gamma A)}{\beta K-\Gamma A}-Jm
=\frac{qy\Gamma A(\beta K+J)}{\beta K(\beta K-\Gamma A)}>0.
\]

Meanwhile (24) puts \(w_L<Jm\), so the required \(w_H>w_L\) follows.

The resource gap is exactly

\[
B-[fb+(1-f)x]
=f\delta+(1-f)(\rho-1)x.
\]

For \(\rho\geq1\), any \(0<f<1\) makes it positive. For \(\rho<1\), positivity is equivalent to

\[
f>\frac{(1-\rho)x}{\delta+(1-\rho)x},
\]

which is (26). The lower bound is strictly below one for every finite admissible \(w_H\). This verifies the economic interpretation in (27): at \(\beta<q\), enough poorer households with rising old-age resources are needed to offset the richer group's declining adult consumption profile.

The preliminary strict preference/finance restrictions themselves are feasible analytically: for any fixed positive \(\alpha,\vartheta,\beta\) and admissible \(q,\phi,\tau\), one may choose

\[
\omega_B>
\max\left\{0,\alpha\Gamma,
\frac{\Gamma A}{\beta}-1-\alpha\right\}.
\]

Then choose \(y>0\), the open endowment intervals, \(f\), arbitrary positive \(\chi,\kappa,r\), and any replacement level in the nonempty interval \(F(p_+)<1/\nu<F(p_-)\). Strict inequalities and continuously varying simple roots give an open set of primitives; this is not an equilibrium obtained only at an equality-tuned parameter point. Any \(\bar H>0\) yields the positive cohort mass \(N=\bar H/D(p^*)\).

## 5. Scope and qualifications to preserve

The proved uniqueness is **price, real household allocations, and cohort mass within the verified interval/regime**. Other regimes are not excluded. Old tenure can tie whenever the optimal home fits the rental cap, and gross mortgage/bond decompositions can tie. The response already appropriately limits its headline uniqueness claim.

The substantive assumptions are the stipulated positive retained-equity requirement with nonpledgeable future \(y\), constrained rentals, the common stock, financeability of the old unrestricted allocation, sufficient rich accumulated wealth **above** \(y\), externally fixed entrant type assignment/endowments, and the public-service tax closure. The additional distributional inequality \(B\geq fb+(1-f)x\) is separately sufficient for the directional equal-weight redistribution claim; that inequality is not needed for the household/tenure construction. The other inequality grouped under Resources, \(\ell<m\), is still used to construct the verified regime and cannot be dropped by invoking the efficiency-gap statement. This audit neither adds assumptions about intermediate paychecks nor questions the approved two-period simplification.

**Required corrections: none in the assigned scope. Unresolved within scope: none.** The planner, its settlement, parental-fertility welfare claims, and the literature comparison remain outside this review's certification.

</details>

<details>
<summary>Earlier rounds and their verification — historical; superseded where the current result differs</summary>

## Interest-serviced mortgage — current assessment

The follow-up repairs the financing problem and supplies an analytical tenure comparison. Under its stated conditions, the dated planner gives young households more housing, and the joint parental planner chooses higher fertility. This is a useful candidate for the illustration. It remains conditional on an income profile with substantial resources arriving later, a particular mortgage payment schedule, and the housing-market and welfare choices below. It is not a completed theory of the policy transition. No specification has been adopted.

The lead checked the planner, fertility and resource restrictions. A separate Astra/max review verified the mortgage cash flows at every payment date, the primitive cutoffs and all tenure deviations. A further narrow check permits the last coupon to be paid from sale proceeds, under the short-interval restriction below.

[Exact follow-up response](oracle_essential_theory_finance_followup_response.md) · [Pro conversation](https://chatgpt.com/c/6aa1e33d-d328-83ea-971e-0541c9e3e359) · [Exact question](../../../docs/prompts/oracle_essential_theory_finance_followup.md)

### The economic change

Mortgage interest is paid during the age; the loan principal is repaid when the house is sold. The same rule applies to young and old. Let \(P\) be the house price, \(q\) the full-age discount factor, \(\phi\) the fraction of the purchase price that can be borrowed, and \(\tau\) the property-tax rate. Rent per unit is \(p=(1+\tau-q)P\).

There is no income between age boundaries. The household must therefore provide for interest payments using current resources. Buying \(h\) at maximum leverage requires, before consumption,

\[
\underbrace{(1-\phi+\tau)Ph}_{\text{equity and property tax}}
+\underbrace{(1-q)\phi Ph}_{\text{present value of interest}}
=\underbrace{ph}_{\text{rental user cost}}
+\underbrace{q(1-\phi)Ph}_{\text{equity retained until sale}}.
\]

The extra liquidity requirement is positive for every \(0<q,\phi<1\). The old condition \(\phi<q\) disappears. This does not make an origination loan limit by itself sufficient: it is the combination of the deposit, interest service and delayed access to later income.

Interest-only mortgages with scheduled interest payments and later principal repayment are an established contract form. That supports the payment schedule, not the model's boundary-only income assumption or its empirical fit. [CFPB explanation](https://www.consumerfinance.gov/ask-cfpb/what-is-an-interest-only-loan-en-101/).

The exact displayed formula assumes every coupon is funded before liquidation. **The result also survives paying the final coupon from the sale proceeds**, provided \(\phi R_{\rm last}<1\), where \(R_{\rm last}\) is the gross interest factor for the final payment interval. In that variant, replace \(\eta=q(1-\phi)\) throughout by \(\eta=q(1-\phi R_{\rm last})\). This is a condition on the last mortgage-payment interval, not on the entire utility age. For equal payment intervals and \(m\) coupons, it is \(\phi q^{-1/m}<1\). A single coupon at the very end reproduces the original problem; regular intervening payments do the work. The supporting review proves feasibility at every earlier payment date as well.

### One complete conditional statement

Households live for two ages. Their lifetime utility is \(u^y+\beta u^o\), with

\[
u^y=\log(c^y-\chi n)+\alpha\log(h^y-\kappa n)+\vartheta\log n,
\qquad
u^o=\log c^o+\alpha\log h^o+\omega_B\log e.
\]

Here \(e\) is the net estate at death; \(\chi\) and \(\kappa\) are goods and space needs per child. All coefficients are positive. A fraction \(f\in(0,1)\) receives \(W\) when young and known, nonpledgeable \(Y\) when old. Call these the deferred-income families. Everyone else receives \(W+qY\) when young. The types have the same present-value resources. They choose tenure, consumption, housing, saving and fertility. Rentals are limited to \(r\) units; any size can be owned from the common divisible stock.

Define

\[
d=1+\tau-q,\quad \eta=q(1-\phi),\quad
K=1+\alpha+\omega_B,\quad J=1+\alpha+\vartheta,
\]
\[
L=\frac{W+qY}{J+\beta K},\qquad
B=f\frac{Y}{K}+(1-f)\frac{\beta L}{q}.
\]

The quantity \(L\) is the liquid young type's unrestricted adult consumption; \(B\) is average old consumption at unrestricted choices. Both are explicit functions of primitives.

**Proposition.** Suppose:

\[
\begin{array}{ll}
\text{Unrestricted choices are financeable:}&
\omega_B d\geq\alpha\eta,\quad
\beta Kd\geq\eta(\alpha+\vartheta);\\[3pt]
\text{Sufficient later resources:}&
L\geq W,\quad B\geq fW+(1-f)L;\\[3pt]
\text{Replacement lies in the rental regime:}&
F(p_+)<1/\nu<F(p_-).
\end{array}
\]

The price bounds and fertility function in the last line are computed from primitives as follows. For any candidate rent \(0<p<W/r\), solve the unique interior quadratic fertility choice

\[
\frac{\vartheta}{n_C(p)}
=\frac{\chi}{W-pr-\chi n_C(p)}
+\frac{\alpha\kappa}{r-\kappa n_C(p)}.
\]

Write \(x_C=W-pr-\chi n_C\), \(s_C=r-\kappa n_C\), and
\(\mathcal M(p)=\alpha x_C/s_C\). Set

\[
\sigma=1-\frac{\beta KW}{qY}>0,\qquad
\mathcal M(p_-)=\left(1+\frac{\sigma\eta}{d}\right)p_-,
\qquad \mathcal M(p_+)=p_+,
\]
\[
F(p)=fn_C(p)+(1-f)\frac{\vartheta L}{\chi+\kappa p}.
\]

The two thresholds are unique positive roots of the explicit quadratic in the response. Thus the test contains neither an assumed equilibrium allocation nor an unobserved multiplier.

Then a positive stationary equilibrium exists with rent \(p^*\in(p_-,p_+)\) solving \(F(p^*)=1/\nu\). Deferred-income young families strictly prefer renting \(r\) to every feasible purchase, including purchases with a different fertility choice. Liquid young families own larger homes. Old households choose unrestricted consumption and housing. Housing clearing determines the cohort size. Prices, real quantities and population are unique within this regime; other regimes are not ruled out. Old tenure can tie, and gross mortgage/saving portfolios need not be unique.

For the equally weighted dated planner, fix individual fertility, young continuation wealth, old net estates and future prices, and choose all current consumption and housing. Define young average adult goods and adult space by

\[
\bar x=fx_C+(1-f)L,\qquad
\bar s=fs_C+(1-f)\frac{\alpha L}{p^*},
\qquad \bar h_o=\frac{\alpha B}{p^*}.
\]

The assumptions imply \(B>\bar x\) and \(\bar h_o>\bar s\). The planner assigns

\[
X=\frac{\bar x+B}{2},\quad S=\frac{\bar s+\bar h_o}{2},
\qquad
c_i^{y,F}=X+\chi n_i,\quad h_i^{y,F}=S+\kappa n_i,
\quad c_j^{o,F}=X,\quad h_j^{o,F}=S.
\]

Consequently,

\[
\boxed{\bar h_y^F-\bar h_y^E
=\frac{\bar h_o-\bar s}{2}>0.}
\]

Every deferred-income young family receives more consumption and housing. If it then chooses fertility at that bundle, it chooses more children. If the dated planner also chooses fertility, valuing only the current parents, its mean fertility satisfies
\(\boxed{n^F>\bar n^E}\).

### Why this is more than redistribution under equal weights

The separate, friction-specific statement is

\[
\frac{u_h^y}{u_c^y}=\frac{\alpha x_C}{s_C}>p^*
=\frac{\alpha c^o}{h^o}=\frac{u_h^o}{u_c^o}.
\]

A small amount of housing transferred from an old household to a deferred-income young family can compensate the old in goods and leave the family better off, preserving estates and continuation wealth. Finance must be relaxed for this allocation. This is a dated efficiency result relative to that authority, not constrained inefficiency for an authority retaining the mortgage limit.

The stronger aggregate direction under equal weights also uses the later-resource restriction. It should not all be attributed to the mortgage. Nor does the compensated transfer by itself sign fertility: that family gives up consumption. Higher fertility follows for the full allocation above, where it receives both more goods and more housing.

### How demanding are the conditions?

There is no restriction on the ordering of \(\beta\) and \(q\), or of \(\phi\) and \(q\). The financing restrictions still require enough desired saving or estate provision to keep liquid young and old choices unconstrained. The old-estate condition is sufficient, not a claim that bequests necessarily cause the misallocation.

The first resource condition is exactly

\[
\frac{Y}{W}\geq\frac{\alpha+\vartheta+\beta K}{q}.
\]

It can require substantial deferred income. The second is

\[
f\left(\frac{Y}{K}-W\right)
\geq(1-f)\left(1-\frac{\beta}{q}\right)L.
\]

When \(\beta\geq q\), the first resource condition already implies the second. When \(\beta<q\), the deferred-income households must offset the other type's declining consumption. Writing \(M=J+\beta K\), the exact additional requirement is

\[
Y\left[\frac{fM}{K}+(1-f)(\beta-q)\right]
\geq W\left[fM+(1-f)(1-\beta/q)\right].
\]

For \(\beta<q\), a finite sufficiently large \(Y\) can satisfy it precisely when
\(f>(q-\beta)K/(J+qK)\). The replacement-price test must also hold; simply increasing \(Y\) at fixed remaining primitives does not guarantee the theorem.

The price condition identifies a region where desired rental housing exceeds the cap, while ownership's equity requirement makes renting preferable. The new supporting-hyperplane proof checks all ownership choices. Unlike the first answer, it does not require every cap-sized home to be literally unaffordable. As \(\phi\) approaches one, the sufficient interval narrows. This is a conditional mechanism, not a claim of large misallocation at every 80% LTV.

These inequalities are jointly possible analytically, without a numerical reference point: the financing and resource restrictions can be met with finite estate and deferred-income parameters; \(p_-<p_+\) then gives a nonempty replacement interval. At fixed demographic conversion \(\nu\), scaling both child costs \(\chi,\kappa\) by the same positive factor leaves the price bounds unchanged and scales both fertility bounds inversely. An open interval of such scales satisfies the replacement test. This establishes nonvacuity, not empirical plausibility.

### What remains a decision, and what remains unproved

The essential choice is whether an illustration of illiquid later wealth is the right first model for the paper. All income already received is usable, but no working-age earnings arrive between housing choices. Mortgage interest must therefore be funded out of the initial resources. A wage-accumulation story needs additional timing and choices.

The analytical result covers two liquidity types with equal lifetime resources, not an arbitrary income–wealth distribution. Household choices and the tenure cutoffs are explicit; the equilibrium rent is characterized by one unique scalar root. The response does not supply a short closed-form expression for that equilibrium rent.

The result also uses equal housing tastes across ages, old secured borrowing, no owner size bound, a common stock able to change tenure, exogenous entrant liquidity types, and estates that do not fund entrants. Property taxes are positive but finance public services rather than household rebates. These remain unadopted changes to the earlier note. Separate physical tenure stocks or an owner cap would need their own feasibility check.

This theorem does not establish that every old household occupies a larger home than every young household. It establishes that the planner gives the young more housing, and each deferred-income family more housing individually. It also does not establish transfers-only constrained inefficiency, a property-tax implementation, or the policy transition for this repaired mortgage. The earlier transition proof belongs to the first mortgage specification; it must not be carried over without checking the dated payment and wealth equations.

### Lead verification of welfare and fertility

For fixed individual continuations, financial settlement is possible with
\(T_i=\Delta c_i+p\Delta h_i\) and
\(\Delta a_i=-qP\Delta h_i^{\rm owned}\), with the opposite changes to intermediary positions. Aggregate current consumption and housing are unchanged, so transfers sum to zero. At the preserved targets \(z_C=Y\), \(z_L>0\), and \(e>0\), allowing mortgage balances up to the full house value is sufficient. For example, choose \(\ell=\max\{Ph-(z-Y),0\}\) and \(b=\max\{z-Y-Ph,0\}\). Both ages obey the same interest and sale accounting.

For joint fertility, let \(\mathcal C=\bar x+B+\chi\bar n\) and
\(\mathcal H=\bar s+\bar h_o+\kappa\bar n\). At a candidate common fertility \(n\),
\[
X(n)=\frac{\mathcal C-\chi n}{2},\qquad
S(n)=\frac{\mathcal H-\kappa n}{2}.
\]
Strict concavity makes equal fertility optimal across otherwise identical parents. The derivative of the reduced planner objective is
\[
\frac{\vartheta}{n}-\frac{\chi}{X(n)}
-\frac{\alpha\kappa}{S(n)},
\]
which is strictly decreasing. The function
\((\chi/x+\alpha\kappa/s)^{-1}\) is concave on positive adult goods and space. Applying Jensen to private fertility choices, and then using
\(X(\bar n)>\bar x\), \(S(\bar n)>\bar s\), makes this derivative positive at competitive mean fertility. Its unique zero is therefore larger. This proves the joint parental result without a negligible child-goods cost assumption.

Capture verification: all 120 displayed/inline math objects retained. Browser and local whitespace-normalized text match exactly: 15,826 characters, FNV-1a 534341554.

<details>
<summary>Independent verification of the repaired mortgage, all tenure choices, and the final coupon at sale</summary>

Lead clarification: the constructive payment schedules below use nonnegative deterministic interest over each interval, so \(D_j\geq q\). Constant interest over the long age is a sufficient schedule. For the closing-coupon variant, “100% balance” includes both principal and the final accrued interest; allowing that total balance is sufficient for the dated planner's financial settlement.

# Independent review: interest-serviced mortgage follow-up

Reviewed the new `oracle_essential_theory_finance_followup_response.md` as an **unadopted stationary proposal**. Scope: mortgage cash flows at both ages, exact feasibility, the primitive price cutoffs, and every tenure deviation. No model run, transition proof, literature review, planner/fertility adjudication, or maintained-source edit.

**Verdict: the stationary finance and tenure proposition passes in the stated model.** I find no algebraic failure in (1)–(14), (F), or the cutoff polynomial. The interpretation must remain an interest-serviced, boundary-income model: its cash requirement covers the whole age's debt service, not only the origination deposit. Equilibrium uniqueness needs the usual real-allocation qualification, and additionally does not identify gross loan/saving positions.

| Object | Verdict | Necessary qualification |
|---|---|---|
| PV interest and within-age solvency | Pass, constructively | Coupons and liquid saving must use the same deterministic discount schedule. All specified coupons are serviced before liquidation; income at the next boundary is unavailable for those coupons. |
| Projection \(a\ge-q\phi Ph\) | Exact | This is a projection of origination credit **plus prefunded future interest**, not an origination LTV redefinition. |
| Both ages and owner FOC | Pass | The old also may originate secured credit and must service interest. Purchase and resale prices here are the same stationary \(P\). |
| (F): old and liquid-young feasibility | Pass; sufficient | The second inequality uses a strict upper bound on the young owner's housing expenditure and need not be necessary. |
| \(n_C(p)\), \(\mathcal M(p)\), cutoff quadratic | Pass | Positive child-goods and child-space costs and all log-domain restrictions are maintained. |
| (R), \(\sigma>0\), (P) | Pass for this scope | (R) is substantive endowment heterogeneity; (P) explicitly selects replacement fertility inside the verified regime. |
| All type-C ownership deviations | Pass, globally | The concavity bound permits different consumption, housing, saving, and fertility. No owner-size minimum is used. |
| Every type L owns \(h_L>r\) | Pass, strictly | An unrestricted optimum is financeable, and its housing exceeds the rental ceiling. |
| Old tenure / uniqueness | Qualification required | If unrestricted old housing is at most \(r\), renting and owning can tie. Gross borrowing and saving generally also tie. The scalar price root and real allocations are unique **within this regime**, not all equilibrium objects or all regimes. |

## Exact cash-flow construction

Let coupon dates be \(j=1,\dots,m\), with discount factors \(D_0=1\), \(D_m=q\), and gross one-step returns \(R_j=D_{j-1}/D_j\). A fixed principal \(\ell\) pays coupon \(I_j=(R_j-1)\ell\), with the final specified coupon paid before the house's sale and principal settlement. Therefore

\[
\sum_{j=1}^m D_j I_j=(1-q)\ell.
\]

To leave liquid wealth \(b\ge0\) immediately before sale, reserve initially

\[
S_0=(1-q)\ell+qb.
\]

After each coupon the saving balance is exactly

\[
S_j=R_jS_{j-1}-I_j
 =\left(1-\frac q{D_j}\right)\ell+\frac q{D_j}b\ge0.
\]

Thus the terminal condition really does imply an implementable, nonnegative saving path at every payment date; there is no unfinanced interest hidden between the two utility ages. At origination,

\[
w=c+\tau Ph+(Ph-\ell)+S_0.
\]

For any proposed \(a\ge-q\phi Ph\), the constructive choice

\[
\ell=\max\{0,-a/q\},\qquad b=\max\{0,a/q\}
\]

satisfies \(0\le\ell\le\phi Ph\), \(b\ge0\), and \(a=q(b-\ell)\). Conversely these inequalities immediately imply the projected bound. Exactly the same proof applies to an old owner, replacing continuation wealth by the net estate and setting later income to zero.

At a binding borrowing constraint, \(\ell=\phi Ph\), \(b=0\), and the cash needed besides consumption is \((1+\tau-q\phi)Ph\). For example, \(q=.5,\phi=.8\) gives a 20% origination deposit **and** a 40%-of-price PV reserve for interest: 60% of price must be supplied from current resources before tax. The net equity at sale remains 20%. Calling this only a 20% down-payment model would obscure the mechanism.

This schedule is coherent as specified. It explicitly holds the chosen house through the utility age and allows sale/resizing and fresh secured borrowing at the boundary for either age. It does not impose an extra no-upgrade rule between the two modeled ages or silently prohibit old credit. It does, however, exclude intervening wage receipts, consumption/housing revisions, and the use of sale proceeds to pay a coupon that the model orders before sale. Letting an interest payment fall into the liquidation settlement would require rechecking the exact cash bound; the displayed result should not be attributed to every possible interest-only payment calendar.

## FOC and unrestricted feasibility

With budget multiplier \(\lambda\) and multiplier \(\mu\) on \(a+q\phi Ph\ge0\),

\[
\lambda=\beta V'(z)/q+\mu,
\qquad
u_h=(1+\tau)P\lambda-P\beta V'(z)-q\phi P\mu.
\]

Consequently

\[
u_h/u_c=p+q(1-\phi)P\mu/\lambda.
\]

For old wealth \(z\), the unrestricted allocation satisfies \(qe=\omega_Bz/K\), \(ph=\alpha z/K\). Its ownership requirement \(qe\ge\eta Ph\) is precisely \(\omega_Bd\ge\alpha\eta\). Thus its global value is indeed \(K\log z\) plus a price constant, including when rental tenure can tie.

For young L, \(qz_L=\beta KL\), and

\[
ph_L=\left(\alpha+\vartheta\frac{\kappa p}{\chi+\kappa p}\right)L
 <(\alpha+\vartheta)L.
\]

The second part of (F) therefore makes this complete unconstrained lifetime allocation financeable; it accounts for interest funding through the exact projection above.

## Primitive roots and global comparisons

Write \(A=W-pr\), \(x=A-\chi n_C\), \(s=r-\kappa n_C\), and

\[
D=\frac{\vartheta}{n_C^2}+\frac{\chi^2}{x^2}+\frac{\alpha\kappa^2}{s^2}>0.
\]

The unique interior fertility choice obeys

\[
n_C'=-\frac{\chi r}{x^2D}<0,\qquad
x'=-r\left(1-\frac{\chi^2}{x^2D}\right)<0,\qquad s'>0.
\]

Hence \(\mathcal M'=\alpha(x's-xs')/s^2<0\). Also \(\mathcal M(0)>0\) and \(\mathcal M(p)\to0\) as \(p\uparrow W/r\). Every equation \(\mathcal M(p)=tp\), \(t>0\), has one root in that interval. Substitution of \(x=tp s/\alpha\) and \(n_C=\vartheta x/(\chi+\kappa tp)\) gives exactly the response's polynomial:

\[
\kappa tr(t+\alpha+\vartheta)p^2+
\left[\chi r\{\alpha+t(1+\vartheta)\}-\kappa t(\alpha+\vartheta)W\right]p
-\alpha\chi W=0.
\]

Its leading coefficient is positive and constant negative, so it has exactly one positive root. Since \(k>1\), \(p_-<p_+\). From \(L\ge W\),

\[
qY\ge(\alpha+\vartheta+\beta K)W>\beta KW,
\]

which verifies \(0<\sigma<1\). Both \(n_C\) and \(n_L\) decrease strictly with rent, so (P) gives one replacement root and a finite positive housing-clearing cohort mass.

For C, rental saving has derivative \(-\zeta<0\), the housing ceiling has multiplier \(\xi>0\), and the fertility derivative is zero. Concavity verifies the global rental optimum. For **any** feasible owner allocation, the common consolidated budget gives exactly

\[
U_O-U_C\le\xi(h-r)-q\zeta(z-Y).
\]

Since \(q(z-Y)\ge\eta Ph\), \(p>p_-\), and \(\zeta\ge\sigma/x\), the right side is at most \(-\xi r<0\). No fixed-fertility restriction is smuggled into this comparison.

For L, the C first-order condition together with \(\mathcal M(p)>p\) yields

\[
n_C<\frac{\vartheta x_C}{\chi+\kappa p}<n_L,
\qquad \frac{\alpha L}{p}>s_C.
\]

Therefore \(h_L=\alpha L/p+\kappa n_L>r\). Strict concavity of the unrestricted lifetime problem makes its financeable ownership allocation strictly better than every rental alternative.

Finally, a fixed allocation determines only \(a=q(b-\ell)\). Any

\[
\ell\in[\max\{0,-a/q\},\phi Ph],\qquad b=a/q+\ell,
\]

implements it. Unless that interval degenerates, gross mortgages and saving are indeterminate because their rates coincide. This is an additional reason to qualify uniqueness as the regime's price, population, and real allocation.

## Narrow settlement corollary: final coupon paid at closing

**Verified as a limited, distinct corollary.** Permit only the last scheduled coupon, \(I_m=(R_m-1)\ell\), to be paid from the sale proceeds. Earlier coupons must still be funded from nonnegative liquid balances. Here \(b\) denotes the liquid balance **net of the final closing coupon**; a negative \(b\) is a closing-account entry, not unsecured borrowing before sale. Actual liquid cash entering the closing is \(b+I_m\ge0\).

The exact feasible projection becomes

\[
b\ge-(R_m-1)\ell,
\qquad
a=q(b-\ell)\ge-D_{m-1}\ell\ge-D_{m-1}\phi Ph.
\]

Sufficiency is constructive: choose

\[
\ell=\max\{0,-a/D_{m-1}\},\qquad b=a/q+\ell.
\]

If \(a<0\), this gives \(b=-(R_m-1)\ell\), with saving after each earlier coupon exactly

\[
S_j=\left(1-\frac{D_{m-1}}{D_j}\right)\ell\ge0,
\qquad j=0,\ldots,m-1.
\]

If \(a\ge0\), choose no mortgage and nonnegative saving. Hence no earlier payment is secretly financed from later sale proceeds. Sale then pays principal and the final coupon together. The equivalent terminal restriction and owner wedge are

\[
z-Y\ge(1-R_m\phi)Ph,
\qquad
\frac{u_h}{u_c}=p+\eta_{\mathrm{close}}P\frac\mu\lambda,
\qquad
\eta_{\mathrm{close}}=q(1-R_m\phi).
\]

The original stationary proof carries through with this replacement in (F), \(k\), and the tenure bounds **provided \(R_m\phi<1\)**, with the same convention at both ages. Minimum current owner funds besides consumption become \((1+\tau-D_{m-1}\phi)Ph=ph+\eta_{\mathrm{close}}Ph\).

This is a restriction on the final coupon interval, not the full utility age. It permits \(\phi\ge q\) when that final interval is sufficiently short. It is not universal over coupon calendars: with a single coupon only at sale, \(D_{m-1}=1\), \(R_m=1/q\), and the restriction reverts to \(\phi<q\). If \(R_m\phi\ge1\), the renter-to-owner replication obstruction returns. References to a “100% balance” in this corollary must count the final coupon together with principal: total closing debt reaches the house value at \(\phi=1/R_m\), not at \(\phi=1\).

</details>

<details>
<summary>Earlier answer: checked results and the financing problem that prompted the follow-up</summary>

The first response contains a coherent conditional allocation and fertility argument. It is not yet a satisfactory replacement for the paper's theory: its mortgage restriction excludes the earlier long-period financing benchmark. One focused follow-up is running on that issue. No new model specification has been adopted.

[Full Pro response](oracle_essential_theory_response.md) · [Pro conversation](https://chatgpt.com/c/6aa1e33d-d328-83ea-971e-0541c9e3e359) · [Exact follow-up](../../../docs/prompts/oracle_essential_theory_finance_followup.md)

## What the proposed model explains

Two types have the same present-value resources. One receives \(W\) while young and a known, nonpledgeable payment \(Y\) when old. The other receives \(W+qY\) while young. Here \(q\) is the price of one unit of next-period goods. The first type therefore has fewer resources available for the down payment, without being poorer over its lifetime.

Both types choose housing, consumption, fertility and saving. A common divisible housing stock can be occupied by owners or renters; renters cannot choose more than \(r\) units. Old households can sell, resize and borrow against housing. The mortgage share is \(\phi\).

The proposed equilibrium has constrained young households renting at \(r\), liquid young households owning larger homes, and unconstrained old housing choices. This assignment is derived under sufficient primitive inequalities. It is not assumed to be every equilibrium of the model.

## The result, separated from its assumptions

Let \(x_i=c_i^y-\chi n_i\) denote consumption beyond children's goods needs, and let \(s_i=h_i^y-\kappa n_i\) denote housing beyond children's space needs. Equal current utility weights and the same housing coefficient \(\alpha\) at both ages imply that the dated planner pools adult goods and adult space:

\[
X=\frac{\bar x+\bar c_o}{2},\qquad
S=\frac{\bar s+\bar h_o}{2}.
\]

The planner chooses all current consumption and housing, initially fixes individual fertility, and preserves young continuation wealth and old net estates. Hence

\[
c_i^{y,F}=X+\chi n_i,\quad h_i^{y,F}=S+\kappa n_i,
\qquad c_j^{o,F}=X,\quad h_j^{o,F}=S.
\]

The equilibrium conditions establish

\[
\bar c_o\geq\bar x,\qquad
\bar h_o>\bar s.
\]

Consequently,

\[
\bar h_y^F-\bar h_y^E
=\frac{\bar h_o-\bar s}{2}>0.
\]

Every constrained young family receives more consumption and housing. The allocation may reduce liquid young families' utility and old households' utility; this is an equally weighted utilitarian comparison.

There is also a separate efficiency argument that does not rely on equal welfare weights. In this equilibrium, constrained renters value housing relative to goods more highly than unconstrained old households:

\[
\frac{\alpha x_C}{s_C}>p
=\frac{\alpha c_o}{h_o}.
\]

A small housing transfer can therefore compensate an old donor in goods while improving the constrained family's utility, keeping continuation wealth and estates fixed. The authority must relax the relevant private financial constraints. This is not a transfers-only policy theorem.

## Why the fertility result is useful

For the proposed Stone–Geary preferences, both additional consumption and additional housing raise the marginal utility of fertility. Each constrained family would therefore choose more children if offered its fixed-fertility planner bundle.

The joint planner result is stronger. Private fertility satisfies

\[
n_i=\vartheta
\left(\frac{\chi}{x_i}+\frac{\alpha\kappa}{s_i}\right)^{-1}.
\]

The reciprocal on the right is concave in \((x_i,s_i)\). Averaging private choices and using the fact that the planner gives the young more average adult goods and space yields

\[
\frac{\vartheta}{\bar n}
\geq \frac{\chi}{\bar x}+\frac{\alpha\kappa}{\bar s}
>
\frac{\chi}{X(\bar n)}+\frac{\alpha\kappa}{S(\bar n)}.
\]

The joint planner's fertility first-order condition then implies \(n^F>\bar n\). This proof does not require the child-goods cost to be negligible. The current parents' utility is the objective; unborn households receive no welfare weight.

## What the primitive conditions do

The source's inequalities have three separate roles:

1. The young with deferred income can afford the largest rental home but cannot finance buying that much housing. The unrestricted rental demand exceeds the rental cap.
2. Later resources and liquid households' resources are large enough to finance unconstrained old choices, make constrained young saving zero, and ensure that old average adult consumption is at least young average adult consumption.
3. Replacement fertility occurs inside the price interval where those household choices hold. A strictly decreasing fertility function determines rent; housing clearing determines the stationary population level.

These are genuine sufficient restrictions on primitives. They remain substantial. In particular, the price-window condition and liquid-young financing condition imply

\[
\frac{q-\phi}{1+\tau-q}(\alpha+\vartheta)>1,
\qquad
\beta(1+\alpha+\omega_B)>1.
\]

The response does not establish empirical mildness. It also does not establish that every old household occupies more total housing than every young household. The theorem compares old housing with young housing net of children's needs. Within the liquid type, old housing is smaller than its own earlier young home when \(\beta<q\).

The separate bound \(\beta<q\) is less fundamental. For the stationary allocation and dated fertility results, the proof extends beyond it if constrained young nonsaving is ensured directly by \(qY/(\beta K)\geq b\), retaining the other finance, resource and price-window conditions. Here \(K=1+\alpha+\omega_B\) and \(b=(q-\phi)W/(1+\tau-\phi)\). This extension has been checked for the core result only. It does not remove the \(\phi<q\) issue.

Uniqueness must refer to prices, real allocations and population within the specified regime. It need not include old tenure labels or portfolios: if an old household's optimal home is below the rental cap, renting and owning can tie.

## The decisive financing problem

The owner budget in Pro's proposed model is

\[
c+(1+\tau)Ph+a=w,\qquad a\geq-\phi Ph,
\qquad z=Y_{\rm next}+a/q+Ph.
\]

Let \(\lambda=u_c>0\), and let \(\mu\geq0\) be the multiplier on its collateral constraint. The housing first-order condition gives the exact identity

\[
\frac{u_h}{u_c}
=p+(q-\phi)P\frac{\mu}{\lambda},
\qquad p=(1+\tau-q)P.
\]

Thus the sign of the financial distortion depends on \(q-\phi\). If \(\phi>q\), a binding collateral constraint makes housing a way to obtain current liquidity and puts the marginal housing valuation below the user cost. This is a feature of these repayment and consumption dates, not a universal claim about mortgages.

The renter argument fails for the same reason. When \(\phi\geq q\), any renter plan can be reproduced by owning with \(a_O=a_R-qPh\), leaving consumption and continuation wealth unchanged and satisfying the mortgage limit. A household with a strictly binding rental cap then benefits from choosing a slightly larger owned home.

Consequently, \(\phi<q\) is not just an inconvenient proof bound for Pro's argument. It is necessary for the desired positive housing distortion in this model. At \(\phi=0.8\), it requires \(q>0.8\), and excludes the prior \(q=0.5\) benchmark. The single follow-up asks whether a conventional repayment structure or a minimum owner size can remove this problem while retaining the intended mechanism and a simple proof.

## What the transition section adds

Within the same strict regime, Pro derives a policy that raises the maximum rental size. An additional goods-versus-space restriction,

\[
(2\alpha+\vartheta)\chi\leq\kappa p_\ell,
\]

signs the stationary response: a higher rental limit raises constrained families' fertility and the terminal population. A lower fertility taste lowers the terminal population. Both destinations have mean replacement fertility.

The local stability calculation and the initial capital-gain boundary pass independent checking. The announcement comparison retains inherited bonds and houses, rather than fixing the initial old's market-valued wealth. For a taste change, their inherited saving choice must also be retained for the first settlement.

The signed path implication is higher cumulative reproduction between the same initial population and the higher terminal population. Neither a higher fertility rate at every date nor an unconditional impact increase is proved. The policy increases access to rental space; it does not implement the result through property taxes or increase constrained young homeownership.

## Model changes that remain author decisions

Pro replaces the morning/afternoon idea with exogenous later resources, permits old borrowing, removes the owner upper bound, makes housing tastes equal across ages, and uses property taxes for additively valued public services rather than household rebates. Its common stock can move freely between tenure categories. These changes must remain visible; none has been integrated into the main deck or manuscript.

The comparison with Coven et al. is broadly correct: their analytical model focuses on housing as an asset and tax capitalization, with old liquidation; Proposition 2 gives cohort-specific welfare effects. Their quantitative model contains richer lifecycle income, an owner-size minimum, separate tenure markets and borrowing against old housing. Their quantitative fiscal system also includes endogenous household transfers; the proposed public-service-only closure should not be attributed to that model. [August 2026 paper, sections 3–4](https://abdouecon.github.io/research/papers/Property_Tax.pdf).

## Verification record

The complete response was captured from the visible conversation, including 142 inline/displayed mathematical objects. Whitespace-normalized text matches the browser extraction exactly (20,980 characters; FNV-1a 360915802). The lead derived the financial-distortion identity and checked the principal planner and fertility inequalities. Two separate Astra/max reviewers checked the stationary household/planner block and the policy/transition block. Their technical reviews are reproduced below as supporting checks.


<details>
<summary>Households, planner and fertility — independent review</summary>

# Independent review of the Pro liquidity proposal

Reviewed September 9, 2026. Scope: stationary households and tenure, primitive conditions (4)–(6), the full dated planner and settlement, and the two dated fertility claims. No transition, literature, calibration, or alternative-model audit. Source: `output/model/simplified_olg_amendments/oracle_essential_theory_response.md`, left unchanged. Mandatory project context was read. This is a review of an unadopted proposal.

## Verdict

| Claim | Verdict | Qualification or smallest repair |
|---|---|---|
| Stationary household quantities and positive replacement equilibrium | **PASS** | Conditional on every stated inequality and the external-finance closure. |
| Liquid young necessarily own a home larger than the rental ceiling | **PASS** | This follows from the price window and \(L\ge b\); the response omits the short proof. |
| Old choices are financeable | **PASS** | Old households may rent or own when their optimal house is at most \(r\). Positive financial estates are not implied. |
| Equilibrium is unique within the regime | **FAIL literally; small wording repair** | Prices, housing/consumption quantities, fertility, and cohort mass are unique. Old tenure labels and associated portfolios need not be unique. |
| Conditions are jointly nonvacuous | **PASS** | They can be strong. The price window forces \(\beta K>1\) and a sharper restriction than \(\phi<q\). For fixed demographic conversion, the replacement bracket also imposes an upper bound on deferred income. |
| Full dated utilitarian planner raises aggregate young housing and each C family's consumption/housing at fixed fertility | **PASS** | Both goods and housing are optimized. Individual continuation wealth and estates are fixed. |
| Financial settlement preserves resources and promised wealth | **PASS** | Intermediaries' matching bond changes must be included. The planner explicitly relaxes private financing limits. |
| C families' fertility rises when offered their fixed-fertility planner bundles | **PASS** | Their consumption and housing are held at the offered bundle and continuation opportunities remain fixed. |
| Joint dated, parents-only planner raises mean fertility | **PASS** | The harmonic-mean Jensen argument has the correct direction; the full planner can be shown to choose common fertility. |
| Ordinary LTV makes the desired financial housing wedge general | **FAIL as a broader interpretation** | This budget's positive wedge requires \(\phi<q\). When \(\phi>q\), constrained ownership has the opposite housing wedge. This is a limitation of the specified mortgage/timing structure, not a general impossibility theorem. |

## 1. Household and tenure checks

Write \(t=(q-\phi)/d>0\), \(A=1+\alpha+\vartheta\), and \(M=A+\beta K=G+\vartheta\). An owner satisfies the consolidated budget and mortgage condition

\[
c+ph+qz=w+qY_{\rm next},\qquad
q(z-Y_{\rm next})\ge(q-\phi)Ph.
\]

For the old, the unconstrained shares are exactly (9). Ownership can implement them iff \(\omega_B\ge\alpha t\). When \(h_o\le r\), renting implements exactly the same allocation with saving \(qe>0\). Ownership instead uses

\[
a_o=\left(\omega_B-\frac{q\alpha}{d}\right)c_o,
\]

which can be negative even when all proposition conditions hold. Thus the old are free to select tenure; their optimal quantities do not depend on that selection. Their value is \(K\log z\) plus a constant, as used in the young problem.

The liquid young's unrestricted shares are correct. The second finance inequality is sufficient because

\[
ph_L=L\left[\alpha+\frac{\vartheta}{1+\chi/(\kappa p)}\right]
<(\alpha+\vartheta)L.
\]

For a C owner, \(c+\delta Ph\le W\) implies \(h<W/(\delta P)<r\) whenever \(p>p_\ell\). Every owner plan is then replicated by a legal renter with financial saving

\[
a_R=q(z-Y)=a_O+qPh\ge(q-\phi)Ph>0.
\]

At the proposed cap choice \(x_C=W-pr-\chi n_C<b\). Condition (5), with \(\beta<q\), gives

\[
f(Y/K-b)\ge(1-f)(1-\beta/q)L>0.
\]

Thus \(Y/K>b>x_C\), and the derivative toward positive saving is strictly negative. The cap binds because \(p<p_u\). Strict concavity in consumption, adult space, fertility and continuation wealth establishes the global rental optimum; replication excludes all owner deviations.

The omitted proof that liquid young own larger homes is short. Since \(\widehat h_C(p_\ell)>r\),

\[
t(\alpha+\vartheta)>1,
\qquad
b=\frac{Wt}{1+t}>\frac W{1+\alpha+\vartheta}.
\]

Consequently \(L\ge b>W/A\), so

\[
h_L=\frac{AL}{W}\widehat h_C(p)>\widehat h_C(p)>r.
\]

The fertility root is unique: \(F_p<0\), both bracket endpoints are in the feasible domain, and (6) places the root strictly between them. Equation (7) then gives the unique positive stationary cohort mass. Equal young/old masses are justified by replacement. External bonds and a residual rental intermediary are maintained closures, so no additional domestic bond-clearing condition has been omitted.

**Tenure-uniqueness repair.** Say “the equilibrium prices, quantities, fertility and population are unique within this regime, with old tenure potentially indeterminate.” The hypotheses do not force old ownership. For example, the exact internal witness

\[
\alpha=\vartheta=W=r=\kappa=1,\quad\chi=1/10,
\quad q=9/10,\ \phi=4/5,\ \tau=1/100,
\quad\beta=1/2,\ \omega_B=3,\ f=9/10,\ Y=14/5
\]

at \(p=3/5\), \(\nu=1/F(p)\), satisfies all inequalities. Here \(p_\ell=11/21<p\), \(\widehat h_C(p)=65/63>r\), and the old housing choices are \(14/15\) and \(16/27\), both below \(r\). Their financial positions satisfy \(a_o=-(57/11)c_o<0\). This witness only checks nonvacuity and refutes strict tenure uniqueness; it is not a recommended parameter family.

## 2. Exact economic force and nonvacuity

Let \(z_\ell=\chi/(\kappa p_\ell)>0\). The first price-window inequality is exactly

\[
\frac{q-\phi}{d}>
\frac{1+(1+\vartheta)z_\ell}{\alpha+\vartheta+\alpha z_\ell}
>\frac1{\alpha+\vartheta}.
\]

It therefore requires

\[
q>\frac{1+\tau+\phi(\alpha+\vartheta)}{1+\alpha+\vartheta},
\qquad
\beta K\ge\frac{q-\phi}{d}(\alpha+\vartheta)>1.
\]

The second implication does not contradict \(\beta/q<1\): a sufficiently large housing or estate weight can make \(\beta K>1\). But the response should not suggest that the impatient region is unrestricted. For \(\phi=.8\), \(q>.8\) is necessary, and the displayed sharper bound may be considerably higher. This is consequential when an age lasts many years.

Conditions (5) impose real restrictions on type composition and deferred income. At fixed other primitives they admit a finite \(Y\) iff

\[
f>f_{\min}:=\frac{(q-\beta)K}{A+qK}<1.
\]

Conditional on that inequality, their exact lower bounds are

\[
Y\ge\frac{Mb-W}{q},\qquad
Y\ge
\frac{fbM+(1-f)(1-\beta/q)W}
{fM/K-(1-f)(q-\beta)}.
\]

Thus there is a nonempty lifecycle-resource region, although small \(f\) cannot always be offset by arbitrarily large \(Y\): liquid types also receive the present value of that larger deferred endowment.

For a **fixed** \(\nu\), replacement adds a finite interval. Define

\[
Y_{\rm rep}(p)=\frac1q\left\{
\frac{M(\chi+\kappa p)}{(1-f)\vartheta}
\left[\frac1\nu-fn_C(p)\right]-W\right\}.
\]

Condition (6) requires \(Y_{\rm rep}(p_\ell)<Y<Y_{\rm rep}(p_u)\). This interval must intersect the lifecycle lower bounds; “sufficiently large deferred income” alone does not prove replacement for a predetermined \(\nu\). Joint nonvacuity is nevertheless exact: after satisfying finance and lifecycle restrictions, the nonempty fertility bracket permits a \(\nu\); for fixed \(\nu\), joint scaling of \(\chi,\kappa\) rescales fertility without changing housing demands or the price window. None of this makes the region empirically mild.

**Structural sign limitation.** If \(\phi\ge q\), any rental plan is replicated by ownership using \(a_O=a_R-qPh\). At a binding rental cap, an owner can then expand housing, so the capped-renter branch is impossible. More generally, for a borrowing-constrained young owner let

\[
\zeta=\frac1x-\frac{\beta K}{qz}\ge0.
\]

The exact housing first-order condition is

\[
\frac{\alpha x}{h-\kappa n}
=p+(q-\phi)Px\zeta.
\]

Hence \(\phi>q\) gives a negative housing wedge: housing collateral helps households obtain current resources. Old households are unrestricted in this region because their positive net estate automatically satisfies the consolidated mortgage constraint. The source's positive young housing wedge cannot be extended to this region by merely admitting C ownership. An equally weighted redistributive planner might still favor the young, but that would be a different argument.

**Small optional generalization.** The core does not otherwise need \(\beta<q\). Replace its use in proving C nonsaving by the direct condition \(qY/(\beta K)\ge b\), and retain finance, both resource inequalities and the price bracket. All stationary allocation and dated fertility proofs continue to hold. This statement does not extend the policy or transition proof.

## 3. Full dated planner and financial settlement

At fixed fertility and individual promised wealth, equal current household weights give common adult consumption and adult space. Strict concavity yields exactly

\[
X=(\bar x+B_o)/2,\qquad S=(\bar s+\bar h_o)/2.
\]

Here the resource condition is strictly stronger than the needed comparison because \(B_o\ge fb+(1-f)L>\bar x\). Also

\[
s_C<\alpha x_C/p<\alpha L/p=s_L,\qquad
\bar s<\alpha\bar x/p<\bar h_o.
\]

Therefore aggregate young housing rises by \((\bar h_o-\bar s)/2\) per young household, and every C household receives \(X>x_C\) and \(S>s_C\), hence more total consumption and housing at its own fixed fertility. Total old housing need not exceed total young housing before the intervention; the relevant comparison is old housing versus young adult space.

The financial settlement is valid. For each household, keep promised wealth fixed and set

\[
\Delta a_i=-qP\Delta h_i^{\rm owned},\qquad
T_i=\Delta c_i+p\Delta h_i.
\]

Let intermediary rental inventory be \(I=\bar H-\sum_i h_i^{\rm owned}\). Its bond position can be written \(a_I=-qPI\). Then

\[
\Delta a_I=qP\sum_i\Delta h_i^{\rm owned}
=-\sum_i\Delta a_i.
\]

Aggregate financial claims and goods are unchanged, and \(\sum_iT_i=0\). The tax base is the same total stock, so public spending remains fixed. Assigning ownership to allocations above \(r\) is feasible because the stock is common and owner housing has no upper bound. The planner relaxes private finance; this is not a transfer-only implementation with the original limits intact.

The separate small compensated housing transfer also has the correct derivative. It establishes a dated fixed-fertility gain under the relaxed financing permissions, not that the full utilitarian allocation is Pareto superior.

## 4. Both dated fertility claims

At a fixed offered bundle, the fertility derivative is strictly decreasing in \(n\), and its derivatives with respect to both \(c\) and \(h\) are positive. Every C family's fixed-fertility planner bundle therefore leads to a strictly larger privately chosen \(n\), holding the bundle and continuation opportunity fixed.

For the joint planner, common fertility is a result, not an extra restriction. Conditional on mean fertility, resources depend only on that mean; equalizing fertility maximizes the sum of \(\vartheta\log n_i\). Consumption and adult space are also equalized. Thus (12) describes the full optimum among current parents.

The function

\[
H(x,s)=\left(\frac\chi x+\frac{\alpha\kappa}s\right)^{-1}
=\frac{xs}{\chi s+\alpha\kappa x}
\]

is concave on positive \((x,s)\). Its Hessian is

\[
\nabla^2H=-\frac{2\chi\alpha\kappa}
{(\alpha\kappa x+\chi s)^3}
\begin{pmatrix}s^2&-xs\\-xs&x^2\end{pmatrix}\preceq0.
\]

Private first-order conditions imply \(\bar n\le\vartheta H(\bar x,\bar s)\). Since \(X(\bar n)>\bar x\) and \(S(\bar n)>\bar s\), the planner's fertility derivative at \(\bar n\) is strictly positive. Its derivative with respect to \(n\) is

\[
-\frac{\vartheta}{n^2}
-\frac{\chi^2}{2X(n)^2}
-\frac{\alpha\kappa^2}{2S(n)^2}<0.
\]

It has a unique interior zero, establishing \(n^F>\bar n\). This part requires no small child-goods-cost restriction. It remains a dated parental welfare comparison, not a dynamic population optimum.

</details>


<details>
<summary>Policy and transition — independent review</summary>

# Independent review: stationary policy and local transition

Reviewed 2026-09-09. Source: `output/model/simplified_olg_amendments/oracle_essential_theory_response.md`, particularly lines 360–415 and Appendices B–C, lines 490–671. Scope is policy, demographic dynamics and the surprise-price boundary only. The allocation/planner theorem and literature were not independently adjudicated here. No numerical examples, sweeps or simulations were used as proof.

## Verdict

**PASS, conditional on the stated strict household regime.** The stationary rental-cap and fertility-taste signs, Appendix B's sufficient inequalities, the demographic Jacobian and Jury tests, and the announcement-price derivative are analytically correct. The latter derivative correctly corresponds to *fixed initial cohort sizes*; the first perturbed future rent satisfies a different initialization from a generic demographic-state perturbation.

**One explicit presentation repair is necessary for a usable transition statement:** a change in the fertility preference changes the liquid young's desired saving. Consequently date-zero old wealth differs from the new stationary old-wealth coefficient even with no capital gain. Equation (14), if used literally with the inherited bond and housing holdings, already includes this correction. It must not be replaced by new-parameter stationary old wealth plus the capital gain alone. The two-dimensional autonomous population map applies from date 1 after a permanent date-zero change; date 0 needs the separate initial balance-sheet equation below.

The evidence supports local feasible convergent paths and the signed cumulative-reproduction comparison. It does not give a signed impact-fertility effect, all-date ordering, global uniqueness or convergence, large reforms crossing regimes, or an endogenous inheritance/type-transmission model.

## 1. Definitions and the stationary rental-access derivative

Write \(\theta=\vartheta\), \(K=1+\alpha+\omega_B\), \(G=1+\alpha+\beta K\), and \(L=(W+qY)/(G+\theta)\). The mean consumption resources of the old are

\[
B_o=fY/K+(1-f)\beta L/q.
\]

Per-young-cohort and per-old-cohort housing demand are

\[
D_y=fr+(1-f)\left(\kappa n_L+\frac{\alpha L}{p}\right),\qquad
D_o=\frac{\alpha B_o}{p},\qquad D=D_y+D_o.
\]

The source's household implicit derivatives, equations (503)–(506), follow directly from differentiating its fertility first-order condition. In particular,

\[
n_{C,p}=-\frac{\chi r}{x_C^2Q},\qquad
n_{C,r}=\frac{\alpha\kappa/s_C^2-\chi p/x_C^2}{Q}.
\]

At fixed price, \(D_r=f\) and \(F_r=fn_{C,r}\). Thus, on the replacement root \(F(p^*,r)=1/\nu\),

\[
p_r^*=\frac{fn_{C,r}}{-F_p},\qquad
\frac{dD^*}{dr}=f+D_pp_r^*.
\]

Therefore \(N_r^*>0\) is equivalent to the source's inequality (15), \((-D_p)n_{C,r}>-F_p\).

### Exact verification of the dimensionless bounds

Use the source's \(z,c,s,m,u,a\), where \(s+m=1\), and define

\[
\widetilde Q=\frac{\theta}{m^2}+\frac{z^2}{c^2}+\frac{\alpha}{s^2},\qquad
B=\frac{\alpha}{s^2}-\frac z{c^2}.
\]

The cap and fertility conditions imply

\[
\frac\theta m=\frac zc+\frac\alpha s,\qquad
\alpha c>s,\qquad
m<\frac\theta{\alpha+\theta},\qquad
c>\frac1{\alpha+\theta}.
\]

The last two follow because the positive term \(z/c\) implies \(\theta s>\alpha m\). The resource condition implies \(a>\alpha(fc+2u)\) because \(b>x_C\) and \(f>0\); the weak bound used by the source is sufficient.

The difference in (15), multiplied by the positive factor \(\kappa p\widetilde Q/r\), is exactly

\[
J=aB-\frac{fz}{c^2}
-\frac{u\theta}{(1+z)^2}
 \left(\frac\theta{m^2}+\frac{z(1+z)}{c^2}\right).
\]

Condition (13) and \(p>p_\ell\) imply \(z<1/(2\alpha+\theta)\). In particular \(B>0\), because \(\alpha c>s\) gives \(\alpha c^2>s^2/\alpha>zs^2\).

After inserting the bound on \(a\), the coefficient of \(f\) is

\[
\alpha cB-z/c^2
>\frac{c(1-\alpha z)-z}{c^2}>0.
\]

For the final inequality, \(1-\alpha z>0\) and

\[
c(1-\alpha z)-z
>\frac{1-z(2\alpha+\theta)}{\alpha+\theta}\geq0.
\]

For the coefficient of \(u\), put \(v=\alpha c/s>1\). Multiplication by \(c^2\) gives

\[
2v^2-\frac{(z+v)^2}{(1+z)^2}
-2\alpha z-\frac{\theta z}{1+z}.
\]

The first two terms form an increasing function of \(v\geq1\), with value 1 at \(v=1\). Hence the displayed expression exceeds

\[
1-z(2\alpha+\theta)\geq0.
\]

Thus \(J>0\), \(p_r^*>0\) and \(N_r^*>0\). The equilibrium cash-poor fertility derivative is also correctly reported:

\[
\frac{dn_C^*}{dr}
=n_{C,r}\frac{(1-f)n_{L,p}}{F_p}>0.
\]

Old consumption resources are fixed in this comparison, so old housing falls with rent. The source's expression for \(pD_y\) increases in both \(p\) and \(r\), whereas \(pD_o=\alpha B_o\) is fixed, proving a larger young housing share.

## 2. Stationary taste sign

Direct differentiation gives

\[
n_{C,\theta}=\frac1{n_CQ}>0,\qquad
n_{L,\theta}=\frac{(W+qY)G}{(G+\theta)^2(\chi+\kappa p)}>0.
\]

Thus \(F_\theta>0\) and \(p_\theta^*>0\). To verify the important inequality \(-D_p>\kappa(-F_p)\), observe that it is equivalent to

\[
a>\frac{fz}{c^2\widetilde Q}.
\]

But

\[
\frac{z}{c^2\widetilde Q}
<\frac{zs^2}{\alpha c^2}<\alpha z<\alpha c,
\]

where the final inequality uses \(z<1/(2\alpha+\theta)<1/(\alpha+\theta)<c\). The lower bound on \(a\) proves the result.

Since \(L_\theta<0\) and \(B_{o,\theta}=(1-f)\beta L_\theta/q<0\),

\[
D_\theta<(1-f)\kappa n_{L,\theta}.
\]

Consequently

\[
\frac{dD^*}{d\theta}
<(1-f)\kappa n_{L,\theta}-\kappa F_\theta
=-\kappa f n_{C,\theta}<0,
\]

so \(N_\theta^*>0\). A permanent taste decline lowers the stationary population, conditional on remaining on this branch.

## 3. Demographic dynamics and stability

With constant post-announcement parameters and all later choices anticipated, the liquid young choose continuation resources

\[
z_{L,t+1}=\beta KL/q
\]

independently of current or future prices. Cash-poor young choose zero financial saving and have old resources \(Y\). Because entrant type shares are fixed independently of parental type, the two cohort sizes suffice after the special initial date.

Linearizing static clearing in log cohort sizes gives

\[
\widehat p_t=A_y\widehat N_t^y+A_o\widehat N_t^o.
\]

The source's demographic matrix is therefore exactly

\[
M=\begin{pmatrix}1-\varepsilon A_y&-\varepsilon A_o\\1&0\end{pmatrix}.
\]

For cash-poor fertility,

\[
\varepsilon_C=\frac{z}{mc^2\widetilde Q}
<\frac{z}{(1+z)c}<1.
\]

Indeed \(mc\widetilde Q=z+\alpha c/s+mz^2/c+\alpha mc/s^2>1+z\); the last bound follows from the inequalities for \(c,z\) above. Liquid fertility has elasticity \(1/(1+z)<1\), and aggregate elasticity is their positive fertility-weighted average.

Also \(-pD_p>D_o\), hence \(0<A_o<1\). Finally,

\[
D_y-D_o-\kappa F
=fs_C+\frac{\alpha((1-f)L-B_o)}p<0.
\]

Together with \(-D_p>\kappa(-F_p)\), this yields \(\varepsilon(A_y-A_o)<1\). The characteristic polynomial is

\[
\lambda^2-(1-\varepsilon A_y)\lambda+\varepsilon A_o.
\]

All three Jury expressions are strictly positive:

\[
\varepsilon(A_y+A_o)>0,\qquad
2-\varepsilon(A_y-A_o)>1,\qquad
1-\varepsilon A_o>0.
\]

Hence the nonlinear autonomous demographic map is locally asymptotically stable, with geometric convergence in a sufficiently small neighborhood.

## 4. The initial balance sheet, including the saving correction

Let \(\theta_-\) describe a pre-shock stationary economy and \(\theta_+\) the permanent new preference. Let \(L_-\) and \(L_+\) be the corresponding expenditure shares. Fix the inherited bonds \(a_{L,-1}\), house \(h_{L,-1}\), and both initial cohort sizes. Then

\[
z_{L,0}=\frac{a_{L,-1}}q+P_0h_{L,-1}
=\frac{\beta KL_-}q+h_{L,-1}(P_0-P_-).
\]

Thus, **relative to the new stationary old resources**, the full discrepancy is

\[
\delta z_{L,0}
=\frac{\beta K}{q}(L_--L_+)
 +h_{L,-1}(P_0-P_-).
\]

The first term is the required one-period saving-state correction. A taste decline raises \(L_+\), so inherited saving alone leaves the initial old liquid type below the new stationary old resources; capital gains are a separate term with their own sign.

The exact initial market equation is

\[
\bar H=N_0^yD_y(p_0;\theta_+,r_+)
+\frac{\alpha N_0^o}{Kp_0}
\left[fY+(1-f)\left(\frac{a_{L,-1}}q+P_0h_{L,-1}\right)\right].
\]

Equivalently, use new-parameter \(D_o\) and add \(\alpha N_0^o(1-f)\delta z_{L,0}/(Kp_0)\). In a derivative with respect to \(\theta_+\), the \(-\beta KL_{+,\theta}/q\) part of this correction cancels the apparent direct change in initial old resources introduced by \(D_o(p_0;\theta_+)\). Omitting it would give incorrect impact derivatives, even though the post-impact stability calculation is unaffected.

At date 1 the new liquid young's optimized continuation resources are already \(\beta KL_+/q\), so no persistent saving state remains for a permanent surprise. For a changing preference sequence, the ordinary old resources instead use \(L(\theta_{t-1})\), not \(L(\theta_t)\). At a later surprise along a baseline transition, replace \(P_-\) by the inherited contract's previously expected date-\(T\) price and retain the actual inherited portfolio. These are direct consequences of equation (14), not new economic assumptions.

## 5. Announcement-price derivative with fixed initial cohorts

For the derivative with respect to the candidate \(p_0\), hold all primitives and inherited cohort sizes fixed. Write \(x_t=dp_t/dp_0\), evaluated at a stationary reference. The initial cohorts do not change, so

\[
x_0=1,\qquad
\frac{d\log N_1^y}{d\log p_0}=-\varepsilon,\qquad
\frac{d\log N_1^o}{d\log p_0}=0,
\qquad x_1=-\varepsilon A_y.
\]

In particular, \(x_1\) is **not** \(1-\varepsilon A_y\). For later dates,

\[
x_t=-\varepsilon(A_y,A_o)M^{t-1}(1,0)'\quad(t\geq1).
\]

Using the bounded/no-bubble price solution and \(\rho=q/(1+\tau)\),

\[
\frac{\partial P_0}{\partial p_0}
=\frac1{1+\tau}\left[1-\rho\varepsilon(A_y,A_o)(I-\rho M)^{-1}(1,0)'\right]
=\frac{1-\rho}{(1+\tau)\Delta},
\]

\[
\Delta=1-\rho+\rho\varepsilon A_y+\rho^2\varepsilon A_o.
\]

This is exactly the source's expression and lies strictly between zero and \(1/(1+\tau)\).

At the same reference, initial static rent feedback from an exogenous \(P_0\) change is

\[
k_P=\frac{\alpha(1-f)h_L}{K(-pD_p)}.
\]

Since \(h_L=(L/p)[\alpha+\theta/(1+z)]\) and

\[
-pD_p>(1-f)\frac Lp\left[\alpha+\frac\theta{(1+z)^2}\right],
\]

\[
k_P<\frac\alpha K\frac{\alpha+\theta/(1+z)}{\alpha+\theta/(1+z)^2}
\leq\frac{\alpha(1+z)}K<1.
\]

Therefore \(1-k_P\partial P_0/\partial p_0>0\): the initial market equation and forward discounted-price equation have a locally unique solution. This is a local implicit-function argument with the inherited *portfolio* fixed; it does not freeze inherited market-valued wealth.

## 6. Conditions and the cumulative comparison

The local continuation requires strict actual mortgage feasibility for old owners and liquid young, strict cash-poor exclusion from owning at the rental cap, strict cap binding and zero saving for cash-poor renters, and positive consumption/estate margins. The weak stationary inequalities (4) alone do not ensure a neighborhood with unchanged old financing status; use the response's stated strict-financing qualification. Along nonstationary paths the relevant owner constraint is

\[
q(z_{t+1}-Y_{t+1})\geq(qP_{t+1}-\phi P_t)h_t,
\]

so it must be checked by continuity from the strict stationary regime rather than replacing \(P_{t+1}\) with \(P_t\). Bounded convergent prices exclude the explosive asset-price solution. Finite types make preservation of uniform strict margins immediate for sufficiently small paths; no global shock radius is supplied.

The cumulative-reproduction identity is exact when baseline and policy share \(N_T^y\): summing their log entry laws telescopes to

\[
\sum_{t=T}^{\infty}\log\frac{F_t^{\rm policy}}{F_t^{\rm baseline}}
=\log\frac{N^{*,\rm policy}}{N^{*,\rm baseline}}.
\]

Local geometric convergence makes this sum well defined. Its positivity follows from the already-verified stationary rental-access sign. This statement includes the common inherited initial old through the announcement boundary, but supplies no separate impact or all-date fertility sign.

</details>

</details>

</details>

</details>

</details>
