# Conventional mortgage finance: a lifetime housing-reallocation theorem

Independent first-round mathematical review, 2026-09-07/08. This is a proposed
analytical result, not an amendment of the author manuscript. It uses the
ordinary mortgage in `docs/model/simplified_olg_conventional_finance_proposal.md`,
including repayment of the one-period loan at gross return \(1/q\). It does
not change that contract. The comparison fixes each household's fertility and
tenure, and therefore fixes every future cohort size.

## Result and scope

There is a useful positive result even when \(q<\phi\), without a small
discount factor or a numerical reference equilibrium. Write
\[
 K=1+\gamma+\omega_B,\quad B=\beta K,\quad
 p=(1-q+q\tau^p)P,\quad L=(1-\phi+q\tau^p)P,
 \quad \delta=p-L=(\phi-q)P>0.
\]
Here \(p\) is the beginning-of-period cost of housing services, and \(L\)
is the cash required per unit of owner housing. A sufficient condition is
\[
 \boxed{\quad B L\ge\delta
 \quad\Longleftrightarrow\quad
 \beta(1+\gamma+\omega_B)(1-\phi+q\tau^p)\ge\phi-q.\quad}
 \tag{1}
\]
At every strictly mortgage-constrained owner optimum with slack size caps and
the old estate floor slack, removing the young financing constraint, **holding
that household's fertility fixed**, then increases its desired young housing
and decreases its old-age housing. This is an exact comparison of lifetime
optima. It is not obtained by attaching an arbitrary positive housing transfer
to an unrelated gain from borrowing.

Under the explicit donor and recipient conditions below, a small positive
mass of such young owners can make this lifetime adjustment. Housing moves
from current old owners to current young owners. Compensating trades next
period absorb the participating households' reduced old-age housing. All
current old households and all future cohorts are weakly better off; selected
young households are strictly better off. All new external borrowing is
settled after two dates, and existing mortgage repayments are preserved.

The conclusion is a local Pareto reallocation with a specified housing
direction. It does not claim that every globally efficient allocation has more
young housing, that a uniform LTV policy is Pareto improving, or that young
housing has a higher current consumption-equivalent marginal value. The last
claim is false here.

## 1. Exact fixed-fertility household comparison

Take a stationary price \(P\) and rebate \(T\). For a selected owner, let
\(w=y^y+b+T\), \(v=y^o+T\), and let \(n\) be its actual equilibrium
fertility. Fix this \(n\) throughout the welfare comparison. Define
\[
 x=c-\chi n,\quad s=h-\kappa n,\quad S=qz,
 \quad a=w-(\chi+L\kappa)n,
 \quad M=w+qv-(\chi+p\kappa)n.
 \tag{2}
\]
The symbol \(a\) in this section is **adult cash after reserving fixed child
requirements**, not financial wealth. Both \(a\) and \(M\) are positive at
a feasible owner allocation. Conditional on slack old housing caps and estate
floor, the old value is \(K\log z+C_O\), so the relevant objective and budgets
are exactly
\[
 \log x+\alpha\log s+B\log S+\text{constant},\qquad
 x+ps+S=M,\qquad x+Ls\le a.
 \tag{3}
\]
The old estate floor is slack when
\[
 \omega_Bp>q\gamma P.
 \tag{4}
\]
In that regime the old choices are
\[
 c^2=z/K,\qquad h^2=\gamma z/(Kp),\qquad
 e=\omega_Bz/(Kq).
 \tag{5}
\]
An equilibrium allocation optimizing over fertility is also the optimum of
(3) conditional on its chosen \(n\). Thus this reduction does not replace the
equilibrium fertility rule.

With complete young finance and no young size cap, write \(D_0=1+\alpha+B\).
The resource-optimal comparison is
\[
 x_F=M/D_0,\qquad s_F=\alpha M/(pD_0),\qquad
 S_F=BM/D_0. \tag{6}
\]
The mortgage constraint strictly restricts (3) exactly when
\[
 \frac{M}{a}>\frac{D_0}{1+\alpha L/p}. \tag{7}
\]
At a strictly constrained optimum let \(\lambda>0\) and \(\mu>0\) be the
multipliers on lifetime resources and current cash. Put
\[
 r=\frac{\mu}{\lambda+\mu}\in(0,1),\qquad
 \rho=p-r\delta,\qquad j=(\lambda+\mu)^{-1}.
\]
The actual allocation is
\[
 x_C=j,\qquad s_C=\frac{\alpha j}{\rho},\qquad
 S_C=\frac{Bj}{1-r},\qquad
 M=j\left(1+\frac{\alpha p}{\rho}+\frac B{1-r}\right).
 \tag{8}
\]
Direct subtraction, using (6), gives
\[
 \boxed{\quad
 s_C-s_F=
 \frac{\alpha j r\{\delta(1-r)-BL\}}
 {pD_0\rho(1-r)}.
 \quad} \tag{9}
\]
Therefore (1) makes \(s_C<s_F\) for every \(r>0\), including equality in
(1). Total housing has the same sign because \(\kappa n\) is fixed.

There is also an exact income condition, useful when (1) fails:
\[
 \boxed{\quad h_F>h_C
 \iff LM>pa
 \iff Lqv>\delta(w-\chi n).
 \quad} \tag{10}
\]
All equivalences in (10) concern a **strictly constrained** allocation of
(3). Binding finance itself still depends on \(B\). In particular, the
absence of \(B\) from the last inequality is not a statement that discounting
does not matter.

One direct verification of (10) is to substitute \(x=a-Ls\) and
\(S=M-a-\delta s\) into (3). Its strictly concave objective has derivative
\[
 -\frac L{a-Ls}+\frac\alpha s
 -\frac{B\delta}{M-a-\delta s}.
\]
At \(s_F\), the derivative has the sign of \(pa-LM\), because the
relaxed cash excess is positive. If \(s_F\) exceeds the cash-feasible upper
bound, then \(s_C<s_F\) directly, and \(LM>pa\) also holds. This covers the
otherwise missing boundary case.

Old resources always decline under the relaxation, even if (1) or (10) fails.
Indeed, (8) implies
\[
 \frac{S_C}{M}
 =\frac B{B+(1-r)(1+\alpha p/\rho)}
 >\frac B{1+\alpha+B}=\frac{S_F}{M},
 \tag{11}
\]
since
\[
 1+\alpha-(1-r)(1+\alpha p/\rho)
 =r(1+\alpha L/\rho)>0.
\]
Equation (5) then makes old consumption, old housing and estate spending fall
in the complete-finance comparison. The household's lifetime welfare still
increases strictly because a genuinely binding constraint was removed.

### Young cap and finite adjustment

The complete-finance \(h_F\) need not fit the owner cap. If the actual young
owner has a strict cap margin, choose a fixed positive step toward (6), small
enough to retain that margin. More explicitly, when \(h_F>h_C\), choose
\[
 0<\theta\le
 \min\left\{\frac12,
 \frac{h_O^{\max}-h_C}{2(h_F-h_C)}\right\}.
 \tag{12}
\]
Use the convex combination of actual and comparison \((x,s,S)\). It respects
the lifetime budget, has greater young housing and lower old resources, and
strictly raises lifetime utility by strict concavity. Its old housing remains
below the actual old housing, so a previously slack old size cap stays slack.
The estate floor remains slack by (4). The whole argument therefore works
with the existing physical owner cap, without enlarging it.

## 2. A full lifetime Pareto construction

Assume a stationary equilibrium has the following groups. Finite income types
and continuously distributed ownership tastes are sufficient; a continuum of
income types is not needed.

1. A positive mass of strictly constrained young owners satisfies (1), or the
   exact condition (10), with strictly slack young and old housing caps and
   old estate floor. Their equilibrium fertility is held fixed.
2. A positive mass \(m_D\) of current old owners has slack owner cap and estate
   floor. Their common baseline bundle is \((c_D,h_D,e_D)\).
3. A disjoint positive mass \(m_R\) of current young owners can serve as old
   housing recipients next period. Their baseline old bundle
   \((c_R,h_R,e_R)\) has strict cap and estate-floor margins. Their young
   allocations are unchanged.

The subscripts \(D,R\) in this section mean donor and recipient, not renter.
Both groups are owners. Their old first-order conditions give
\(\gamma c_D/h_D=\gamma c_R/h_R=p\). Ownership taste does not change because
no household changes tenure. A positive owner group can be split into separate
selected and recipient subsets, so a special extra income type is unnecessary.

For each selected young household, let a finite step from (12) have increments
\(\Delta c,\Delta h,\Delta c^2,\Delta k,\Delta e\), where
\(\Delta h>0\) and \(\Delta k<0\). Let its utility gain be \(g>0\), and its
new adult consumption before a fee be \(x_A>0\). The lifetime budget implies
\[
 \Delta c+p\Delta h+q\Delta c^2+qp\Delta k+q^2\Delta e=0.
 \tag{13}
\]
Deduct current consumption
\[
 \eta=x_A(1-e^{-g/2})>0. \tag{14}
\]
This leaves exactly \(g/2\) of the lifetime gain and releases \(\eta\)
units of current goods per selected household.

Select a mass \(\epsilon>0\), leaving the recipient group untouched. At the
current date, give each selected household \(\Delta h\) extra housing and
take \(s_D=\epsilon\Delta h/m_D\) from each old donor. Keep each donor's
estate fixed and give it extra current consumption
\[
 D(s_D)=c_D\left[(1-s_D/h_D)^{-\gamma}-1\right]. \tag{15}
\]
This preserves its utility **exactly**, not merely to first order. Reducing
its housing relaxes its estate floor.

Next period, selected households use \(\Delta k<0\) less old housing. Give
each old recipient \(s_R=\epsilon(-\Delta k)/m_R\) extra housing, keep its
estate fixed, and collect consumption
\[
 R(s_R)=c_R\left[1-(1+s_R/h_R)^{-\gamma}\right]. \tag{16}
\]
Its old utility and hence its lifetime utility are exactly unchanged. Choose
\(\epsilon\) to retain its strict size-cap and estate-floor margins.

Housing clears at both affected dates. At the first date housing moves from
old to young owners. At the second, the selected households' smaller old homes
are absorbed by other old owners. Aggregate owner housing and renter housing
are unchanged at each date. The rental intermediary's stock and portfolio are
unchanged. No separate tenure-specific housing stock is introduced.

### Explicit goods bound

Let \(k=-\Delta k>0\). Choose \(s_D\le h_D/2\), and define
\[
 A_D=\frac{c_D\gamma(\gamma+1)2^{\gamma+2}}{h_D^2},\qquad
 A_R=\frac{c_R\gamma(\gamma+1)}{h_R^2},\qquad
 C=\frac{A_D(\Delta h)^2}{2m_D}
   +q\frac{A_Rk^2}{2m_R}.
 \tag{17}
\]
Taylor's theorem with the displayed global derivative bounds gives
\[
 m_DD(s_D)-p\epsilon\Delta h
 +q\{p\epsilon k-m_RR(s_R)\}\le C\epsilon^2. \tag{18}
\]
Combining (13) and the consumption fee, the present value of all additional
goods and terminal estate spending is therefore at most
\[
 C\epsilon^2-\eta\epsilon. \tag{19}
\]
It is strictly negative for a sufficiently small explicitly bounded positive
mass. For example, take \(\epsilon\) below the available selected mass and
below each of
\[
 \frac{\eta}{2C},\qquad
 \frac{m_Dh_D}{2\Delta h},\qquad
 \frac{m_R}{2k}
 \min\{h_O^{\max}-h_R,\ e_R/P-h_R\}.
 \tag{20}
\]
All numbers in (20) are positive under the stated conditions. This is an
analytical feasibility proof with an explicit mass bound, not existence in an
unspecified numerical neighborhood.

### External borrowing, repayment and estate endpoints

Let the current intervention date be \(t\). The additional dated goods uses,
after internal title payments cancel, are
\[
 G_0=\epsilon(\Delta c-\eta)+m_DD(s_D),\qquad
 G_1=\epsilon\Delta c^2-m_RR(s_R),\qquad
 G_2=\epsilon\Delta e.
 \tag{21}
\]
Existing old estates, existing mortgage repayments, and all other dated
external claims are held fixed. New external financing at \(t\) covers
\(G_0\); the next date balance is \(G_0/q+G_1\). At the estate date the final
balance is
\[
 \frac{G_0}{q^2}+\frac{G_1}{q}+G_2
 =\frac{G_0+qG_1+q^2G_2}{q^2}\le0. \tag{22}
\]
Thus new lending is repaid at the stipulated return. A negative final balance
is a surplus, which may be returned without hurting anyone. There is no
rollover beyond the selected households' estate date and no Ponzi argument.

Internal title trades can be settled at the original prices while preserving
each pre-existing private mortgage note. Current old donors replace the title
removed from their fixed estate with a bond. Selected owners' additional
bridge funding and later payments are part of (21)-(22). Next period the
offsetting old title trades preserve aggregate old death sales. The stock of
housing presented to future entrants and all of their allocations are
unchanged after the intervention.

Estates remain warm-glow expenditure. They are not assigned to entrants, used
as their initial wealth, or counted as a new resource twice. Selected
households may choose smaller estates as part of their strictly preferred
lifetime plan; this utility cost is already included in \(g\). Donors' and
recipients' estates stay exactly fixed. Property-tax receipts and common
rebates stay fixed because total housing is unchanged each date.

## 3. Primitive restrictions and stationary equilibrium coverage

The central restriction (1) is entirely in primitives and does not depend on
the equilibrium house price. At zero property tax it is
\[
 \beta(1+\gamma+\omega_B)\ge\frac{\phi-q}{1-\phi}. \tag{23}
\]
It asks that the weight on future utility be large enough relative to the
net extraction of future income permitted per unit of young cash. It does not
require a small \(\beta\); greater patience makes this sufficient housing
direction easier to obtain. It should nevertheless be reported as a real
economic parameter restriction, not assumed automatically.

Finance must also actually restrict a positive owner group. With uncapped
comparison choices, put \(E=1+\alpha+\vartheta\), \(D=E+B\),
\(\ell=L/p\). The exact primitive-at-prices binding condition already in the
proposal is
\[
 \frac{w+qv}{w}>
 \frac D{1+\alpha\ell+
  \vartheta(\chi+L\kappa)/(\chi+p\kappa)}. \tag{24}
\]
A stronger but price-uniform sufficient income condition is
\[
 \frac{w+qv}{w}>\frac D{1+\ell(\alpha+\vartheta)}, \tag{25}
\]
because \((\chi+L\kappa)/(\chi+p\kappa)\ge\ell\). At zero tax, \(w,v\)
are themselves primitive endowments, and both (23) and (25) are independent of
\(P\). For positive tax, include the equilibrium common rebate in \(w,v\);
it cannot silently be treated as a fixed independent endowment.

For a constrained owner, \(qz=qv-\delta h<qv\). Thus the sufficient
old-cap condition
\[
 \frac{\gamma v}{Kp}<h_O^{\max} \tag{26}
\]
ensures the actual old cap is slack and, by (11), also the comparison cap.
The young cap is slack if \(w/L<h_O^{\max}\). One may impose these conditions
for one positive-mass income type rather than the entire population. The
owner taste distribution has full support, so any type with finite tenure
values has a positive owner share; that group supplies selected households,
current old donors, and next-period old recipients in a stationary economy.

These inequalities are genuine restrictions on endowments, patience, finance
and physical size limits. They replace a bare endogenous MRS-gap assumption.
They do not, by themselves, prove that an arbitrary independently chosen
housing stock generates a stationary equilibrium satisfying (26) and the
young cap margins. Completing the fully primitive stationary-equilibrium
existence statement still requires an analytical housing-market construction
or price bounds. This review does not label an assumed equilibrium as a
proved existence result. The welfare theorem is complete conditional on that
explicit equilibrium coverage; an independent equilibrium-construction task
can combine with it.

## 4. Exact obstruction to a current housing-value proof

At a strictly constrained young owner with a slack housing cap,
\[
 \frac{u_h}{u_c}=\rho=p-r\delta<p. \tag{27}
\]
An old owner with slack floor and cap has MRS \(p\), and old borrowing or size
constraints weakly raise the relevant MRS above \(p\). Therefore a compensated
transfer of current housing from old owners to young owners, holding all of
the young recipients' future allocations fixed, cannot be justified by a
positive current MRS gap in the original \(q<\phi\) model.

The lifetime theorem does not contradict (27). It changes the selected
household's old resources and housing, respects the associated intertemporal
cost, and exploits the actual complete-finance lifetime optimum. The young
household trades some later consumption, housing and estate expenditure for
its preferred current bundle.

If (1) fails, there are strictly constrained owner configurations in which
finance relief lowers young housing. Equations (9)-(10) characterize the
failure exactly. It is therefore incorrect to assert a universal old-to-young
allocation direction from a binding ordinary mortgage constraint alone.

## 5. Renting and the discrete access margin

Renting can be quantitatively relevant and limited by \(h_R^{\max}\), but a
renter-to-owner switch is unnecessary for the preceding theorem. Avoiding it
also avoids changing the household's ownership taste contribution or
reoptimizing its fertility as a hidden part of the welfare comparison.

There is a useful exact obstruction to a simpler access proof. Suppose a
renter's old bundle satisfies the owner estate floor and size cap. At constant
prices with \(q<\phi\), every renter lifetime bundle can be implemented as
ownership with the same real allocation: if the renter has financial wealth
\(a'_R\ge0\), let \(a'_O=a'_R-Ph\). Then
\[
 qa'_O+\phi Ph=qa'_R+(\phi-q)Ph\ge0. \tag{28}
\]
Consequently a rental cap does not automatically create an otherwise
unfinanceable ownership bundle at the same size. A housing increase financed
solely by reducing current consumption at service cost \(p\), keeping future
resources fixed, remains privately finance-feasible when its initial point
has \(z\ge v\). A finite ownership-taste penalty cannot simply be omitted.

A near-indifferent renter can gain from a genuinely relaxed owner lifetime
plan, but if its fertility is held fixed at the renter choice, that plan must
be compared with the privately optimized owner value, which used a potentially
different fertility. The standard "near the tenure cutoff" argument alone
does not prove this comparison. A rigorous discrete-access theorem needs
additional restrictions or a separate fixed-fertility tenure equilibrium.
The owner theorem above avoids this outstanding step while retaining positive
renting, heterogeneity and the physical rental cap in the economy.

## Verification and remaining obligation

The fixed-fertility formulas follow directly from the two budgets and three
first-order conditions. The lead independently derived (9) and (11). The
expanded endogenous-fertility comparison was also checked symbolically, but
is not needed for the welfare theorem and is omitted to keep the core result
fixed-population. No model solve, numerical reference equilibrium, paper edit,
or change to the mortgage contract was performed.

The remaining substantive obligation for a fully primitive illustrative
equilibrium theorem is stationary market existence with the explicit cap
margins and a positive constrained owner group. The next useful independent
task is that bounded market construction, rather than changing mortgage
timing or seeking another current MRS-gap proof.
