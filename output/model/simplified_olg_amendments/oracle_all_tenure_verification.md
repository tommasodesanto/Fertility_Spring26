# Independent analytical verification of the overnight theory argument

This is supporting work, not an additional reader or adopted manuscript. The current short statement is in oracle_essential_theory_assessment.md. Original worker notation is retained below; e sometimes denotes an auxiliary spending index and must not replace estate notation in the paper.


---

## Household choices, income interval, and private-fertility counterexample

# All-tenure household characterization and a simpler housing-wedge theorem

Bounded original derivation for the September 10 theory discussion. This is a
research scratch deliverable, not adopted paper prose. The stipulated two-period
budgets, positive common retirement income, and all ownership deviations are
retained. No replacement-fertility condition or numerical reference point is
used in the main household theorem. The separate stationary witness below
explicitly includes demographic replacement.

## Main result

**There is no need to prove that affected young households rent.** At a given
rent, the incomes with a strictly positive housing distortion form one explicit
interval. Poorer members of that interval choose capped rentals; richer members
choose constrained ownership. The tenure cutoff inside the interval need not be
computed for the allocation theorem.

Write

\[
d=1+\tau-q,\qquad p=dP,\qquad
\Gamma=\frac{q(1-\phi)}d>0,\qquad
K=1+\alpha+\omega_B,\qquad \delta=\beta K,
\]
\[
Y=qy,\qquad A=\alpha+\vartheta,\qquad J=1+A,\qquad
L(p)=\chi+\kappa p,\qquad
A(p)=\alpha+\frac{\vartheta\kappa p}{L(p)}.
\]

Assume the exact retirement implementability condition

\[
\boxed{\omega_B\geq\alpha\Gamma.}
\]

Then an old household of resources \(z>0\) can implement its unrestricted
allocation by ownership:

\[
c_o=z/K,\quad h_o=\alpha z/(Kp),\quad qe=\omega_Bz/K,
\qquad V_o(z)=K\log z+C(p).
\]

Let

\[
x_r(p)=\frac{pr}{A(p)},\qquad
w_r(p)=Jx_r(p)+\max\{\delta x_r(p)-Y,0\},
\]
\[
w_O(p)=
\begin{cases}
\displaystyle \frac{Y[J+\Gamma A(p)]}{\delta-\Gamma A(p)},
&\delta>\Gamma A(p),\\[1ex]
+\infty,&\delta\leq\Gamma A(p).
\end{cases}
\]

Here \(w_r\) is the income at which the optimal uncapped rental first reaches
the rental ceiling. The threshold \(w_O\) is the income at which the fully
unrestricted lifetime allocation first satisfies ownership finance.

**Proposition.** At fixed \(p>0\), for every young income \(w>0\), every globally
optimal tenure plan satisfies

\[
M(w,p):=\frac{u_h}{u_c}
=\frac{\alpha(c-\chi n)}{h-\kappa n}\geq p.
\]

Its inequality is strict if and only if

\[
\boxed{w_r(p)<w<w_O(p).}
\]

At a tenure tie inside that interval, both optimal plans have a strict wedge.
The interval is nonempty precisely when

\[
\boxed{\delta pr<A(p)[Y+\Gamma pr].}
\]

If it is nonempty, there is a unique income cutoff \(w_T(p)\) strictly inside
\((w_r,w_O)\), with the upper endpoint interpreted as infinity when necessary.
The household rents below \(w_T\) and owns above it. Renting is capped on
\((w_r,w_T)\); ownership is financially constrained on \((w_T,w_O)\).
The cutoff \(w_T\) is a strictly monotone scalar value-comparison root, usually
not an elementary formula. The wedge result does not depend on that root.

If \(w_r\geq w_O\), the tenure union implements the optimal uncapped rental
allocation for every income, so no young household has a housing wedge. Tenure
can tie where that allocation fits both tenure constraints.

### An even simpler sufficient income condition

Let \(m=Y/\delta\). If

\[
\boxed{w<Jm,\qquad
\frac{w}{J}\left(\frac\alpha p+
\frac{\vartheta\kappa}{\chi+\kappa p}\right)>r,}
\]

then the household has a strict wedge under whichever tenure it chooses. The
first restriction says the young would like to borrow against retirement
resources at an uncapped rental. The second says that rental would exceed the
ceiling. There is no lower price cutoff and no rich-household finance condition.
In particular, cheap ownership does not remove the result: the constrained
owner is allowed to replace the capped renter as the affected household.

## 1. Unified household problem and exact finance condition

Use adult goods \(x=c-\chi n\), adult space \(s=h-\kappa n\), and retirement
expenditure \(t=qz\). A young household maximizes, up to constants,

\[
\log x+\alpha\log s+\vartheta\log n+\delta\log t
\]

subject to

\[
x+ps+L(p)n+t=w+Y.
\]

Renting requires

\[
t\geq Y,\qquad s+\kappa n\leq r.
\]

Ownership requires

\[
t\geq Y+\Gamma p(s+\kappa n).
\]

This is exactly the given \(a\geq-q\phi Ph\) restriction: neither \(Y\) nor
its retirement payoff can be counted as accumulated housing equity. The owner
constraint implies \(t>Y\) at every feasible positive home. With endogenous
continuous fertility, both tenure problems are feasible at every \(w>0\).
Each tenure problem is strictly concave over a convex set, and has a unique
real allocation. The tenure union is the union of these two feasible sets.

## 2. Renting: all four constraint combinations

### Uncapped rental solution

If the rental ceiling is removed but nonnegative saving is retained, then

\[
x_U=\min\left\{\frac wJ,\frac{w+Y}{J+\delta}\right\},
\quad n_U=\frac{\vartheta x_U}{L(p)},
\quad h_U=\frac{A(p)x_U}{p},
\quad t_U=\max\{Y,\delta x_U\}.
\]

Saving is zero at \(w\leq Jm\), positive at \(w>Jm\). The ceiling is slack
exactly at \(w<w_r(p)\), with the zero-shadow boundary at equality.

### Capped rental solution

For \(w>w_r(p)\), set \(h=r\). Define \(n_m\in(0,r/\kappa)\) as the unique
solution of

\[
\frac{\vartheta}{n_m}=\frac\chi m+
\frac{\alpha\kappa}{r-\kappa n_m}.
\]

The saving threshold conditional on the cap is

\[
w_s^{\rm cap}(p)=pr+m+\chi n_m.
\]

Below this threshold use \(g=1\), \(C=w-pr\), and \(t=Y\); above it use
\(g=1+\delta\), \(C=w+Y-pr\), and \(t=\delta x\). In both cases fertility is
the smaller admissible root of

\[
\chi\kappa(\vartheta+g+\alpha)n^2
-\{\kappa C(\vartheta+\alpha)+\chi r(\vartheta+g)\}n
+\vartheta Cr=0,
\]

and

\[
x=\frac{C-\chi n}{g},\qquad s=r-\kappa n.
\]

These formulas apply only when the cap is optimal. If \(x_r<m\), the sequence
as income rises is uncapped/no saving, capped/no saving, capped/positive saving.
If \(x_r>m\), the sequence is uncapped/no saving, uncapped/positive saving,
capped/positive saving. At \(x_r=m\), the two thresholds coincide. Thus no
constraint combination is silently excluded.

The rental MRS is \(p\) if the cap is slack and strictly exceeds \(p\) if it
binds with positive shadow value. The saving constraint by itself does not
produce a housing MRS wedge.

## 3. Ownership: unrestricted solution or one unique cubic root

The unrestricted lifetime solution is

\[
x_F=\frac{w+Y}{J+\delta},\quad
n_F=\frac{\vartheta x_F}{L(p)},\quad
ph_F=A(p)x_F,\quad t_F=\delta x_F.
\]

Its exact ownership feasibility test is

\[
[\delta-\Gamma A(p)]x_F\geq Y,
\]

which gives exactly \(w\geq w_O(p)\). No sufficient bound replacing \(A(p)\)
by \(A\) is necessary.

When this test fails, let \(u\in(0,1)\) be the borrowing multiplier divided by
the present-value resource multiplier. Set \(b=1+\Gamma u\in(1,1+\Gamma)\).
Then

\[
M=pb,\qquad
s=\frac{\alpha x}{pb},\qquad
n=\frac{\vartheta x}{\chi+\kappa pb},\qquad
t=\frac{\delta\Gamma x}{1+\Gamma-b}.
\]

Put \(v=\kappa p\), \(L_0=1+\Gamma\), and

\[
D(b)=\frac\alpha b+\frac{\vartheta v}{\chi+vb},
\quad
T(b)=1+\frac{L_0\alpha}{b}
+\frac{\vartheta(\chi+L_0v)}{\chi+vb}.
\]

The unique constrained-owner root solves

\[
\boxed{\frac{\delta\Gamma}{L_0-b}-\Gamma D(b)
-\frac YwT(b)=0,\qquad 1<b<L_0.}
\]

The left side is strictly increasing in \(b\), is negative at 1 precisely
when the unrestricted finance test fails, and tends to positive infinity at
\(L_0\). Therefore this equation identifies the owner globally, not merely a
stationary candidate. Recover adult consumption from either exact identity

\[
x=\frac w{T(b)}
=\frac{w+Y(L_0-b)/\Gamma}{J+\delta}.
\]

For an explicit polynomial, multiply the boxed equation by the positive
quantity \(wb(\chi+vb)(L_0-b)\). The result is the cubic

\[
\begin{aligned}
0={}&\delta\Gamma w b(\chi+vb)
-\Gamma w[\alpha\chi+v(\alpha+\vartheta)b](L_0-b)\\
&-Y\{vb^2+[\chi(1+\vartheta)+L_0v(\alpha+\vartheta)]b
+L_0\alpha\chi\}(L_0-b).
\end{aligned}
\]

Only its one root in \((1,L_0)\) is admissible. The resulting loan constraint
binds, and \(p<M<(1+\Gamma)p\). As income rises within this owner regime,
\(b\) falls strictly, while \(x,s,n,h,t\) rise strictly. Adult space is **not**
monotone over the entire tenure union: capped renters have \(s=r-\kappa n\),
which falls as their fertility rises.

## 4. Why the income interval is exact and why there is one tenure cutoff

First, every ownership plan is feasible in the hypothetical uncapped rental
problem, because ownership implies \(t\geq Y\). Its value is therefore no
greater than the uncapped rental value. If \(w\leq w_r\), that optimal rental
is actually available, so it is globally optimal. Its housing MRS is \(p\).

Second, when \(w\geq w_O\), the unrestricted lifetime optimum is implementable
by ownership. It is globally optimal. When its house exceeds \(r\), it strictly
dominates renting; when it fits within \(r\), the two tenure labels can tie.
Its housing MRS is \(p\).

Third, at \(w_r<w<w_O\), the rental optimum has a positive cap multiplier,
while the owner optimum has a positive loan multiplier. Hence both candidates
have \(M>p\), whichever has greater utility. This covers **every** saving,
fertility and ownership deviation. It is not an assumption about the desired
MRS ordering.

The nonempty-interval test follows because \(x_r\) is the adult-consumption
level at rental-cap onset, and the owner-finance threshold in these units is
\(Y/[\delta-\Gamma A(p)]>m\) whenever finite. Thus \(w_r<w_O\) is exactly
\(\delta pr<A(p)(Y+\Gamma pr)\).

For completeness, the within-interval tenure comparison is strictly
single-crossing in income. Let \(V_O,V_R\) denote optimized tenure values.
The envelope theorem gives

\[
\frac{d}{dw}(V_O-V_R)=\frac1{x_O}-\frac1{x_R}.
\]

For every constrained owner, \(T(b)>J\) and the second formula for \(x\)
above imply

\[
x_O<\min\{w/J,(w+Y)/(J+\delta)\}=x_U.
\]

For a capped renter, let \(b_R=M_R/p>1\). Its first-order conditions imply

\[
x_R=\begin{cases}
[w+(b_R-1)pr]/J,&t=Y,\\
[w+Y+(b_R-1)pr]/(J+\delta),&t>Y.
\end{cases}
\]

In either case \(x_R>x_U\). An unrestricted owner has \(x_O=x_U\), so
\(V_O-V_R\) is strictly increasing on all \(w>w_r\). It is negative at
\(w_r\) if \(w_r<w_O\), and positive at finite \(w_O\). If \(w_O=\infty\),
ownership's value grows like \((J+\delta)\log w\), whereas capped rental value
grows like \((1+\delta)\log w\); therefore the difference eventually becomes
positive. This proves the unique \(w_T\).

At \(w_T\), the owner house must strictly exceed \(r\). Otherwise renting
would implement the owner plan and allow a beneficial reduction in forced
saving. Hence housing jumps upward at a tenure switch, although fertility's
jump need not have a universal sign.

## 5. Fixed-fertility version, including owner infeasibility

Fix \(n>0\), initially with \(\kappa n<r\) and \(w>L(p)n\). Let
\(B_n=1+\alpha\), \(R=r-\kappa n\), and \(C=w-L(p)n\). The uncapped rental is

\[
x_U=\min\{C/B_n,(C+Y)/(B_n+\delta)\},\quad
s_U=\alpha x_U/p,\quad t_U=\max\{Y,\delta x_U\}.
\]

Define

\[
x_r=pR/\alpha,\qquad
w_r=L(p)n+B_nx_r+\max\{\delta x_r-Y,0\}.
\]

Above \(w_r\), the renter takes \(h=r\). With
\(C_r=w-pr-\chi n\), its adult consumption is
\(x_R=\min\{C_r,(C_r+Y)/(1+\delta)\}\), and saving is positive exactly when
\(C_r>m\).

The unrestricted owner's adult consumption is
\((w+Y-L(p)n)/(B_n+\delta)\). Define

\[
x_O^{\rm threshold}=
\frac{Y+\Gamma p\kappa n}{\delta-\Gamma\alpha},\quad
w_O=L(p)n+(B_n+\delta)x_O^{\rm threshold}-Y
\]

when \(\delta>\Gamma\alpha\), and set \(w_O=\infty\) otherwise. Again the
strict-wedge interval is exactly \((w_r,w_O)\); an infeasible owner candidate
simply has value negative infinity. The interval is nonempty exactly when

\[
\delta p(r-\kappa n)<\alpha[Y+\Gamma pr].
\]

Ownership is feasible only if

\[
v_0=w-[\chi+(1+\Gamma)p\kappa]n>0.
\]

Below \(w_O\), but when feasible, put \(Y_0=Y+\Gamma p\kappa n\). Adult
housing expenditure \(e=ps\) is the unique positive root of the quadratic

\[
\Gamma(1+\Gamma)(B_n+\delta)e^2
+\{(1+\Gamma)B_nY_0-\Gamma(\alpha+\delta)v_0\}e
-\alpha v_0Y_0=0.
\]

Then \(x=v_0-(1+\Gamma)e\), \(t=Y_0+\Gamma e\), and \(h=e/p+\kappa n\).
The admissible root lies in \((0,v_0/(1+\Gamma))\).

If fixed fertility has \(\kappa n\geq r\), renting is infeasible; the owner
formulas and finance threshold still apply. If neither tenure's current
resource condition holds, that imposed fertility is infeasible. Endogenous
continuous fertility avoids this low-income feasibility issue by adjusting
\(n\); it does not borrow against the pension.

## 6. Optional dated-equilibrium checks, not an adopted closure

This subsection takes old households' dated purchasing powers \(z_o\geq y\)
as fixed, with positive, exogenous cohort masses \(N_y,N_o\) and finite resource
means. It is an optional current-market-clearing check under the stipulated
continuation prices, not an adopted OLG equilibrium closure or an existence
proof for an entire transition. It does not close the population through
fertility. Write

\[
W_y=N_y\mathbb E[w],\qquad Z_o=N_o\mathbb E[z_o],\qquad H>0.
\]

The old demand is \(\alpha Z_o/(Kp)\). Every young global optimum satisfies

\[
h(w,p)\leq h_U(w,p)\leq\frac{A(p)w}{Jp}
\leq\frac{Aw}{Jp}.
\]

The first bound follows from the explicit owner formulas and the rental cap:
constrained owners choose less housing than the uncapped rental, while capped
renters take less than its desired home. Therefore every market-clearing rent
satisfies the elementary bounds

\[
\boxed{p_{\min}:=\frac{\alpha Z_o}{KH}\leq p\leq
\frac{AW_y/J+\alpha Z_o/K}{H}=:p_{\max}.}
\]

An equilibrium exists for every \(H>0\), allowing a continuum of households
of an exactly indifferent type to split between its two optimal tenure plans.
This uses only optimal choices at an exact tie. The aggregate housing-demand
correspondence is nonempty, compact, upper hemicontinuous, and interval-valued.
Old demand alone exceeds \(H\) for sufficiently small \(p\), and the displayed
upper bound is below \(H\) for sufficiently large \(p\). A connected price
interval cannot pass between these two signs without its demand correspondence
containing \(H\). With atomless young income heterogeneity, the unique income
tenure threshold instead makes aggregate demand continuous, so the ordinary
intermediate-value theorem suffices.

No global price uniqueness is asserted. If one insists that all households in
an atomic type choose the same pure tenure even at indifference, the demand
jump can skip the stock and existence is not guaranteed.

Now suppose a positive-mass young group has incomes in

\[
0<w_a\leq w\leq w_b<\frac{Jqy}{\beta K}.
\]

The following simple stock-versus-rental-ceiling condition suffices:

\[
\boxed{H>\frac r{w_a}
\left(\frac A\alpha W_y+\frac JK Z_o\right).}
\]

At every equilibrium it implies \(p\leq p_{\max}<\alpha w_a/(Jr)\).
Consequently every member of the group wants an uncapped no-saving rental
larger than \(r\), and cannot implement the unrestricted lifetime allocation
by ownership. Every member therefore has \(M>p\), independently of tenure.

A weaker but still explicit sufficient test replaces the boxed stock bound by
\(h_U(w_a,p_{\max})>r\). A fully sharp support test based on the price bounds
is \(w_a>w_r(p_{\max})\) and \(w_b<w_O(p_{\min})\), since both thresholds
increase with \(p\). These are optional refinements; the stronger boxed bound
keeps the economic content transparent.

The dated result therefore needs no rich-owner condition, no \(\delta>\Gamma A\)
condition, no endogenous cohort mass, and no replacement-price bracket. It does
need a positive group with relatively low current resources, a rental ceiling
that matters relative to the common stock, and \(\Gamma>0\).

### If inherited old titles must be revalued at the clearing price

Fixed \(z_o\) is a substantive dated initial-state convention. If instead
\(z_o(p)=\ell_o+(p/d)h_o^-\), retaining that actual title valuation gives

\[
H_{\rm eff}=H-\frac{\alpha H^-}{Kd},\qquad
D_o(p)=\frac{\alpha L_o}{Kp}+\frac{\alpha H^-}{Kd},
\]

where \(L_o=\int\ell_o\) and \(H^-=\int h_o^-\). Provided
\(\ell_o>0\) and \(H_{\rm eff}>0\), the same existence proof and all price
bounds hold with \(H\) replaced by \(H_{\rm eff}\) and \(Z_o\) by \(L_o\).
If \(H^-\leq H\), \(Kd>\alpha\) is a simple sufficient capitalization
condition. Negative nonhousing positions require a separate solvency/domain
check. The fixed-\(z\) proof must not silently be represented as proving this
different price-dependent inherited-wealth closure.

## 7. An actual stationary-equilibrium nonvacuity family for the new theorem

The earlier Pro construction can supply an **ancillary existence witness** for
the general all-tenure theorem, including the new resource and fertility
conditions. The general theorem need not repeat that witness's forced tenure
pattern or demographic closure. The following establishes a fully explicit
primitive family, not an unspecified neighborhood of a numerical equilibrium.

Fix positive \(\alpha,\vartheta,\kappa,r,y\), the stipulated finance parameters,
and preferences such that

\[
\omega_B>\alpha\Gamma,\qquad \delta>\Gamma A.
\]

Write \(\rho=\beta/q\), \(m=Y/\delta\). Choose

\[
0<w_L<\frac{Jm\min\{1,\rho\}}
{1+\Gamma(1-\rho)_+},
\qquad
w_H>\max\left\{
\frac{Y(J+\Gamma A)}{\delta-\Gamma A},\frac{Jy}{K}
\right\}.
\]

Let

\[
\ell=w_L/J,\quad
k=\frac{1+\Gamma}{1+\Gamma\ell/m},\quad b=k\ell,
\quad x=\frac{w_H+Y}{J+\delta}.
\]

The low-income restriction implies both
\(w_L<Jy/K\) and \(b<y/K\). Thus the numbers

\[
f_C=\frac{w_H-Jy/K}{w_H-w_L},\qquad
f_B=\begin{cases}
\displaystyle\frac{(1-\rho)x}{y/K-b+(1-\rho)x},&\rho<1,\\
0,&\rho\geq1
\end{cases}
\]

are strictly below 1. Choose \(\max\{f_C,f_B\}<f<1\). This yields

\[
\boxed{\frac yK>\frac{fw_L+(1-f)w_H}{J}}
\]

and also preserves the stronger old Pro resource bound
\(f y/K+(1-f)\beta x/q>fb+(1-f)x\).

Finally, choose the positive child-goods coefficient within the explicit range

\[
\boxed{0<\chi<\frac{\kappa w_L}{r(1+\alpha+k)}.}
\]

Let \(n_L(p)\) be Pro's capped, no-saving rental candidate, and
\(M_L(p)=\alpha x_L(p)/s_L(p)\). Its two elementary price cutoffs satisfy
\(M_L(p_-)=kp_-\), \(M_L(p_+)=p_+\). At
\(p_0=\alpha\chi/\kappa\), the candidate is feasible and

\[
x_L(p_0)>w_L-(1+\alpha)\chi r/\kappa,\qquad s_L(p_0)<r.
\]

Consequently the explicit \(\chi\) bound gives

\[
M_L(p_0)>
\frac\alpha r[w_L-(1+\alpha)\chi r/\kappa]
>\frac{k\alpha\chi}{\kappa}=kp_0.
\]

Because \(M_L(p)-kp\) is strictly decreasing, this proves
\(p_0<p_-\). Hence **every** price in the previously verified regime obeys

\[
\boxed{\kappa p>\alpha\chi.}
\]

Define the prior aggregate fertility schedule

\[
F(p)=fn_L(p)+(1-f)\frac{\vartheta x}{\chi+\kappa p}.
\]

Choose the primitive household-entry conversion in the nonempty open interval
\(F(p_+)<1/\nu<F(p_-)\), as in the already verified stationary witness.
For any stock \(H>0\), its unique price in that regime and cohort mass
\(N=H/D(p)\) give the actual stationary equilibrium, including old wealth
generated by the preceding young choices. It has a positive mass of capped
poor renters with a strict wedge, and satisfies both boxed new conditions.

This witness is stronger than is necessary. It proves compatibility with an
actual equilibrium while allowing the new main statement to quantify over
all existing competitive equilibria and all young tenure regimes. It does not
relabel the optional fixed-old-wealth temporary equilibrium as a stationary
OLG equilibrium.

## 8. What is and is not established

The household theorem is exact. The stock inequality is a conservative,
explicit sufficient primitive condition for a positive mass of distorted young
at every dated equilibrium under its stated old-wealth convention. It replaces
the earlier narrow capped-poor-renter/rich-unconstrained-owner construction with
the actual union of all tenure regimes.

The lead owns the dated planner and fertility results. A strict MRS gap alone
does not guarantee that every affected household receives more goods or chooses
more fertility under an equally weighted planner. Additional resource and
child-cost conditions belong in that separate theorem. There is no claim about
policy implementation, a transition, or replacement population dynamics here.

Verification performed: direct Kuhn–Tucker derivation, positivity and endpoint
checks for the unique owner root, exact rental quadratics, envelope
single-crossing proof, and comparison with the stipulated retirement and owner
budgets. Symbolic substitution verified that the endogenous-owner cubic, capped
renter quadratic, and fixed-fertility owner quadratic equal their corresponding
first-order equations after multiplying only by the stated positive
denominators. No model run, numerical sweep, manuscript edit, or compiled
artifact.

## 9. New bounded scope: private fertility after the fixed-fertility allocation

This section addresses the separate September 10 follow-up. The planner first
holds each original fertility \(n_i\) fixed and assigns total bundles
\(c_i^F=X+\chi n_i\), \(h_i^F=S+\kappa n_i\). Parents then separately
rechoose fertility, holding those **total** bundles and their continuation
resources fixed. This is not the joint planner choosing fertility.

**Result: mean private fertility need not rise. It can strictly fall even under
the stronger primitive condition \(y/K\geq\bar w/J\), with a positive mass of
strict housing distortions, unrestricted old choices, and
\(\kappa p>\alpha\chi\).** The counterexample below has explicit primitive
bounds, includes every tenure deviation, and is an actual stationary
equilibrium. No simulation or numerical root supports the result.

### 9.1 Why the cutoff argument cannot sign the mean

Write \(\Delta_i=n_i^{\rm new}-n_i\). The new private first-order condition is

\[
n_i+\Delta_i=f(\Delta_i),\qquad
f(\Delta)=\vartheta\left[
\frac{\chi}{X-\chi\Delta}
+\frac{\alpha\kappa}{S-\kappa\Delta}
\right]^{-1}.
\]

The function \(f\) is decreasing and concave. If
\(X/\chi\ne S/\kappa\), it is strictly concave. Hence the inverse relation
\(n=f(\Delta)-\Delta\) makes \(\Delta(n)\), and therefore
\(n^{\rm new}(n)=n+\Delta(n)\), strictly concave. The individual sign is
indeed determined by

\[
n_i<n_c:=\vartheta
\left[\frac\chi X+\frac{\alpha\kappa}S\right]^{-1}.
\]

But \(n_c>\bar n\) does not sign the mean response of a concave function.
For example, at the undistorted resource boundary \(B=\bar e\), heterogeneous
young have \(x_i=v_i=e_i\), while the fixed-fertility planner has
\(X=\bar e\), \(S=\alpha\bar e/p\). If
\(\kappa p>\alpha\chi\), the private-response function is strictly concave
and fixes the competitive mean at that mean. Strict Jensen then gives a mean
private decline. The following finite primitive range preserves the decline
while introducing a strict housing distortion.

### 9.2 Explicit stationary counterfamily, including the strong income floor

Choose any

\[
0<\epsilon\leq\frac1{10000},\qquad
r=\frac52-\epsilon.
\]

Define the explicit radical

\[
d(\epsilon)=
\frac{\sqrt{49/4+10\epsilon+7\epsilon^2}-7/2-\epsilon}{3},
\qquad a_c=\frac13-\frac d4.
\]

Here \(a_c\) is a common coefficient for children's goods and space needs,
not a financial asset. Set the primitives

\[
\alpha=\vartheta=1,\quad \chi=\kappa=a_c,\quad
q=\frac12,\quad\beta=\frac18,\quad
\tau=\frac12,\quad\phi=\frac12,\quad\omega_B=2,
\]
\[
y=8,\qquad w_L=3,\qquad w_H=9,\qquad f=\frac12,
\qquad \nu=\frac12.
\]

Then \(K=4\), \(J=3\), \(\delta=\beta K=1/2\), \(Y=qy=4\),
\(d_{\rm price}=1+\tau-q=1\), and \(\Gamma=1/4\). The proposed stationary
prices are \(p=P=2\), and both adult cohorts have unit mass. Set the common
housing stock to

\[
H=\frac83-\frac\epsilon2.
\]

To keep the algebra readable, write \(m=a_c n\) for the goods and space
required by a parent's children. Utility differs by a constant from
\(\log(c-m)+\log(h-m)+\log m\). Actual fertility is always
\(n=m/a_c\). This common cost scaling ensures exact replacement at the fixed,
conventional \(\nu=1/2\); it does not change payment timing or any tenure
comparison.

For fixed total goods and housing, the optimal child-needs quantity is

\[
N(c,h)=\frac{c+h-\sqrt{c^2-ch+h^2}}3.
\tag{C1}
\]

The candidate young choices are

\[
\begin{array}{c|ccccc}
&c&h&m&x=c-m&s=h-m\\\hline
L&4/3&5/6&1/3&1&1/2\\
H&4+2\epsilon&5/2-\epsilon&1-d&3+2\epsilon+d&3/2-\epsilon+d.
\end{array}
\tag{C2}
\]

Both young types rent and save zero. Their retirement resources are exactly
\(z_L=z_H=y=8\). The rich fertility in (C2) is exactly (C1) evaluated at its
displayed total bundle. Rationalizing the radical gives

\[
d=\frac{\epsilon(1+2\epsilon)}
{\sqrt{49/4+10\epsilon+7\epsilon^2}+7/2+\epsilon},
\qquad 0<d<\epsilon.
\tag{C3}
\]

#### Global household verification

The poor allocation is the unconstrained rental optimum with zero saving:
\(x=w_L/J=1\), \(m=x/(1+p)=1/3\), and \(h=5/6<r\).
Its saving threshold is \(Y/\delta=8>x\). Every ownership plan belongs to the
uncapped rental feasible set, while this unique optimum violates owner finance
because \(t=Y\) and \(h>0\). Thus the poor strictly prefer renting to every
ownership deviation.

For the rich, the uncapped no-saving home is \(5/2>r\), so the rental cap
binds strictly. The capped choice has \(3<x_H<4<8=Y/\delta\), confirming zero
saving. Define the exact supporting multipliers

\[
\xi=\frac1{s_H}-\frac2{x_H},\qquad
\zeta=\frac1{x_H}-\frac18.
\]

Since \(0<d<\epsilon\),

\[
0<\xi=\frac{4\epsilon-d}{x_Hs_H}<\epsilon,
\qquad \zeta>\frac18,
\qquad \Gamma p\zeta>\frac1{16}>\epsilon.
\tag{C4}
\]

The inequality for \(\xi\) uses \(x_H>3\) and
\(s_H>3/2-\epsilon>4/3\). Concavity of complete lifetime utility, including
fertility and retirement wealth, gives for **any** alternative plan

\[
U-U_H\leq\xi(h-r)-\zeta(t-Y).
\]

An owner must have \(t-Y\geq\Gamma ph\), so its value difference is at most

\[
(\xi-\Gamma p\zeta)h-\xi r<0.
\]

This excludes every ownership, saving and fertility deviation, not just one
candidate purchase. The rich rental has a strictly positive housing wedge.

Old choices are unrestricted because
\(\omega_B=2>\alpha\Gamma=1/4\). Every old household chooses
\(c_o=2\), \(h_o=1\), and \(qe=4\), hence \(e=8\). These are feasible
under the original serviced-interest budget and collateral rule. Old tenure
may tie; its real allocation does not.

#### Actual stationarity and the theorem's assumptions

The young current-expenditure indices are \(e_L=w_L/J=1\), \(e_H=w_H/J=3\).
Consequently

\[
\boxed{B=2=\bar e=\frac yK=\frac{\bar w}{J},\qquad
\kappa p=2a_c>\alpha\chi=a_c.}
\tag{C5}
\]

A positive mass, one half of the young, has the strict housing wedge in (C4).
The stock clears exactly:

\[
\frac12\left(\frac56+\frac52-\epsilon\right)+1
=\frac83-\frac\epsilon2=H.
\]

Moreover \(\bar m=2/3-d/2=2a_c\), so actual mean fertility is exactly
\(\bar n=2\). With \(\nu=1/2\), unit young and old cohort masses satisfy
replacement. Each previous young type saved zero, so the old wealth \(z=8\)
is precisely generated by the prior cohort's choices, including actual title
positions; it is not imposed as a separate fixed-wealth market closure.
Entrant endowments, external finance, common divisible stock, and the
public-services fiscal closure are the same maintained objects as in the
checked Pro model. Thus this is an actual stationary equilibrium for every
primitive point in the stated finite \(\epsilon\) interval.

### 9.3 The mean private response is strictly negative throughout that range

The full fixed-fertility planner assigns common adult goods and adult space

\[
X=2+\frac\epsilon2+\frac d4,
\qquad S=1-\frac\epsilon4+\frac d4.
\]

It then gives parent \(i\) the total bundle \((X+m_i,S+m_i)\).
After private reoptimization, its child-needs quantity is
\(m_i^{\rm new}=N(X+m_i,S+m_i)\), and actual fertility is
\(n_i^{\rm new}=m_i^{\rm new}/a_c\).

At the algebraic limiting endpoint \(\epsilon=0\),

\[
m_L^{\rm new}(0)=\frac{11-\sqrt{37}}9,
\qquad
m_H^{\rm new}(0)=\frac{5-\sqrt7}3.
\]

Thus the mean **child-needs** loss at that endpoint is exactly

\[
G_0=\frac23-
\frac12\left(\frac{11-\sqrt{37}}9+
\frac{5-\sqrt7}3\right)
=\frac{\sqrt{37}+3\sqrt7-14}{18}
>\frac1{1000}.
\tag{C6}
\]

For an elementary rational verification of the last inequality, use
\(\sqrt{37}>60827/10000\) and \(\sqrt7>26457/10000\), each verified by
squaring. These lower bounds actually imply \(G_0>11/10000\).

The positive-\(\epsilon\) conclusion does **not** rely on an unspecified
continuity neighborhood. From the first-order condition defining (C1),

\[
0<N_c<1,\quad0<N_h<1,\quad N_c+N_h<1.
\]

So changing total goods and housing changes \(N\) by at most the largest
absolute coordinate change. Using \(0<d<\epsilon\), the low parent's assigned
bundle moves by at most \(3\epsilon/4\) in that norm; the high parent's by at
most \(\epsilon\). Their mean optimized child-needs quantity therefore moves
by at most \(7\epsilon/8\). Original mean child needs moves by
\(d/2<\epsilon/2\). It follows that

\[
\bar m^{\rm new}(\epsilon)-\bar m(\epsilon)
<-\frac1{1000}+\frac{11}{8}\epsilon
\leq-\frac{69}{80000}<0.
\tag{C7}
\]

Dividing by the positive \(a_c\) proves

\[
\boxed{\bar n^{\rm new}<\bar n=2,}
\]

throughout the explicit primitive range. Meanwhile the conditions of the
separate sharp joint-planner theorem hold strictly where required, so that
joint planner chooses \(n^F>2\). These are different exercises with opposite
mean-fertility signs in the same equilibrium. The private low-income parent
raises fertility; the high-income parent's reduction is larger in the mean.

### 9.4 A genuine all-tenure low-income individual result survives

Global fertility sorting through a tenure switch is not needed for a useful
individual statement. Under \(\kappa p\geq\alpha\chi\), private optimality
and \(z_i\geq y\) imply

\[
n_i\leq\frac{\vartheta e_i}{\chi+\kappa p}
\leq\frac{\vartheta w_i}{J(\chi+\kappa p)}.
\tag{C8}
\]

Let \(\mu=N_o/N_y>0\), and write adult-space value as \(V=pS/\alpha\).
The fixed-fertility planner has
\(X=(\bar x+\mu B)/(1+\mu)\),
\(V=(\bar v+\mu B)/(1+\mu)\). Because both young averages are positive,
\(X,V>\mu B/(1+\mu)\). Therefore the individual fertility cutoff obeys

\[
n_c=\vartheta\left(\frac\chi X+\frac{\kappa p}V\right)^{-1}
>\frac{\vartheta\mu B}{(1+\mu)(\chi+\kappa p)}
\geq\frac{\vartheta\mu y}{(1+\mu)K(\chi+\kappa p)}.
\]

Combining this with (C8) proves the primitive sufficient condition

\[
\boxed{w_i\leq\frac{J\mu y}{(1+\mu)K}
\quad\Longrightarrow\quad n_i^{\rm new}>n_i.}
\tag{C9}
\]

For equal cohorts this is \(w_i\leq Jy/(2K)\). It allows any initial tenure,
uses the full fixed-fertility planner bundle, and needs no global income–fertility
ranking. It is conservative and does not identify every household whose
fertility rises. The more general individual test remains \(n_i<n_c\), and the
sharp planner theorem guarantees it for any parent with \(n_i\leq\bar n\).

Within each tenure branch, fertility rises with income by the household
formulas above. A claim that this ordering also holds at every tenure switch
is **not established here** and is not used in (C9). No mean private-rechoice
claim should be substituted for the proved joint-planner result.

Verification for the new scope: direct quadratic root, exact collateral
supporting-hyperplane inequality covering every ownership deviation, rational
radical bounds, derivative bounds valid on the full stated interval, actual
stationary wealth/housing/replacement accounting, and a separate proof of the
primitive individual cutoff. No numerical optimization, model simulation,
reference calibration, manuscript change, or policy/transition extension.


---

## Full allocation, joint fertility, and independent household review

# Sharper fertility and full-planner result in the unchanged two-period model

Bounded analytical assignment, September 10, 2026. This note owns no source,
deck, or manuscript changes. It uses the serviced-interest financing rule,
common stock, rental ceiling, positive common retirement income, and full dated
planner in the checked positive-retirement-income response. No simulation,
transition calculation, or lifecycle change was used.

## Main finding

The condition that the old have at least as much **adult consumption** as the
young is stronger than necessary. A useful replacement is that old adult
consumption cover a young household's **current expenditure divided by its
sum of utility coefficients**, on average. If children use housing at least as
intensively as adults, this weaker comparison delivers both more young housing
and higher mean fertility, even when the planner takes goods from the young.

In the original capped-poor-renter/unrestricted-rich-owner regime, replace

\[
B\geq f b+(1-f)x
\]

by

\[
\boxed{B\geq f\ell+(1-f)x,\qquad \kappa p\geq\alpha\chi.}
\tag{1}
\]

Here \(\ell=w_L/J<b\), with \(J=1+\alpha+\vartheta\), and all other objects
have their existing definitions. The first inequality is now precisely the
simple primitive bound that suffices for the housing direction. The second
says that children's housing-to-goods spending ratio, \(\kappa p/\chi\), is
at least the unrestricted adult ratio \(\alpha\). It is not an assumption
that both consumption and housing rise.

More broadly, a theorem covering **any young tenure mix** follows from

\[
\boxed{\frac yK\geq\frac{\bar w}{J},\qquad
\kappa p\geq\alpha\chi,\qquad K=1+\alpha+\omega_B.}
\tag{2}
\]

It also requires unrestricted old housing–estate choices, as ensured by the
existing \(\omega_B\geq\alpha\Gamma\), and a positive mass of young with
a strict housing wedge, which the separate household-threshold work identifies.
No ordering of \(\beta\) and \(q\) is required for (2). The income bound in
(2) is robust and simple, but is **not uniformly weaker** than the old
regime-specific condition: it exchanges the regime-specific resource calculation
for a sufficient restriction on average income. Section 6 gives a less
conservative alternative that still allows every tenure mix.

## 1. Objects established by private optimality

For a young household define adult goods and the rental value of adult space:

\[
x_i=c_i-\chi n_i,
\qquad v_i=\frac{p(h_i-\kappa n_i)}{\alpha},
\qquad t_i=\frac{x_i}{v_i}.
\]

All are positive. The ratio \(t_i\) is its housing MRS divided by rent.
An unconstrained housing choice has \(t_i=1\). A rental ceiling or a binding
owner-equity constraint gives \(t_i\geq1\); strictness is the economically
relevant constraint. This uses the branch first-order conditions and does not
select the household's tenure in advance. In particular,

\[
\frac{\alpha}{h_i-\kappa n_i}
=\frac p{x_i}+\xi_i
\]

for a renter, where \(\xi_i\geq0\) is its rental-cap multiplier. For an owner,
the same derivative includes \(\Gamma p\mu_i\), where \(\mu_i\geq0\) is
the multiplier on \(q(z_i-y)\geq\Gamma p h_i\). Hence both menus imply
\(t_i\geq1\).

Let

\[
C=\chi+\kappa p,\qquad \lambda=\frac{\chi}{\kappa p},\qquad
g(x,v)=\left(\frac\chi x+\frac{\kappa p}{v}\right)^{-1}.
\]

The individual fertility first-order condition is exactly

\[
n_i=\vartheta g(x_i,v_i)
=\frac{\vartheta v_i t_i}{\kappa p(\lambda+t_i)}.
\tag{3}
\]

Define a current-expenditure index

\[
e_i=\frac{c_i+p h_i}{J}
=\frac{x_i+\alpha v_i+C n_i}{J}
=\frac{w_i-q(z_i-y)}{J}.
\tag{4}
\]

The auxiliary unrestricted bundle with \(x=v=e_i\) and
\(n=\vartheta e_i/C\) has the same total current expenditure as household
\(i\). It need not be privately feasible. It is used only for the planner
comparison.

Substituting (3) gives two useful identities:

\[
e_i-v_i
=\frac{v_i(t_i-1)}J
\left[1+\frac{\vartheta\lambda}{\lambda+t_i}\right],
\tag{5}
\]

\[
x_i-e_i
=\frac{v_i(t_i-1)}J
\left[\alpha+\frac{\vartheta t_i}{\lambda+t_i}\right].
\tag{6}
\]

Thus \(v_i\leq e_i\leq x_i\), with both inequalities strict at a strict
housing wedge. Moreover,

\[
\boxed{
\frac{e_i}{C}-g(x_i,v_i)
=\frac{v_i(t_i-1)(t_i-\alpha\lambda)}
{J C(\lambda+t_i)}.}
\tag{7}
\]

Under \(\alpha\lambda\leq1\), the private housing distortion reduces
fertility relative to that unrestricted bundle with the **same current
expenditure**. This is the structural step that replaces the earlier assertion
that the planner must give the young more goods. It follows from household
optimality and a relative child-cost restriction, rather than from an assumed
ordering of the planner's consumption and housing.

Another identity, requiring no cost restriction, is

\[
Cg(x_i,v_i)-v_i
=\frac{\lambda v_i(t_i-1)}{\lambda+t_i}\geq0.
\tag{8}
\]

The key comparison is therefore

\[
v_i\leq Cg(x_i,v_i)\leq e_i\leq x_i,
\tag{9}
\]

where the middle inequality requires \(\kappa p\geq\alpha\chi\).

## 2. Fixed-fertility full planner: housing can rise while goods fall

Let overbars average across young households, let \(B\) be average old
consumption, and write \(\mu=N_o/N_y>0\) for the ratio of old to young
households. The maintained stationary comparison has \(\mu=1\); allowing a
different predetermined ratio costs nothing in this dated argument.

Old households have housing \(\alpha c_i^o/p\), so average old adult-space
value is also \(B\). At the competitive fertility choices, the full planner
equalizes adult goods and adult space across all currently living households:

\[
X=\frac{\bar x+\mu B}{1+\mu},\qquad
V=\frac{\bar v+\mu B}{1+\mu},\qquad S=\frac\alpha p V.
\tag{10}
\]

Young household \(i\) receives total goods \(X+\chi n_i\) and total housing
\(S+\kappa n_i\). Old households receive \(X,S\). Current resource totals,
individual young continuation resources, individual old estates, and
continuation prices are preserved, as in the existing dated benchmark.

If

\[
B\geq\bar e,
\tag{11}
\]

and a positive young mass has \(t_i>1\), (5) gives \(B>\bar v\). Hence

\[
\bar h_y^F-\bar h_y^E
=\frac\alpha p\frac\mu{1+\mu}(B-\bar v)>0.
\tag{12}
\]

At the boundary \(B=\bar e\), (6) instead gives

\[
X-\bar x
=-\frac\mu{1+\mu}(\bar x-\bar e)<0.
\tag{13}
\]

Thus the same comparison explicitly permits a strict **decline** in average
young consumption. Because fertility is fixed in this first comparison, the
change in total young consumption equals the change in adult consumption.

## 3. Joint parental fertility: a short complete proof

The function \(g\) is increasing, homogeneous of degree one, and concave.
For completeness,

\[
d^2g
=-\frac{2\chi\kappa p\,(v\,dx-x\,dv)^2}
{(\chi v+\kappa p x)^3}\leq0.
\tag{14}
\]

By (7), the cost condition and (11) imply

\[
\frac BC\geq\frac{\bar e}{C}
>\overline{g(x_i,v_i)}=\frac{\bar n}{\vartheta}.
\tag{15}
\]

The strict inequality follows from a positive mass of strict housing wedges.
Two applications of concavity now give

\[
\begin{aligned}
g(X,V)
&\geq\frac{g(\bar x,\bar v)+\mu g(B,B)}{1+\mu}\\
&\geq\frac{\overline{g(x_i,v_i)}+\mu B/C}{1+\mu}
>\frac{\bar n}{\vartheta}.
\end{aligned}
\tag{16}
\]

The full planner choosing fertility as well as consumption and housing assigns
the same fertility to every parent. Normalize resource totals by young mass,
calling them \(\mathcal C,\mathcal H\). Its reduced objective is

\[
(1+\mu)\log X(n)+(1+\mu)\alpha\log S(n)
+\vartheta\log n+\text{constants},
\]

\[
X(n)=\frac{\mathcal C-\chi n}{1+\mu},\qquad
S(n)=\frac{\mathcal H-\kappa n}{1+\mu}.
\]

Its derivative is

\[
\frac{\vartheta}{n}-\frac\chi{X(n)}
-\frac{\alpha\kappa}{S(n)}.
\tag{17}
\]

This derivative strictly decreases from positive infinity to negative infinity
over the feasible interval. At competitive mean fertility, its adult bundles
are exactly (10), and (16) makes (17) strictly positive. Therefore

\[
\boxed{n^F>\bar n^E.}
\tag{18}
\]

No utility weight on unborn people is used. The objective values existing
parents' own fertility preferences.

The joint fertility change strengthens the housing result. Relative to the
fixed-fertility comparison, mean young housing rises by a further

\[
\frac{\mu\kappa}{1+\mu}(n^F-\bar n^E)>0.
\]

These results concern a dated allocation. They do not establish a stationary
population comparison, policy implementation, or a transition.

### A sharper intermediate condition

The sufficient test

\[
\boxed{B\geq C\bar n^E/\vartheta}
\tag{19}
\]

alone gives both (12) and (18), **without** imposing
\(\kappa p\geq\alpha\chi\). Equation (8) gives the housing direction.
For fertility use (16), now with weak final inequality. Strict concavity of
\(g\) between \((\bar x,\bar v)\) and \((B,B)\) makes it strict, because
\(\bar x>\bar v\) at a positive mass of wedges. This is an informative
structural reduction, but it depends on equilibrium mean fertility and is not
the recommended primitive paper condition. Equations (2) or (1), with the
relative child-cost restriction, are transparent sufficient bounds for (19).

## 4. Why the simple income condition works for every tenure mix

Under the maintained financial contracts every feasible young household has
\(z_i\geq y\). Therefore (4) gives \(e_i\leq w_i/J\). Every old household
with total resources at least \(y\), and the unrestricted old allocation,
has \(c_i^o\geq y/K\). Consequently

\[
\frac yK\geq\frac{\bar w}{J}
\quad\Longrightarrow\quad B\geq\bar e.
\]

This argument neither labels the poor as renters nor requires the rich to
save or own without a binding financial constraint. It applies to any actual
mixture, once the household theorem establishes a positive mass of
\(t_i>1\). It also allows inherited old wealth to differ across households,
provided their resources remain at least the common retirement payment.

The condition has economic content: it requires the common retirement payment
to deliver enough adult consumption, after old housing and estate spending,
relative to young current endowments. Written in income units it is
\(y/\bar w\geq K/J\). Positive retirement income alone does not imply it.

For a dated economy with predetermined \(N_o>0\) and total housing stock
\(\bar H\), old housing demand and clearing imply

\[
p>\frac{\alpha N_o y}{K\bar H}.
\]

Thus \(\kappa N_o y/(K\bar H)\geq\chi\) is an optional price-free sufficient
bound for the child-cost condition. In the stationary construction with an
endogenous cohort mass, \(N_o\) is endogenous: do not call this a fully
primitive stationary condition. One may instead use the explicit lower rent
cutoff supplied by the household construction.

## 5. The sharper result in the original verified mixed regime

In that regime the poor save zero and have \(e_L=w_L/J=\ell\). The rich
choose their unrestricted bundle with \(e_H=x\). Thus (11) becomes exactly
the first inequality of (1). The private identities imply

\[
v_L<\ell<x_L<x,
\qquad
\bar v<f\ell+(1-f)x<\bar x.
\]

This both explains why \(b\) is unnecessary and proves that goods losses are
admissible. The result is strict even at equality in the new resource test.

The poorer parent's **individual** fertility also rises when offered its
fixed-fertility planner bundle. For equal cohort masses define

\[
(X_L^\circ,V_L^\circ)
=\left(\frac{x_L+\ell}{2},\frac{v_L+\ell}{2}\right).
\]

Concavity and (7) give
\(g(X_L^\circ,V_L^\circ)>g(x_L,v_L)\). Since \(x>x_L\), the aggregate
planner bundle at \(B\geq f\ell+(1-f)x\) dominates this pair coordinatewise:

\[
X\geq fX_L^\circ+(1-f)x\geq X_L^\circ,
\quad
V\geq fV_L^\circ+(1-f)x\geq V_L^\circ.
\]

Thus \(g(X,V)>n_L/\vartheta\). The parent's fertility derivative at its old
choice is positive, and strict concavity gives a higher private fertility
choice at the assigned total goods and housing. Its consumption need not rise.

For an arbitrary tenure mix, the aggregate theorem implies this conditional
individual result for every parent with \(n_i\leq\bar n^E\). Calling those
parents the poorer type additionally requires a proved fertility–income
ranking. Global tenure switches should not be assumed to preserve that ranking.
The original verified mixed regime supplies it; the general aggregate theorem
does not need it.

### A fully primitive cost bound in this regime

Let \(p_-\) be the original explicit cutoff where \(M(p_-)=kp_-\), with
\(k=(1+\Gamma)/(1+\Gamma\ell/m)\). Since equilibrium rent exceeds \(p_-\),
it suffices to impose \(\kappa p_-\geq\alpha\chi\). If eliminating even the
cutoff notation is useful, substitution into the known quadratic shows this is
equivalent to

\[
\boxed{
\frac{\kappa w_L}{\chi r}
\geq
\frac{\alpha k(k+A)+\alpha+k(1+\vartheta)}{kA+1},
\qquad A=\alpha+\vartheta.}
\tag{20}
\]

The quadratic has one positive root, is negative before that root, and positive
after it; evaluating it at \(p=\alpha\chi/\kappa\) proves the equivalence.
The compact price condition is economically clearer. Equation (20) is supplied
for a purely primitive verification, not as a recommendation to clutter the
paper statement.

## 6. A less conservative primitive income test, with arbitrary actual tenure

This section uses a stationary matching of old and young type distributions.
It is optional; the income floor in Section 4 does not need that matching.

Set \(\rho=\beta/q\) and define the unconstrained lifetime adult consumption
and two explicit reference allocations:

\[
x_i^u=\frac{w_i+qy}{J+\beta K},\qquad
b_i^0=\max\left\{\frac yK,\rho x_i^u\right\},\qquad
e_i^0=\min\left\{\frac{w_i}{J},x_i^u\right\}.
\tag{21}
\]

They obey \(Je_i^0+qK b_i^0=w_i+qy\). At any actual young optimum,
the continuation-wealth KKT condition implies

\[
z_i/K\geq\rho x_i\geq\rho e_i,
\qquad z_i/K\geq y/K.
\]

Combining this with \(Je_i+qz_i=w_i+qy\) gives

\[
z_i/K\geq b_i^0,
\qquad e_i\leq e_i^0.
\tag{22}
\]

Therefore the primitive condition

\[
\boxed{\overline{b_i^0}\geq\overline{e_i^0}}
\tag{23}
\]

implies \(B\geq\bar e\), for any actual tenure configuration. This is
strictly more permissive than using only the income floor whenever private
wealth accumulation supplies some of the required old resources. It does not
require a rich household's unconstrained ownership plan to be feasible.

When \(\beta\geq q\), (23) is automatic. As \(\beta\downarrow0\), (23)
reduces to \(y/K\geq\bar w/J\). Thus there is no universal positive lower
bound on patience once the income profile supplies the needed old resources.

In the original regime with \(w_L<Jm<w_H\), the reference objects reduce to
\(e_L^0=\ell,e_H^0=x,b_L^0=y/K,b_H^0=\rho x\), reproducing (1).
When \(D=y/K-\ell>0\) and \(\beta<q\), that resource test is

\[
fD\geq(1-f)(1-\beta/q)x.
\]

Equivalently, holding the income-ordering branch fixed,

\[
\beta\geq
q\frac{(1-f)(w_H+qy)-fJD}
{(1-f)(w_H+qy)+fqKD}.
\tag{24}
\]

If the right side is negative it imposes no additional restriction. This is a
sufficient patience threshold for a rent-uniform result, not a necessary
threshold at a particular strongly distorted equilibrium. The separate
\(\beta K>\Gamma A\) in the old construction is needed for its unrestricted
rich-owner branch; it is not intrinsic to the all-tenure planner theorem.

## 7. Sharpness and the remaining obstacle

The new relative child-cost bound is sharp for a uniform comparison at the
resource boundary and arbitrarily small housing wedges. With one affected
type of current-expenditure index \(e\), hold its expenditure fixed and let
\(t=1+\epsilon\). Equation (7) shows that its fertility is below the
unrestricted same-expenditure level to first order if
\(\alpha\lambda<1\), and above it if \(\alpha\lambda>1\). At \(B=e\),
the planner moves partway back toward the unrestricted adult bundle. A
first-order expansion at \(t=1\) therefore makes the planner/private
fertility difference negative when \(\alpha\lambda>1\). The concavity gain
is only second order at that point. Thus children being sufficiently housing
intensive is doing real work when a goods loss is permitted.

Nor can a statement based only on positive retirement income and binding
constraints determine the full equal-weight allocation. In the original
verified regime, as rent approaches \(p_+\) from below,
\(x_L,v_L\to\ell\), and both young averages approach
\(T=f\ell+(1-f)x\). If \(B<T\), the full planner gives the young less
housing for rents sufficiently close to this endpoint. At the endpoint itself
it also lowers fertility: all young have \(t_i=1\),
\(\bar n=\vartheta T/C\), while the planner's adult bundles at that mean are
\((T+B)/2\). The replacement condition can select rents arbitrarily close
to the endpoint. A positive but arbitrarily weak constrained mass therefore
does not fix this failure.

The exact unresolved obstacle is now narrow. A clean **aggregate** theorem is
proved with the primitive income and relative child-cost restrictions, for any
young tenure mix. The household work must independently establish a positive
mass with a strict MRS wedge. A claim that **every poorer parent** privately
raises fertility under arbitrary tenure switches additionally needs an
income–fertility ranking or another dominance argument. It is proved here
only in the original checked mixed regime, or for any parent already known to
have fertility below the competitive mean. No general poor-type ranking has
been assumed or certified.

## Verification performed

The argument was checked directly against the private fertility condition,
the consolidated serviced-interest budget, the full dated resource constraints,
and the old allocation \(c^o=z/K,ph^o=\alpha z/K\). Identities (5)–(8),
the quadratic substitution in (20), and the affine inequality leading to (24)
were independently reduced algebraically. A symbolic manipulation was used
only to explore and check an intermediate one-type formula; no numerical
example or simulation supports the theorem. Startup context and the current
checked response/assessment were read. No protected author files were touched.

## Independent cross-review of the all-tenure interval and stationary witness

Second bounded assignment: independently audited
`tmp/theory_designs/overnight_all_tenures.md`, including its completed Section 7.
This cross-review does not repeat or certify the fertility theorem above;
another reviewer owns that independent check.

**Verdict: PASS for the exact income interval and the stated stationary
compatibility family.** No mathematical correction is needed in those claims.
One opening sentence should be qualified: the assertion that no
replacement-fertility condition is used is true of the main household theorem,
but the ancillary stationary witness in Section 7 explicitly and appropriately
uses the replacement bracket. Its inherited old wealth is endogenous to the
preceding cohort's actual choices; it does not rely on the optional fixed-old-
wealth dated closure in Section 6.

### A. Exact affected-income interval, including every optimal tenure plan

For this cross-review use the other note's definitions
\(\delta=\beta K\), \(Y=qy\),
\(A_p=\alpha+\vartheta\kappa p/(\chi+\kappa p)\),
\(J=1+\alpha+\vartheta\), and \(m=Y/\delta\).
The old implementability condition \(\omega_B\geq\alpha\Gamma\) correctly
gives \(V_o(z)=K\log z+\text{constant}\), so the young objectives and feasible
sets in the other note are exactly the stipulated model.

1. Remove only the rental ceiling, retaining \(t=qz\geq Y\). The unique
   solution has
   \[
   x_U=\min\{w/J,(w+Y)/(J+\delta)\},\quad
   h_U=A_px_U/p,\quad t_U=\max\{Y,\delta x_U\}.
   \]
   The inverse of the strictly increasing map \(w\mapsto x_U\) is
   \(w=Jx+\max\{\delta x-Y,0\}\). Thus
   \[
   x_r=pr/A_p,\qquad
   w_r=Jx_r+\max\{\delta x_r-Y,0\}
   \]
   is the exact cap-onset income, including the saving boundary. At equality,
   the cap has zero shadow value. Strictly above it, the constrained rental
   optimum has a strict housing wedge.

2. The fully unrestricted lifetime solution has
   \(x_F=(w+Y)/(J+\delta)\), \(ph_F=A_px_F\), and \(t_F=\delta x_F\).
   Its owner constraint is exactly
   \[
   (\delta-\Gamma A_p)x_F\geq Y.
   \]
   It is feasible if and only if
   \[
   w\geq w_O=\frac{Y(J+\Gamma A_p)}{\delta-\Gamma A_p}
   \]
   when \(\delta>\Gamma A_p\), and never feasible when
   \(\delta\leq\Gamma A_p\). Accordingly, \(w_O=+\infty\) in the latter
   case. Equality has zero owner-constraint multiplier. Below this threshold,
   the owner optimum has a strictly positive multiplier and strict housing
   wedge. Endogenous fertility makes the owner feasible at every \(w>0\):
   choose sufficiently small positive goods, space, and fertility and a
   retirement balance just above the positive equity requirement.

3. Every owner plan is feasible in the hypothetical uncapped-rental problem,
   since its equity restriction implies \(t>Y\). When \(w\leq w_r\), the
   unique uncapped-rental optimum is actually available as a rental. Strict
   concavity means that another globally optimal tenure plan must reproduce
   the same real allocation; a distinct distorted plan cannot tie it. When
   \(w\geq w_O\), the analogous argument applies to the unique fully
   unrestricted lifetime optimum. Thus **every** globally optimal plan has
   MRS equal to \(p\) outside the open interval, including the two endpoints.
   Within \((w_r,w_O)\), both possible tenure optima have strict wedges, so
   every global maximizer does too, even at an exact tenure tie.

4. When finite, the owner threshold in adult-consumption units is
   \(x_O=Y/(\delta-\Gamma A_p)>m\). Its income is the same increasing inverse
   map evaluated at \(x_O\). Therefore
   \[
   w_r<w_O
   \iff x_r<x_O
   \iff \delta pr<A_p(Y+\Gamma pr).
   \]
   If \(\delta\leq\Gamma A_p\), the last inequality always holds, in agreement
   with a finite \(w_r\) and infinite \(w_O\). If \(w_r\geq w_O\), the
   no-wedge rental and ownership regions cover every income. The exact
   interval claim is therefore correct in all parameter cases, without
   imposing the desired tenure or MRS pattern in its assumptions.

The supplementary uniqueness of the within-interval tenure cutoff also
passes. For a constrained owner, the other note's root \(b\in(1,1+\Gamma)\)
gives \(T(b)>J\) and
\[
x_O=w/T(b)
=\frac{w+Y(1+\Gamma-b)/\Gamma}{J+\delta}<x_U.
\]
For a capped renter the displayed identities imply \(x_R>x_U\).
Consequently the optimized owner-minus-renter value is strictly increasing
in income by the envelope theorem. It is negative at \(w_r\), positive at
finite \(w_O\), and eventually positive when \(w_O=\infty\): owner value
has coefficient \(J+\delta\) on \(\log w\), while capped renter value has
coefficient \(1+\delta\). Both housing and fertility remain bounded at a cap,
which explains the latter coefficient. The unique scalar crossing and the
strict wedge under both plans at the tie follow.

### B. Primitive stationary compatibility, rather than a fixed-wealth closure

The following checks concern only the actual stationary family in Section 7.

* The low-income bound implies \(\ell=w_L/J<m\) and
  \(b=(1+\Gamma)\ell/(1+\Gamma\ell/m)<y/K\). For
  \(\rho=\beta/q<1\), the latter inequality is exactly
  \(\ell[1+\Gamma(1-\rho)]<y/K\); for \(\rho\geq1\), it follows from
  \(b<m\leq y/K\). Hence \(w_L<Jy/K<w_H\), so \(f_C\in(0,1)\).
  The positive gap \(y/K-b\) likewise makes \(f_B<1\). Choosing
  \(\max\{f_C,f_B\}<f<1\) is therefore possible and establishes both the
  new strict income condition and the retained stronger Pro resource bound.

* The richer-income restriction uses the upper bound
  \(A_p<\alpha+\vartheta\). It therefore makes the **exact** richer financing
  condition hold at every rent in the witness interval. It includes the
  essential subtraction of \(Y=qy\); no retirement income is double-counted
  as accumulated equity.

* The explicit restriction
  \(0<\chi<\kappa w_L/[r(1+\alpha+k)]\) is a nonempty positive interval.
  At \(p_0=\alpha\chi/\kappa\), it ensures the candidate renter is feasible
  and
  \[
  x_L(p_0)>w_L-(1+\alpha)\chi r/\kappa>k\chi r/\kappa,
  \quad s_L(p_0)<r.
  \]
  Therefore \(M_L(p_0)>kp_0\). The candidate fertility first-order condition
  gives \(n_L'(p)<0\), \(x_L'(p)<0\), and \(s_L'(p)>0\), so
  \(M_L(p)-kp\) strictly decreases. Its unique zero satisfies
  \(p_0<p_-<p_+\). Every rent in the witness interval thus satisfies the
  strict child-cost inequality \(\kappa p>\alpha\chi\), using a finite
  primitive restriction rather than a numerical continuity argument.

* The actual aggregate fertility schedule is
  \[
  F(p)=fn_L(p)+(1-f)\frac{\vartheta x}{\chi+\kappa p}.
  \]
  Both components strictly decrease, and \(p_-<p_+\). Choosing the primitive
  entry conversion from
  \(F(p_+)<1/\nu<F(p_-)\) consequently produces exactly one price in the
  verified regime with \(\nu F(p)=1\). The root is implicit but uniquely
  defined; it is not a claimed closed form or an unspecified local existence
  neighborhood. The choice of \(\nu\) is an explicit nonvacuity construction,
  not proof that an arbitrary externally fixed \(\nu\) passes the bracket.

* Most importantly, the old resources used at that root are actually inherited
  from the previous stationary young choices. The poor rent with zero saving,
  so \(z_L=y\). The rich choose
  \[
  a_H=\delta x-Y-qPh_H,
  \qquad z_H=y+a_H/q+Ph_H=\delta x/q.
  \]
  At the same stationary price these are exactly the resources of the next
  old cohort. Old consumption therefore averages
  \(B=fy/K+(1-f)\beta x/q\), and old housing averages \(\alpha B/p\).
  With young housing \(fr+(1-f)h_H\), the stock-clearing cohort mass is
  \[
  N=\frac{H}{fr+(1-f)h_H+\alpha B/p}>0.
  \]
  The replacement equation makes the next young mass equal to this same
  \(N\). The entry-type fractions remain \(f,1-f\) by the maintained entry
  rule, while old wealth is generated by the preceding young. Estate choice
  is feasible by \(\omega_B>\alpha\Gamma\). Thus the household choices,
  old wealth, stationary age/type distribution, and common housing stock
  close jointly. Section 6's fixed-\(z_o\) price bounds do not enter this
  argument.

This witness deliberately retains the older **stronger** resource bound. It
therefore proves compatibility and nonempty interior for the new theorem,
but does not itself exhibit the newly admitted case in which young goods
consumption falls. That distinction is not a defect in the requested
compatibility result and should remain explicit in any synthesis.


---

## Necessity and independent full-theorem/counterexample reviews

# Housing direction: necessary restrictions and the weakest useful results

Independent bounded mathematical pass, September 10, 2026. Proposal assessment;
no model, deck, manuscript, or external review was changed. This note holds the
two-period model and serviced-interest contract fixed. The numerical fractions
below are exact analytical counterexamples, not simulations or calibration.

## Judgment

Finance and rental restrictions establish a **relative-price distortion**. They
do not, by themselves, establish the direction of the full equally weighted
planner's allocation, even for a poorer household that is strictly constrained.
An economic restriction on current resources or available space is indispensable.

The simplest honest housing theorem need not identify which young households
rent or own. If retirement choices are unrestricted, define each young family's
current expenditure divided by its utility-weight sum as

\[
e_i=\frac{w_i-q(z_i-y)}{J},\qquad J=1+\alpha+\vartheta.
\]

Then a mean-old-consumption condition \(B\geq\mathbb E[e_i]\), together with a
positive mass of distorted young choices, gives the young more aggregate
housing. It is strictly weaker than \(B\geq\mathbb E[x_i]\), where
\(x_i=c_i-\chi n_i\) is adult goods consumption. A transparent primitive
sufficient condition is

\[
\boxed{\frac yK\geq\frac{\mathbb E[w_i]}J},\qquad
K=1+\alpha+\omega_B,
\]

provided every old household enters with retirement income \(y\) plus
nonnegative net assets. This condition makes a substantive claim about the
resource distribution; it is not generated by finance. It can be conservative
because it ignores old accumulated wealth and young saving. An alternative
simple stationary sufficient restriction is \(\beta\geq q\), with the same
unrestricted-old condition. Neither restriction is necessary.

An individual housing claim has a different, weaker boundary than the
age-wide claim. Some poor households can gain while the young lose housing in
aggregate. Conversely, an aggregate gain does not automatically imply that
every affected young household gains.

## 1. Exact planner geometry: no household regime is needed

Write

\[
x_i=c_i^y-\chi n_i>0,\qquad s_i=h_i^y-\kappa n_i>0.
\]

The planner preserves each young household's continuation wealth, each old
household's estate, and initially each individual's fertility. It chooses all
current goods and housing with equal weight on each living household. Let
\(N_y,N_o>0\) be the two cohort masses. Let \(\bar x,\bar s\) be young means,
and \(B,\bar h_o\) old mean consumption and housing. Strict concavity gives

\[
X=\frac{N_y\bar x+N_o B}{N_y+N_o},\qquad
S=\frac{N_y\bar s+N_o\bar h_o}{N_y+N_o}.
\]

Every young household receives \((X+\chi n_i,S+\kappa n_i)\), and every old
household receives \((X,S)\). Hence the **necessary and sufficient** boundaries
are

\[
\boxed{\Delta H_y>0\iff\bar h_o>\bar s},
\qquad
\Delta H_y=\frac{N_yN_o}{N_y+N_o}(\bar h_o-\bar s),
\]

\[
\boxed{\Delta h_i^y>0\iff S>s_i},\qquad
\boxed{\Delta c_i^y>0\iff X>x_i}.
\]

An old home need not exceed every young home. The aggregate comparison uses
old housing versus **young housing after children's space needs**, averaged
within each cohort. Unequal cohort masses change the magnitude but not the
aggregate sign.

These formulas also prove why willingness to pay for housing is insufficient:
the private wedge is a ratio \(\alpha x_i/s_i\), while the full planner's
individual direction depends on the level \(s_i\) relative to a resource mean.

## 2. All-tenure housing inequalities

Set \(p=(1+\tau-q)P\), \(\Gamma=q(1-\phi)/(1+\tau-q)>0\), and
\(\delta=\beta K\). Suppose

\[
\omega_B\geq\alpha\Gamma.
\]

Then unrestricted old demand is feasible as ownership, so

\[
V^o(z)=K\log z+\text{constant},\quad
c^o=z/K,\quad h^o=\alpha z/(Kp).
\]

Thus \(\bar h_o=\alpha B/p\). This old-finance restriction is sufficient,
not necessary: unrestricted realized old choices could also be implemented by
uncapped rental. Some guarantee of undistorted old demand is required for the
following simple formulas; a young-side restriction alone does not supply it.

For any optimal young tenure, define

\[
t_i=\frac{\alpha x_i}{p s_i}\geq1.
\]

The inequality follows from the housing first-order condition. The rental
ceiling or owner equity requirement adds a nonnegative housing multiplier.
There is no need to identify the global tenure winner. A strictly distorted
household has \(t_i>1\).

At an interior private fertility choice,

\[
n_i=\frac{\vartheta x_i}{\chi+t_i\kappa p}.
\]

Both tenure budgets imply

\[
w_i-q(z_i-y)=x_i+ps_i+(\chi+\kappa p)n_i.
\]

Let \(d_i=ps_i/\alpha=x_i/t_i\). Direct substitution gives the useful
sandwich

\[
\boxed{d_i\leq e_i\leq x_i},\qquad
e_i=\frac{w_i-q(z_i-y)}J\leq\frac{w_i}J.
\]

Both inequalities in the first sandwich are strict when \(t_i>1\). To verify
the first without assuming a housing ordering,

\[
\frac{J e_i}{d_i}
=t_i+\alpha+
\frac{\vartheta t_i(\chi+\kappa p)}{\chi+t_i\kappa p}
\geq1+\alpha+\vartheta=J.
\]

To verify the second, use \(ps_i\leq\alpha x_i\) and
\((\chi+\kappa p)n_i\leq\vartheta x_i\). Finance gives \(z_i\geq y\), so
the last upper bound follows from the budget.

Consequently the exact aggregate housing condition is

\[
\boxed{B>\mathbb E[d_i]
=\bar x-\mathbb E\!\left[x_i\left(1-\frac1{t_i}\right)\right]}.
\]

The old can have **less** mean consumption than the young have adult
consumption, by as much as the mean housing distortion on the right. This is
the precise answer to how much weaker the space condition is.

There are three nested tests:

\[
B>\mathbb E[d_i]\quad\text{(exact housing boundary)},
\]
\[
B\geq\mathbb E[e_i]\quad\text{(simple sufficient test with a strict wedge)},
\]
\[
B\geq\mathbb E[x_i]\quad\text{(stronger goods-and-space test)}.
\]

The gaps can be measured exactly:

\[
x_i-e_i=
\frac{x_i(t_i-1)}J
\left[\frac\alpha{t_i}+
\frac{\vartheta\kappa p}{\chi+t_i\kappa p}\right],
\]
\[
e_i-d_i=
\frac{x_i(t_i-1)}{Jt_i}
\left[1+\frac{\vartheta\chi}{\chi+t_i\kappa p}\right].
\]

The primitive bound \(y/K\geq\mathbb E[w_i]/J\) implies the middle test if
old wealth is at least \(y\). The ratio

\[
K/J=(1+\alpha+\omega_B)/(1+\alpha+\vartheta)
\]

adjusts retirement versus young resources for the utility shares allocated to
estates and children. If \(\omega_B=\vartheta\), the primitive condition is
simply retirement income at least mean young endowment. It is not a claim
that retirement income is empirically that high. Also, the primitive test is
not uniformly weaker than every existing primitive restriction: replacing
actual old wealth by its income floor and ignoring young saving can tighten it.

At an arbitrary inherited date, \(z_o\geq y\) must be stated as an inherited
net-resource condition. It is automatic in the maintained stationary
constant-price budgets, but a price decline after mortgage origination can
make net housing equity negative. No transition claim is being made here.

## 3. A second simple route: stationary smoothing

For an optimal renter, the saving first-order condition gives
\(z_i\geq\delta x_i/q\). For an owner, let \(\mu_i\geq0\) multiply
\(q(z_i-y)-\Gamma ph_i\geq0\). The retirement-wealth condition is

\[
\delta/z_i=q(1/x_i-\mu_i)\leq q/x_i,
\]

so exactly the same inequality follows. Hence, when the old are the stationary
predecessors of these young choices,

\[
B=\mathbb E[z_i/K]\geq(\beta/q)\bar x.
\]

Therefore \(\beta\geq q\), unrestricted old demand, and a positive young
housing distortion imply aggregate young housing gains for **any** young
tenure mix and endowment distribution. There is no requirement to support a
particular richer-owner/poorer-renter configuration.

The economics is simple: preferences and the intertemporal price do not
intrinsically tilt adult consumption toward youth. When \(\beta<q\),
unrestricted rich households naturally consume less when old; that can
overwhelm the housing distortion. The exact boundary above allows many
\(\beta<q\) economies. This is an alternative sufficient theorem, not a
proposed new restriction to impose silently.

## 4. Individual gains and the existing mixed regime

With equal cohort masses, household \(i\) gains housing exactly when

\[
B+\mathbb E[d_j]>2d_i.
\]

For any strictly distorted young choice, the budget-and-fertility calculation
also gives a primitive upper bound

\[
\boxed{s_i<\frac{\alpha w_i}{Jp}}.
\]

Thus a transparent sufficient individual test is

\[
\boxed{S\geq\frac{\alpha w_i}{Jp}}.
\]

This is useful for a genuinely low-endowment household. It compares the
available adult-space average with an upper bound on that household's private
space, rather than assuming the desired private/planner ordering. Given stock
and predetermined cohort masses,

\[
S=\frac{\bar H-\kappa N_y\bar n}{N_y+N_o}.
\]

For fixed exogenous fertility instead of an endogenous private optimum, use
the corresponding bound
\(s_i\leq\alpha[w_i-(\chi+p\kappa)n_i]/[p(1+\alpha)]\), not the
\(J\)-bound. The planner may hold endogenous equilibrium fertility fixed; the
\(J\)-bound remains valid in that comparison because the original private
fertility first-order condition still describes the initial allocation.

In the existing Pro regime, let \(L\) be capped renters and \(H\) unrestricted
owners, and retain its notation \(\ell=w_L/J\), \(x=x_H\), and low share
\(f\). The proof already establishes \(x>\ell\). The new bound gives

\[
d_L<\ell<x=d_H.
\]

Consequently:

* Aggregate housing gains have the exact boundary
  \(B>fd_L+(1-f)x\). The **simpler sufficient condition** is
  \(B\geq f\ell+(1-f)x\). The previous \(b>\ell\) bound is unnecessary for
  housing alone.
* Poorer individual housing gains have the exact boundary
  \(B+(1-f)x>(2-f)d_L\). A sufficient primitive form is
  \(B+(1-f)x\geq(2-f)\ell\). This is strictly weaker than the aggregate
  sufficient condition because \(x>\ell\).
* A still simpler, somewhat stronger individual condition is \(B\geq\ell\).
  It does not require \(B\) to exceed mean young adult consumption.
* These housing-only conditions do not by themselves establish goods gains or
  fertility gains. Those require their own exact tests.

Do not claim that adult space is generally monotone in endowment across all
young regimes. Among capped renters, total housing is fixed while chosen
fertility rises with income, so adult space can fall. Individual claims should
use the exact boundary or a proved bound, not an unproved rank argument.

## 5. Exact primitive counterexample: even the affected poorer young can lose

All coefficients are positive, the financed share is 80%, the property tax is
positive, and fertility is privately chosen. Set

\[
q=\tfrac12,\quad\phi=\tfrac45,\quad\tau=\tfrac1{100},\quad
\beta=\tfrac14,\quad\alpha=\omega_B=1,\quad
\vartheta=\tfrac15,\quad\chi=\kappa=\tfrac1{10},\quad y=\tfrac95.
\]

Then \(K=3\), \(J=11/5\), \(\Gamma=10/51\), and \(m=6/5\). Choose
\(p=1\), \(P=100/51\), and

\[
r=\frac{1106}{1005},\quad
w_L=\frac{44441}{20100},\quad w_H=\frac{441}{100},\quad f=\frac9{10}.
\]

The exact household choices are

| Object | Poorer young | Richer young |
|---|---:|---:|
| Adult goods \(x\) | \(101/100\) | \(9/5\) |
| Adult space \(s\) | \(1\) | \(9/5\) |
| Fertility \(n\) | \(202/201\) | \(9/5\) |
| Total housing \(h\) | \(1106/1005=r\) | \(99/50\) |
| Total goods \(c\) | \(22321/20100\) | \(99/50\) |
| Retirement wealth \(z\) | \(9/5=y\) | \(27/10\) |
| Tenure | Strict capped renter | Unrestricted owner |

Both fertility first-order conditions hold exactly. The poorer type has
\(t_L=101/100>1\), \(x_L<m\), and strictly prefers renting to **every**
ownership plan. To verify the latter, the concavity bound in the response
applies whenever

\[
t_L-1<\Gamma(1-x_L/m).
\]

Here it reads

\[
\frac1{100}<\frac{19}{612}.
\]

The inequality proves all tenure, saving, and fertility deviations are covered,
not merely a local rental optimum. Old unrestricted ownership is feasible
because \(1>10/51\). The rich exact equity test is

\[
q(z_H-y)=\frac9{20}>\frac{33}{85}=\Gamma p h_H.
\]

Even the response's stronger price-independent rich-finance test holds:
\((\beta K-\Gamma A)x_H=63/68>9/10=qy\), where \(A=6/5\).

Old consumption and housing are \(3/5\) for the poor type and \(9/10\) for
the rich type. Thus

\[
\bar s=\frac{27}{25},\quad \bar x=\frac{1089}{1000},\quad
B=\bar h_o=\frac{63}{100}.
\]

At equal cohort mass one, the full planner gives

\[
S=\frac{171}{200}=0.855,
\qquad X=\frac{1719}{2000}=0.8595.
\]

Hence the **strictly constrained poorer young loses** \(29/200\) units of
housing and \(301/2000\) units of goods. The aggregate young housing loss is
\(9/40\). This is not a knife-edge counterexample: every relevant inequality
has a strict margin.

It is a complete stationary construction. Set

\[
\nu=\frac{3350}{3633},\qquad \bar H=\frac{30459}{16750},\qquad N_y=N_o=1.
\]

Then \(\nu\bar n=1\) and housing clears exactly. Competitive external finance
gives the stated \(p,P\). Taxes fund the same additive public service as the
response. No new timing, taste, or asset restriction is used. The constructed
price lies inside the strict capped-renter interval; equivalently, selecting
the displayed replacement factor implements its strictly interior fertility
root. This example obeys the existing regime's finance and saving restrictions
but violates its additional resource condition.

The indispensable economic issue is \(\beta/q=1/2\): both types enter old
age with lower adult consumption than the poor young type's adult-space
equivalent. Full utilitarian redistribution therefore sends resources toward
old households despite the poor young housing wedge.

## 6. Exact counterexample separating individual from age-wide housing gains

Keep the same coefficients, poorer endowment, rental ceiling, and \(p=1\).
Change only

\[
f=\tfrac12,\qquad x_H=3,\qquad w_H=\tfrac{159}{20}.
\]

The rich unrestricted optimum has \(n_H=s_H=3\), \(h_H=c_H=33/10\),
\(z_H=9/2\). Its finance constraint is strictly slack. Then

\[
\bar s=2,\quad\bar x=\frac{401}{200},\quad B=\frac{21}{20},
\quad S=\frac{61}{40},\quad X=\frac{611}{400}.
\]

The poor household gains \(21/40\) units of housing and \(207/400\) units
of goods, while the young in aggregate **lose** \(19/40\) units of housing.
Stationarity and housing clearing hold at cohort mass one with

\[
\nu=\frac{402}{805},\qquad \bar H=\frac{6533}{2010}.
\]

This establishes that the individual claim is substantively broader than the
aggregate age-direction claim. It is not a reason to abandon the latter; it
specifies which further resource inequality the latter needs.

## Recommended presentation scope

Use the all-tenure distortion result to establish the friction. State the
dated planner's exact aggregate space boundary, then one transparent resource
sufficient condition such as \(y/K\geq\mathbb E[w]/J\). State individual
gains through their own proved low-resource cutoff. Do not carry over the
previous goods-gain restriction merely because it makes a different fertility
proof convenient. If desired, retain \(\beta\geq q\) as an alternative simple
stationary benchmark, clearly labeled rather than adopted by default.

The central limitation is economic rather than algebraic: an equal-weight
planner redistributes goods and space jointly. An assertion about who faces a
housing wedge cannot determine that planner's direction unless it is paired
with information about where current resources and adult space are located.

## 8. Independent certification of the joint-fertility theorem

Second, distinct bounded scope assigned September 10, 2026 at 04:10 UTC:
independently cross-check `tmp/theory_designs/overnight_fertility_sharp.md`.
**Verdict: PASS for the general full dated planner, aggregate housing, and
joint mean-fertility theorem.** No consequential algebraic correction is
needed. The old resource condition and child-cost condition remain substantive
assumptions, and the stated limits on individual conclusions are correct.

This review independently checked identities (5)–(9), the full planner
first-order conditions and resources, Jensen strictness for every positive
cohort ratio, the primitive income reduction, and the optional sharper
mean-fertility and stationary-resource tests. It did not reopen private global
tenure selection or independently certify the optional cutoff polynomial (20);
those are separate household-verification objects.

### 8.1 The comparison chain is exact, including its strictness

To avoid confusing the current-expenditure index with estates, continue to use
\(e_i=(c_i+ph_i)/J\) only in this section. Let

\[
a=\kappa p>0,\quad C=\chi+a,\quad
\lambda=\chi/a,\quad v_i=ps_i/\alpha,\quad t_i=x_i/v_i\geq1.
\]

The independent expansion starts from

\[
g(x,v)=\frac{xv}{\chi v+a x},\qquad
n_i=\vartheta g(x_i,v_i),\qquad
\frac{e_i}{v_i}
=\frac{t_i+\alpha+\vartheta t_i(1+\lambda)/(\lambda+t_i)}J.
\]

Subtracting \(v_i\), \(x_i\), and \(Cg(x_i,v_i)\) yields, respectively,

\[
e_i-v_i=\frac{v_i(t_i-1)}J
\left(1+\frac{\vartheta\lambda}{\lambda+t_i}\right),
\]
\[
x_i-e_i=\frac{v_i(t_i-1)}J
\left(\alpha+\frac{\vartheta t_i}{\lambda+t_i}\right),
\]
\[
e_i-Cg(x_i,v_i)
=\frac{v_i(t_i-1)(t_i-\alpha\lambda)}{J(\lambda+t_i)},
\]
\[
Cg(x_i,v_i)-v_i
=\frac{\lambda v_i(t_i-1)}{\lambda+t_i}.
\]

Every displayed identity matches the fertility note. Under
\(\kappa p\geq\alpha\chi\), equivalently \(\alpha\lambda\leq1\),

\[
\boxed{v_i\leq Cn_i/\vartheta\leq e_i\leq x_i}.
\]

Every inequality is strict for \(t_i>1\). In particular, when the cost
restriction holds at equality, the numerator of the middle difference is
\(v_i(t_i-1)^2\), so the strict result does not require a strict child-cost
inequality. Conversely, a constraint that binds with a zero multiplier does
not supply strictness. The theorem properly needs a positive mass with a
strict housing MRS wedge.

### 8.2 Arbitrary cohort masses and Jensen are valid

Let \(\mu=N_o/N_y>0\), and let \(\bar g\) average \(g(x_i,v_i)\) over
young households. The independently checked Hessian is

\[
d^2g=-\frac{2\chi a(v\,dx-x\,dv)^2}
{(\chi v+a x)^3}\leq0.
\]

The function is increasing and concave, with its only linear directions
given by rays through the origin. The full fixed-fertility planner has

\[
(X,V)=\frac{(\bar x,\bar v)+\mu(B,B)}{1+\mu}.
\]

Under \(B\geq\bar e\), the chain in Section 8.1 and positive distorted mass
give \(B/C>\bar g=\bar n/\vartheta\). Therefore

\[
g(X,V)\geq
\frac{g(\bar x,\bar v)+\mu B/C}{1+\mu}
\geq\frac{\bar g+\mu B/C}{1+\mu}
>\bar g.
\]

All population factors are correct; neither equal cohort masses nor equality
of the individual young bundles is used. This proof does not require
\(X>\bar x\). Indeed, at \(B=\bar e\), strict wedges imply
\(X<\bar x\), exactly as reported in the note.

The optional sharper condition

\[
B\geq C\bar n/\vartheta
\]

also passes and dispenses with the relative child-cost restriction. First,
\(\bar v<C\bar n/\vartheta\leq B\) proves the housing sign. For fertility,
the first Jensen inequality becomes **strict**: positive wedges imply
\(\bar x>\bar v\), so \((\bar x,\bar v)\) is not proportional to
\((B,B)\). The Hessian formula verifies strict concavity along that segment.
It follows that \(g(X,V)>\bar n/\vartheta\), even when the sharper resource
condition holds at equality. This test is an equilibrium-object test, not a
fully primitive restriction.

### 8.3 This is the full joint optimum, not a housing-only improvement

Fixed individual young continuation targets and fixed individual old estates
enter the dated welfare objective as constants. The planner still chooses
every current household's goods and housing. Its goods first-order conditions
make every adult surplus equal to \(X\); its housing conditions make every
adult space surplus equal to \(S\). Each parent's fertility condition then
sets the same unique \(n\). Heterogeneous continuation targets create no
remaining parent-specific term in that condition.

Normalizing current resources by young mass gives

\[
X(n)=\frac{\mathcal C-\chi n}{1+\mu},\qquad
S(n)=\frac{\mathcal H-\kappa n}{1+\mu}.
\]

The full reduced objective is

\[
(1+\mu)\log X(n)+(1+\mu)\alpha\log S(n)
+\vartheta\log n+\text{constants}.
\]

Its derivative and second derivative are

\[
D(n)=\frac{\vartheta}{n}-\frac\chi{X(n)}
-\frac{\alpha\kappa}{S(n)},
\]
\[
D'(n)=-\frac{\vartheta}{n^2}
-\frac{\chi^2}{(1+\mu)X(n)^2}
-\frac{\alpha\kappa^2}{(1+\mu)S(n)^2}<0.
\]

The derivative goes from positive to negative infinity on the feasible
interval. At \(n=\bar n^E\), the adult bundles are exactly the
fixed-fertility bundles above. Hence \(g(X,V)>\bar n^E/\vartheta\) means
\(D(\bar n^E)>0\), proving the unique joint optimum satisfies
\(n^F>\bar n^E\).

Its mean young housing is

\[
h_y(n)=S(n)+\kappa n
=\frac{\mathcal H}{1+\mu}+\frac{\mu\kappa}{1+\mu}n.
\]

Thus joint fertility adds exactly
\(\mu\kappa(n^F-\bar n^E)/(1+\mu)>0\) to the already positive
fixed-fertility housing gain. No unallocated extra goods or space appear.

The objective counts the fixed set of existing parents and old households.
The \(\vartheta\log n\) term is the existing parent's fertility preference;
there is no added welfare weight, resource endowment, or utility term for
unborn people. This remains a dated parental-welfare comparison, not a
population or policy implementation theorem.

### 8.4 Financial settlement preserves the promised continuation objects

The feasibility claim also passes under the stated permission to relax
finance. For any reassigned owned home, keeping a young target \(z_i\geq y\)
requires net financial position

\[
a_i'=q(z_i-y)-qPh_i'.
\]

Allowing the financed share to reach one makes this feasible because
\(a_i'\geq-qPh_i'\). For an old household with fixed positive estate
\(E_i\), use \(a_i'=qE_i-qPh_i'\), with the same feasibility argument.
Current transfers are \(T_i=\Delta c_i+p\Delta h_i\). They sum to zero
because the full planner preserves aggregate current consumption and the
common housing stock. Intermediary positions adjust consistently. The tax
base remains \(P\bar H\), so the maintained public-service expenditure is
unchanged. This settlement does not require lifting the fixed estate or
continuation targets, but it does use the expressly allowed relaxation of
finance. It is not a transfers-only claim with the original credit limit.

### 8.5 The primitive reduction and its qualifications are correct

Every actual young budget has \(z_i\geq y\), so
\(e_i\leq w_i/J\). With unrestricted old demand and inherited old resources
at least \(y\), \(B\geq y/K\). Hence

\[
\boxed{y/K\geq\bar w/J\ \Longrightarrow\ B\geq\bar e}
\]

uses no ordering between \(\beta\) and \(q\), no prescribed tenure mix, and
no ranking of young fertility. It does require the stated old resource floor
at an arbitrary inherited date. The old-finance bound
\(\omega_B\geq\alpha\Gamma\) is a valid sufficient way to make old
unrestricted demand attainable, including at equality.

The relative child-cost condition still contains equilibrium rent. The
optional price-free sufficient bound in the note is correct for predetermined
old mass and stock: positive young housing and market clearing imply

\[
p>\frac{\alpha N_o y}{K\bar H}.
\]

Therefore \(\kappa N_o y/(K\bar H)\geq\chi\) implies
\(\kappa p>\alpha\chi\). In a stationary construction with endogenous
cohort mass it is not an entirely primitive condition, exactly as the note
acknowledges.

The optional stationary refinement (21)–(23) also passes. With
\(x_i^u=(w_i+qy)/(J+\beta K)\), set

\[
b_i^0=\max\{y/K,(\beta/q)x_i^u\},\qquad
e_i^0=\min\{w_i/J,x_i^u\}.
\]

The continuation KKT condition gives
\(z_i/K\geq(\beta/q)x_i\geq(\beta/q)e_i\). Combining this with
\(Je_i+qz_i=w_i+qy\) proves
\(z_i/K\geq b_i^0\) and \(e_i\leq e_i^0\). Thus
\(\overline{b_i^0}\geq\overline{e_i^0}\) is a valid less conservative
primitive route when current old and young type distributions are matched in
stationarity. It must not be imported into an arbitrary inherited-date claim
without that matching. No hidden unrestricted-rich-owner assumption is used.

Finally, no claim that all poor parents gain goods or that fertility is
globally increasing in endowment appears in the general proof. The
conditional individual implication is precisely limited: because
\(g(X,V)>\bar n^E/\vartheta\), a parent with \(n_i\leq\bar n^E\) has a
positive fertility derivative at its fixed-fertility planner bundle. Calling
that parent the poorer type requires the separate ranking argument that the
note correctly declines to assume under arbitrary tenure changes.

**Certification scope:** the general housing and joint mean-fertility theorem
passes as stated, with its resource, relative child-cost, old-choice, and
positive-wedge assumptions retained. Verification was analytical; no model
solve, simulation, browser action, or build was used in this second pass.

## 9. Independent check of the finite private-rechoice counterfamily

Third, distinct bounded scope, September 10, 2026: verify Section 9.2–9.3 of
`tmp/theory_designs/overnight_all_tenures.md`. **Verdict: PASS throughout
\(0<\epsilon\leq1/10000\). No consequential correction found.** The family
is an actual stationary equilibrium of the maintained model. At its full
fixed-fertility planner bundles, subsequent **private** rechoice lowers mean
fertility despite the strong primitive income condition. This does not
contradict the independently certified **joint** planner result.

The following checks were completed independently, without a simulation or
numerical equilibrium root.

1. **Exact child choice and radical.** With the common positive child-cost
   coefficient \(a_c\), putting \(m=a_cn\) changes young utility only by a
   constant. Its first-order condition is
   \(3m^2-2(c+h)m+ch=0\), and the unique feasible root is
   \(N(c,h)=[c+h-\sqrt{c^2-ch+h^2}]/3\). At
   \((c_H,h_H)=(4+2\epsilon,5/2-\epsilon)\), this is exactly \(1-d\)
   with the stated radical. Rationalization gives (C3) exactly; its positive
   numerator and denominator prove \(0<d<\epsilon\). Therefore
   \(a_c=1/3-d/4>0\) throughout the stated range.

2. **Both young rental optima are global.** The poor choice has
   \((x,s,m)=(1,1/2,1/3)\), fits the rental cap, and uniquely maximizes
   utility over the whole uncapped rental set with zero saving. Every owner
   plan is in that larger set, while this unique optimum violates owner
   equity finance at \(qz=Y=4\); every owner deviation is strictly worse.
   The rich capped choice has \(3<x_H<4<Y/\delta=8\), so saving zero is
   optimal. Its exact multipliers satisfy
   \(0<\xi<\epsilon\), \(\zeta>1/8\), and
   \(\Gamma p\zeta>1/16>\epsilon\). The supporting-hyperplane bound
   \(U-U_H\leq\xi(h-r)-\zeta(qz-Y)\) follows directly from the common
   budget and includes changing fertility, saving, and housing. Owner finance
   makes it strictly negative. Rental deviations are also excluded by the
   same bound and strict concavity. There is no unexamined tenure deviation.

3. **Old choices and inherited wealth are genuine.** Here
   \(K=4\), \(J=3\), \(\Gamma=1/4\), and old unrestricted demand is
   feasible because \(\omega_B=2>1/4\). Every old household enters with
   \(z=8\) because the actual previous young saved zero and rented; the
   old resources are not separately imposed. Old choices
   \((c_o,h_o,E_o)=(2,1,8)\) satisfy the consolidated budget
   \(2+2+\tfrac12 8=8\). They also satisfy the original owner budget with
   net financial position \(a_o=3\), producing estate
   \(a_o/q+Ph_o=8\). Rental implementation gives the same optimum. Estate
   payments remain external to entrant endowments, as maintained in the model.

4. **Replacement, market clearing, and the resource floor are exact.** Mean
   child needs are \(\bar m=2/3-d/2=2a_c\), so \(\bar n=2\) exactly and
   \(\nu=1/2\) gives unit cohort replacement. Housing use is
   \(\tfrac12(5/6+5/2-\epsilon)+1=8/3-\epsilon/2=H\).
   Prices satisfy \(p=(1+\tau-q)P=2\). Finally,
   \(B=\bar e=y/K=\bar w/J=2\) and
   \(\kappa p=2a_c>\alpha\chi=a_c\). Half the young have a strict
   housing wedge. Every premise of the sharp joint-planner theorem is met.

5. **The full planner bundles and finite error bound check.** Current resource
   pooling gives exactly
   \(X=2+\epsilon/2+d/4\),
   \(S=1-\epsilon/4+d/4\). The private rechoice bundles are correctly
   \((X+m_i,S+m_i)\), with each original \(m_i\) retained in the assigned
   totals before rechoice. The limiting mean loss is exactly
   \(G_0=(\sqrt{37}+3\sqrt7-14)/18\).
   The rational lower bounds in (C6) are valid: their squared numerators are
   \(60827^2=3699923929<37\cdot10^8\) and
   \(26457^2=699972849<7\cdot10^8\). They give
   \(G_0>11/10000>1/1000\).
   Implicit differentiation of the private first-order condition gives
   \(0<N_c,N_h\) and
   \(N_c+N_h=[(c-N)^{-2}+(h-N)^{-2}]/
   [N^{-2}+(c-N)^{-2}+(h-N)^{-2}]<1\).
   Thus the stated sup-norm Lipschitz bound is valid on the positive domain.
   The poor and rich assigned bundles move by less than
   \(3\epsilon/4\) and \(\epsilon\), respectively. Their mean optimized
   child needs move by at most \(7\epsilon/8\), and the old mean child
   needs fall by \(d/2<\epsilon/2\). Consequently
   \[
   \bar m^{\mathrm{new}}-\bar m
   <-\frac1{1000}+\frac{11}{8}\epsilon
   \leq-\frac{69}{80000}<0.
   \]
   Dividing by the same positive \(a_c\) at each primitive point proves
   \(\bar n^{\mathrm{new}}<2\) on the entire stated interval. The argument
   uses a finite explicit bound, not an unspecified continuity neighborhood.

The conceptual distinction is also correct. A positive fertility derivative
at the competitive mean does not sign the mean of separate nonlinear private
responses to heterogeneous assigned total bundles. In this family the poor
parent's gain is outweighed by the rich parent's reduction. Jointly choosing
fertility and redistributing goods and housing again is a different feasible
optimization and still produces \(n^F>2\). No mean private-rechoice claim
should replace that joint-planner conclusion.
