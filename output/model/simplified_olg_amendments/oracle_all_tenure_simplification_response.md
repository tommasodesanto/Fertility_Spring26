The all-tenure simplification is correct, including the weak child-intensity inequality. The proposed mean-private-fertility conclusion is false, even under the stronger retirement-income floor. Your counterfamily establishes that failure.

The useful result is therefore: the full planner reallocates housing toward the young in aggregate; the joint parents-only planner chooses higher fertility; and a compensated housing transfer raises an affected parent’s fertility. These conclusions do not imply that every constrained parent benefits from the full redistribution, or that mean private fertility rises after parents receive its bundles.

Proposition: housing allocation without selecting young tenure

Retain the frozen model and its proposed external-finance, independent-estate, and tax-financed-public-service closure. Write

\[\delta=\beta K,\qquad Y=qy,\qquad C=\chi+\kappa p,\]

and define

\[A(p)=\alpha+\frac{\vartheta\kappa p}{C},
\qquad x_r=\frac{pr}{A(p)}.\]

The two income thresholds are

\[w_r(p)=Jx_r+\max\{\delta x_r-Y,0\},\]

\[w_O(p)=
\begin{cases}
\displaystyle\frac{Y[J+\Gamma A(p)]}{\delta-\Gamma A(p)},
&\delta>\Gamma A(p),\\[6pt]
+\infty,&\delta\leq\Gamma A(p).
\end{cases}\]

Let
\[\mathcal D(p)=(w_r(p),w_O(p))\]
, interpreted as empty when its endpoints are reversed.

Proposition. At a stationary competitive equilibrium with price
\[p\]
, suppose

\[\boxed{
\omega_B\geq\alpha\Gamma,\qquad
\frac yK\geq\frac{\mathbb E w}{J},\qquad
\Pr\{w\in\mathcal D(p)\}>0.
}
\tag{A}\]

Then the full equally weighted dated planner, initially fixing competitive fertility, assigns more aggregate housing to the young.

If additionally

\[\boxed{\kappa p\geq\alpha\chi,}
\tag{F}\]

the joint dated planner choosing consumption, housing, and parental fertility chooses higher mean fertility. A sufficiently small compensated housing transfer to any financially/rental-distorted young household also raises that parent’s private fertility.

After receiving the fixed-fertility full-planner bundles, every parent whose initial fertility is at or below competitive mean fertility raises its private fertility. There is, however, no general increase in mean private fertility, nor an individual guarantee for every constrained household.

The resource condition can be sharpened to

\[B\geq\mathbb E\widehat x,\qquad
B=\mathbb E c^o,\qquad
\widehat x_i=\frac{c_i+ph_i}{J}.
\tag{R}\]

Under this sharper condition, the allocation and fertility arguments apply to any dated comparison with old/young mass ratio
\[\mu>0\]
, using the old households’ actual resources.

The price-dependent conditions identify the affected income range; they do not assume a constraint multiplier or preselect tenure. An analytical equilibrium construction appears below.

Proof

1. The affected-income interval is exact

Under
\[\omega_B\geq\alpha\Gamma\]
, unrestricted old choices are feasible at every positive wealth level:

\[c^o=\frac zK,\qquad ph^o=\frac{\alpha z}{K},
\qquad qe=\frac{\omega_Bz}{K}.\]

Thus
\[V^o(z)=K\log z+\text{constant}\]
, including the weak-inequality case.

Consider a hypothetical uncapped rental menu that still requires
\[z\geq y\]
. Its unique optimum has

\[x^0=\min\left\{\frac wJ,\frac{w+Y}{J+\delta}\right\},
\qquad
n^0=\frac{\vartheta x^0}{C},
\qquad
h^0=\frac{A(p)x^0}{p}.
\tag{1}\]

Ownership is a subset of this menu because it additionally requires

\[q(z-y)\geq\Gamma ph.\]

The optimum (1) fits within the rental ceiling precisely when
\[w\leq w_r\]
. The fully unrestricted lifetime optimum becomes ownership-feasible precisely when
\[w\geq w_O\]
. Whenever either optimum is available, strict concavity rules out a different, distorted allocation tying it.

For the remaining incomes
\[w_r<w<w_O\]
, the rental optimum has a strictly positive ceiling multiplier, while the owner optimum has a strictly positive financing multiplier. On either branch,

\[\frac{u_h}{u_c}\geq p,\]

strictly in this interval. Consequently,

\[\boxed{
\text{Every globally optimal choice has }u_h/u_c>p
\iff w\in\mathcal D(p).
}
\tag{2}\]

This includes both choices at a tenure tie. No rent–own value crossing needs to be solved.

The interval is nonempty exactly when

\[\delta pr<A(p)[Y+\Gamma pr].\]

The packet’s simpler sufficient test is also correct:

\[w<\frac{JY}{\delta},
\qquad
\frac wJ\left(\frac{\alpha}{p}+\frac{\vartheta\kappa}{C}\right)>r.\]

Very poor households can be undistorted because a small rental suffices.

2. The spending identity signs the housing reallocation

For each competitive young choice define

\[x_i=c_i-\chi n_i,\qquad
v_i=\frac{p(h_i-\kappa n_i)}{\alpha},\qquad
t_i=\frac{x_i}{v_i}\geq1,\]

and

\[g(x,v)=\left(\frac{\chi}{x}+\frac{\kappa p}{v}\right)^{-1}.\]

The fertility first-order condition gives
\[n_i=\vartheta g(x_i,v_i)\]
. Consolidated spending gives

\[J\widehat x_i
=x_i+\alpha v_i+C\vartheta g(x_i,v_i)
=w_i-q(z_i-y).
\tag{3}\]

All four identities in the clarification check. The two decisive ones, with
\[\lambda=\chi/(\kappa p)\]
, are

\[\widehat x-v
=\frac{v(t-1)}J
\left(1+\frac{\vartheta\lambda}{\lambda+t}\right),
\tag{4}\]

\[\widehat x-Cg(x,v)
=\frac{v(t-1)(t-\alpha\lambda)}
       {J(\lambda+t)}.
\tag{5}\]

In particular,
\[v\leq\widehat x\leq x\]
, with strict inequalities at distorted choices. No child-intensity restriction is needed for this ordering.

The full dated benchmark chooses both consumption and housing, includes the old, and preserves young continuation resources and old estates.  Its fixed-fertility solution is

\[X=\frac{\bar x+\mu B}{1+\mu},
\qquad
V=\frac{\bar v+\mu B}{1+\mu},
\qquad S=\frac{\alpha V}{p},
\tag{6}\]

\[(c_i^{y,F},h_i^{y,F})=(X+\chi n_i,S+\kappa n_i),
\qquad
(c_j^{o,F},h_j^{o,F})=(X,S).\]

Positive mass in
\[\mathcal D(p)\]
, together with
\[B\geq\overline{\widehat x}\]
, implies
\[B>\bar v\]
. Hence

\[\boxed{
\bar h_y^F-\bar h_y^E
=\frac{\alpha\mu}{p(1+\mu)}(B-\bar v)>0.
}
\tag{7}\]

The primitive retirement floor suffices at stationarity because

\[B=\frac{\mathbb E z}{K}\geq\frac yK
\geq\frac{\mathbb E w}{J}
\geq\overline{\widehat x}.\]

It must not be substituted automatically for actual inherited old resources away from stationarity.

Importantly, young mean consumption changes by

\[\bar c_y^F-\bar c_y^E
=\frac{\mu}{1+\mu}(B-\bar x).
\tag{8}\]

At
\[B=\overline{\widehat x}\]
, this is strictly negative. The theorem does not obtain its fertility result by assuming a consumption gain.

Financial settlement uses

\[T_i=\Delta c_i+p\Delta h_i,\qquad \sum_iT_i=0,\]

and
\[\Delta a_i=-qP\Delta h_i^{\rm owned}\]
, with matching intermediary positions. Young
\[z_i\]
 and old estates remain fixed. Allowing 100%-balance finance suffices because
\[z_i\geq y\]
; no additional goods or net claims are created.

3. Joint parental fertility

Under (F),
\[\alpha\lambda\leq1\]
, so (5) yields

\[v_i\leq Cg(x_i,v_i)\leq\widehat x_i\leq x_i,\]

with the relevant inequalities strict whenever
\[t_i>1\]
, including at equality in (F). Consequently,

\[\frac BC>\frac{\bar n}{\vartheta}.
\tag{9}\]

The function
\[g\]
 is concave:

\[d^2g
=-\frac{2\chi\kappa p\,(v\,dx-x\,dv)^2}
        {(\chi v+\kappa p x)^3}\leq0.\]

Applying Jensen directly to the young allocations and the old point
\[(B,B)\]
,

\[g(X,V)
\geq
\frac{\mathbb E g(x_i,v_i)+\mu B/C}{1+\mu}
>
\frac{\bar n}{\vartheta}.
\tag{10}\]

For the joint planner, let
\[\mathcal C,\mathcal H\]
 denote fixed current resources per young household. Identical parental preferences imply common chosen fertility, with

\[X(n)=\frac{\mathcal C-\chi n}{1+\mu},
\qquad
S(n)=\frac{\mathcal H-\kappa n}{1+\mu}.\]

The reduced objective has derivative

\[D(n)=\frac{\vartheta}{n}
-\frac{\chi}{X(n)}-\frac{\alpha\kappa}{S(n)},\]

which is strictly decreasing. Equation (10) gives
\[D(\bar n)>0\]
, proving

\[\boxed{n^F>\bar n.}\]

Young aggregate housing rises further by

\[\frac{\mu\kappa}{1+\mu}(n^F-\bar n).\]

This objective values existing parents only; it is not a new stationary-population comparison.

The three fertility experiments must remain distinct

Compensated local reallocation. Give a distorted young household
\[\epsilon\]
 additional housing. An old household’s exact compensation costs

\[c_o\left[
\left(\frac{h_o}{h_o-\epsilon}\right)^\alpha-1
\right]
=p\epsilon+O(\epsilon^2).\]

Taking these goods from the young recipient gives it a strictly positive first-order utility gain,

\[\left(\frac{\alpha}{h_i-\kappa n_i}-\frac p{x_i}\right)\epsilon>0.\]

Its private fertility response has the sign of

\[\frac{\alpha\kappa}{(h_i-\kappa n_i)^2}
-\frac{\chi p}{x_i^2},\]

so a positive first-order response occurs exactly when

\[\kappa p\,t_i^2>\alpha\chi.
\tag{11}\]

Condition (F) suffices for every strictly distorted recipient. This is the friction-specific Pareto improvement, independently of equal welfare weights.

Private rechoice after the full allocation. At assigned total bundles
\[(X+\chi n_i,S+\kappa n_i)\]
, parent
\[i\]
 increases fertility exactly when

\[n_i<n^\dagger,\qquad
n^\dagger=\frac{\vartheta}{\chi/X+\alpha\kappa/S}
=\vartheta g(X,V).
\tag{12}\]

The joint proof establishes
\[n^\dagger>\bar n\]
. Thus every initially below-mean-fertility parent gains fertility, but this cutoff is not the common fertility parents privately choose.

There is also a short income-based sufficient condition. At stationarity, under the primitive retirement floor,

\[\boxed{w_i\leq\frac{Jy}{2K}
\quad\Longrightarrow\quad
n_i^{\rm private,new}>n_i.}
\tag{13}\]

Indeed,

\[n_i\leq\frac{\vartheta w_i}{JC},
\qquad
n^\dagger>\frac{\vartheta B}{2C}
\geq\frac{\vartheta y}{2KC}.\]

This is conservative, but requires no forced-tenure regime.

Joint reoptimization. Equation (10) proves the higher mean fertility of the planner who reallocates consumption and housing as fertility changes. It does not prove higher mean private rechoice at previously assigned bundles.

The counterfamily is valid

Use the parameters in the clarification:

\[\alpha=\vartheta=\chi=\kappa=1,\quad
(q,\beta,\tau,\phi)=\left(\tfrac12,\tfrac18,\tfrac12,\tfrac12\right),\]

\[\omega_B=2,\quad y=8,\quad p=P=2,\]

with equally weighted incomes
\[3,9\]
, and
\[r=5/2-\epsilon\]
.

Both young types rent and save zero. The low type has
\[n_L=1/3\]
. The high type has consumption
\[4+2\epsilon\]
, and its fertility satisfies

\[1-\epsilon<n_H<1
\qquad(0<\epsilon\leq10^{-4}).
\tag{14}\]

The first-order condition evaluated at these two endpoints verifies this bracket.

Ownership deviations are excluded globally. For the high type put

\[\xi=\frac1{r-n_H}-\frac2{4+2\epsilon-n_H},
\qquad
\zeta=\frac1{4+2\epsilon-n_H}-\frac18.\]

The bounds (14) give

\[0<\xi<2\epsilon,\qquad
\Gamma p\zeta>\frac1{16}.\]

Concavity therefore gives, for every owner alternative,

\[U_O-U_R
\leq(\xi-\Gamma p\zeta)h-\xi r<0.
\tag{15}\]

The low type’s unrestricted rental allocation already dominates all alternatives.

Old resources are
\[8\]
, so

\[B=2=\overline{\widehat x}
=\frac yK=\frac{\bar w}{J}.\]

At
\[\epsilon=0\]
, the full planner has
\[X=2,S=1\]
. Private rechoice gives

\[n_L^{\rm new}=\frac{11-\sqrt{37}}9,\qquad
n_H^{\rm new}=\frac{5-\sqrt7}3.\]

Thus mean private fertility falls below
\[2/3\]
 by

\[\Delta_0
=\frac{\sqrt{37}+3\sqrt7-14}{18}
>\frac1{1000}.
\tag{16}\]

A simpler perturbation bound than the proposed one suffices. Set
\[u=1-n_H\in(0,\epsilon)\]
. Then

\[X=2+\frac\epsilon2+\frac u4,\qquad
S=1-\frac\epsilon4+\frac u4.\]

For these preferences, the private fertility function

\[\mathcal N(c,h)=\frac{c+h-\sqrt{c^2-ch+h^2}}3\]

has positive derivatives whose sum is below one. The mean rechoice changes by at most
\[7\epsilon/8\]
; competitive mean fertility changes by at most
\[\epsilon/2\]
. Therefore

\[\overline n_{\rm private,new}(\epsilon)-\bar n_E(\epsilon)
\leq-\Delta_0+\frac{11}{8}\epsilon<0\]

throughout the stated positive-
\[\epsilon\]
 interval.

Housing clears with per-unit-cohort stock
\[8/3-\epsilon/2\]
, and
\[\nu=1/\bar n_E(\epsilon)\]
 supplies stationary replacement. This is an actual stationary counterfamily. Its poorer parent raises fertility; its cap-distorted richer parent lowers fertility after the full redistribution. The joint planner nevertheless raises mean fertility.

Analytical nonvacuity, including young consumption losses

A short all-renter construction shows that the proposition also has a nonempty region with a strict primitive retirement floor and falling young mean consumption. These construction restrictions are not additional assumptions of the proposition.

Choose arbitrary positive heterogeneous incomes and probabilities, write

\[\ell_i=\frac{w_i}{J},\qquad B_0=\bar\ell,\qquad y_0=KB_0,\]

and impose

\[\omega_B>\alpha\Gamma,\qquad
\frac{\beta}{q}<\frac{\bar w}{w_H},\qquad
0<\chi<\frac{\kappa w_H}{(1+\alpha)r}.
\tag{17}\]

Then
\[m_0=qy_0/\delta>\ell_H\]
. Define

\[k=\frac{1+\Gamma}{1+\Gamma\ell_H/m_0}>1.\]

Let
\[p_i^+\]
 solve the uncapped, nonsaving rental demand equation

\[\ell_i\left(\frac{\alpha}{p}+\frac{\vartheta\kappa}{C}\right)=r.\]

For the high type’s capped, nonsaving candidate, let
\[M_H(p)\]
 be its housing MRS, and define
\[p_H^-\]
 by
\[M_H(p_H^-)=kp_H^-\]
.

These thresholds satisfy

\[p_L^+<p_H^+,\qquad p_H^-<p_H^+,\qquad
\frac{\alpha\chi}{\kappa}<p_H^+.\]

The last inequality follows by evaluating high-type uncapped demand at
\[p=\alpha\chi/\kappa\]
, where it equals
\[\kappa w_H/[(1+\alpha)\chi]>r\]
.

Consequently the interval

\[\left(
\max\{p_L^+,p_H^-,\alpha\chi/\kappa\},\ p_H^+
\right)
\tag{18}\]

is nonempty. Throughout it, the low type rents below the ceiling; the high type rents at the ceiling. Both save zero. For the high type, the same support calculation as (15) applies:
\[x_H<k\ell_H<m_0\]
 and
\[M_H/p<k\]
 imply
\[\xi<\Gamma p\zeta\]
.

Mean fertility
\[F(p)\]
 is strictly decreasing on (18). Choosing replacement
\[1/\nu\]
 strictly between its endpoint values gives a unique price
\[p^*\]
 within this regime. At that price,

\[\bar x>B_0,\]

because the capped high type has
\[x_H>\ell_H\]
.

Now choose retirement income from the nonempty interval

\[KB_0<y<K\bar x.
\tag{19}\]

Young choices remain unchanged: higher retirement income strengthens the zero-saving and rental-dominance inequalities. Thus
\[p^*\]
 remains unchanged, while

\[\frac{\bar w}{J}<B=\frac yK<\bar x.\]

The primitive floor holds strictly, the planner increases young housing and joint fertility, and young mean consumption falls. Housing clears at

\[N=\frac{\bar H}
{f h_L+(1-f)r+\alpha y/(Kp^*)}.\]

All inherited old wealth is generated by the preceding young choices: it is exactly
\[y\]
.

The substantive conditions are therefore four: old households can attain their unrestricted housing–estate mix; retirement resources provide the stated aggregate spending cushion; some incomes fall in the derived affected interval; and children are sufficiently housing-intensive for the fertility result. The resource floor is conservative and is not uniformly weaker than the previous regime-specific bound.

The efficiency gap identifies distorted households. Equal-weight redistribution need not identify those same households as its recipients. Keeping that distinction explicit is what makes this shorter theorem both useful and honest.
