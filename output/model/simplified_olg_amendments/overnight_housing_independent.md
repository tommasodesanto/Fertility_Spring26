# Independent housing check

September 9, 2026 UTC. Read-only model review; no live Pro output consulted. The dated planner and competitive household model remain unchanged. This note establishes a cap-valid primitive theorem and an exact competitive-equilibrium counterexample to the weaker capped extension.

## 1. A primitive theorem with all physical caps retained

Consider a zero-property-tax reference equilibrium. This is an explicit parameter subcase, not an amendment. Write current purchase resources as \(w_i=y_i^y+b_i>0\) and old income as \(v_i=y_i^o\). Define
\[
a=\alpha+\vartheta,\quad E=1+a,\quad K=1+\gamma+\omega_B,\quad
m=\min\{1,(1-\phi)/(1-q)\},
\]
\[
\Gamma=\min\{\gamma,(\gamma+\omega_B)(1-q)\},\qquad
\delta=\frac{a}{qE}(1/m-1).
\]
The coefficient \(\Gamma\) allows either old-owner estate regime. The number \(m\) compares minimum current housing cash cost with its service cost.

**Proposition.** Assume \(\alpha\ge\gamma\), a positive stationary competitive equilibrium, and positive unused aggregate housing capacity, \(\bar H<2N\int H_d\,dQ\). The last condition concerns the reference equilibrium; Section 2 supplies a primitive certificate. If
\[
\boxed{\quad
\frac{v_i}{w_i}>
\delta+K\max\left\{\frac{a}{Em\Gamma},\,\frac{\beta}{q}\right\}
\quad\text{for every endowment type},\quad} \tag{P}
\]
then every young financing multiplier is positive, matched old households occupy weakly more total housing than young households, and the full dated consumption–housing planner strictly increases aggregate young housing. Competitive and planner caps may bind. Neither side of \(\beta/q=1\) is imposed.

This is a strong future-income/current-liquidity restriction, not a claim of empirical plausibility. The first term inside the maximum suffices for the allocation direction; the second certifies strictly restrictive finance. For the latter conclusion it is enough that its income inequality hold on a positive-mass group.

**Proof.** Put \(p=(1-q)P\), \(L_R=p\), \(L_O=(1-\phi)P\), \(x=c-\chi n\), and \(s=h-\kappa n\). The young first-order conditions imply, with \(\rho=\alpha x/s\),
\[
\rho\ge mp,\qquad \vartheta x=n(\chi+\kappa\rho),\qquad
ax=\rho h+\chi n.
\]
Thus \(ac=\rho h+E\chi n>mph\). Since \(c+L_dh\le w\) and \(L_d\ge mp\),
\[
ph<\frac{aw}{Em}. \tag{1}
\]
The lifetime budget and financing inequality give
\[
qz=w+qv-c-ph\ge qv-(1-m)ph,\qquad z\ge v-\delta w. \tag{2}
\]
For old housing, strict concavity and the two unrestricted-cap estate solutions imply
\[
h_d^o=\min\{H_d,\Gamma_dz/(Kp)\},\qquad
\Gamma_R=\gamma,\quad\Gamma_O=\Gamma. \tag{3}
\]
Clipping is legitimate because the old objective, after optimizing goods and estates conditional on housing, is concave. Its Euler identity also gives \(c^o\ge z/K\), including binding physical caps.

If old housing is uncapped, (P), (1)–(3) imply
\[
h_d^o>\frac{aw}{Emp}>h_d^y.
\]
If old housing is capped, \(h_d^o=H_d\ge h_d^y\). If a young financing multiplier were zero, its goods condition would imply
\[
x=q c^o/\beta\ge qz/(\beta K)>w,
\]
contradicting current cash feasibility. Hence finance is strictly restrictive.

The fixed-fertility planner has
\[
h_i^{y,F}=\min\{H_{d_i},\kappa n_i+\alpha/\lambda\},\quad
h_i^{o,F}=\min\{H_{d_i},\gamma/\lambda\}.
\]
Consequently \(h_i^{y,F}\ge h_i^{o,F}\), strictly whenever old housing is below its cap. Positive unused aggregate capacity excludes all old being capped. Therefore
\[
H_y^F>\bar H/2\ge H_y^{eq}.
\]
Consumption is optimized too; separability makes the housing solution independent of that optimization.

**Individual and local meanings.** Any young household below its retained cap with \(s_i\le\bar s\) gains housing: market old-total ordering implies \(\alpha/\lambda>\bar s\). This set need not include every constrained household or have positive mass without further restrictions. Every uncapped young household has a strictly beneficial paired housing transfer before adjustment costs, since \(h_i^o\ge h_i^y\) and \(\alpha\ge\gamma\). These statements are distinct from the aggregate result.

## 2. Analytical compatibility, including binding planner caps

The conditions are jointly nonempty for every finite \(\beta>0\). Choose a bounded, nondegenerate distribution of \(w\), with lower bound \(\underline w>0\). Choose \(\chi>0\) small and \(H_R>\kappa/\nu\) large enough that
\[
\vartheta\nu>\chi(E+\beta K)/\underline w+
\frac{\alpha\kappa}{H_R-\kappa/\nu}.
\tag{4}
\]
At sufficiently low prices all conditional young choices have housing at least \(H_R\); the young Euler bound \(x\ge w/(E+\beta K)\) then puts fertility above replacement. At high prices cash feasibility puts fertility below replacement. Continuity supplies a positive stationary price; housing clearing determines \(N\). This establishes existence, not global uniqueness.

Choose bounded old incomes high enough to satisfy (P). No upper bound on old income conflicts with (4). Fully capped young choices have fertility above replacement under (4), so some young choices must be uncapped in equilibrium; aggregate capacity is therefore not exhausted.

Binding caps can be economically active in this construction. From (1),
\[
p<p_{\max}:=\nu a\bar w/(Em\kappa).
\]
Impose also \(v_i-\delta w_i>Kp_{\max}H_O/\Gamma\). Every competitive old household then occupies its tenure cap. Choose \(\phi\ge q\), \(\omega_B(1-q)>q\gamma\), ownership-taste location zero, and
\[
H_O>3H_R+2\kappa/\nu.
\]
An owner can attain the renter optimum, so its ownership probability is at least one half; logistic tastes still give both tenures positive probability. Hence mean old capped housing exceeds \(2H_R+\kappa/\nu\). If the planner had \(\alpha/\lambda\le H_R\), its total housing per cohort mass would be at most \(2H_R+\kappa/\nu\), a contradiction. Thus all planner young renters are capped. This is an analytical parameter construction, not a numerical point plus continuity.

## 3. Exact counterexample to the weaker capped extension

Positive financing multipliers, \(\beta=q\), positive old financial estates, and paired young marginal housing utility exceeding old marginal utility do **not** suffice once planner caps bind.

Set
\[
q=\phi=\beta=\tfrac12,\quad
\alpha=\gamma=\vartheta=\chi=\kappa=1,\quad\omega_B=2,\quad
p=1,\ P=2,\quad H_R=\tfrac{29}{20},\ H_O=\tfrac{17}{4}.
\]
There are two endowment types, each of mass one in each age:
\[
(w_L,v_L)=(27/10,4),\qquad (w_H,v_H)=(15/2,52/5).
\]
For either low-type tenure, and for high-type ownership, the competitive choices are:

| Type/tenure | \(x\) | \(n\) | \(h^y\) | \(c^o=h^o\) | Financing multiplier |
|---|---:|---:|---:|---:|---:|
| Low, either | \(9/10\) | \(9/20\) | \(27/20\) | \(1\) | \(1/9\) |
| High, owner | \(5/2\) | \(5/4\) | \(15/4\) | \(13/5\) | \(1/65\) |

Old owner estates are \(e=4c^o>Ph^o\); all these caps are slack.

The high renter instead has both housing caps binding and finance slack. Its fertility and adult consumption are
\[
n_R=\frac{1045-\sqrt{652501}}{360}\in(13/20,7/10),\qquad
x_R=\frac{421}{100}-\frac25 n_R.
\]
Its old consumption is \(x_R\), and its estate is \(4x_R\). These solve
\[
\tfrac52x_R+n_R=\tfrac{421}{40},\qquad
1/n_R=1/x_R+1/(H_R-n_R).
\]
Cash spending is \(x_R+n_R+H_R<15/2\); both cap multipliers are positive. Thus every displayed choice satisfies the original household optimality conditions.

The high type's ownership value advantage is the explicit constant
\[
\begin{aligned}
\Delta={}&2\log(5/2)+\log(5/4)+\log(13/5)+\log(52/5)\\
&-\{(3/2)\log x_R+\log(H_R-n_R)+\log n_R
+(1/2)\log H_R+\log(4x_R)\}.
\end{aligned}
\]
It is strictly positive: the renter optimum is feasible for ownership, while the unique owner optimum violates the rental caps. The low type's advantage is zero. For any \(\varepsilon\in(0,1/2)\), set
\[
\bar\xi=-\Delta/2,\qquad
\sigma_\xi=\frac{\Delta}{2\log((1-\varepsilon)/\varepsilon)}.
\]
The endogenous owner probabilities are exactly \(\varepsilon\) for the low type and \(1-\varepsilon\) for the high type. Both tenures occur at both endowments; each tenure's total mass is one.

Choose
\[
\nu=\frac{2}{9/20+(1-\varepsilon)5/4+\varepsilon n_R},\qquad
\bar H=\frac{87}{10}-\frac{69}{20}\varepsilon .
\]
Replacement fertility and housing clearing hold with \(N=2\), giving an exact positive stationary equilibrium. Current earnings and entry wealth can each be half of the specified \(w\).

Competitive young housing is
\[
H_y^{eq}=51/10-(23/10)\varepsilon .
\]
The full planner caps both renter ages, leaves both owner ages uncapped, and has
\[
1/\lambda=91/40-(53/40)\varepsilon.
\]
These regimes hold throughout \(0<\varepsilon<1/2\). Planner young housing is
\[
H_y^F=199/40-(17/8)\varepsilon,
\qquad
\boxed{H_y^F-H_y^{eq}=(-5+7\varepsilon)/40<0.}
\]
Every paired young marginal housing utility exceeds its old counterpart's. Finance is strictly restrictive except for high renters. The decline is therefore a genuine competitive-equilibrium cap obstruction, not merely an arbitrary feasible allocation.

## 4. The weaker adult-space result and next decision

Replacing \(a\) by \(\gamma\) only in the housing part of (P) gives
\[
v_i/w_i>\delta+K\gamma/(Em\Gamma).
\]
This derives \(h_i^o>(\gamma/\alpha)s_i\) when old housing is uncapped; capped old also satisfy it if \(\alpha\ge\gamma\). If the planner is uncapped,
\[
H_y^F-H_y^{eq}=\frac{N}{\alpha+\gamma}
(\alpha\bar h^o-\gamma\bar s)>0.
\]
Section 3 shows why this reasoning cannot simply retain binding planner caps.

A sharper distributional result is available in the explicit \(\phi=q,\tau=0\), competitive-cap-slack benchmark. Put
\[
\bar x=\int\min\{w/E,(w+qv)/(E+\beta K)\}\,dF.
\]
Tenure probabilities are constant across endowments; let
\(\bar\Gamma=(1-\pi)\gamma+\pi\Gamma\). With an uncapped planner, the exact primitive housing criterion is
\[
\overline{w+qv}>(E+qK\gamma/\bar\Gamma)\bar x.
\]
It measures whether current borrowing shortfalls offset patience and estate effects, without imposing either side of \(\beta/q=1\).

**Proposed follow-up:** verify the exact capped counterexample above, then seek one economically meaningful restriction on the joint distribution of current cash, future income, and retained caps that weakens (P) while excluding its type/tenure sorting mechanism. Do not attempt another proof based only on paired marginal-utility ordering. Empirical plausibility of the strong income floor and the funded tax transition remain unassessed.
