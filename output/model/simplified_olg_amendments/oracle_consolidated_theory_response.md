# Pro consolidated review

Captured from the completed rendered response on September 9, 2026. UI: 6 Pro; worked for 72m 6s.



**The housing and joint-fertility results can be completed, including binding planner caps. There is also a genuine local tax-transition result—but the tax does not generally implement the dated planner or give young households larger homes.**

**Theory section and complete proof appendix — PDF**
**Editable, self-contained LaTeX source**

The proposed paper text is on **pages 2–7**. The full 28-page document includes the proofs, analytical counterexamples, existence audit, separate general-preferences investigation, and six-slide outline.

## 1. The principal conclusions

**Housing.** The relevant comparison is old housing against young **adult space**, not old housing against young total housing. A single primitive condition combining patience, current cash, old income, and estate preferences generates the required comparison in an explicitly solved heterogeneous benchmark. It permits either \(\beta<q\) or \(\beta>q\). A separate construction establishes compatibility with strictly binding rental caps and strictly positive borrowing multipliers.

**Fertility.** The joint dated planner raises average fertility under the two resource comparisons below and a capacity condition substantially weaker than requiring its solution to be uncapped. However, **average private fertility after the fixed-fertility allocation can fall**. I give an analytical counterfamily in which that happens despite strictly binding young finance and an increase in aggregate young housing.

**Transition.** In an explicit local regime, the equilibrium—not an imposed fertility path—has a unique convergent response to a fertility-taste decline and a subsequent unexpected permanent property-tax increase. The tax raises fertility on impact and raises terminal adult population. The proof derives the stability conditions. Its limits are material: competitive physical caps and old financial-estate floors remain slack locally, and no global convergence theorem is claimed.

# 2. Housing: adult space is the useful comparison

I retain the specified one-date criterion: equal weights on current households’ remaining utilities, fixed fertility initially, fixed incumbent continuation opportunities and old net estates, and optimization over **all current consumption and housing**. This is not a compensated-Pareto comparison or a stationary lifetime comparison that leaves out the initial old.
Pasted text

Write reference adult consumption and space as

\(x_i=c_i^{y,\mathrm{eq}}-\chi n_i,
\qquad
s_i=h_i^{y,\mathrm{eq}}-\kappa n_i.\)

Bars denote averages under the matched stationary distribution \(Q\).

The full fixed-fertility planner has the solution

\(x^F=c^{o,F}
=\frac{\bar x+\bar c^{o,\mathrm{eq}}}{2},\)
\(h_i^{y,F}
=\min\left\{H_{d_i},\kappa n_i+\frac{\alpha}{\lambda_F}\right\},
\qquad
h_i^{o,F}
=\min\left\{H_{d_i},\frac{\gamma}{\lambda_F}\right\},\)

where \(\lambda_F\) clears housing.

## The uncapped comparison is exact

When planner housing caps are slack,

\(\boxed{
H_Y^F-H_Y^{\mathrm{eq}}
=
\frac{N}{\alpha+\gamma}
\left(\alpha\bar h^{o,\mathrm{eq}}-\gamma\bar s\right).
}\)

Thus the necessary and sufficient comparison is

\(\boxed{
\bar h^{o,\mathrm{eq}}>
\frac{\gamma}{\alpha}\bar s.
}\)

There is no need to establish
\(\bar h^{o,\mathrm{eq}}>\bar h^{y,\mathrm{eq}}\).
Children already account for part of young housing.

## A cap-valid extension

The same aggregate direction holds if

\(\bar h^{o,\mathrm{eq}}>
\frac{\gamma}{\alpha}\bar s,\)

and

\(\boxed{
H_{d_i}-\kappa n_i\ge\bar s
\quad Q\text{-almost everywhere},
\qquad
Q\{H_{d_i}-\kappa n_i>\bar s\}>0.
}
\tag{CF}\)

This says that every retained housing maximum can accommodate that household’s fixed child space plus the reference mean adult space, with some strict room. **Planner caps may nevertheless bind.**

The proof is short. Put \(t=\alpha/\lambda_F\). At \(t=\bar s\), condition (CF) makes aggregate young adult space exactly \(\bar s\), while desired old housing is at most \((\gamma/\alpha)\bar s\). Together these demand less housing than the reference allocation contains. Housing clearing therefore requires \(t_F>\bar s\). The positive mass with additional capacity then makes aggregate young housing strictly larger.

The individual comparison is also exact:

\(h_i^{y,F}-h_i^{y,\mathrm{eq}}
=
\min\left\{
H_{d_i}-h_i^{y,\mathrm{eq}},
\frac{\alpha}{\lambda_F}-s_i
\right\}.\)

Consequently, every reference-uncapped young household with \(s_i\le\bar s\) gains housing. A fully capped young household cannot receive strictly more under the retained menu.

Without a capacity restriction, the comparison between the two means is insufficient. The appendix gives both an exact capped allocation test and a counterexample showing how a small residual rental capacity can redirect redistribution toward the old.

# 3. Generating the housing gap from equilibrium primitives

The maintained budgets reduce at stationarity to

\(c+ph+qz=w+qv,\qquad c+L_dh\le w,\)

where

\(w=y^y+b+T,\qquad v=y^o+T,\)
\(p=(1-q+q\tau^p)P,\qquad
L_R=p,\qquad L_O=(1-\phi+q\tau^p)P.\)

Current income is included in \(w\); old owners retain the ability to sell and resize.
Pasted text

Define

\(E=1+\alpha+\vartheta,\qquad
K=1+\gamma+\omega_B,\qquad
D=E+\beta K.\)

When old housing is uncapped, the two estate regimes give

\(c_i^o=\frac{z_i}{K},
\qquad
ph_i^{o,d}=\Gamma_dc_i^o,\)

with

\(\Gamma_R=\gamma,\qquad
\Gamma_O=
\min\left\{
\gamma,\,
\frac{(\gamma+\omega_B)(1-q+q\tau^p)}
{1+q\tau^p}
\right\}.\)

The old financial-estate floor therefore matters even though the old can sell their houses. When it binds, it lowers old housing demand relative to the unrestricted expenditure-share formula.

## A single combined restriction

Consider the explicitly solved benchmark

\(\tau^p=0,\qquad \phi=q,\)

and verify that its computed competitive choices are below their finite tenure caps. Define directly from endowments

\(x_i=\min\left\{
\frac{w_i}{E},
\frac{w_i+qv_i}{D}
\right\},
\qquad
c_i^o=\frac{w_i+qv_i-Ex_i}{qK}.\)

Replacement fertility determines the service price:

\(p=\frac{\nu\vartheta\bar x-\chi}{\kappa}>0.\)

Conditional young choices coincide across tenures, and the old tenure-value difference is constant across endowments. Hence the logistic choice rule gives a constant owner share \(\pi\in(0,1)\), calculated explicitly in the appendix. Put

\(\bar\Gamma=(1-\pi)\gamma+\pi\Gamma_O.\)

The useful primitive condition is

\(\boxed{
\frac{
\displaystyle\int x_i
\max\left\{\frac{\beta}{q},\frac{Ev_i}{Kw_i}\right\}\,dF
}{\bar x}
>
\frac{\gamma}{\bar\Gamma}.
}
\tag{P}\)

Also require

\(F\{qEv_i>\beta Kw_i\}>0\)

to identify a positive mass for whom finance is strictly restrictive.

The derivation is exact:

\(\boxed{
\frac{c_i^o}{x_i}
=
\max\left\{\frac{\beta}{q},\frac{Ev_i}{Kw_i}\right\},
\qquad
\mu_i>0
\Longleftrightarrow
qEv_i>\beta Kw_i.
}\)

Moreover,

\(\bar s=\frac{\alpha\bar x}{p},
\qquad
\bar h^{o,\mathrm{eq}}
=\frac{\bar\Gamma\,\bar c^o}{p}.\)

Thus (P) implies both

\(\bar c^o>\bar x,
\qquad
\bar h^o>\frac{\gamma}{\alpha}\bar s.\)

The competitive cap checks also imply (CF) in this benchmark, without requiring planner caps to be slack.

**Economic interpretation.** Old resources relative to current adult consumption are high because the household is patient, because old income is high relative to current cash, or both. A binding old-estate floor raises the threshold through \(\gamma/\bar\Gamma\). This is why neither side of \(\beta/q=1\) needs to be imposed separately.

An analytical nonempty family is immediate. Take genuinely heterogeneous bounded \(w>0\), set \(v=aw\), and choose

\(a>
\max\left\{
\frac{\beta K}{qE},
\frac{\gamma K}{\bar\Gamma E}
\right\},
\qquad
0<\chi<\frac{\nu\vartheta\bar w}{E}.\)

Choose finite caps above the computed competitive demands. Every young household then strictly wants additional borrowing.

A positive-mass group that is **both constrained and an individual recipient** is identified by

\(\boxed{
F\{w_i<E\bar x,\ qEv_i>\beta Kw_i\}>0.
}\)

These households gain consumption as well as housing.

The inactive competitive rental cap is a limitation of this particularly simple closed-form benchmark. Appendix C gives a different exact construction with heterogeneous cash endowments, strictly constrained young households, **strictly capped high-cash renters**, and uncapped young owners. It satisfies the resource and planner-cap conditions without a numerical reference point or continuity certificate.

## Why young borrowing constraints alone are insufficient

There is a clean stationary counterfamily. Set

\(\tau^p=0,\quad \phi=q,\quad \beta=q,
\quad \omega_B(1-q)<q\gamma,\)

so \(\bar\Gamma<\gamma\). Let \(w\) be nondegenerate and \(v=aw\), with

\(\frac KE<a<\frac{\gamma K}{E\bar\Gamma}.\)

Every young household has a strictly positive borrowing multiplier, but

\(H_Y^F-H_Y^{\mathrm{eq}}
=
\frac{N\alpha\bar w}{p(\alpha+\gamma)}
\left(\frac{\bar\Gamma a}{K}-\frac{\gamma}{E}\right)
<0.\)

Relaxing the old financial-estate restriction can benefit old housing demand enough to reverse the intended direction. This is not mechanical attachment to an inherited home.

The welfare interpretation remains utilitarian. At \(\beta=q\), with slack old-estate and physical constraints, the paired marginal housing gap reduces to the young financing wedge. With all finance slack, heterogeneous redistribution can improve this welfare criterion while leaving aggregate young housing unchanged.

# 4. Fertility: private adjustment and the joint optimum differ

## Private fertility at a changed bundle

At a given gross bundle, private fertility satisfies

\(\frac{\vartheta}{n}
=
\frac{\chi}{c-\chi n}
+
\frac{\alpha\kappa}{h-\kappa n}.\)

Its derivative is

\(\boxed{
dn=
\frac{\chi\,dc/x^2+\alpha\kappa\,dh/s^2}
{\vartheta/n^2+\chi^2/x^2+\alpha\kappa^2/s^2}.
}\)

More housing therefore raises fertility **holding consumption fixed**. For a housing gain \(\Delta h\ge0\) accompanied by a consumption loss \(\delta\ge0\), the exact finite-change condition is

\(\boxed{
n_1\ge n_0
\Longleftrightarrow
\delta\le
\frac{\alpha\kappa x^2\Delta h}
{\chi s(s+\Delta h)+\alpha\kappa x\Delta h},
}\)

provided the new bundle accommodates the initial fertility.

At the **actual fixed-fertility planner allocation**, the appropriate test is

\(\boxed{
n_i^S>n_i
\Longleftrightarrow
\frac{\chi}{x^F}
+\frac{\alpha\kappa}{s_i^F}
<
\frac{\vartheta}{n_i},
\qquad
s_i^F=h_i^{y,F}-\kappa n_i.
}\)

The constrained, low-resource recipients identified above satisfy this test. Other households may not.

### Average sequential fertility can fall

This is more than a missing proof. The appendix establishes a counterfamily with

\(\alpha=\gamma,\quad \beta=q,\quad \phi=q,\quad \tau^p=0,\)

a slack old-estate floor, heterogeneous \(w_i\), and

\(v_i=\frac KE(1+\delta)w_i,\qquad \delta>0.\)

Every young household is strictly constrained. The fixed-fertility planner increases aggregate young housing, and the joint planner increases average fertility, yet private fertility after the fixed-fertility assignment can fall on average.

The argument uses the concavity of the private fertility function \(\mathfrak n(c,h)\), not a numerical example. With \(X=\bar x\), \(S=\bar s\), \(m=\bar n\), define

\(g(n)=\mathfrak n(X+\chi n,S+\kappa n).\)

When \(p\kappa\ne\alpha\chi\), this function is strictly concave along the displayed bundle line. Thus

\(G:=\int g(n_i)\,dF<g(m)=m.\)

For the explicit nonempty interval

\(0<\delta<2(m/G-1),\)

the sequential allocation satisfies \(\bar n^S<m\).

So the paper should not substitute “reallocate at fixed fertility and then let households adjust” for the joint-planner theorem.

## A joint-fertility theorem with binding caps

Here is the stronger result.

Suppose

\(\bar c^{o,\mathrm{eq}}\ge\bar x,
\qquad
\bar h^{o,\mathrm{eq}}\ge\frac{\gamma}{\alpha}\bar s,\)

with at least one inequality strict, and

\(\boxed{H_R>\bar h^{y,\mathrm{eq}}.}
\tag{CJ}\)

Then the joint dated planner satisfies

\(\boxed{\bar n^J>\bar n,\qquad H_Y^J>H_Y^{\mathrm{eq}}.}\)

**Its housing caps need not be slack.**

The capacity requirement concerns the **reference mean young home**, not the largest home the joint planner would choose.

### Proof

The joint optimum has common adult consumption \(x^J\), and within each tenure \(d\),

\(n_d^J=\frac{\vartheta}{\chi/x^J+\kappa r_d},
\qquad
s_d^J=\frac{\alpha}{r_d},
\qquad
r_d=\lambda_H+\eta_d.\)

First, multiply the reference fertility condition by **\(n_i^2\)** and integrate:

\(\vartheta\bar n
=
\chi\int\frac{n_i^2}{x_i}\,dQ
+\alpha\kappa\int\frac{n_i^2}{s_i}\,dQ.\)

Cauchy–Schwarz gives

\(\boxed{
\frac{\vartheta}{\bar n}
\ge
\frac{\chi}{\bar x}
+\frac{\alpha\kappa}{\bar s}.
}\)

Suppose \(\bar n^J\le\bar n\). Goods clearing then gives \(x^J\ge\bar x\). Some tenure must have \(n_d^J\le\bar n\); its fertility condition and the preceding inequality imply \(s_d^J\le\bar s\). Hence

\(h_d^{y,J}\le\bar s+\kappa\bar n
=\bar h^{y,\mathrm{eq}}<H_R,\)

so this tenure is uncapped.

All uncapped tenures have the same fertility, and capped tenures have lower fertility because their \(r_d\) is higher. Therefore every tenure has \(n_d^J\le\bar n\). Repeating the argument makes **every young tenure uncapped**, with adult space no greater than \(\bar s\).

Old mean housing is then at most \((\gamma/\alpha)\bar s\). This contradicts housing clearing when the old-housing resource comparison is strict. If only the consumption comparison is strict, \(x^J>\bar x\) forces strictly smaller young adult space and yields the same contradiction. Thus \(\bar n^J>\bar n\).

For housing, if any young cap binds, all joint-planner young homes are at least \(H_R>\bar h^{y,\mathrm{eq}}\). Otherwise, common adult space and housing clearing imply the aggregate gain directly.

This result retains heterogeneous reference households throughout. It does not replace them by their mean bundle.

# 5. A funded tax transition, rather than a chosen fertility path

The policy exercise uses the specified permanent property-tax increase, equal rebates to young and old decision units, and identical inherited states at intervention. It is unexpected at \(t_p\) and known thereafter. Private borrowing restrictions remain.
Pasted text

## A proved local regime

The local theorem starts from a zero-tax stationary benchmark with \(\phi=q\), strictly slack competitive housing caps, and a slack old financial-estate floor. Define

\(W=\int w_i\,dF,\qquad V=\int v_i\,dF,\)

and require

\(qEv_i>\beta Kw_i\)

uniformly on the compact endowment support. This **derives** strictly binding young finance.

From the explicit stationary solution, define

\(\epsilon=\frac{p\kappa}{\chi+p\kappa},
\qquad
A=\alpha+\vartheta\epsilon,
\qquad
\zeta=\frac{\bar h^o}{\bar h^y}
=\frac{\gamma EV}{KWA},
\qquad
\pi=\operatorname{logit}^{-1}(\bar\xi/\sigma_\xi).\)

One analytical sufficient region is

\(\boxed{
\begin{gathered}
V\ge W,\qquad
\frac{\alpha}{1+\alpha}<\epsilon<q,\\
1<\zeta<\frac{2E\epsilon}{A}-1,\\
\omega_B(1-q)>q\gamma,\qquad
\frac{\pi\gamma}{K}
<
\min\{(1-q)^2,q\zeta\}.
\end{gathered}
}
\tag{T}\)

These restrictions have identifiable roles. The lower bound on \(\epsilon\) is \(p\kappa>\alpha\chi\): child space costs are sufficiently important relative to goods costs. The upper bound on \(\zeta\) ensures tax capitalization lowers the terminal purchase price. The last restriction bounds the feedback from old housing exposure.

Under (T):

A sufficiently small exogenous fertility-taste decline lowers fertility initially and starts a convergent transition to a lower stationary adult population.

An unexpected sufficiently small permanent tax increase introduced along that same local transition raises fertility relative to baseline on impact.

The policy converges to a higher terminal cohort mass \(N_1>N_0\), although both terminal fertility rates equal \(1/\nu\).

The appendix constructs a nonempty primitive family satisfying all these restrictions. Its admissible interval

\(0<\beta<\frac{qA\zeta}{\gamma}\)

contains values on both sides of \(q\) when \(\alpha\ge\gamma\).

### What establishes convergence?

In this regime, the exact nonlinear equilibrium has inherited state

\((Y_t,O_t,P_{t-1},B_{t-1}),\)

where \(B_{t-1}\) is the preceding cohort’s mean owner titles, and \(P_t\) is the jump variable.

The proof derives the household demands, rebate consistency, and housing-clearing derivatives from the budgets. Its characteristic cubic has **one root above one and two roots strictly inside the unit circle**; the full state system adds two zero stable roots. The stable manifold projects invertibly onto inherited states. This establishes a unique local bounded equilibrium path and geometric convergence.

Thus stability is proved in the stated parameter region, not assumed. The result does not cover a large transition that crosses physical-cap or estate regimes.

## The terminal effect is explicit

Write \(\tau=\tau^p\), and define

\(B_0=\frac{(\alpha+\vartheta)W}{E}
+\frac{\gamma V}{K}-\frac{\chi}{\nu},
\qquad
b_0=\frac{\alpha+\vartheta}{E}+\frac{\gamma}{K}.\)

The stationary tax equilibrium in this regime is

\(T(\tau)=
\frac{q\tau B_0}{2(1-q)+q\tau(2-b_0)},\)
\(p(\tau)=
\frac{\nu\vartheta(W+T)/E-\chi}{\kappa},
\qquad
N(\tau)=
\frac{\bar H\,p(\tau)}{B_0+b_0T(\tau)}.\)

Since \(T'(\tau)>0\),

\(\boxed{
\frac{\partial N}{\partial T}
=
\frac{\bar H}{\kappa(B_0+b_0T)^2}
\left[
\frac{\nu\vartheta\gamma}{EK}(V-W)
+\chi\left(\frac{\alpha}{E}+\frac{\gamma}{K}\right)
\right]>0
}\)

when \(V\ge W\).

This is not a universal tax sign. For example, take \(v_i=aw_i\) with

\(\frac{\beta K}{qE}<a<1.\)

Young finance is still strictly binding, but sufficiently small positive \(\chi\), specifically

\(\chi<
\frac{\nu\vartheta\gamma W(1-a)}
{\alpha K+\gamma E},\)

makes the stationary population derivative negative.

## The important housing qualification

In the positive tax-result regime, stationary **aggregate young housing increases**, but the mean home per young household satisfies

\(\bar h^y
=
\frac{\kappa}{\nu}
+\frac{\alpha(W+T)}{Ep},\)

and therefore

\(\boxed{
\frac{d\bar h^y}{dT}
=
-\frac{\alpha\chi}{E\kappa p^2}<0.
}\)

More young households occupy a larger aggregate share of the stock, but a smaller home on average.

The tax also does not have an unconditional impact home-size sign under the theorem’s assumptions. The appendix gives its exact derivative separately. The positive fertility response therefore must not be described as necessarily operating through larger homes.

## Population and welfare

With common initial young mass,

\(\frac{Y_T^P}{Y_T^B}
=
\prod_{t=t_p}^{T-1}
\frac{\bar n_t^P}{\bar n_t^B}.\)

The proved local convergence gives

\(\boxed{
\sum_{t=t_p}^{\infty}
\log\frac{\bar n_t^P}{\bar n_t^B}
=
\log\frac{N_1}{N_0}>0.
}\)

Fertility need not be higher at every subsequent date. Terminal adult population is \(2N\).

A temporary reform returning to the same primitives and the same unique local terminal equilibrium has no permanent population-level effect: its cumulative fertility gap must eventually be offset.

For intervention welfare, Appendix H includes current young continuation-price effects, initial-old title losses, and the estate-floor terms. It also derives a nonempty sufficient condition for a local tax welfare improvement, based on an explicit marginal-utility-of-cash moment. **Without that additional sign condition, higher fertility or higher terminal population is not a welfare theorem.**

# 6. General preferences: what survives without shifts

The separate investigation begins with gross utility

\(U^y(c,h,n),\qquad U^o(c^o,h^o,e),\)

rather than assuming linear child-needs offsets. This follows the requested distinction between the maintained model and the exploratory branch.
Pasted text

The competitive wedge survives:

\(\boxed{
U_h^y-U_h^o
=
(\beta/q-1)pm
+L_d\mu+\eta^y-P\rho-\eta^o,
}\)

where \(m=U_c^o\) and \(\rho\) is the old-owner estate-floor multiplier. No logarithm is required.

At an attained interior fertility optimum,

\(\boxed{
dn=-\frac{U_{nc}^y\,dc+U_{nh}^y\,dh}{U_{nn}^y}.
}\)

Positive resource effects require additional cross-partial restrictions. Concavity alone does not supply them.

For example,

\(U^y(c,h,n)=\sqrt c+\sqrt{h+n}+A_n\sqrt n-kn\)

is strictly concave and increasing in goods and housing, with a unique interior fertility optimum. But \(U_{nh}^y<0\), so more housing lowers fertility.

For the intermediate class

\(U^y=a(c,n)+b(h,n)+v(n),\)

positive \(a_{cn}\) and \(b_{hn}\) give conditional complementarity. A common cardinal housing component for old households, \(b(h^o,0)\), yields a useful planner age ordering. The appendix supplies additional restrictions that turn it into an equilibrium-to-planner theorem, together with a nonshift example satisfying those restrictions.

The mapping back is precise: the maintained shifts deliver positive resource complementarity; equal adult-space allocation does not itself require logarithms when age housing weights match. However, the explicit cash-income shares, constant adult-space ratio with unequal weights, and the particular heterogeneous fertility aggregation argument use additional functional structure.

# 7. Checks and the paper’s remaining choices

The supplied stationary existence certificate is valid, including its rebate bound and low- and high-price arguments. Its low-price fertility restriction contains \(D=E+\beta K\), so it does impose a finite upper bound on \(\beta\) at fixed other primitives. That is a feature of this **sufficient certificate**, not a necessary finite-life existence restriction. The appendix checks compatibility with the new constructions and distinguishes conditional household uniqueness, planner uniqueness, within-regime stationary uniqueness, and local transition uniqueness.
Pasted text

The literature supports the ingredients, not the entire theorem. Coven et al.’s checked August 1, 2026 version analyzes capitalization and intergenerational housing allocation; its two-period benchmark does not have this equal-cash-rebate parental-fertility problem. Van Doornik et al.’s December 2024 model incorporates fertility, education, and space requirements, but not the maintained mortgage and persistent-tenure architecture. Neither is used as authority for the new welfare or transition propositions. [Banco Central do Brasil+3Abdou Econ+3Abdou Econ+3](https://abdouecon.github.io/research/papers/Property_Tax.pdf)

The six-slide outline follows **environment → equilibrium and resource gap → housing allocation → fertility → funded reform → demographic transition**. Its housing illustration uses marginal utilities consistent with the dated objective. Its two-panel transition illustration marks the proved impact and endpoint movements while leaving the intervening shape schematic; it does not pass off chosen curves as a solved path.

The author’s substantive choice is whether the restricted \(\phi=q\), locally uncapped tax regime is acceptable as the market-policy illustration. The dated housing and joint-fertility argument does not require a new model amendment. The stronger claim that the tax gives young households larger homes, improves living-household welfare without an additional condition, and does so along a large transition with binding rental caps remains outside what is proved.

**The clean paper claim is therefore:** financial and estate restrictions can generate too little young adult space under the specified dated welfare criterion; reallocating resources jointly can increase parental fertility; and a funded housing tax can alter the demographic destination through an equilibrium transition—but these are distinct propositions, not automatic steps of one unrestricted welfare theorem.
