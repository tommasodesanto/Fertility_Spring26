# Pro response — housing restrictions and fertility through parents

Source: https://chatgpt.com/c/6aa06dde-50d0-83ea-ab19-bf66d86b0354

Third response, following attachment Pasted text(20260908-235319).txt. Completed in 22 minutes 19 seconds. Captured September 8 Eastern (September 9 UTC), from visible text and mathematical source labels. Display delimiters reconstructed for readability. This is an external review, not an adopted revision to the model, note or slides.

# Housing allocation and fertility valued through parents

**The half-stock lemma is correct. A savings-based restriction can generate the required market housing ordering without the earlier lower-price certificate. Fertility, however, requires a separate argument: neither more housing nor higher parental welfare alone establishes its direction.**

Throughout, retain the specified one-date benchmark: existing adults, retained tenure, total current goods and housing, existing obligations, incumbent continuation opportunities, and old estate payments. The planner relaxes both young financing and old nonnegative financial saving; the competitive economy does not. Pasted text

Write \(H_d=h_d^{\max}\), \(d\in\{R,O\}\), and normalize resources by the reference cohort mass:

\[
\bar c=\frac{C^{eq}}N,\qquad
\bar h=\frac{\bar H}N.
\]
These are totals across both ages **per cohort mass**, not per living household. Reference young quantities are

\[
x_i=c_i^{y,eq}-\chi n_i,\qquad
s_i=h_i^{y,eq}-\kappa n_i,\qquad
\bar n=\frac1\nu.
\]
Superscripts \(F,S,J\) below denote the fixed-fertility planner, subsequent conditional fertility choice, and joint planner.

## 1. Housing: a theorem allowing binding physical caps

Define

\[
\begin{gathered}
E=1+\alpha+\vartheta,\qquad K=1+\gamma+\omega_B,\\
d_p=1-q+q\tau^p,\qquad a_p=1+q\tau^p,\\
\Gamma_R=\gamma,\qquad
\Gamma_O=\min\left\{\gamma,\frac{(\gamma+\omega_B)d_p}{a_p}\right\},
\qquad \Gamma=\Gamma_O.
\end{gathered}
\]
A transparent sufficient preference restriction is

\[
\boxed{\alpha\ge\gamma,\qquad
       \frac{\beta\Gamma}{q}\ge\alpha+\vartheta.}
\tag{H}
\]

### Housing proposition

Take any positive stationary competitive equilibrium satisfying (H) and

\[
\bar H<2N\int H_{d_i}\,dQ.
\]
**Individual comparison.** Every young household in

\[
\boxed{\mathcal I=\{i:h_i^{y,eq}<H_{d_i},\ s_i\le\bar s\}}
\tag{1}
\]
receives strictly more housing at the full fixed-fertility optimum. A positive-mass conclusion requires \(Q(\mathcal I)>0\); an explicit primitive recipient group is supplied below.

**Aggregate comparison.** Equilibrium generates

\[
h_i^{o,eq}\ge h_i^{y,eq}
\quad\text{for every paired type,}
\]
and the full optimum satisfies

\[
\boxed{H_y^F>\frac{\bar H}{2}\ge H_y^{eq}.}
\tag{2}
\]
No income ratio, slack old estate restriction, or absence of binding housing caps is assumed in this proposition. It is conditional on a positive equilibrium, not a general existence theorem.

### Proof

The full fixed-fertility solution is

\[
\begin{aligned}
c_i^{y,F}&=\chi n_i+x^F,&c_i^{o,F}&=x^F,
&x^F&=\frac{\bar x+\overline{c^{o,eq}}}{2},\\
h_i^{y,F}&=\min\{H_{d_i},\kappa n_i+\alpha/\lambda_F\},&
h_i^{o,F}&=\min\{H_{d_i},\gamma/\lambda_F\},
\end{aligned}
\tag{3}
\]
where \(\lambda_F>0\) clears housing.

To establish the market ordering, let

\[
p=d_pP,\qquad L_R=p,\qquad L_O=(1-\phi+q\tau^p)P.
\]
Let \(\Lambda_i,\mu_i,\eta_i\) be the lifetime-budget, young-finance, and young-cap multipliers. Household optimality gives

\[
\begin{gathered}
\Lambda_i=\frac{\beta}{q c_i^{o,eq}},\qquad
\frac1{x_i}=\Lambda_i+\mu_i,\\
s_i=\frac{\alpha}{p\Lambda_i+L_{d_i}\mu_i+\eta_i},\\
n_i=\frac{\vartheta}
{\chi(\Lambda_i+\mu_i)+\kappa(p\Lambda_i+L_{d_i}\mu_i+\eta_i)}.
\end{gathered}
\]
Consequently,

\[
h_i^{y,eq}
\le
\frac{\displaystyle\alpha+
          \frac{\vartheta p\kappa}{\chi+p\kappa}}
     {p\Lambda_i}
<
\frac{\alpha+\vartheta}{p\Lambda_i}.
\tag{4}
\]
If old housing is uncapped, the two estate regimes give

\[
h_i^{o,eq}
=\frac{\Gamma_{d_i}c_i^{o,eq}}p
=\frac{\beta\Gamma_{d_i}}{q p\Lambda_i}.
\]
Condition (H) therefore implies the ordering. If old housing is capped, it equals the young household’s retained cap and the ordering follows immediately.

Inside the planner allocation, \(\alpha\ge\gamma\) and \(n_i>0\) imply \(h_i^{y,F}\ge h_i^{o,F}\), strictly wherever old housing is uncapped. If every old household were capped, every young household would also be capped, contradicting the stock restriction. This proves (2).

Finally, put

\[
S=\bar h-\kappa\bar n=\bar s+\overline{h^{o,eq}}.
\]
Market ordering gives \(S\ge2\bar s+\kappa\bar n\), while (3) gives

\[
S\le\frac{\alpha+\gamma}{\lambda_F}.
\]
Hence \(\alpha/\lambda_F>\bar s\). Both arguments of the minimum in (3) exceed reference housing for every \(i\in\mathcal I\). ∎

### Interpretation and sharper restrictions

The old-income and home-sale proceeds are already included in old resources. Nothing in this proof prevents old owners from selling or downsizing. The smaller coefficient \(\Gamma_O\) captures their inability to borrow against the estate **after** selling and choosing their new home. Pasted text

The restriction \(\alpha\ge\gamma\) compares adult-space weights relative to the common consumption-log coefficient. Children weaken the restrictions needed for particular conclusions:

\[
\frac{\alpha}{s_i}>\frac{\gamma}{h_i^{o,eq}}
\quad\Longleftrightarrow\quad
\frac{\alpha}{\gamma}>
\frac{h_i^{y,eq}-\kappa n_i}{h_i^{o,eq}}.
\]
Thus, with \(h_i^{o,eq}\ge h_i^{y,eq}\), the sufficient threshold is
\(\alpha/\gamma>1-\kappa n_i/h_i^{y,eq}\), not necessarily one.

Separately, when the fixed-fertility planner is uncapped,

\[
H_y^F>\bar H/2
\quad\Longleftrightarrow\quad
\boxed{\frac{\alpha}{\gamma}>
1-\frac{2\kappa\bar n}{\bar h}.}
\]
With binding caps, the exact minimum formulas must be used; a sufficient typewise comparison is \(\alpha+\lambda_F\kappa n_i\ge\gamma\), with strictness on a positive-mass uncapped-old set. These weakenings do not automatically preserve every individual conclusion above.

Condition (H) can be demanding. When the estate restriction binds and \(d_p\) is small, \(\Gamma_O\) can be much smaller than \(\gamma\). The required two-age patience can then be high. No calibration compatibility follows without a mapping from the quantitative model.

Saving and old income can substitute for one another. Define

\[
A(p)=\alpha+\frac{\vartheta p\kappa}{\chi+p\kappa}.
\]
For renters, and for owners in the explicit mortgage subcase \(\phi=q\), the sharper condition

\[
\boxed{
\Gamma_d\max\left\{\frac{\beta}{q},
                  \frac{E v}{K w}\right\}\ge A(p)
}
\tag{5}
\]
suffices for \(h^{o,eq}\ge h^{y,eq}\), including caps. With both caps slack, it is exact. The patience part follows from (4). For the income part, current expenditure \(B=c+ph\le w\) implies

\[
h^y\le\min\{H_d,A(p)w/(Ep)\},\qquad z\ge v,
\]
which yields the result. This is considerably more informative than bounding housing by the entire cash endowment.

### An explicit primitive recipient group

For a solved illustration, take **\(\phi=q,\ \tau^p=0\)** and require the following computed market demands to lie below their finite tenure caps. Define

\[
\begin{gathered}
w_{0i}=y_i^y+b_i,\qquad v_{0i}=y_i^o,\qquad D=E+\beta K,\\
x_i=\min\left\{\frac{w_{0i}}E,\frac{w_{0i}+qv_{0i}}D\right\},\\
c_i^{o,eq}=\frac{w_{0i}+qv_{0i}-Ex_i}{qK},\qquad
p=\frac{\nu\vartheta\bar x-\chi}{\kappa}>0,\\
n_i=\frac{\vartheta x_i}{\chi+p\kappa},\qquad
h_i^{y,eq}=\frac{A(p)x_i}{p},\qquad
h_i^{o,eq}=\frac{\Gamma_{d_i}c_i^{o,eq}}p.
\end{gathered}
\tag{6}
\]
These equations solve the stationary subcase, with \(N\) determined by housing clearing.

Conditional young choices coincide across tenures. The old value difference is constant across endowments, so the original logistic taste gives a constant owner share \(\pi\in(0,1)\). Explicitly,

\[
\pi=\operatorname{logit}^{-1}
\left(\frac{\bar\xi+\beta\Delta}{\sigma_\xi}\right),
\]
where \(\Delta=0\) on the unrestricted-estate branch; otherwise, writing \(J=\gamma+\omega_B\),

\[
\Delta=J\log J+\gamma\log\frac{1-q}{\gamma}
                    +\omega_B\log\frac q{\omega_B}.
\]
Let \(\bar\Gamma=(1-\pi)\gamma+\pi\Gamma_O\). The **aggregate**, rather than all-type, market-ordering restriction is exactly

\[
\boxed{\bar\Gamma\,\overline{c^{o,eq}}\ge A(p)\bar x.}
\tag{7}
\]
Under (7) and \(\alpha\ge\gamma\), every \(x_i\le\bar x\) gains both consumption and housing in the fixed-fertility planner, even if planner caps bind. A primitive positive-mass group with **strictly binding young finance** is

\[
\boxed{
F\left\{w_{0i}<E\bar x,\quad
qv_{0i}>\frac{\beta K}{E}w_{0i}\right\}>0.
}
\tag{8}
\]
Thus the relevant group is cash-poor relative to the cross-section but sufficiently future-resource-rich. Heterogeneity and both tenures remain.

These are redistributive utilitarian results. Under (H), patience itself contributes to the age gap. Only in the diagnostic case \(\beta=q,\Gamma_d=\gamma\), with slack physical caps, does the paired marginal housing gap reduce to \(\mu_iL_d\).

Finally, your correction to the previous certificate is right: it implied

\[
\beta<\frac{\nu\vartheta\bar w_0/\chi-E}{K}.
\]
Calling that certificate “beta-unrestricted” was incorrect. The new general ordering proof does not use it; the solved subcase has its own explicit \(p>0\) feasibility requirement. Pasted text

## 2. Fertility chosen privately within the assigned bundle

Let \(\mathfrak n(c,h)\) denote the parent’s optimal fertility when the current bundle and continuation opportunities are fixed. It is the unique zero, on \(0<n<\min\{c/\chi,h/\kappa\}\), of

\[
f(n;c,h)=\frac{\vartheta}{n}
-\frac{\chi}{c-\chi n}
-\frac{\alpha\kappa}{h-\kappa n}.
\]
This is exactly the maintained parental fertility condition, not a valuation of unborn utility. Pasted text

Since \(f_n<0\), differentiation confirms

\[
\boxed{
dn=
\frac{(\chi/x^2)\,dc+(\alpha\kappa/s^2)\,dh}
{\vartheta/n^2+\chi^2/x^2+\alpha\kappa^2/s^2}.
}
\tag{9}
\]
For a finite change accommodating baseline fertility \(n_0\),

\[
\boxed{
n_1\ge n_0
\Longleftrightarrow
\frac{\chi}{c_1-\chi n_0}
+\frac{\alpha\kappa}{h_1-\kappa n_0}
\le\frac{\vartheta}{n_0}.
}
\tag{10}
\]
If the new bundle cannot accommodate \(n_0\), fertility falls.

In particular, let housing rise by \(\Delta h>0\), while consumption falls by \(\delta\ge0\). Evaluating \(x,s\) at the original bundle,

\[
\boxed{
n_1\ge n_0
\Longleftrightarrow
\delta\le
\frac{\alpha\kappa x^2\Delta h}
{\chi s(s+\Delta h)+\alpha\kappa x\Delta h}.
}
\tag{11}
\]
Strict inequality gives strictly higher fertility.

Locally, the allowable consumption loss per unit of housing is

\[
-\frac{dc}{dh}<\frac{\alpha\kappa x^2}{\chi s^2}.
\]
For the illustrative exchange \(dc=-p\,dh\), at an uncapped, financially unrestricted reference bundle \(s=\alpha x/p\), this becomes

\[
\boxed{p\kappa>\alpha\chi.}
\]
The space component of child costs must be sufficiently important. This is a bundle-composition test, not an assertion that a particular tax or credit policy implements the exchange.

### Applying the test to the actual fixed-fertility optimum

Put \(s_i^F=h_i^{y,F}-\kappa n_i\). Then

\[
\boxed{
n_i^S>n_i
\Longleftrightarrow
\frac{\chi}{x^F}+\frac{\alpha\kappa}{s_i^F}
<\frac{\vartheta}{n_i}.
}
\tag{12}
\]
The low-resource recipients identified in (6)–(8) gain both components and therefore have strictly higher conditional fertility. Other young households may face a consumption–space tradeoff.

Crucially, \(x^F,s_i^F\) evaluate the new bundle **at original fertility**. After fertility responds, adult goods and space are

\[
x_i^S=x^F-\chi(n_i^S-n_i),\qquad
s_i^S=s_i^F-\kappa(n_i^S-n_i).
\]
An aggregate result needs to include losers. One cap-valid sufficient test is

\[
\boxed{
\int\mathfrak n(x^F,s_i^F)\,dQ
>\frac{1+\alpha}{E}\bar n
\quad\Longrightarrow\quad
\bar n^S>\bar n.
}
\tag{13}
\]
To prove it, \(\mathfrak n\) is homogeneous and concave: its upper contour sets are convex by (10), and homogeneity turns this into concavity. Hence it is superadditive, and

\[
n_i^S
=\mathfrak n(x^F+\chi n_i,s_i^F+\kappa n_i)
\ge\mathfrak n(x^F,s_i^F)+\frac{\vartheta}{E}n_i.
\]
Integrating proves (13). This is a sufficient allocation test, not a primitive equilibrium restriction.

No market-policy sign has been claimed. The direct assignments are financed by the maintained balanced transfers

\[
t_k=\Delta c_k+u_t\Delta h_k,
\]
with the previously specified financial adjustments preserving incumbent future resources and estates. A cash grant or credit reform instead requires household reoptimization, fiscal funding, price clearing, and—when tenure is free—the ownership-selection contribution to average fertility.

## 3. The same planner chooses fertility, valuing parents only

Use \(x_i=c_i^y-\chi n_i\) and \(s_i=h_i^y-\kappa n_i\). The normalized problem is

\[
\max\int\left[
\log x_i+\alpha\log s_i+\vartheta\log n_i
+\log c_i^o+\gamma\log h_i^o
\right]dQ
\]
subject to

\[
\begin{aligned}
\int(x_i+\chi n_i+c_i^o)dQ&=\bar c,\\
\int(s_i+\kappa n_i+h_i^o)dQ&=\bar h,\\
s_i+\kappa n_i&\le H_{d_i},\qquad h_i^o\le H_{d_i},
\end{aligned}
\tag{14}
\]
and positive log arguments. This merely rewrites total \(c,h\); it does not add child costs a second time. No replacement-fertility constraint is imposed. Pasted text

The objective is strictly concave in these variables and the constraints are affine. Averaging within retained tenure preserves feasibility and raises utility, reducing existence and characterization to a finite-dimensional problem. The positive reference supplies feasibility.

Let \(\lambda_C>0,\lambda_H>0\) be the goods and housing resource multipliers. Let

\[
r_d=\lambda_H+\eta_d,\qquad \eta_d\ge0
\]
include the young tenure-cap multiplier. The unique optimum satisfies

\[
\boxed{
\begin{aligned}
x_d^J=c_d^{o,J}&=\frac1{\lambda_C},\\
n_d^J&=\frac{\vartheta}{\chi\lambda_C+\kappa r_d},\\
h_d^{y,J}&=\frac{\alpha}{r_d}+\kappa n_d^J,\\
h_d^{o,J}&=\min\{H_d,\gamma/\lambda_H\}.
\end{aligned}}
\tag{15}
\]
For an uncapped young tenure, \(r_d=\lambda_H\). Otherwise \(r_d>\lambda_H\) is the unique solution of

\[
\frac{\alpha}{r_d}
+\frac{\kappa\vartheta}{\chi\lambda_C+\kappa r_d}=H_d.
\]
Writing \(\pi_d=Q(d)\), the resource equations determining the multipliers are

\[
\frac2{\lambda_C}+\chi\sum_d\pi_dn_d^J=\bar c,\qquad
\sum_d\pi_d(h_d^{y,J}+h_d^{o,J})=\bar h.
\tag{16}
\]
Thus individual fertility rises exactly for \(n_i<n_{d_i}^J\); individual housing rises exactly for \(h_i^{y,eq}<h_{d_i}^{y,J}\). The fixed-fertility recipient set need not remain unchanged.

### Useful cap-aware sufficient inequalities

The optimality conditions imply

\[
\lambda_C\le U_C:=\frac{2+\vartheta}{\bar c},\qquad
r_d\le U_d:=
\max\left\{\frac{\alpha+\gamma+\vartheta}{\bar h},
           \frac{\alpha+\vartheta}{H_d}\right\}.
\]
Indeed, \(\chi\lambda_C n_d^J<\vartheta\); summing the goods and housing optimality identities gives the resource bounds, while a binding young cap implies \(r_d<(\alpha+\vartheta)/H_d\).

Consequently,

\[
\underline n_d=
\frac{\vartheta}{\chi U_C+\kappa U_d},\qquad
\underline h_d=\frac{\alpha}{U_d}+\kappa\underline n_d
\]
are lower bounds on joint-planner fertility and young housing. Therefore

\[
\boxed{
n_i<\underline n_{d_i}\Rightarrow n_i^J>n_i,\qquad
h_i^{y,eq}<\underline h_{d_i}\Rightarrow h_i^{y,J}>h_i^{y,eq},
}
\tag{17}
\]
and, separately,

\[
\boxed{\sum_d\pi_d\underline n_d>\frac1\nu
\Rightarrow\bar n^J>\bar n.}
\tag{18}
\]
These are conservative **resource restrictions**, not primitive equilibrium inequalities. In the solved benchmark (6), the resource totals are explicit functions of primitives.

Importantly, under \(\alpha\ge\gamma\), (15) also gives

\[
H_y^J>\bar H/2.
\]
If young housing is capped it is at least old housing; otherwise \(r_d=\lambda_H\) and positive child space makes it larger. Thus the housing theorem’s aggregate gain survives free fertility **even with caps**, independently of whether average fertility rises.

### A sharper fertility theorem when the joint optimum is uncapped

Then fertility is common across all young households and is the unique solution of

\[
\boxed{
\frac{\vartheta}{n^J}
=\frac{2\chi}{\bar c-\chi n^J}
+\frac{\kappa(\alpha+\gamma)}{\bar h-\kappa n^J}.
}
\tag{19}
\]
The exact aggregate comparison is

\[
\boxed{
n^J>\bar n
\Longleftrightarrow
\frac{\vartheta}{\bar n}>
\frac{2\chi}{\bar c-\chi\bar n}
+\frac{\kappa(\alpha+\gamma)}{\bar h-\kappa\bar n}.
}
\tag{20}
\]
**Under (H), this inequality follows from household optimality**, provided the solution of (19) respects the finite caps. The Euler condition gives \(\overline{c^{o,eq}}>\bar x\), and the housing theorem gives

\[
\overline{h^{o,eq}}\ge\bar s+\kappa\bar n>
\frac{\gamma}{\alpha}\bar s.
\]
Moreover, multiplying the reference fertility condition by \(n_i\), integrating, and applying Cauchy–Schwarz yields

\[
\frac{\vartheta}{\bar n}
\ge\frac{\chi}{\bar x}+\frac{\alpha\kappa}{\bar s}.
\]
The right side of (20) is strictly smaller than this last expression. Hence \(n^J>\bar n\), and every initially below-average-fertility parent gains fertility.

If **both** planner solutions are uncapped,

\[
H_y^J-H_y^F
=\frac{N\gamma\kappa}{\alpha+\gamma}(n^J-\bar n)>0.
\tag{21}
\]
Concavity of \(\mathfrak n\), applied to the sequential bundles, also gives \(\bar n^S<n^J\) in this regime. There is no corresponding unrestricted componentwise ordering between \(J\) and \(S\). Generally, only their welfare ordering is automatic:

\[
\mathcal W^J\ge\mathcal W^S\ge\mathcal W^F.
\]
Free fertility need not rise without the restrictions. In the uncapped, young-finance-slack, unrestricted-estate subcase,
\(c_i^{o,eq}=(\beta/q)x_i\) and \(h_i^{o,eq}=\gamma c_i^{o,eq}/p\). Setting \(r=\beta/q\), the fertility test has the sign of

\[
(r-1)\left[\frac{\chi}{1+r}
+\frac{p\kappa\gamma}{\alpha+\gamma r}\right].
\]
Thus \(r<1\) gives \(n^J<\bar n\), including with heterogeneous incomes. This is an age-weighting effect, not a value assigned to unborn people.

## Bridge to the transition

The conditional tests (9)–(12) apply at each date using the same exogenous \(\vartheta_t\) on both paths; they require no convergence. The stationary housing-pairing shortcut does not automatically apply when age distributions differ. A policy path additionally requires fiscal funding, prices, financial obligations, housing clearing, and tenure selection. From common initial cohorts, date-by-date fertility ordering implies cohort ordering through

\[
\frac{Y_t^{policy}}{Y_t^{baseline}}
=\prod_{s=t_0}^{t-1}
\frac{\bar n_s^{policy}}{\bar n_s^{baseline}},
\qquad O_{t+1}=Y_t.
\]
It does not prove that a policy produces that fertility ordering. Positive stationary endpoints both have \(\bar n=1/\nu\); population-level differences can reflect fertility accumulated along the path. Today’s parent-only exercise neither imposes replacement fertility after reallocation nor fixes the future allocations of an unchanged set of unborn households.
