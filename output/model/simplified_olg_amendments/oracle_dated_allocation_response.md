# Pro response — one-date allocation, individual and aggregate comparisons

Source: https://chatgpt.com/c/6aa06dde-50d0-83ea-ab19-bf66d86b0354

Completed response: “Worked for 23m 41s”. Captured 2026-09-08 from visible browser text and mathematical source labels. Display delimiters reconstructed for readability. This is an external review, not an accepted change to the model or paper.

## Result

**A full-optimum individual result is possible, but “the young household is borrowing-constrained” is not sufficient.** The additional restriction must compare its current cash with the resources available across the cross-section. Below, that comparison is expressed entirely through parameters and the endowment distribution—not equilibrium housing choices or assumed marginal-utility gaps.

The aggregate conclusion requires a separate condition. In fact, there is an analytical family in which **every young household is financially constrained, yet the full dated planner reduces aggregate young housing**.

## 1. The dated planner and financial feasibility

Write \(H_d=h_d^{\max}\), for tenure \(d\in\{R,O\}\). Let \(Q\) be the common stationary distribution of complete types and retained tenure, and let \(N\) be each age-group’s reference mass. Throughout, individual fertility, tenure, future real allocations, and estate payments remain fixed. This is the one-date, equally weighted remaining-utility comparison specified in the follow-up. Pasted text

Let \(C^{\mathrm{eq}}\) be total reference current consumption expenditure. The planner solves

\[
\max N\int\left[
\log(c_i^y-\chi n_i)+\alpha\log(h_i^y-\kappa n_i)
+\log c_i^o+\gamma\log h_i^o
\right]dQ
\]
subject to

\[
N\int(c_i^y+c_i^o)dQ=C^{\mathrm{eq}},
\qquad
N\int(h_i^y+h_i^o)dQ=\bar H,
\]
the positive log arguments, and

\[
h_i^y,h_i^o\le H_{d_i}.
\]
Consumption and housing separate. Consequently,

\[
c_i^{y,SP}=\chi n_i+x^*,\qquad
c_i^{o,SP}=x^*,\qquad
x^*=\frac12\left(\frac{C^{\mathrm{eq}}}{N}-\frac{\chi}{\nu}\right),
\]
and

\[
\boxed{
h_i^{y,SP}=\min\{H_{d_i},\kappa n_i+\alpha/\lambda\},
\qquad
h_i^{o,SP}=\min\{H_{d_i},\gamma/\lambda\}.
}
\tag{1}
\]
Here \(\lambda>0\) clears housing when the stock is below total retained capacity. These are the **full dated optimum’s** housing choices, not a consumption-fixed approximation.

### Financial implementation

For any such reallocation, hold the price path fixed and use the proposed owner adjustments

\[
\Delta a_i'=-P_{t+1}\Delta h_i^y,
\qquad
\Delta a_j^e=-qP_{t+1}\Delta h_j^o.
\]
They preserve, respectively,

\[
\Delta(a_i'+P_{t+1}h_i^y)=0,
\qquad
\Delta(q^{-1}a_j^e+P_{t+1}h_j^o)=0.
\]
The required current transfer to each household is exactly

\[
t_k=\Delta c_k+u_t\Delta h_k.
\]
Thus \(\int t_k\,dk=0\). These identities follow from the original purchase, mortgage-repayment, and estate budgets. Pasted text

Including rental-intermediary financing, the change in aggregate net financial payoffs is

\[
-P_{t+1}\Delta H^{\mathrm{own}}
-P_{t+1}\Delta H^{\mathrm{rent}}=0.
\]
Existing obligations need not be written down; offsetting financial positions implement the changes.

**No choice between internal and external estate recipients is needed:** each estate payment is unchanged. The old owner’s retained physical cap is \(H_O\), not \(\min\{H_O,e/P_{t+1}\}\). The latter would incorrectly reimpose nonnegative old financial saving, which this planner is authorized to relax. Pasted text

## 2. One proposition with primitive sufficient conditions

The following conditions are conservative certificates, but they allow binding physical caps and both old-owner estate regimes.

Define

\[
\begin{gathered}
w_0=y^y+b,\qquad v_0=y^o,\qquad
a=\alpha+\vartheta,\qquad E=1+a,\\
K=1+\gamma+\omega_B,\qquad D=E+\beta K,\\
d_p=1-q+q\tau^p,\qquad d_L=1-\phi+q\tau^p,\\
\ell_R=1,\qquad \ell_O=d_L/d_p,\\
\Gamma=\min\left\{\gamma,\frac{(\gamma+\omega_B)d_p}{1+q\tau^p}\right\}.
\end{gathered}
\tag{2}
\]
Assume bounded endowments, \(w_0\ge\underline w>0\), \(v_0\ge0\), and

\[
\boxed{0\le\tau^p<2,\qquad \phi\ge q.}
\tag{3}
\]
The second restriction gives \(0<\ell_d\le1\). There is **no restriction on \(\beta>0\)** and no requirement that \(\alpha\ge\gamma\).

The following are primitive upper bounds on the rebate and service price:

\[
\bar M_0=\int(w_0+qv_0)dF,\qquad
\bar T=\frac{\tau^p\bar M_0}{(1-q)(2-\tau^p)},\qquad
\bar p=\frac{\nu[\bar M_0+(1+q)\bar T]}{\kappa}.
\tag{4}
\]
Choose a primitive lower price certificate \(\underline p>0\) as described in the short appendix below.

Define the **cross-sectional adult-space resource bound**

\[
\boxed{
\mathcal R_-=
\int\left[
\alpha\min\left\{\frac{w_0}{D},\frac{\underline p H_R}{a}\right\}
+
\min\left\{\underline p H_R,
\frac{\beta\Gamma(w_0+qv_0)}{qD}\right\}
\right]dF.
}
\tag{5}
\]
This uses clipped moments of current and lifetime resources; it does not require high old income for every type.

Finally, put

\[
B_d=
1+\alpha\ell_d+
\vartheta\frac{\chi+\ell_d\bar p\kappa}{\chi+\bar p\kappa},
\qquad
k_d=\frac{D}{B_d}-1.
\tag{6}
\]

### Proposition

Consider any positive stationary competitive equilibrium of the maintained model.

**Individual conclusion.** For tenure \(d\), define the endowment set \(\mathcal S_d\) by the two inequalities

\[
\boxed{
\begin{aligned}
w_0+\bar T
&<
E\ell_d
\min\left\{
\frac{\mathcal R_-}{\alpha+\gamma},
\frac{\underline p H_d}{a}
\right\},\\
qv_0
&>
k_dw_0+(k_d-q)_+\bar T.
\end{aligned}
}
\tag{I}
\]
Every current young household whose endowments belong to \(\mathcal S_d\) and whose equilibrium tenure is \(d\):

- is strictly below its market housing cap;

- has a strictly positive multiplier on its young financing restriction;

- receives strictly more housing at the full dated optimum:

\[
\boxed{h_i^{y,SP}>h_i^{y,\mathrm{eq}}.}
\]
If \(F(\mathcal S_d)>0\) for either tenure, these households have positive \(Q\)-mass. Logistic tastes give each feasible tenure strictly positive probability at every endowment type.

**Separate aggregate conclusion.** Let \(\widehat w=w_0+\bar T\), and define

\[
G_d(w_0)=
\min\left\{
\left[\underline p H_d-\frac{a\widehat w}{E\ell_d}\right]_+,\;
\alpha\left[
\frac{\mathcal R_-}{\alpha+\gamma}
-\frac{\widehat w}{E\ell_d}
\right]
\right\}.
\tag{7}
\]
If the additional distributional inequality

\[
\boxed{
\mathcal G\equiv
\int\min\{G_R(w_0),G_O(w_0)\}\,dF>0,
}
\tag{A}
\]
holds, then

\[
\boxed{
H_Y^{SP}-H_Y^{\mathrm{eq}}
\ge\frac{N}{p}\mathcal G>0.
}
\tag{8}
\]
Condition (A) allows some young households to lose housing: their potentially negative contributions are explicitly subtracted. It does not infer aggregation from the individual conclusion.

The conclusions hold in **every positive reference equilibrium** satisfying these primitive restrictions. This is not an equilibrium-existence or uniqueness theorem.

## 3. Proof

### Market bounds, including binding caps

Fix an equilibrium and a conditional tenure. Write

\[
w=w_0+T,\qquad v=v_0+T,\qquad M=w+qv,
\]
and let \(\Lambda,\mu,\eta_y\) be the lifetime-budget, young-finance, and young-cap multipliers. With \(x=c-\chi n\) and \(s=h-\kappa n\),

\[
\frac1x=\Lambda+\mu,\qquad
\frac{\alpha}{s}=\Lambda p+\mu\ell_dp+\eta_y,\qquad
\frac{\vartheta}{n}=\frac{\chi}{x}+\frac{\alpha\kappa}{s}.
\tag{9}
\]
Multiplying the household first-order conditions by quantities and summing gives

\[
\Lambda M+\mu w\le D.
\tag{10}
\]
The omitted terms are nonnegative physical-cap contributions; the homogeneous estate restriction contributes zero. In particular,

\[
x\ge w/D,\qquad \Lambda\le D/M.
\]
Because \(\ell_d\le1\), the effective current housing price

\[
\rho=\frac{\alpha x}{s}
\]
satisfies \(\rho\ge\ell_dp\). Combining this with the fertility condition and the cash constraint gives

\[
\boxed{
ps\le\frac{\alpha w}{E\ell_d},
\qquad
ph<\frac{aw}{E\ell_d}.
}
\tag{11}
\]
For example, the fertility condition implies

\[
ax=\rho h+\chi n,
\]
so \(\ell_dph<ac\), which proves the second inequality. For the first, \(x\ge\ell_dps/\alpha\) and

\[
n\ge\frac{\vartheta\ell_dps}
{\alpha(\chi+\ell_dp\kappa)};
\]
substitution into \(c+\ell_dph\le w\) proves the claim.

These bounds also establish a lower resource bound. If young housing is uncapped, \(\rho\le p\), hence \(ps\ge\alpha w/D\). If it is capped, the fertility condition gives \(s>\alpha H_d/a\). Therefore

\[
ps\ge\alpha\min\{w/D,pH_d/a\}.
\tag{12}
\]
For old housing, define \(\Gamma_R=\gamma\) and \(\Gamma_O=\Gamma\). In either uncapped old-estate regime,

\[
c^o=z/K,\qquad ph^o=\Gamma_d z/K.
\]
The minimum defining \(\Gamma\) explicitly includes zero old financial saving; it is not an assumption of a slack estate restriction. Pasted text

Let \(m=1/c^o\). Since \(\Lambda=\beta m/q\), (10) gives

\[
m\le\frac{qD}{\beta M}.
\]
If old housing is uncapped, \(ph^o=\Gamma_d/m\); if capped, \(ph^o=pH_d\). Consequently,

\[
ph^o\ge
\min\left\{pH_d,\frac{\beta\Gamma_dM}{qD}\right\}.
\tag{13}
\]
Using \(p\ge\underline p\), \(w\ge w_0\), \(M\ge w_0+qv_0\), and taking the lower envelope across tenures in (12)–(13), we obtain

\[
\boxed{
p\left(\frac{\bar H}{N}-\frac{\kappa}{\nu}\right)
=p\int(s_i+h_i^o)dQ\ge\mathcal R_-.
}
\tag{14}
\]

### Comparing the same household with the full optimum

From the capped planner formula (1),

\[
\frac{\bar H}{N}-\frac{\kappa}{\nu}
\le\frac{\alpha+\gamma}{\lambda}.
\]
Thus

\[
\frac{p\alpha}{\lambda}
\ge\frac{\alpha\mathcal R_-}{\alpha+\gamma}.
\tag{15}
\]
The first inequality in (I), together with (11), proves both

\[
h_i^{y,\mathrm{eq}}<H_d,
\qquad
ps_i^{\mathrm{eq}}
<
\frac{\alpha\mathcal R_-}{\alpha+\gamma}
\le\frac{p\alpha}{\lambda}.
\]
Both arguments of the minimum defining \(h_i^{y,SP}\) therefore exceed \(h_i^{y,\mathrm{eq}}\). This proves the individual full-optimum comparison.

### Establishing binding finance

Suppose instead that \(\mu=0\) for one of these households. Its young cap is already known to be slack. Equation (9) then implies

\[
s=\frac{\alpha x}{p},\qquad
n=\frac{\vartheta x}{\chi+p\kappa},
\]
and its cash expenditure is

\[
c+\ell_dph
=
x\left[
1+\alpha\ell_d+
\vartheta\frac{\chi+\ell_dp\kappa}{\chi+p\kappa}
\right].
\]
The bracket decreases with \(p\), so it is at least \(B_d\). Moreover, (10) with \(\mu=0\) gives \(x\ge M/D\). Cash feasibility therefore requires

\[
w\ge B_dM/D,
\quad\text{equivalently}\quad qv\le k_dw.
\]
But the second inequality in (I) guarantees \(qv>k_dw\) for every \(T\in[0,\bar T]\). Contradiction. Hence \(\mu>0\).

Notice that this argument **does not require the old housing cap or old estate restriction to be slack**.

### Aggregation with losses included

Exactly,

\[
p(h_i^{y,SP}-h_i^{y,\mathrm{eq}})
=
\min\left\{
p(H_d-h_i^{y,\mathrm{eq}}),\;
\frac{p\alpha}{\lambda}-ps_i^{\mathrm{eq}}
\right\}.
\]
Equations (11) and (15) bound this below by \(G_d(w_{0i})\). Integrating, and using \(G_{d_i}\ge\min\{G_R,G_O\}\), proves (8). ∎

## 4. Interpretation and a precise obstruction

The first inequality in (I) identifies **cash-poor households relative to the economy’s available adult space**, while also guaranteeing room below their retained tenure cap. The second establishes that future resources are sufficiently large relative to current cash to make finance strictly restrictive.

For renters, \(B_R=E\), so the finance test without taxes is simply

\[
\frac{v_0}{w_0}>\frac{\beta K}{qE}.
\]
At \(\beta=q\), this is \(K/E\), rather than the conservative \(K/\gamma\) needed to force old housing above young housing type by type. **These tests establish different facts:** the lower threshold establishes binding finance; the cash and distributional conditions establish the full housing comparison.

Binding finance alone cannot replace those additional conditions. Here is an analytical counterfamily.

Take the explicitly tax-free subcase

\[
\phi=q,\qquad \beta=q,\qquad
\omega_B(1-q)<q\gamma,
\]
with finite caps large enough to be inactive. Then

\[
\Gamma=(\gamma+\omega_B)(1-q)<\gamma.
\]
Conditional tenure-value differences are constant across endowments, so logistic tastes can generate any constant owner share \(\pi\in(0,1)\). Put

\[
\bar\Gamma=(1-\pi)\gamma+\pi\Gamma,
\]
and choose

\[
\boxed{\frac KE<r<\frac{\gamma K}{E\bar\Gamma}.}
\tag{16}
\]
Let current cash be genuinely heterogeneous and let \(v_0=rw_0\).

Every young household has strictly binding finance, because \(r>K/E\). Nevertheless,

\[
ps_i=\frac{\alpha w_{0i}}E,
\qquad
p\bar h^o=\frac{\bar\Gamma r}{K}\,\overline{w_0}.
\]
The uncapped full dated optimum consequently satisfies

\[
\boxed{
H_Y^{SP}-H_Y^{\mathrm{eq}}
=
\frac{N\alpha\overline{w_0}}{p(\alpha+\gamma)}
\left(\frac{\bar\Gamma r}{K}-\frac{\gamma}{E}\right)<0.
}
\tag{17}
\]
This is compatible with \(\alpha\ge\gamma\). It is also a stationary family: choose positive

\[
\chi<\frac{\vartheta\nu\overline{w_0}}E,
\qquad
p=\frac{\vartheta\nu\overline{w_0}/E-\chi}{\kappa},
\]
and determine \(N\) from housing clearing. No numerical example or neighborhood argument is involved.

The obstruction is substantive: **old owners’ inability to borrow against their estate changes their housing demand sufficiently that relaxing finance need not move housing toward the young.**

Conversely, equal-weight redistribution can produce individual housing gains without any financing distortion. In the uncapped, slack-finance, positive-financial-estate benchmark with \(\beta=q\),

\[
h_i^{y,SP}-h_i^{y,\mathrm{eq}}
=\frac{\alpha}{p}(\bar x-x_i),
\qquad
H_Y^{SP}-H_Y^{\mathrm{eq}}=0.
\]
Thus the proposition identifies households that **both are constrained and gain housing**. It does not attribute their entire gain exclusively to the borrowing constraint.

## Technical appendix: the primitive price certificate

The bounds \(T\le\bar T\) and \(p<\bar p\) follow by aggregating the original lifetime budget and using replacement fertility; the packet provides the corresponding rebate and price bounds. Pasted text

A sharper lower bound than one based only on minimum wealth can use the distribution. For each \(w_0\), let \(r_c(w_0)\in(0,H_R/\kappa)\) solve the elementary quadratic equation

\[
\frac{\vartheta}{r_c}
=
\frac{\chi D}{w_0}
+\frac{\alpha\kappa}{H_R-\kappa r_c}.
\]
Choose any \(\underline p>0\) satisfying

\[
\boxed{
\int
\min\left\{
\frac{\vartheta w_0}{D(\chi+\underline p\kappa)},
r_c(w_0)
\right\}dF
>\frac1\nu.
}
\tag{18}
\]
To verify it, uncapped young households satisfy

\[
n\ge\frac{\vartheta w_0}{D(\chi+p\kappa)};
\]
capped young households satisfy \(n\ge r_c(w_0)\), because \(x\ge w_0/D\) and their cap is at least \(H_R\). If \(p\le\underline p\), (18) would imply fertility above replacement, contradicting stationarity.

All quantities in (I) and (A) are therefore functions of primitives and the chosen scalar certificate. Failure of this sufficient certificate—or of (A)—makes these bounds inconclusive; it does not establish the opposite allocation.

**The resulting scope is precise:** the individual theorem identifies a positive-mass group without assuming its constraint pattern; the aggregate theorem separately accounts for possible young losers; and both retain the original mortgage, positive child costs, persistent tenure, and physical caps.
