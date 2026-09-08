# Pro response — September 8, 2026

Source: https://chatgpt.com/c/6aa06dde-50d0-83ea-ab19-bf66d86b0354

Captured from the completed browser response. Formatting reconstructed from visible text and equation labels; claims remain advisory until independently checked.

## Verdict

**Use the equally weighted lifetime welfare of a stationary cohort, with population and aggregate external wealth fixed at their reference values.** This preserves the warm-glow household model, but requires an explicit completion of its accounting.

Under the completion below, every positive competitive steady state fails the full planner’s **consumption** optimality condition when \(q<1\). That failure therefore cannot be attributed entirely to borrowing restrictions or housing segmentation.

The stronger housing conclusion is conditional. I establish conditions under which the **fully optimized allocation**, including consumption, estates, and tenure, gives more housing to the young. I also give an analytical family in which it gives them **less**.

## 1. Welfare: one recommendation and one comparator

Let \(Q\) denote the distribution of the complete type \(i=(y_i^y,b_i,y_i^o,\xi_i)\). Fix \(n_i=n_i^{\rm eq}\) for every type, including when its assigned tenure changes. Write \(d_i=1\) for ownership.

### Recommended criterion: stationary cohort lifetime welfare

Use

\[
\boxed{
\mathcal W_C
=N\int\left[
\log(c_i-\chi n_i)+\alpha\log(h_i-\kappa n_i)
+\vartheta\log n_i+\xi_i d_i
+\beta\bigl(\log c_i^2+\gamma\log h_i^2+\omega_B\log e_i\bigr)
\right]\,dQ .
}
\tag{1}
\]

This gives each household equal weight on utility measured at entry. It retains the note’s preferences, including warm-glow estates—not descendants’ continuation utility. With fixed \(N\), dividing by \(N\), or by the living population \(2N\), changes only normalization. Pasted text

The literature does **not** select these equal weights uniquely. Golosov–Jones–Tertilt, §3.4, Result 1 and equation (1), consider positive weights on potential people; Result 2 and equation (2) instead consider weights on initial agents, with a uniqueness requirement. These characterize different efficiency-supporting problems, not a uniquely mandated utilitarian objective. [Tertilt Lab](https://tertilt.vwl.uni-mannheim.de/research/optimality_Econometrica.pdf)

Nor does a preference for dynastic conventions make this household a dynasty. In a recursive dynasty, social weights across founders multiply utilities that already incorporate descendants through **private altruism weights**. Here \(\beta\) weights the household’s own old age. The recursive preferences and inherited-bequest budgets explicitly present in Golosov–Jones–Tertilt, §4, pp. 1059–1060, are absent. I do not recommend adding them. [Tertilt Lab](https://tertilt.vwl.uni-mannheim.de/research/optimality_Econometrica.pdf)

Farhi–Werning make the distinction particularly clear: their supplied manuscript’s §2, equations (1) and (5)–(7), separates private dynastic discounting from society’s additional direct concern for descendants. That argument does not convert warm glow into altruism. Equation references here are to that verified manuscript; the supplied Becker–Barro NBER copy was inaccessible. [Massachusetts Institute of Technology](https://web.mit.edu/iwerning/Public/inequality_social_screen_old.pdf)

### Comparator: discounted welfare over successive cohorts

For an actual path, use

\[
\mathcal W_D
=\rho N\int u_{i,0}^{o,\mathrm{initial}}\,dQ
+\sum_{t=0}^{\infty}\delta^tN
\int\left[u_{i,t}^y+\xi_i d_{i,t}+\beta u_{i,t+1}^o\right]dQ,
\qquad 0<\delta<1.
\tag{2}
\]

Here \(\delta\) discounts successive cohort lifetimes; \(\rho\) weights the initial old. Neither is determined by \(\beta\).

At date \(t\ge1\), young utility has coefficient \(\delta^t\), and current old utility has coefficient \(\beta\delta^{t-1}\). Thus the relative age weight in dated first-order conditions is

\[
\boxed{\frac{\text{old weight}}{\text{young weight}}=\frac{\beta}{\delta}.}
\tag{3}
\]

More generally, cohort weights \(a_t\) give \(\beta a_{t-1}/a_t\). Extending geometric birth weights backward gives \(\rho=\beta/\delta\); other initial-old treatment must be stated.

Restricting (2) to stationary allocations produces an old-flow coefficient

\[
\beta+\rho(1-\delta)
\]

after normalization—not generally \(\beta\). Moreover, the unrestricted dynamic planner has inherited-state and asset-accumulation conditions. With interior external saving, its resource multipliers satisfy \(q\lambda_t=\lambda_{t+1}\); constant positive consumption with geometric cohort weights requires \(\delta=q\). A stationary maximization therefore need not describe its limiting allocation.

Equal weights on every cohort correspond to \(\delta=1\), requiring an appropriate infinite-horizon comparison convention. They are not the note’s remaining-lifetime objective, whose stationary old-flow coefficient is \(1+\beta\). Pasted text

## 2. The accounting completion

I recommend the following explicit additions.

**Outside transfers.** Entrant wealth \(b_i\) is an exogenous outside remittance. Estates are delivered to outside recipients in goods, after liquidation, at the end of old age. The house itself remains domestic: selling its title is not a resource loss. Welfare covers the modeled domestic households, whose warm glow values the estate payment.

**Domestic rental intermediaries.** Rental titles are held by zero-equity domestic intermediaries. At date \(t\), prepaid rent and borrowing against resale finance acquisition and prepaid taxes:

\[
u_tH_t^R+qP_{t+1}H_t^R
=(1+q\tau^p)P_tH_t^R.
\tag{4}
\]

At the next date, sale proceeds repay the face-value debt \(P_{t+1}H_t^R\). There are no omitted foreign landlords or endowed intermediary equity. Titles outside living owners are held by these intermediaries or domestic estate executors pending liquidation.

**Real estates and unrestricted internal finance.** The planner chooses an actual goods payment \(e_i\), not an index increased by changing house prices. Use \(P=P^{\rm eq}\) only for settlement bookkeeping:

\[
a_i^e=q(e_i-Pd_i h_i^2).
\tag{5}
\]

The direct planner may make this negative. It relaxes young renters’ borrowing prohibition, owners’ mortgage limit, and old owners’ nonnegative financial-saving restriction. The last restriction is financial, not physical, exactly as the note’s estate equations imply. Pasted text

Define aggregate external bond payoffs **before estate remittances**. At stationarity,

\[
H^R=N\int(1-d_i)(h_i+h_i^2)\,dQ,
\qquad
B^{\rm ext}
=N\int(a_i'+q^{-1}a_i^e)\,dQ-PH^R.
\tag{6}
\]

The last term is intermediary debt. Compute \(B^*=B^{\rm ext,eq}\) from the reference equilibrium and hold it fixed.

With \(E_t\) denoting estates actually remitted at date \(t\), consolidation gives

\[
C_t+qB_{t+1}^{\rm ext}+E_t
=Y_t^g+I_t+B_t^{\rm ext}.
\tag{7}
\]

Consequently,

\[
\boxed{C+E=Y^g+I+(1-q)B^*\equiv\Omega,}
\tag{8}
\]

where

\[
C=N\int(c_i+c_i^2)dQ,\quad
E=N\int e_i\,dQ,\quad
Y^g=N\int(y_i^y+y_i^o)dQ,\quad I=N\int b_i\,dQ.
\]

Child goods are already in \(c_i\). Fully rebated property taxes cancel. Existing-house purchases cancel.

The timing of the fixed wealth stock matters: fixing \(B-E\) instead would yield a different restriction, \(C+qE=Y^g+I+(1-q)(B-E)\). These are not interchangeable definitions of “the same resources.”

## 3. The full stationary planner

Maximize (1) over measurable

\[
(c_i,c_i^2,h_i,h_i^2,e_i,d_i),\qquad d_i\in\{0,1\},
\]

subject to

\[
\begin{aligned}
N\int(c_i+c_i^2+e_i)dQ&=\Omega,\\
N\int(h_i+h_i^2)dQ&=\bar H,\\
c_i>\chi n_i,\qquad h_i&>\kappa n_i,\\
c_i^2,h_i^2,e_i&>0,\\
h_i,h_i^2&\le h_{d_i}^{\max},
\end{aligned}
\tag{9}
\]

with the **same \(d_i\)** at both ages.

Fixed objects include \(Q,n_i,N,\bar H,B^*\), all incomes and preference parameters. In particular, \(\bar n=1/\nu\), but \(N\) is not chosen. The financial positions and balanced transfers implementing (9) are reconstructed below; they are not omitted allocation controls.

This is a full-information, redistributive planner. Private budgets are not additional feasibility restrictions after individual finance is relaxed.

## 4. Main proposition

Define

\[
\mathcal A=\frac{\Omega}{N}-\chi\bar n,\qquad
S=\bar H-\kappa N\bar n,
\]

and

\[
\widehat s=\frac{\alpha S}{N(\alpha+\beta\gamma)},
\qquad
\widehat h^2=\frac{\beta\gamma S}{N(\alpha+\beta\gamma)}.
\tag{10}
\]

**Proposition.** Assume a positive reference steady state with finite resources and welfare.

**Full-planner consumption and estates.** Every optimum has

\[
\boxed{
x_i^*=x^*=\frac{\mathcal A}{1+\beta+\beta\omega_B},
\qquad c_i^{2*}=\beta x^*,
\qquad e_i^*=\beta\omega_Bx^*,
}
\tag{11}
\]

where \(x_i=c_i-\chi n_i\). Every competitive steady state fails these optimality conditions when \(q<1\).

**Full-optimum housing.** Suppose the candidate

\[
d_i^*=\mathbf1\{\xi_i\ge0\},\qquad
h_i^*=\kappa n_i+\widehat s,\qquad h_i^{2*}=\widehat h^2
\tag{12}
\]

respects its tenure-specific caps almost everywhere. Then it is the full planner’s housing and tenure solution, and

\[
\boxed{
H_Y^*=\kappa N\bar n+
\frac{\alpha}{\alpha+\beta\gamma}
(\bar H-\kappa N\bar n).
}
\tag{13}
\]

Furthermore, if every competitive old housing cap is slack and

\[
\boxed{
(\gamma+\omega_B)(1-q+q\tau^p)
>q\gamma(1+q\tau^p),
}
\tag{14}
\]

then \(H_Y^*>H_Y^{\rm eq}\). No restriction on \(\beta>0\), or assumption of binding young finance, is needed.

### Proof

The goods and housing problems separate because fertility is fixed and \(e_i\) is a real payment independent of physical housing. The goods first-order conditions are

\[
\frac1{x_i}=\frac{\beta}{c_i^2}
=\frac{\beta\omega_B}{e_i}=\lambda_G.
\]

Their resource constraint gives (11).

At the competitive allocation, let \(m_i=1/c_i^2\). The old value function’s derivative with respect to resources is \(m_i\), including when old housing or financial-saving restrictions bind. The reduced young problem therefore implies

\[
\boxed{\frac1{x_i}=\frac{\beta}{q}m_i+\mu_i,\qquad \mu_i\ge0.}
\tag{15}
\]

For renters, the reduced cash requirement uses \(L_R=p\); for owners, \(L_O=L\).

Hence

\[
\frac1{x_i}-\frac{\beta}{c_i^2}
=\beta(q^{-1}-1)m_i+\mu_i>0.
\tag{16}
\]

Increase young consumption by \(\varepsilon\), decrease its stationary old counterpart’s consumption by \(\varepsilon\), and leave housing, estates, and financial positions unchanged. Balanced current transfers implement this variation. It preserves (9), and its welfare derivative is (16). Selecting a positive-measure subset with uniform positive margins gives a finite improvement. Thus the competitive allocation cannot solve the full planner.

For housing, first remove the caps. Strict concavity gives common adult young space \(\widehat s\) and common old housing \(\widehat h^2\). Independently, the ownership taste is maximized by \(d_i=\mathbf1\{\xi_i\ge0\}\). If this unconstrained solution satisfies the original caps, it solves the original full problem, proving (12)–(13).

For the aggregate comparison, competitive housing conditions give

\[
\frac{\alpha}{s_i}
=\frac{\beta}{q}pm_i+\mu_iL_{d_i}+\eta_i^y.
\tag{17}
\]

For uncapped old renters, \(\gamma/h_i^2=pm_i\). For uncapped old owners, the two financial-estate regimes derived in the appendix imply that (14) gives

\[
\frac{\gamma}{h_i^2}<\frac pq m_i.
\]

Thus \(\alpha/s_i>\beta\gamma/h_i^2\) for every type. Finally,

\[
\boxed{
H_Y^*-H_Y^{\rm eq}
=\frac{N}{\alpha+\beta\gamma}
\int(\alpha h_i^2-\beta\gamma s_i)\,dQ>0.
}
\tag{18}
\]

This proves an ordering of the **full optimum**, not merely an improving direction. ∎

The cap-feasibility condition in (12) is substantive. A stronger, easily checked sufficient condition is

\[
h_R^{\max}\ge
\max\{\kappa\,\operatorname*{ess\,sup}_i n_i+\widehat s,\widehat h^2\}.
\]

With active planner caps, tenure and housing must instead be optimized jointly through the capped menus. Equation (18) cannot then be used with the uncapped formula.

## 5. What the result means—and does not mean

For the proposed matched-owner variation,

\[
\Delta h=\varepsilon,\quad \Delta h^2=-\varepsilon,\quad
\Delta a'=-P\varepsilon,\quad \Delta a^e=qP\varepsilon,
\tag{19}
\]

with transfers \(+p\varepsilon\) and \(-p\varepsilon\), entering-old resources and estates remain unchanged. Moreover,

\[
\Delta B^{\rm ext}
=\Delta a'+q^{-1}\Delta a^e=0.
\]

Thus the candidate variation is feasible under this completion.

On the uncapped positive-financial-estate branch, its derivative is indeed

\[
\frac{\alpha}{s}-\beta\frac{\gamma}{h^2}
=\beta(q^{-1}-1)pm+\mu L.
\tag{20}
\]

The first term survives without a mortgage distortion. The appendix shows why an additional negative term appears when old financial saving is zero.

A same-resource frictionless diagnostic makes the distinction transparent: whenever financial restrictions and housing caps are inactive, competitive first-order conditions imply

\[
H_Y^{F}
=\kappa N\bar n+
\frac{\alpha q}{\alpha q+\beta\gamma}S
<
\kappa N\bar n+\frac{\alpha}{\alpha+\beta\gamma}S.
\tag{21}
\]

The stationary cohort criterion already produces an age-allocation gap. Equal weights additionally generate ordinary redistribution across heterogeneous households. Neither effect should be relabeled a mortgage inefficiency.

These are failures to maximize (1). They are **not** Pareto-improvement or constrained-policy results. In particular, the initial old do not receive the young-age benefit of a permanent reallocation. A common inherited state also fixes outstanding estate payments and asset/title distributions—not merely aggregate \(B^*\). Those restrictions belong in a transition comparison.

## Essential appendix

### A. Financial reconstruction of the full feasible allocation

Given any allocation satisfying (9), set \(a_i^e\) by (5), and choose young financial positions with aggregate

\[
N\int a_i'dQ
=B^*-E+P(\bar H-H_Y^O),
\qquad
H_Y^O=N\int d_i h_i\,dQ.
\tag{22}
\]

This is exactly (6). Distribute that aggregate across types arbitrarily; individual financing restrictions have been relaxed.

Define age/type transfers as the residuals in the original household budgets, retaining the reference rebate \(T\). Summing those residuals, including the intermediary accounts, gives

\[
\int N(t_i^y+t_i^o)dQ
=C+E-\bigl[Y^g+I+(1-q)B^*\bigr]=0.
\]

Thus every allocation in (9) has a stationary financial implementation with balanced domestic transfers. Negative \(a_i^e\) represents estate-secured borrowing; the delivered estate remains \(e_i>0\). Changing the bookkeeping price alone creates neither goods nor estate utility.

### B. Both old-owner regimes, and a binding-finance test

Put

\[
d_p=1-q+q\tau^p,\qquad a_p=1+q\tau^p,\qquad J=\gamma+\omega_B.
\]

For an uncapped old owner,

\[
c^2=\frac zK,\qquad h^2=g_Oz,
\]

where

\[
g_O=
\begin{cases}
\dfrac{\gamma}{Kp},&
\omega_Bd_p\ge q\gamma,\\[6pt]
\dfrac{J}{Ka_pP},&
\omega_Bd_p<q\gamma.
\end{cases}
\tag{23}
\]

In the first case \(e=\omega_Bz/(Kq)\); in the second, \(e=Ph^2\) and \(a^e=0\).

Let \(\zeta\) be the multiplier on \(e-Ph^2\ge0\). On the regular zero-financial-estate branch,

\[
\zeta=m\,\frac{q\gamma-\omega_Bd_p}{J}>0,
\qquad
\frac{\gamma}{h^2}=pm+P\zeta.
\]

Consequently, with both housing caps slack,

\[
\boxed{
\frac{\alpha}{s}-\beta\frac{\gamma}{h^2}
=\beta Pm\left(\frac{d_p}{q}-\frac{\gamma a_p}{J}\right)+\mu L.
}
\tag{24}
\]

A strictly positive mortgage multiplier alone need not dominate the negative first term. An active old housing cap adds the further term \(-\beta\eta^o\).

**Finance need not be assumed binding.** Let

\[
M=w+qv,\qquad D=1+\alpha+\vartheta+\beta K.
\]

The sufficient cap bound

\[
h_O^{\max}>\max\{M/p,\;g_OM/q\}
\tag{25}
\]

makes both owner caps slack in the actual problem and its cash-unconstrained comparison. In both old-estate regimes, the uncapped continuation value is \(K\log z+\text{constant}\). Removing only the young cash restriction therefore gives cash expenditure

\[
\mathcal C^0
=\frac{M}{D}
\left[
1+\alpha\frac Lp+
\vartheta\frac{\chi+L\kappa}{\chi+p\kappa}
\right].
\tag{26}
\]

If \(\mathcal C^0>w\), strict concavity implies \(\mu>0\). This uses the original, endogenous-fertility household problem to certify its reference choice; fertility is frozen only afterward.

For a simpler income test, set

\[
\ell=\frac{1-\phi+q\tau^p}{d_p},\qquad
C_{\min}=1+\alpha\ell+\vartheta\min\{1,\ell\}.
\]

Then

\[
\frac{qv}{w}>\frac{D}{C_{\min}}-1
\tag{27}
\]

suffices for binding finance. This explains why the note’s income test extends to the zero-financial-estate branch, while its old-housing cap coefficient must change to \(g_O\). Pasted text

Even when (14) fails, a sufficient income condition for a **positive housing derivative** on that branch is

\[
\frac vw>
\frac{(\phi/q-1)_+}{1-\phi+q\tau^p}
+\frac{\beta\gamma K a_p}
{\alpha J(1-\phi+q\tau^p)}.
\tag{28}
\]

Indeed, the cash constraint gives \(s<w/L\), while mortgage repayment gives

\[
z\ge v-\frac{(\phi/q-1)_+w}{1-\phi+q\tau^p}.
\]

Substituting these bounds into \(\alpha/s-\beta\gamma/(g_Oz)\) proves (28). Thus a borrowing-related local result can be generated by income restrictions rather than assumed marginal-value gaps.

### C. An analytical reversal of the aggregate housing conclusion

Suppose

\[
Jd_p<q\gamma a_p,
\tag{29}
\]

young financial constraints are slack, and competitive housing caps are inactive. Owners then have zero financial estate. Conditional young choices are identical across tenures because the old value functions differ only by a constant. Hence ownership has a common probability \(\pi\), despite heterogeneous endowments, and

\[
s_i=\frac{\alpha M_i}{Dp},\qquad
h_i^{2R}=\frac{\beta\gamma M_i}{Dqp},\qquad
h_i^{2O}=\frac{\beta J M_i}{Dqa_pP}.
\tag{30}
\]

When the planner candidate (12) is feasible, the sign of \(H_Y^*-H_Y^{\rm eq}\) is the sign of

\[
\frac{\gamma(1-q)}{d_p}
+\pi\left(\frac{J}{a_p}-\frac{\gamma}{d_p}\right).
\]

Therefore,

\[
\boxed{
\pi>
\frac{\gamma(1-q)a_p}{\gamma a_p-Jd_p}
\quad\Longrightarrow\quad
H_Y^*<H_Y^{\rm eq}.
}
\tag{31}
\]

Under (29), the ownership threshold lies strictly between zero and one.

This is a nonempty analytical family, not a numerical example. For instance, take \(\tau^p=0\), \(q>1/2\), and

\[
0<\omega_B<\frac{\gamma(2q-1)}{1-q}.
\]

Any \(\beta>0\) is permitted. Choose bounded heterogeneous endowments satisfying

\[
\phi\ge q,\qquad
\frac{qv_i}{w_i}<\frac{\beta K}{1+\alpha+\vartheta},
\]

which makes young finance slack in both tenures. Choose sufficiently large finite caps and a sufficiently high finite logistic location to satisfy (31), retaining both tenures with positive mass. Positive child costs can be rescaled together to obtain \(\nu\bar n=1\), without changing goods or space used by children. Population is then taken from this reference equilibrium and held fixed in the planner comparison.

The reversal reflects the old owner’s inability to borrow against housing: relaxing that restriction can make the full planner shift housing **toward old age**.

## Decisions that remain

**Normative:** approve equal weights on stationary cohort lifetime utility. A transition analysis separately requires social generation discounting and an initial-old weight; dynastic literature does not determine them here.

**Accounting:** approve outside entrant remittances and estate recipients, domestic zero-equity rental intermediaries, real estate payments, and the pre-remittance definition of fixed external wealth. Internal inheritance would require a different, explicitly closed recipient account.

**Household behavior:** no change is needed for the recommended benchmark. Standard mortgages, heterogeneous endowments, positive child costs, tenure persistence, and warm glow remain. The planner relaxes financial restrictions—not physical caps. Genuine dynastic altruism would be a separate model change, not a welfare-weight convention.
