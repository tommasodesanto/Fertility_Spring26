# Positive child costs and stationary population

## Result and scope

There is a sufficient condition for the original mixed-tenure model with genuinely positive child goods costs, homogeneous entrants and zero property tax. It compares an upper bound on the extra mean young housing needed for replacement with a lower bound on old housing released. The condition contains household ratios and primitives, not an equilibrium derivative or a borrowing multiplier. It covers an exact economy with \(\chi=3/20\), owner share \(11/21\), and every positive ownership-taste scale. A separate rational example with owner share \(33/83\) proves that credit can instead lower stationary population, even when old owners occupy more than capped renters and every original household restriction is strict.

The economic sentence is: **credit raises population when old-owner downsizing releases more space than young households need to maintain replacement fertility.**

These are stationary results. All original utilities, values, budgets, mortgage repayment and estate dating remain fixed. No transition, planner or new instrument is introduced. The exact positive-cost result sets \(\tau^p=0\); positive tax is covered only locally by continuity around a strict example, with no quantified tax neighborhood.

## Maintained original conditions

Use the original homogeneous entrant income and liquid wealth \((y,b)\), the fixed rental cap \(a=h_R^{\max}\), and positive stationary replacement fertility \(1/\nu\). Write \(h=h^O\), \(h_2=h^{2,O}\), \(n_O,n_R\), \(x_m=c_m-\chi n_m\), and \(s_m=h_m-\kappa n_m\). Maintain
\[
 a<h<h_O^{\max},\quad 0<h_2<h,\quad e_O>P h_2,
 \quad a'_O,a'_R>0,
\]
\[
 \alpha x_O/s_O>u,\quad\alpha x_R/s_R>u,
 \quad \beta\gamma x_R/(qa)>u,
\]
with positive goods, adult space, fertility and estates. Thus owners are strictly constrained by their down payment, renters are strictly capped at both ages, and old owner retention and estate bounds are slack. In addition assume
\[
 n_O>n_R,\qquad h_2\ge a,\qquad0<\pi<1.
\]
All weights, \(\chi,\kappa,\nu\) are positive, \(0<q,\phi<1\), and the original logistic parameters \((\bar\xi,\sigma_\xi)\), \(\sigma_\xi>0\), are held fixed during each comparison.

Source equations: [original household problems](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:104), [mixed proof conditional budgets and FOCs](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/mixed_transition_proof.md:50), and [stationary population](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:364). The earlier zero-cost result is [transition extensions, section 2](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/transition_extensions.md:33).

## Conditional fertility and tenure derivatives

Put
\[
 \ell=1-q,\quad d=b/(1-\phi),\quad w=y+b,\quad
 \rho_O=1+\beta(1+\gamma+\omega_B),\quad
 \rho_R=1+\beta(1+\omega_B).
\]
The original stationary budgets and FOCs are
\[
 P=d/h,\quad u=\ell P,\quad
 \rho_Ox_O+\chi n_O=w-\ell d,\quad
 \rho_Rx_R+\chi n_R=w-(1+q)ua,
\]
\[
 \vartheta=\frac{\chi n_m}{x_m}
 +\frac{\alpha\kappa n_m}{s_m},\qquad
 h_2=\frac{\beta\gamma x_O}{qu}.
\]
Define the child goods and space ratios and a positive curvature quantity:
\[
 g_m=\frac{\chi n_m}{x_m},\quad t_m=\frac{\kappa n_m}{s_m},\quad
 K_m=g_m+\alpha t_m(1+t_m)+g_m^2/\rho_m>0.
\]
For \(R_m=\rho_mx_m+\chi n_m\), differentiating these two original equations gives the useful identity
\[
 d\log n_m=
 \frac{\alpha t_m(1+t_m)}{K_m}\,d\log h_m
 +\frac{g_m}{\rho_m K_m}\frac{dR_m}{x_m}
 +\frac{d\vartheta}{K_m}.
\]
For \(H=d\log h\) and \(D=d\log d\), write
\[
 m=\ell d/x_O,\quad v=(1+q)ua/x_R,\qquad
 \eta=\frac{\alpha t_O(1+t_O)}{K_O},\quad
 \lambda_O=\frac{g_Om}{\rho_OK_O},\quad
 \lambda_R=\frac{g_Rv}{\rho_RK_R}.
\]
Then
\[
 d\log n_O=\eta H-\lambda_O D+d\vartheta/K_O,
 \qquad
 d\log n_R=\lambda_R(H-D)+d\vartheta/K_R. \tag{1}
\]
Here \(\eta\) is the owner fertility response to housing at fixed lifetime resources. The two \(\lambda\)'s measure fertility responses to the housing-payment change. They are explicit functions of the displayed original household ratios.

The original value envelope, including the owner's later rental-service cost, gives
\[
 d\Delta=A H-FD+L\,d\vartheta,\qquad \Delta=W^O-W^R,
\]
\[
 A=\alpha h/s_O+\beta\gamma-v,\quad
 F=m+\beta\gamma-v,\quad L=\log(n_O/n_R)>0. \tag{2}
\]
The old renter cap and owner retention imply \(F>0\), exactly as in the zero-cost proof. Strict owner purchase gives \(A-F=\alpha h/s_O-m>0\). Thus
\[
 A>F>0
\]
for general positive child costs, without a new preference restriction. Logistic tenure choice gives \(d\pi=k\,d\Delta\), where \(k=\pi(1-\pi)/\sigma_\xi>0\).

## Replacement and a transparent credit condition

Let
\[
 U=\pi n_O\eta+(1-\pi)n_R\lambda_R,\quad
 V=\pi n_O\lambda_O+(1-\pi)n_R\lambda_R,\quad
 W=\pi n_O/K_O+(1-\pi)n_R/K_R.
\]
Differentiating replacement and substituting (1)–(2) gives the complete scalar stationary equation
\[
 [U+(n_O-n_R)kA]H
 =[V+(n_O-n_R)kF]D-[W+(n_O-n_R)kL]d\vartheta. \tag{3}
\]
Its coefficient is strictly positive. Thus the stationary root is locally unique. Mean fertility is strictly increasing in \(h\), hence decreasing in \(P=d/h\), throughout a connected branch with the maintained fertility gap. There is at most one root on such a branch.

The owner-housing elasticity with respect to credit expenditure is a weighted average of two transparent ratios:
\[
 E=\frac{d\log h}{d\log d}
 =\frac{U E_0+(n_O-n_R)kA E_1}{U+(n_O-n_R)kA},
 \qquad E_0=V/U,\quad E_1=F/A\in(0,1). \tag{4}
\]
\(E_0\) is the size response needed for replacement at a fixed tenure share; \(E_1\) holds the ownership value difference fixed. The first credit restriction is
\[
 E_0\le E_1. \tag{C1}
\]
It implies \(0<E\le E_1<1\), so the stationary price rises and the owner share weakly falls.

For a bound on total housing, define only household ratios:
\[
 C=h_2/h\in(0,1),\quad \delta=a/h,\quad r=n_R/n_O<1,\quad
 j=g_O/\rho_O,\quad
 \varepsilon_R=\frac{1-\pi}{\pi}r\lambda_R,
\]
\[
 \Gamma=1-j\eta<1,\quad
 D_0=1+m/\rho_O-j\lambda_O>1,\quad
 \Lambda=\frac{1-\delta}{1-r}.
\]
The inequality \(D_0>1\) follows from \(g_O^2/(\rho_OK_O)<1\). The symbol \(D_0\) is a coefficient, distinct from the differential \(D\). Let \([z]_+=\max\{z,0\}\). The following bounds measure housing changes per \(\pi h\) and per unit increase in \(\log d\). A sufficient second restriction is
\[
 \boxed{\quad
 \underbrace{\Lambda(\lambda_O+\varepsilon_R)
  +[1-\Lambda(\eta+\varepsilon_R)]_+E_1}_{\text{upper bound on extra young housing}}
 <
 \underbrace{C[D_0-[\Gamma]_+E_1]}_{\text{lower bound on old housing released}}.
 \quad} \tag{C2}
\]
Conditions (C1)–(C2) contain no taste scale or equilibrium derivative. Together with the maintained original conditions, they imply
\[
 P_\phi^*>0,\qquad N_{{\rm hh},\phi}^*>0.
\]
They are sufficient bounds, not necessary conditions or a renamed determinant sign.

**Proof of the bound.** A dot denotes \(d/d\log d\). The owner budget gives
\[
 \dot{\log h_2}=\Gamma E-D_0.
\]
Replacement gives
\[
 \dot\pi/\pi=
 \frac{\lambda_O+\varepsilon_R-(\eta+\varepsilon_R)E}{1-r}\le0.
\]
Consequently, mean housing at the two ages satisfies the exact decomposition
\[
 \frac{\dot{\bar h}^{Y}}{\pi h}
 =\Lambda(\lambda_O+\varepsilon_R)
  +[1-\Lambda(\eta+\varepsilon_R)]E,
\]
\[
 \frac{\dot{\bar h}^{O}}{\pi h}
 =C(\Gamma E-D_0)
 +(C-\delta)\frac{\lambda_O+\varepsilon_R-(\eta+\varepsilon_R)E}{1-r}.
\]
Since \(C\ge\delta\), the composition term in old housing is nonpositive. Bounding \(E\) by \(E_1\) gives precisely the upper and lower bounds in (C2). Thus \(S=\bar h^Y+\bar h^O\) falls and \(N_{\rm hh}^*=2\bar H/S\) rises. This separates within-tenure housing, goods-induced fertility changes, and tenure changes.

## A sufficient fertility-weight condition

Write \(H_\vartheta=\partial\log h/\partial\vartheta\). Equation (3) immediately gives \(H_\vartheta<0\) and hence \(P_\vartheta>0\). The owner-share response has no assigned sign. For population, one additional sufficient allocation-ratio condition is
\[
 \mathcal T\equiv1+C\Gamma
 -\frac{1+C-2\delta}{1-r}(\eta+\varepsilon_R)\ge0. \tag{C3}
\]
This bounds the fertility response from enlarging owner housing relative to the fertility difference between tenures. To verify the sign, put
\[
 \Omega=1/K_O+\frac{1-\pi}{\pi}\frac{r}{K_R}>0.
\]
Directly differentiating total housing and eliminating the tenure change with replacement yields
\[
 \frac{S_\vartheta}{\pi h}
 =\mathcal T H_\vartheta-Cj/K_O
 -\frac{1+C-2\delta}{1-r}\Omega<0. \tag{5}
\]
Thus (C3) implies \(N_{{\rm hh},\vartheta}^*>0\). A fall in \(\vartheta\) lowers stationary price and population. At zero child cost, \(\mathcal T=\delta(1-C)/(1-\delta)>0\), so this condition retains the earlier limiting result. It is sufficient and need not be necessary; an unrestricted positive-cost fertility-weight theorem is not proved here.

## Exact positive-cost family

Keep the original positive-cost allocations from [mixed proof, section 6](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/mixed_transition_proof.md:200), but set tax to zero and specify income and liquid wealth accordingly:
\[
 q=\tfrac12,\ \phi=\tfrac45,\ b=\tfrac{9717}{46250},\ y=\tfrac{779829}{370000},
 \quad\alpha=\beta=\omega_B=\tfrac25,\ \gamma=\tfrac3{10},
\]
\[
 \chi=\tfrac3{20},\ \kappa=\tfrac12,\ \vartheta=\tfrac{141}{400},\ \nu=2,
 \quad a=\tfrac14,\ h_O^{\max}=2,\ \tau^p=0,\ \bar H=\tfrac{68104}{68019}.
\]
There is an exact stationary allocation
\[
 P=\tfrac{9717}{9250},\quad Y=O=1,\quad\pi=\tfrac{11}{21},
\]
\[
 (x_O,h,n_O)=(1,1,\tfrac34),\quad
 (x_R,a,n_R)=(\tfrac{99}{74},\tfrac14,\tfrac9{40}),\quad
 h_2=\tfrac{1480}{3239}>a.
\]
For each \(\sigma_\xi>0\), choose the primitive location once as
\(\bar\xi=\sigma_\xi\log(11/10)-W^O+W^R\) evaluated at these allocations. Hold both taste parameters fixed in every comparison. This is a family of primitive economies, not a policy rule that reselects tastes.

All three conditions hold with strict margins, computed as exact rational numbers in the companion check:
\[
 E_0=0.09284084\ldots< E_1=0.81272542\ldots,
\]
\[
 \text{young upper bound}=0.22429667\ldots
 <0.24523919\ldots=\text{old lower bound},\quad
 \mathcal T=0.39070122\ldots>0.
\]
The (C2) margin is exactly
\(12376832748143069677/590990492246573075475>0\).
Owner purchase, retention, estate, physical-cap, and both saving margins are strictly positive; the smallest displayed estate margin is \(4/25\). The original two renter cap inequalities also hold strictly. This establishes the signs for every positive taste scale in a family with material renting and child costs, not just continuity from \(\chi=0\).

At \(\sigma_\xi=1\), the original stationary equations give
\[
 N_\phi=2.119588329\ldots,\quad P_\phi=3.966238469\ldots,
 \qquad N_\vartheta=4.723781115\ldots,\quad P_\vartheta=3.068095745\ldots.
\]

At this same point, \(\bar h^Y_\phi=0.315866236\ldots\) and \(\bar h^O_\phi=-1.376984773\ldots\): old housing falls by more than young housing rises.

## Exact counterexample with material renting

Take
\[
 q=\phi=\tfrac45,\ b=\tfrac15,\ y=\tfrac{199}{62},\quad
 \alpha=\tfrac25,\ \beta=\tfrac1{10},\ \gamma=\tfrac15,
 \ \omega_B=\tfrac{139}{155},
\]
\[
 \chi=2,\ \kappa=\tfrac12,\ \vartheta=\tfrac{12}{5},\ \nu=2,
 \quad a=\tfrac1{10},\ h_O^{\max}=2,\ \sigma_\xi=1,
 \quad\tau^p=0,\ \bar H=\tfrac{377}{664}.
\]
Specify \(\bar\xi=\log(33/50)-W^O+W^R\) once at the following exact allocation:
\[
 P=Y=O=1,\quad\pi=\tfrac{33}{83},\quad
 (x_O,h,n_O)=(1,1,1),\quad
 (x_R,a,n_R)=(\tfrac{51}{20},\tfrac1{10},\tfrac{17}{100}).
\]
Here \(h_2=1/8>a\), \(c_O^2=1/8\), \(e_O=139/992\),
\(a'_O=13/62\), and \(a'_R=1549/3100\). Owner purchase, retention, estate, and physical-cap margins are respectively
\(3/5,7/8,15/992,1\). Young and old renter cap margins are \(339/5\) and \(7/16\). All goods, space and fertility are positive. Original household feasibility is therefore strict, including the estate and saving restrictions.

Exact rational differentiation gives
\[
 N_\phi=-\frac{4083925586718912}{2831234115298355}
 =-1.442454216\ldots<0,
\]
while \(P_\phi=3.528796622\ldots>0\). In this economy
\(E_0=0.413452913\ldots>E_1=0.255474453\ldots\): the credit comparison raises the owner share. The conditional fertility of both tenures falls, so replacement requires more households in the higher-fertility owner tenure. Larger young homes and more ownership use more housing than old downsizing releases. Here \(\bar h^Y_\phi=0.618608011\ldots\) but \(\bar h^O_\phi=-0.209116114\ldots\). The reversal occurs with the same binding constraints, material renting, and old-owner housing above the rental cap.

The same counterexample satisfies (C3), with \(\mathcal T=0.951897157\ldots\) and \(N_\vartheta=2.992338647\ldots>0\). Its reversal concerns the credit comparison.

## Verification and limits

The companion script is [/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/verify_simplified_olg_positive_costs.py](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/verify_simplified_olg_positive_costs.py). It uses exact fractions for the conditions and credit derivatives, then imports the [original helper](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/verify_simplified_olg_local_transition.py:56) without executing its main. Run with `PYTHONDONTWRITEBYTECODE=1`. It verifies the original dated budgets, fertility FOCs, all branch restrictions, and separately optimizes both original household problems in both examples.

Maximum original budget errors are below \(7\times10^{-16}\), FOC errors below \(2\times10^{-15}\), and optimization choice discrepancies below \(7\times10^{-7}\). Central differences with step \(10^{-6}\), keeping both taste parameters fixed, match the analytical population and price derivatives within \(9\times10^{-10}\). The counterexample sign and the positive family's conditions use exact rational arithmetic; the numerical roots are supporting checks.

The sufficient conditions are conservative and need not describe every positive-population response. They require homogeneous entrants, the original strict branches, \(n_O>n_R\), and zero tax for the direct proof. They provide explicit positive-child-cost examples, but not a certified large policy-change interval. Signs extend to finite comparisons along an interval that preserves the conditions. Strict signs and regularity persist under sufficiently small positive taxes and compact entrant heterogeneity with uniform original margins; their size is not quantified. An unrestricted positive-cost population theorem is false for credit, as the counterexample shows. The weakest sufficient fertility-weight condition remains open. This supporting proof does not amend the main paper, slides, quantitative model, or author decisions.
