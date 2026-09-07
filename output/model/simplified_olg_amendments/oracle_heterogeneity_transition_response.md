Verdict

A useful theorem can be proved with genuine income–wealth heterogeneity, positive property taxes, mixed tenure, and renters who remain below replacement. The appropriate equilibrium claim is:

Explicit aggregate fertility brackets guarantee at least one level-stationary equilibrium with a positive-mass compensated housing improvement. Every stationary equilibrium inside the bracket has that improvement; equilibria outside it are not ruled out.

The key improvement over the previous argument is to separate young recipients from old donors. Young recipients need sufficiently scarce liquid wealth relative to their demand for current housing. Old donors need sufficiently much inherited housing relative to their old-age resources. These can be different income–wealth groups, with different masses. There is no reason to impose the earlier two-sided liquidity restriction on every household.

For transitions, I obtain a dated theorem using actual inherited assets and titles, including positive taxes, and a separate sufficient condition that propagates donor eligibility from earlier decisions. I do not establish existence or convergence of the complete post-decline/post-reform equilibrium paths. There is also a concrete obstruction to an unrestricted global transition theorem: sufficiently large fertility-preference declines are incompatible with an infinite equilibrium path under the model’s fixed stock, finite housing caps, and exact housing-clearing requirement.

Throughout, entrant wealth \(b\) remains exogenous. Estates are preserved in the welfare comparison but do not determine the next entrant distribution, exactly as specified in the packet.

1. The two household results that make the argument work

Write

$$ A=1+\gamma+\omega_B,\qquad M=1+\beta A. $$

At fixed fertility and estates, the relevant consumption-valued housing benefits are

$$ MV^Y=\frac{\alpha x}{h-\kappa n}, \qquad MV^O=\frac{\gamma c^2}{h^2}, \qquad x=c-\chi n. $$

The young expression is a current-flow marginal value. A positive down-payment multiplier alone would not establish the required comparison: part of that multiplier might reflect the value of carrying a larger house into old age.

1.1 Old owners: an exact solution and an exact donor cutoff

At date \(t\), define

$$ u_t=(1+q\tau^p)P_t-qP_{t+1}>0, \qquad Z=a+P_tH+T_t, $$

and

$$ v_t= \max\left\{ u_t,\, \frac{\gamma}{\gamma+\omega_B}(1+q\tau^p)P_t \right\}. $$

The old owner’s solution is

$$ h^2=\min\left\{H,\frac{\gamma Z}{Av_t}\right\}, \tag{1} $$ $$ c^2= \min\left\{ \frac{Z-u_th^2}{1+\omega_B}, \ Z-(1+q\tau^p)P_th^2 \right\}, \qquad e=\frac{Z-u_th^2-c^2}{q}. \tag{2} $$

These formulas solve the packet’s actual old-owner problem, including the estate floor \(e\ge P_{t+1}h^2\).

In particular,

$$ Z<\frac{Av_tH}{\gamma} \quad\Longrightarrow\quad h^2<H,\qquad c^2=\frac ZA,\qquad MV^O=v_t. \tag{3} $$

Thus the donor value need not equal user cost. When the financial-estate constraint binds, \(v_t>u_t\).

There is a useful stronger characterization. For any cutoff \(C\ge v_t\), define

$$ \mathcal Z_t(C)= \max\left\{ u_t+\frac{1+\omega_B}{\gamma}C,\, (1+q\tau^p)P_t+\frac C\gamma \right\}. $$

Then

$$ \boxed{ MV^O\le C \quad\Longleftrightarrow\quad \frac ZH\le\mathcal Z_t(C). } \tag{4} $$

This also admits donors whose retention constraint binds. At \(C=v_t\),

$$ \mathcal Z_t(v_t)=\frac{Av_t}{\gamma}. $$

Derivation. Without the retention cap, the estate-slack solution has housing cost \(u_t\). With a binding estate floor, housing and the required estate jointly cost \((1+q\tau^p)P_t\), with combined utility weight \(\gamma+\omega_B\). Both branches give \(c^2=Z/A\), producing (1). Conditional on retained housing, optimal consumption is (2). When retention binds, \(MV^O=\gamma c^2/H\); substituting (2) gives (4).

1.2 Young owners: a lower bound that does not assume future retention is slack

The rebate-inclusive transformation is

$$ z=a'+P_{t+1}h+T_{t+1}, \qquad W_t=y+b+T_t+qT_{t+1}. $$

The young budget becomes

$$ x+\chi n+qz+u_th=W_t, \tag{5} $$

while actual gross bond holdings are

$$ qa'+\phi_tP_th = qz-qT_{t+1}+(\phi_tP_t-qP_{t+1})h. \tag{6} $$

This distinction matters when taxes and future rebates are positive.

Temporarily omit the gross-bond constraint. Let \(\mathcal V(z,h)\) denote the future old-owner value, expressed in total resources and inherited housing. Its joint homogeneity gives

$$ z\mathcal V_z+h\mathcal V_H=A, \qquad \mathcal V_H\ge0. $$

The saving first-order condition therefore implies

$$ qz=\beta x(A-h\mathcal V_H)\le\beta Ax. \tag{7} $$

Also, \(\mathcal V_z=1/c^2\), so \(c^2=\beta x/q\) and \(qz>\beta x\).

For \(m,E,H>0\), define \(\mathcal N_m(E,H)\) by

$$ mx+\chi n=E,\qquad \frac{\vartheta}{n} =\frac{\chi}{x}+\frac{\alpha\kappa}{H-\kappa n}. \tag{8} $$

Its explicit solution is

$$ \mathcal N_m(E,H)= \frac{2\vartheta EH} {K_m+\sqrt{K_m^2- 4\chi\kappa(m+\alpha+\vartheta)\vartheta EH}}, $$

where

$$ K_m=\kappa E(\alpha+\vartheta)+\chi H(m+\vartheta). \tag{9} $$

It increases strictly in \(E\) and \(H\).

At fixed current housing \(h=H\), put

$$ E=W_t-u_tH>0,\qquad n^{\mathrm{bench}}=\mathcal N_M(E,H),\qquad x^{\mathrm{bench}}=\frac{E-\chi n^{\mathrm{bench}}}{M}. $$

Then the actual optimum of the gross-bond-relaxed problem satisfies

$$ x\ge x^{\mathrm{bench}},\qquad n\ge n^{\mathrm{bench}}, $$

and hence

$$ \boxed{ MV^Y\ge \frac{\alpha x^{\mathrm{bench}}} {H-\kappa n^{\mathrm{bench}}}. } \tag{10} $$

To prove this, (5)–(7) imply \(E\le Mx+\chi n\). At fixed \(H\), the fertility first-order condition makes \(n\) strictly increasing in \(x\). The benchmark is the unique point satisfying equality.

This is the main tightening: the benchmark supplies a lower bound even when the young owner’s future retention constraint binds. Recipient eligibility therefore does not require a lower liquidity bound that guarantees future downsizing.

2. Main proposition: heterogeneous stationary existence, including positive taxes

I first state all the restrictions. They involve primitives, a specified price interval, and elementary functions evaluated at that interval—not unknown equilibrium multipliers or an assumed equilibrium constraint pattern.

Let

$$ \rho=h_R^{\max},\qquad \bar h=h_O^{\max}, $$

and impose the original positive preference and child-cost parameters,

$$ 0<q<1,\quad q\le\phi<1,\quad \tau^p\ge0,\quad 0<\rho<\bar h<\infty,\quad \bar H>0. $$

Ownership tastes have finite location \(\bar\xi\) and finite positive scale \(\sigma_\xi\). Let \(F\) be a probability distribution with compact support in \((0,\infty)^2\), with nondegenerate income and wealth heterogeneity.

Set

$$ w_i=y_i+b_i,\qquad B_i=\frac{b_i}{1-\phi}, $$ $$ d=1-q+q\tau^p,\qquad v=\max\left\{d,\frac{\gamma}{\gamma+\omega_B}(1+q\tau^p)\right\}, \qquad D=M+\alpha+\vartheta. $$

Choose primitive price endpoints

$$ 0<P_{\mathrm{low}}<P_{\mathrm{high}}<\infty. $$

Define

$$ T_{\mathrm{high}}=q\tau^pP_{\mathrm{high}}\bar h, \qquad w_i^{\mathrm{high}}=w_i+(1+q)T_{\mathrm{high}}, \qquad E_i=w_i-dB_i. \tag{11} $$

For types with \(E_i>0\), define

$$ h_i^{\mathrm{dp}}(P)=\frac{B_i}{P}, $$ $$ n_i^{\mathrm{bench}}(P)= \mathcal N_M\!\left(E_i,\frac{B_i}{P}\right), \qquad x_i^{\mathrm{bench}}(P)= \frac{E_i-\chi n_i^{\mathrm{bench}}(P)}M. \tag{12} $$
Explicit value bounds for the fertility bracket

Let \(\mathcal R(P,W)\) be the fully uncapped renter value, with both age-specific housing caps and the saving lower bound removed. It is evaluated at

$$ x_f=\frac WD,\qquad n_f=\frac{\vartheta x_f}{\chi+\kappa dP},\qquad h_f=\frac{\alpha x_f}{dP}+\kappa n_f, $$ $$ (c_f^2,h_f^2,e_f)= \left( \frac{\beta x_f}{q}, \frac{\beta\gamma x_f}{qdP}, \frac{\beta\omega_Bx_f}{q^2} \right). \tag{13} $$

Thus \(\mathcal R(P,W)=u^y(x_f+\chi n_f,h_f,n_f)+\beta u^o(c_f^2,h_f^2,e_f)\).

This is an upper bound, not an assertion about actual renter policies.

Let \(\mathcal V_P(Z,H)\) be the old-owner value from (1)–(2) at constant price \(P\). Define the feasible zero-rebate owner test

$$ \mathcal U_i^{\mathrm{test}}(P)= u^y\!\left( x_i^{\mathrm{bench}}+\chi n_i^{\mathrm{bench}}, h_i^{\mathrm{dp}},n_i^{\mathrm{bench}} \right) + \beta\mathcal V_P\!\left( \frac{\beta A x_i^{\mathrm{bench}}}{q},h_i^{\mathrm{dp}} \right). \tag{14} $$

All quantities on the right are evaluated at \(P\).

Writing \(\Lambda(a)=1/(1+e^{-a})\), put

$$ p_i^{\mathrm{low}}= \Lambda\!\left( \frac{ \mathcal U_i^{\mathrm{test}}(P_{\mathrm{low}}) -\mathcal R(P_{\mathrm{low}},w_i^{\mathrm{high}}) +\bar\xi }{\sigma_\xi} \right). \tag{15} $$

Finally, define the primitive housing floor

$$ h_i^{\mathrm{floor}}(P)= \min\left\{ \bar h,\frac{B_i}{P}, \frac{w_i}{D} \left[ \frac{\alpha}{dP} +\frac{\kappa\vartheta}{\chi+\kappa dP} \right] \right\}. \tag{16} $$
Proposition 1

Suppose there are measurable sets \(S_Y,S_O\) with positive \(F\)-mass satisfying the following conditions.

Financial feasibility. For every type in \(S_Y\cup S_O\),

$$ \boxed{ E_i>0,\qquad \frac{\beta E_i}{M+\vartheta}>qT_{\mathrm{high}}. } \tag{S} $$

Young recipients. For every \(i\in S_Y\),

$$ \boxed{ \frac{B_i}{P_{\mathrm{low}}}<\bar h,\qquad \frac{\alpha x_i^{\mathrm{bench}}(P_{\mathrm{low}})} {B_i-\kappa P_{\mathrm{low}} n_i^{\mathrm{bench}}(P_{\mathrm{low}})} >v. } \tag{Y} $$

Old donors. For every \(i\in S_O\),

$$ \boxed{ \frac{\beta\gamma}{qM} \left[ \frac{w_i^{\mathrm{high}}} {P_{\mathrm{low}}h_i^{\mathrm{floor}}(P_{\mathrm{low}})} -d \right]<v. } \tag{O} $$

Aggregate replacement bracket.

$$ \boxed{ \int_{S_Y} p_i^{\mathrm{low}} n_i^{\mathrm{bench}}(P_{\mathrm{low}})\,dF(i) >\frac1\nu, } \tag{B-low} $$

and

$$ \boxed{ \frac{\vartheta}{\kappa(\alpha+\vartheta)} \int \max\left\{ \rho,\, \min\left(\bar h,\frac{B_i}{P_{\mathrm{high}}}\right) \right\}\,dF(i) <\frac1\nu. } \tag{B-high} $$

Then:

There exists a level-stationary equilibrium with

$$ P\in(P_{\mathrm{low}},P_{\mathrm{high}}). $$

Every level-stationary equilibrium in this interval admits a positive-mass compensated Pareto improvement of the specified kind.

In every such equilibrium, selected young owners have a strictly binding down-payment constraint, a slack physical cap, and positive gross bond holdings. Selected old owners retain strictly less than their inherited housing and have \(MV^O=vP\). Their financial estate may be zero. No constraint pattern is imposed on renters, on unselected owners, or on the future retention decisions of selected young recipients.

The stationary population level is determined by housing clearing; it is not prescribed independently.

Proof
A. Gross saving and young housing choices

At any candidate stationary price and rebate,

$$ P\in[P_{\mathrm{low}},P_{\mathrm{high}}], \qquad 0\le T\le T_{\mathrm{high}}, $$

a relaxed owner has effective resources \(W_i=w_i+(1+q)T\).

Because \(Ph\le B_i\), (5)–(7) and the fertility first-order condition imply

$$ x>\frac{W_i-dPh}{M+\vartheta} \ge\frac{E_i}{M+\vartheta}. $$

Actual gross bonds satisfy

$$ qa'+\phi Ph=qz-qT+(\phi-q)Ph >\frac{\beta E_i}{M+\vartheta}-qT_{\mathrm{high}}>0. $$

Thus, for selected types, omitting the gross-bond constraint changes nothing.

Now consider the owner benchmark that drops future retention. Its optimized value as a function of current housing is concave. Its derivative has the sign of \(MV^{Y,\mathrm{bench}}-dP\). The actual owner has at least this current-flow marginal value, by (10), and an additional nonnegative benefit from relaxing future retention.

At \(h=B_i/P\), the normalized benchmark value

$$ \frac{MV_i^{Y,\mathrm{bench}}}{P} = \frac{\alpha x_i^{\mathrm{bench}}(P)} {B_i-\kappa Pn_i^{\mathrm{bench}}(P)} \tag{17} $$

increases with \(P\). Indeed, \(n_i^{\mathrm{bench}}\) decreases with \(P\), while its housing share
\(\kappa Pn_i^{\mathrm{bench}}/B_i\) increases; consequently the numerator increases and the denominator decreases. Increasing the rebate also increases the benchmark marginal value.

Condition (Y) therefore implies, throughout the interval,

$$ h_i^Y=\frac{B_i}{P}<\bar h, \qquad MV_i^Y>vP\ge dP. \tag{18} $$

The down-payment multiplier is strictly positive: its first-order condition has the strictly positive current-flow term \((MV_i^Y-dP)/x_i\), plus a nonnegative retention term.

B. Donor eligibility follows from their earlier optimal choices

The same housing argument establishes, for any selected donor type,

$$ h_i^Y\ge h_i^{\mathrm{floor}}(P). \tag{19} $$

This does not require that its down payment bind.

Combining \(qz\le\beta Ax\) with the budget gives

$$ z\le\frac{\beta A}{qM}(W_i-dPh_i^Y). \tag{20} $$

Moreover, \(P h_i^{\mathrm{floor}}(P)\) is nondecreasing in \(P\). Hence

$$ \frac{\gamma z}{APh_i^Y} \le \frac{\beta\gamma}{qM} \left[ \frac{w_i^{\mathrm{high}}} {P_{\mathrm{low}}h_i^{\mathrm{floor}}(P_{\mathrm{low}})} -d \right] <v. $$

At stationarity these are precisely the resources and titles inherited by the corresponding old owners. Equation (3) therefore establishes their slack retention and \(MV_i^O=vP\).

C. The bracket is aggregate, including endogenous tenure

At \(P_{\mathrm{low}}\), selected owners choose their down-payment-limited house and have fertility at least \(n_i^{\mathrm{bench}}(P_{\mathrm{low}})\).

Their ownership probability is at least \(p_i^{\mathrm{low}}\). The owner test in (14) is feasible at zero rebates because \(\phi\ge q\); positive current and future rebates only expand its feasible consumption opportunities. Conversely, \(\mathcal R(P_{\mathrm{low}},w_i^{\mathrm{high}})\) is an upper bound on actual renter utility. The logistic rule therefore gives the claimed probability bound.

Thus (B-low) makes aggregate fertility exceed replacement, uniformly over the candidate rebates. Notice that it counts only selected owners. Renters and other owners contribute additional positive fertility.

For every household, irrespective of its financial constraints,

$$ \frac{\vartheta}{n} =\frac{\chi}{x}+\frac{\alpha\kappa}{h-\kappa n} >\frac{\alpha\kappa}{h-\kappa n}, $$

so

$$ n< \frac{\vartheta h}{\kappa(\alpha+\vartheta)}. \tag{21} $$

The physical and down-payment caps then make aggregate fertility strictly below replacement at \(P_{\mathrm{high}}\), by (B-high).

D. Simultaneous fertility and rebate balance

Let \(s=T/P\). For any candidate \((P,s)\), solve the stationary household problems and let \(h_{\mathrm{tot}}(P,s)\) denote average young-plus-old housing per entrant.

Conditional policies, values, and logistic probabilities are continuous on the compact candidate rectangle. Also,

$$ 0<h_{\mathrm{tot}}(P,s)<2\bar h. $$

The rebate equation, after housing clearing, is

$$ s=\frac{q\tau^p}{2}h_{\mathrm{tot}}(P,s). \tag{22} $$

This follows directly from the model’s government budget; goods and bonds impose no additional domestic clearing equation.

For \(\tau^p>0\), consider the continuous self-map of

$$ [P_{\mathrm{low}},P_{\mathrm{high}}] \times[0,q\tau^p\bar h] $$

given by

$$ (P,s)\longmapsto \left( \operatorname{proj}_{[P_{\mathrm{low}},P_{\mathrm{high}}]} \bigl[P+\bar n(P,s)-1/\nu\bigr], \ \frac{q\tau^p}{2}h_{\mathrm{tot}}(P,s) \right). $$

Brouwer gives a fixed point. The strict fertility signs exclude both price boundaries, so the fixed point satisfies aggregate replacement and (22). For zero taxes, the ordinary intermediate value theorem suffices.

Set

$$ N=\frac{\bar H}{h_{\mathrm{tot}}(P,s)},\qquad Y=O=N, $$

and take the old distribution induced by the young choices. All stationary equilibrium conditions follow.

This establishes existence inside the bracket, not exclusion of equilibria outside it.

E. A positive-mass, financially settled improvement

Finite logistic tastes give positive owner mass in both selected sets. Choose equal positive submasses—not necessarily identical types—from young recipients and old donors. Restrict further to uniform strict margins when necessary.

Transfer \(\epsilon>0\) units from an old donor to a young recipient. Preserve fertility and the donor’s estate. The exact compensation is

$$ D_j(\epsilon)= c_j^2\left[ \left(\frac{h_j^2}{h_j^2-\epsilon}\right)^\gamma-1 \right], \qquad D_j'(0)=vP. $$

The young recipient gains

$$ \Delta U_i= \log\left(1-\frac{D_j(\epsilon)}{x_i}\right) +\alpha\log\left(1+\frac{\epsilon}{h_i-\kappa n_i}\right), $$

whose derivative at zero is

$$ \frac{MV_i^Y-vP}{x_i}>0. $$

A common sufficiently small transfer therefore yields a genuine positive-mass Pareto improvement.

For settlement, set

$$ L(\epsilon)=D_j(\epsilon)-u\epsilon. $$

The old seller’s sale proceeds, compensation, and released tax reserve finance its consumption compensation and a bond costing \(qP\epsilon\). The young buyer borrows exactly that amount, sells the additional title next period, and repays it. Estates, future real allocations, total tax receipts, rebates, and existing creditors’ repayments remain unchanged. This is the packet’s settlement, with its positive-tax terms retained.

\(\square\)

3. Nonemptiness: fixed \(\nu\), positive taxes, heterogeneous entrants, substantial renting

The preceding conditions are not merely formally compatible. Here is an analytic family satisfying them.

Fix \(q,\alpha,\gamma,\omega_B,\vartheta,\chi,\kappa,\nu>0\), with \(q<1\), and fix \(q\le\phi<1\). Do not adjust \(\nu\).

Choose desired tenure-share bounds \(p,\eta>0\) with \(p+\eta<1\), and choose

$$ n_{\mathrm{target}}>\frac{1}{p\nu}. $$

Set \(P_{\mathrm{low}}=1\), and choose

$$ \rho< \frac{\kappa(\alpha+\vartheta)}{\vartheta\nu}, \qquad B_{\mathrm{low}}> \frac{\kappa(\alpha+\vartheta)n_{\mathrm{target}}}{\vartheta}, \qquad B_{\mathrm{high}}>B_{\mathrm{low}}, \qquad \bar h>B_{\mathrm{high}}. \tag{23} $$

The first inequality makes renter fertility below replacement at every price.

Define

$$ C_f= \frac{\vartheta}{n_{\mathrm{target}}} -\frac{\alpha\kappa} {B_{\mathrm{low}}-\kappa n_{\mathrm{target}}}>0. $$

Choose \(W_{\mathrm{low}}\) so that \(E_{\mathrm{low}}=W_{\mathrm{low}}-B_{\mathrm{high}}\) satisfies

$$ E_{\mathrm{low}}> \max\left\{ \frac{(1+q)(2+\vartheta)B_{\mathrm{high}}}{\alpha}, \ \chi n_{\mathrm{target}}+\frac{2\chi}{C_f} \right\}, \tag{24} $$

and choose any finite \(W_{\mathrm{high}}>W_{\mathrm{low}}\).

Take any genuinely two-dimensional distribution on

$$ w\in[W_{\mathrm{low}},W_{\mathrm{high}}], \qquad B\in[B_{\mathrm{low}},B_{\mathrm{high}}], $$

with

$$ b=(1-\phi)B,\qquad y=w-b>0. $$

These are explicit support bands, not a neighborhood of a numerical reference.

Choose

$$ P_{\mathrm{high}}> \max\left\{ 1,\frac{\vartheta B_{\mathrm{high}}\nu} {\kappa(\alpha+\vartheta)} \right\}. $$

Writing

$$ v_0=\max\left\{1-q,\frac{\gamma}{\gamma+\omega_B}\right\}, $$

choose a strictly positive \(\beta\) satisfying

$$ \beta< \min\left\{ \frac1A,\, \frac{qv_0B_{\mathrm{low}}} {\gamma(W_{\mathrm{high}}+1+q)} \right\}, \tag{25} $$

and then choose a strictly positive tax rate satisfying

$$ \tau^p< \min\left\{ 1,\, \frac1{qP_{\mathrm{high}}\bar h},\, \frac{\beta E_{\mathrm{low}}} {q^2P_{\mathrm{high}}\bar h(2+\vartheta)} \right\}. \tag{26} $$

These inequalities establish (S), (Y), and (O) for the whole support. They also imply

$$ n_i^{\mathrm{bench}}(1)>n_{\mathrm{target}}, $$

and establish (B-high).

It remains to choose finite tastes. This does not require an infinite-scale or infinite-location limit. On the specified compact price–support rectangle, the displayed owner test and unrestricted renter bound give a finite lower bound \(L\) on \(W^O-W^R\). A feasible capped renter test and the same unrestricted bound give a finite upper bound \(U\).

Both bounds are primitive functions. For example, the renter test chooses

$$ h=\min\{\rho,h_f(P;w)\}, $$

uses the split \(\mathcal N_M(w-dPh,h)\), and then solves the capped old renter problem:

$$ h^2=\min\left\{\rho,\frac{\gamma Z}{Au}\right\}, \quad c^2=\frac{Z-uh^2}{1+\omega_B}, \quad e=\frac{\omega_Bc^2}{q}. $$

Thus \(L,U\) can be taken as extrema of explicit continuous functions on the specified compact rectangle.

Let \(\ell(p)=\log[p/(1-p)]\). Choose finite

$$ \sigma_\xi> \frac{U-L}{\ell(1-\eta)-\ell(p)} $$

and

$$ \sigma_\xi\ell(p)-L < \bar\xi < \sigma_\xi\ell(1-\eta)-U. \tag{27} $$

Ownership probabilities then lie between \(p\) and \(1-\eta\). In particular,

$$ \int p_i^{\mathrm{low}}n_i^{\mathrm{bench}}(1)\,dF >p\,n_{\mathrm{target}}>\frac1\nu. $$

This proves nonemptiness with positive child goods costs, positive taxes, nondegenerate \(F\), a prescribed positive renter share, and renters who never individually attain replacement.

The small-\(\beta\) construction is a nonemptiness device, not a requirement of Proposition 1.

4. Transitions: the dated result uses inherited households, not stationary replacements
Proposition 2: a dated improvement on an equilibrium path

Consider any date \(t\) on an equilibrium path, including a date after an unexpected fertility or finance shock. Use the actual state entering that date and the continuation prices and rebates faced by current households.

No restriction \(\phi_t\ge q\) is imposed here.

Choose \(C_t\ge v_t\). Suppose a positive-\(F\)-mass set of current young types satisfies, with

$$ H_i^{\mathrm{dp}}=\frac{b_i}{(1-\phi_t)P_t}, \qquad W_{it}=y_i+b_i+T_t+qT_{t+1}, $$ $$ E_{it}=W_{it}-u_tH_i^{\mathrm{dp}}>0, $$

and the benchmark from (8)–(9) using \(\vartheta_t\),

$$ \boxed{ H_i^{\mathrm{dp}}<\bar h,\qquad \frac{\alpha x_{it}^{\mathrm{bench}}} {H_i^{\mathrm{dp}}-\kappa n_{it}^{\mathrm{bench}}} >C_t, } \tag{28} $$

together with

$$ \boxed{ \min\left\{ \beta A x_{it}^{\mathrm{bench}} +(\phi_tP_t-qP_{t+1})H_i^{\mathrm{dp}}, \quad \beta x_{it}^{\mathrm{bench}} +(\phi_tP_t+q^2\tau^pP_{t+1})H_i^{\mathrm{dp}} \right\} >qT_{t+1}. } \tag{29} $$

Suppose a positive-\(G_t\)-mass set of actual old owners satisfies

$$ \boxed{ Z_j=a_j+P_tH_j+T_t>0, \qquad \frac{Z_j}{H_j}\le\mathcal Z_t(C_t). } \tag{30} $$

Then the selected young owners have a strictly binding down payment, a slack physical cap, and positive gross bonds. A positive-mass compensated Pareto improvement exists at date \(t\), preserving the rest of the real path and all the specified financial obligations.

Proof. The benchmark housing argument establishes the young housing choice before imposing gross-bond feasibility.

If future retention is slack, \(qz=\beta Ax\), giving the first branch of (29). If future retention binds, the future old budget and estate floor imply

$$ z\ge c^2+(1+q\tau^p)P_{t+1}H_i^{\mathrm{dp}}, \qquad c^2=\frac{\beta x}{q}. $$

Substitution into actual gross bonds gives the second branch. Since \(x\ge x^{\mathrm{bench}}\), (29) verifies feasibility in either case.

Equation (4) establishes the donor bound. The compensation and settlement proof then applies unchanged. \(\square\)

This is not a proposition that assumes a marginal-value gap or a multiplier pattern. It derives recipient behavior from current resources and prices, and donor behavior from an explicit inherited-resource ratio. But it is conditional on an equilibrium path, not a construction of that path.

At an unexpected policy date \(T\), the two paths have the same entering assets, titles, tenure, cohort counts, and existing contracts. Their old resource difference is

$$ Z_j^{\mathrm{policy}}-Z_j^{\mathrm{baseline}} = (P_T^{\mathrm{policy}}-P_T^{\mathrm{baseline}})H_j + (T_T^{\mathrm{policy}}-T_T^{\mathrm{baseline}}). \tag{31} $$

Current \(\phi_T\) does not retroactively change the old household’s mortgage contract. This respects the model’s specified treatment of reforms.

A sufficient condition that propagates donor eligibility

The inherited-state restriction can partly be replaced by restrictions on earlier household decisions. The following zero-tax extension is deliberately more conservative than Proposition 2.

Suppose the initial old cohort was generated by a stationary equilibrium. Subsequently let

$$ 0<\vartheta_t\le\vartheta_{\mathrm{high}}, \qquad \phi_{\mathrm{low}}\le\phi_t\le\phi_{\mathrm{high}}<1. $$

Define

$$ B_i^{\mathrm{low}}=\frac{b_i}{1-\phi_{\mathrm{low}}}, \quad B_i^{\mathrm{high}}=\frac{b_i}{1-\phi_{\mathrm{high}}}, \quad E_i^{\mathrm{low}}=w_i-B_i^{\mathrm{high}}. $$

For the selected source and recipient types, require

$$ E_i^{\mathrm{low}}>0,\qquad \frac{\beta A E_i^{\mathrm{low}}}{M+\vartheta_{\mathrm{high}}}>b_i. \tag{32} $$

Assume explicitly that the equilibrium prices satisfy

$$ P_t\ge P_{\mathrm{floor}}>0,\qquad \frac{P_t}{P_{t-1}}\ge\lambda_{\mathrm{low}}>0, \tag{33} $$

including the initial shock date. These are price-envelope assumptions, not conclusions about shocks.

For recipients, require

$$ \frac{B_i^{\mathrm{high}}}{P_{\mathrm{floor}}}<\bar h, \qquad \frac{\alpha E_i^{\mathrm{low}}} {(M+\vartheta_{\mathrm{high}})B_i^{\mathrm{high}}}>1. \tag{34} $$

For donor source types, define

$$ K_i^{\mathrm{low}}= \min\left\{ P_{\mathrm{floor}}\bar h,\, B_i^{\mathrm{low}},\, \frac{\alpha w_i}{M+\alpha+\vartheta_{\mathrm{high}}} \right\}, $$

and require

$$ \boxed{ \frac{\beta\gamma w_i} {qM\lambda_{\mathrm{low}}K_i^{\mathrm{low}}} +\frac{\gamma}{A} < \frac{\gamma}{\gamma+\omega_B}. } \tag{35} $$

Then eligible recipient and donor masses exist at every date along any equilibrium path satisfying these conditions.

The important carry-forward calculation is as follows. Condition (32) verifies positive gross bonds even without an upper bound on expected appreciation beyond positive user cost. Earlier young choices satisfy

$$ P_{t-1}H_i\ge K_i^{\mathrm{low}}, \qquad \widehat z_i\le\frac{\beta A w_i}{qM}, $$

where \(\widehat z_i=a_i'+\widehat P_tH_i\) uses the price expected when the position was chosen. Actual old resources are

$$ Z_{it}=\widehat z_i+(P_t-\widehat P_t)H_i. $$

Since \(\widehat P_t>0\),

$$ \frac{\gamma Z_{it}}{AP_tH_i} < \frac{\beta\gamma w_i} {qM\lambda_{\mathrm{low}}K_i^{\mathrm{low}}} +\frac{\gamma}{A}. $$

Condition (35) therefore establishes slack old retention, even allowing unexpected price revaluation. Condition (34) gives \(MV^Y>P_t\ge v_t\).

This extension removes the need to assume the desired old distribution at each date. It does not remove (33). Establishing those price bounds—and existence of the post-shock paths—from specified initial states and shocks remains an unproved equilibrium step.

A genuine obstruction to unrestricted global transitions

Equation (21) and the owner physical cap imply, at every date,

$$ \nu\bar n_t< \frac{\nu\vartheta_t\bar h} {\kappa(\alpha+\vartheta_t)}. $$

Suppose a permanent preference decline makes

$$ g_{\max}:= \frac{\nu\vartheta_{\mathrm{new}}\bar h} {\kappa(\alpha+\vartheta_{\mathrm{new}})} <1. \tag{36} $$

Then

$$ Y_t\le g_{\max}^{\,t}Y_0,\qquad O_t=Y_{t-1}. $$

But housing clearing requires

$$ \bar H\le (Y_t+O_t)\bar h. \tag{37} $$

The right side converges to zero, contradicting the fixed positive stock after finitely many dates.

Therefore no infinite equilibrium path of the stated model exists after such a shock, regardless of subsequent changes in \(\phi\). Allowing vacant housing would remove this particular contradiction, but that would be a substantive model extension.

5. Interpretation, the saved reference, and the remaining limits

The restrictions have distinct jobs.

The mechanism is (Y), or its dated counterpart (28): insufficient liquid wealth makes current young-owner housing small relative to adult consumption and child space needs. The benchmark proves a current-flow housing gap without assuming future downsizing.

Donor eligibility is separate: (O) guarantees enough inherited housing relative to old resources. The sharper dated cutoff (30) can admit donors who retain their entire inherited house. Neither condition requires equal young and old type masses.

The aggregate brackets are existence devices: they do not drive the compensated gain. The low bracket explicitly permits—and the nonemptiness family ensures—renters who remain below replacement. Owners offset them in the aggregate. No renter cap, young or old, is assumed slack.

The saving bounds are sufficient, not necessary. The dated two-branch test is appreciably less conservative than the uniform stationary or path-envelope conditions. None of these bounds should be interpreted as necessary empirical signatures of the mechanism.

The saved positive-tax reference is not excluded by the revised household test

Using the supplied reference’s rebate and prices, the dated benchmark gives exactly

$$ E=\frac{717}{400},\qquad n^{\mathrm{bench}}=\frac34,\qquad x^{\mathrm{bench}}=1. $$

The old-owner solution has

$$ Z=\frac{34}{25},\qquad H=1,\qquad h^2=\frac{1480}{3239}<1. $$

Consequently,

$$ MV^Y=\frac{16}{25}, \qquad MV^O=v=\frac{9717}{18500}, $$

with exact gap

$$ MV^Y-MV^O=\frac{2123}{18500}>0. $$

Gross bonds are positive as well. This is a check of the taxed dated criterion, not a numerical existence argument or an application of a zero-tax theorem. The reference and its below-replacement renter ceiling are supplied in the packet.

The unrestricted claim is still false with heterogeneity

For completeness, heterogeneity does not invalidate the efficient class. Set taxes to zero and take a nondegenerate compact \(F\) for which the fully relaxed allocations satisfy all private constraints strictly. With

$$ u_*=\frac{\nu\vartheta\,\mathbb E_F[w]/D-\chi}{\kappa}>0, $$

those allocations have

$$ x_i=\frac{w_i}{D},\qquad n_i=\frac{\vartheta x_i}{\chi+\kappa u_*}, $$ $$ h_i=\frac{\alpha x_i}{u_*}+\kappa n_i, \qquad c_i^2=\frac{\beta x_i}{q}, \qquad h_i^2=\frac{\beta\gamma x_i}{qu_*}. $$

For example, require \( (1-q)\omega_B>q\gamma\), \(\phi\ge q\), and

$$ h_i^2<h_i<\rho,\qquad P_*h_i<B_i $$

throughout the support. Aggregate replacement holds by construction, while every young and old household has marginal housing value \(u_*\). Concavity rules out a consumption- and housing-preserving Pareto improvement. This extends the packet’s slack-constraint counterargument without collapsing \(F\).

What is not established

The dated improvements cannot simply be added across dates. Nor do they prove that a credit reform raises welfare, fertility, or the eventual population level. At two level-stationary endpoints, aggregate fertility is \(1/\nu\) at both; their population levels may differ because their equilibrium housing allocations differ.

Even the household fertility response to compensated extra housing is ambiguous. Holding the post-transfer consumption–housing perturbation fixed and allowing fertility to reoptimize gives

$$ \frac{dn}{d\epsilon} = \frac{ \alpha\kappa/s^2-\chi MV_j^O/x^2 }{ \vartheta/n^2+\chi^2/x^2+\alpha\kappa^2/s^2 }. $$

A housing-welfare gap does not determine this sign. Changing fertility also changes the future demographic path, so it is outside the path-preserving comparison.

The most informative next quantitative diagnostic is therefore one dated joint eligibility calculation: the mass of young owners with physical headroom and high current-flow housing value, alongside the actual old distribution of resources per inherited unit, \(Z/H\), and fixed-estate housing values. A material share of capped renters alone does not establish this comparison, particularly because the planner also respects their cap. The larger model’s moving costs, discrete owner sizes, additional ages, and estate behavior require their own welfare calculation rather than a direct insertion of its demographic normalization into these inequalities.

The completed result is thus heterogeneous, positive-tax stationary existence plus a derived dated transition test. The missing global step is control of the equilibrium price–population path, not the household constraint pattern or the compensated reallocation.