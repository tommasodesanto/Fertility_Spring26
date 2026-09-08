# Pro response: high-beta allocation theorem

Source: rendered accessibility text from the existing Theorem Construction chat, captured in this task. Pro worked for 30m 26s. Formula strings are preserved; inline prose may be split across lines.

Small 
\beta
 is not necessary for the requested reallocation.  In the unchanged positive-tax model, a sufficient condition can instead be placed on  income–liquid-wealth differences across households . The theorem below works at every fixed finite 
\beta>0
, including 
\beta=1
. It does not require old donors to downsize voluntarily, small bequest motives, or negligible child costs.
There is, however, a genuine patience-related obstacle: a constrained young owner can value  current  housing less than its own old-age counterpart. The purchase constraint restricts future housing access as well as current consumption. Comparing different income–wealth groups avoids making that within-household comparison carry the entire argument.
The remaining equilibrium dependence is explicit: the theorem uses the equilibrium price 
P
 and rebate 
T
. It is not a price-free characterization of an arbitrary given 
F
. An analytical construction below verifies that its conditions are jointly compatible with  all  stationary equilibrium equations at fixed 
\beta
.
1. Main proposition
Maintain a level-stationary equilibrium, 
Y=O=N>0
, of the original model. All preference weights and both child costs remain positive; 
q,\phi\in(0,1)
; 
\tau^p>0
; housing caps and logistic taste parameters are finite. Let 
F(y,b)
 have compact, genuinely heterogeneous support in income and liquid wealth. Entrant wealth remains exogenous—not replenished by estates. The original closing rule, net-
a'
 accounting, retention constraint and estate floor are unchanged. 
Pasted text +1
Define
\begin{aligned} D&=1+\alpha+\vartheta+\beta(1+\gamma+\omega_B),\\ d&=1-q+q\tau^p,\qquad u=dP,\\ v&=\max\left\{d,\frac{\gamma(1+q\tau^p)}{\gamma+\omega_B}\right\}, \end{aligned}
and, for each entrant type,
B_i=\frac{b_i}{1-\phi},\qquad W_i=y_i+b_i+(1+q)T,\qquad \ell_i=\frac{W_i}{B_i}.
Here 
B_i/P
 is the largest house permitted by the down payment; 
W_i
 is lifetime resources including rebates; and 
\ell_i
 measures resources relative to purchasing capacity.
Proposition: cross-group old-to-young reallocation
Suppose there are positive-
F
-mass sets 
S_Y
 and 
S_D
, representing young recipients and the entrant types that generated old donors.
For every 
k\in S_Y\cup S_D
, require
\boxed{ B_k<P h_O^{\max}, \qquad qT+(q-\phi)_+B_k < \frac{\beta(1+\omega_B)}{D}W_k. } \tag{F}
For every 
i\in S_Y
 and 
j\in S_D
, require
\boxed{ \ell_j>\frac{dD}{\alpha}, \qquad \frac{\alpha(\ell_i-d)}{D-\alpha} > \max\left\{ v,\, \frac{\beta\gamma[\ell_j-(1+q)d]} {q[1+\beta(1+\omega_B)]} \right\}. } \tag{R}
Then the equilibrium admits a  positive-mass compensated Pareto improvement that transfers housing from old owners to young owners , holding fertility, tenure, estates, existing obligations and the rest of the real path fixed.
The selected types’ young-owner decisions have a strictly binding down-payment constraint, a slack physical cap, and positive gross bond holdings. These properties are conclusions, not assumed constraint patterns.  Neither old retention nor old financial-estate slackness is required.
What the restrictions mean
Condition (F) supplies physical headroom and verifies that the household wants positive  gross  bonds after accounting for future rebates and mortgage repayment. It does not require 
a'>0
.
Condition (R) identifies two groups through their endowments. The first inequality is a sufficient way to establish the donor source types’ original purchases, so their inherited houses are known rather than freely specified. The second says that recipients have sufficiently greater resources relative to down-payment purchasing capacity.
Holding total resources approximately fixed, lower 
b
 means a smaller financed house. Alternatively, increasing income at fixed 
b
 raises consumption and fertility demand without increasing the down-payment ceiling. These are the useful dimensions of heterogeneity—not simply differences in total wealth.
Equilibrium-state versus primitive restrictions.  Only 
P
 and 
T
 are unidentified equilibrium objects in these tests. In particular, 
T
 inside 
\ell_i
 is not a primitive. Housing clearing, aggregate replacement and
T=\frac{q\tau^pP\bar H}{2N}
continue to hold. No price endpoint is being substituted for equilibrium. 
Pasted text
2. Proof
Step 1: verify saving without assuming an old-age regime
Temporarily remove only the young owner’s gross-bond restriction. Set
z=a'+Ph+T.
Combining the two age-specific budgets gives the equivalent lifetime problem
\max\; u^y(x+\chi n,h,n)+\beta u^o(c^2,h^2,e)
subject to
x+\chi n+dPh+qc^2+qdPh^2+q^2e=W, \tag{1}
h\le B/P,\quad h\le h_O^{\max},\quad h^2\le h,\quad e\ge Ph^2.
These are the original retention and estate restrictions, not replacements for them. 
Pasted text
The consumption first-order conditions and the estate constraint imply
c^2=\frac{\beta x}{q}, \qquad q^2e\ge\beta\omega_Bx. \tag{2}
Consequently,
qz=qc^2+qdPh^2+q^2e \ge\beta(1+\omega_B)x. \tag{3}
The objective is logarithmically homogeneous, with total weight 
D
. All restrictions except the two current housing ceilings are homogeneous in choices. Its scaling identity therefore gives
D=\frac Wx+\text{nonnegative housing-cap terms}, \qquad\text{so}\qquad x\ge\frac WD. \tag{4}
Actual gross bonds—not net 
a'
—are
qa'+\phi Ph=qz-qT+(\phi-q)Ph.
Using 
Ph\le B
, (3)–(4), and (F),
qa'+\phi Ph \ge \frac{\beta(1+\omega_B)}D W-qT-(q-\phi)_+B >0. \tag{5}
Thus the relaxed solution satisfies the omitted restriction strictly. It is the original household optimum, and the Euler relation in (2) is valid.
Step 2: establish the down payment and bound current housing demand
At any solution,
MV^Y=\frac{\alpha x}{h-\kappa n} > \frac{\alpha x}{h} \ge \frac{\alpha W}{Dh} \ge \frac{\alpha\ell}{D}P. \tag{6}
The donor inequality in (R) makes this exceed 
dP
. The recipient inequality does so as well: its right side is at least 
v\ge d
, hence it implies 
\ell_i>dD/\alpha
.
An additional unit of current housing also weakly relaxes future retention. Therefore, when 
MV^Y>dP
, the household cannot optimally stop below both current housing ceilings. By (F),
h=\frac BP<h_O^{\max}. \tag{7}
The down-payment multiplier is strictly positive because current housing benefits alone exceed user cost.
A sharper consumption bound is useful for the welfare comparison. Write the old value in resources and inherited housing as 
\mathcal V(z,h)
. Homogeneity and the saving first-order condition give
z\mathcal V_z+h\mathcal V_h=1+\gamma+\omega_B, \qquad qz=\beta x\bigl(1+\gamma+\omega_B-h\mathcal V_h\bigr) \le\beta(1+\gamma+\omega_B)x. \tag{8}
Here 
\mathcal V_h\ge0
; future retention need not be slack.
The fertility first-order condition is
\frac{\vartheta}{n} = \frac{\chi}{x}+\frac{\alpha\kappa}{h-\kappa n},
so 
\chi n<\vartheta x
. Combining this with 
h=B/P
 and the budget,
W-dB=x+\chi n+qz<(D-\alpha)x.
Thus
\boxed{ MV_i^Y> P\,\frac{\alpha(\ell_i-d)}{D-\alpha}. } \tag{9}
This lower bound permits arbitrary positive child goods and space costs and either future-retention regime.
Step 3: bound old donors, including those retaining their entire house
Let an old owner have resources 
Z=a+PH+T
, and write 
A=1+\gamma+\omega_B
. Solving the old problem gives
h^2=\min\left\{H,\frac{\gamma Z}{AvP}\right\}, \tag{10}
with consumption
c^2= \min\left\{ \frac{Z-dPh^2}{1+\omega_B}, \ Z-(1+q\tau^p)Ph^2 \right\}. \tag{11}
To see the two branches, without the retention cap the estate-slack solution spends at housing user cost 
dP
. When the estate floor binds, housing and its required estate cost 
(1+q\tau^p)P
 and carry combined utility weight 
\gamma+\omega_B
. Both branches have 
c^2=Z/A
, producing (10). Conditional consumption is (11).
It follows that
\boxed{ MV^O=\max\left\{vP,\frac{\gamma c^2}{H}\right\}. } \tag{12}
With slack retention, the first term applies. With binding retention, 
h^2=H
, and the second is at least the first.
For an old donor generated by type 
j
’s stationary owner plan,
H_j=B_j/P,\qquad c_j^2=\beta x_j/q.
If retention is slack, 
MV_j^O=vP
.
If retention binds, the old budget and (2) give
qz_j\ge\beta(1+\omega_B)x_j+qdB_j.
Hence
W_j > [1+\beta(1+\omega_B)]x_j+(1+q)dB_j,
and
MV_j^O = \frac{\beta\gamma x_j}{qH_j} < P\, \frac{\beta\gamma[\ell_j-(1+q)d]} {q[1+\beta(1+\omega_B)]}.
In either case,
\boxed{ MV_j^O\le P\max\left\{ v,\, \frac{\beta\gamma[\ell_j-(1+q)d]} {q[1+\beta(1+\omega_B)]} \right\}. } \tag{13}
Equations (9), (13), and (R) establish the required strict valuation gap.
Crucially, the old donor is the actual old-age realization of an earlier stationary household decision. No independent old wealth distribution has been chosen.
Step 4: compensate the old and settle the transaction
Finite logistic tastes imply positive owner mass in both selected groups. Select equal positive  submasses , restricting them further where necessary to obtain uniform strict margins.
Give a young recipient 
\epsilon
 units taken from an old donor. Preserve the donor’s estate and compensate it by
D_j(\epsilon) = c_j^2 \left[ \left(\frac{h_j^2}{h_j^2-\epsilon}\right)^\gamma-1 \right].
This keeps 
u^o
 exactly unchanged, and
D_j'(0)=MV_j^O.
The young recipient, with fertility fixed, gains
\Delta U_i = \log\left(1-\frac{D_j(\epsilon)}{x_i}\right) + \alpha\log\left(1+\frac{\epsilon}{h_i-\kappa n_i}\right).
Its derivative at zero is
\Delta U_i'(0) = \frac{MV_i^Y-MV_j^O}{x_i}>0.
A common sufficiently small positive transfer respects physical headroom and benefits a positive mass of young households.
For financial settlement, set
L(\epsilon)=D_j(\epsilon)-dP\epsilon.
The seller’s title proceeds, supplementary payment and released tax reserve, after consumption compensation, leave
P\epsilon+L(\epsilon)+q\tau^pP\epsilon-D_j(\epsilon) =qP\epsilon.
The seller buys a bond with these funds. The buyer’s additional borrowing is exactly 
qP\epsilon
. Next period the buyer sells the additional title for 
P\epsilon
 and repays that bond; the donor’s estate receives the same amount in place of the lost title. Existing obligations, estates, subsequent real allocations and aggregate tax receipts remain unchanged. This is the packet’s settlement, including its positive-tax accounting. 
Pasted text
\square
3. What high 
\beta
 changes
The admissible region has no upper bound on 
\beta
Whenever 
\ell_j>(1+q)d
, the donor-retention part of (R) can be written
\frac{\ell_i-d}{\ell_j-(1+q)d} > \frac{\beta\gamma}{q\alpha} \frac{D-\alpha}{1+\beta(1+\omega_B)}. \tag{14}
This is the explicit tradeoff: greater patience can require greater separation in resources per unit of purchasing capacity. It does  not  require 
\beta
 to be small.
For example, at 
\beta=1
, the recipient restriction is
\ell_i> d+ \frac{2+\gamma+\omega_B+\vartheta}{\alpha} \max\left\{ v,\, \frac{\gamma[\ell_j-(1+q)d]}{q(2+\omega_B)} \right\}. \tag{15}
All thresholds are finite with the other preference weights held fixed. The relevant intertemporal comparison is 
\beta/q
, because the model’s Euler equation is 
c^2/x=\beta/q
; both objects refer to the same young-to-old interval.
These bounds are sufficient, not necessary. For particular parameters they may require substantial income–wealth separation.  Finite admissibility is not evidence of empirical plausibility.  But nothing must approach an extreme as 
\beta\uparrow1
; the compatibility construction below works with fixed positive preference weights and a fixed positive tax rate.
The within-type obstacle is real
For a stationary owner with positive gross saving, (12) gives the exact identity
MV_i^O = \max\left\{ vP,\frac{\beta\gamma x_i}{qh_i} \right\}. \tag{16}
Meanwhile, the fertility first-order condition implies
MV_i^Y = \frac{x_i}{h_i} \left( \alpha+\vartheta-\frac{\chi n_i}{x_i} \right) < \frac{x_i}{h_i}(\alpha+\vartheta). \tag{17}
Therefore,
\boxed{ \beta\gamma\ge q(\alpha+\vartheta) \quad\Longrightarrow\quad MV_i^Y<MV_i^O. } \tag{18}
This is an actual impossibility result for the  same-type  comparison—not for cross-group reallocation. The positive-tax compatibility construction below also works under (18), so it supplies equilibria with binding down payments and this reversed within-type ordering.
The mechanism is the prohibition on old-age upsizing. A young purchase creates both current housing services and future housing access. A positive purchase-constraint multiplier does not distinguish those benefits. The original model expressly imposes that restriction. 
Pasted text
4. Analytical compatibility at fixed 
\beta
, with all markets respected
Here is a construction of admissible primitives, not a freely chosen price for an already specified economy.
Fix  any finite 
\beta>0
, positive 
\alpha,\gamma,\omega_B,\vartheta,\tau^p
, 
q,\phi\in(0,1)
, the demographic parameter 
\nu>0
, housing stock 
\bar H>0
, and any finite logistic taste parameters. No preference weight or tax rate is sent to zero.
Choose 
0<W_{\mathrm{low}}<W_{\mathrm{high}}<\infty
 and any 
g>0
. Choose a donor resource-ratio interval with lower endpoint satisfying
\ell_{D,\mathrm{low}}> \max\left\{ \frac{dD}{\alpha},\ \frac{W_{\mathrm{high}}}{W_{\mathrm{low}}} \bigl[(1-\phi)+(1+q)q\tau^p\bigr],\ \frac{D W_{\mathrm{high}}} {\beta(1+\omega_B)W_{\mathrm{low}}} \bigl[q^2\tau^p+(q-\phi)_+\bigr] \right\}. \tag{19}
Choose any finite 
\ell_{D,\mathrm{high}}>\ell_{D,\mathrm{low}}
.
Choose a recipient interval whose lower endpoint exceeds both 
\ell_{D,\mathrm{high}}
 and
d+\frac{D-\alpha}{\alpha} \max\left\{ v,\, \frac{\beta\gamma[\ell_{D,\mathrm{high}}-(1+q)d]} {q[1+\beta(1+\omega_B)]} \right\}. \tag{20}
Its upper endpoint can be any larger finite number.
Put a positive density on each of the two rectangles in 
(W,\ell)
, using the resource interval and these ratio intervals. Set
B=W/\ell,
and choose finite caps satisfying
0<h_R^{\max}<\frac{W_{\mathrm{low}}}{\ell_{Y,\mathrm{high}}}, \qquad h_O^{\max}>\frac{W_{\mathrm{high}}}{\ell_{D,\mathrm{low}}}. \tag{21}
Now construct allocations at the candidate price 
P=1
.
First solve the household problems in lifetime resources 
W
, with auxiliary child costs 
\kappa=1,\chi=g
, temporarily omitting young financial lower bounds. These strictly concave problems have well-defined solutions. The preceding proof establishes 
h^Y_O=B
. Use the prescribed finite tastes to determine tenure probabilities.
Let 
\bar h_{\mathrm{life}}
 be mean young-plus-old housing, and let 
\bar n_{\mathrm{aux}}>0
 be mean auxiliary fertility. Since owners retain no more than they purchase and both renter ages face their cap,
0<\bar h_{\mathrm{life}} \le \frac{2W_{\mathrm{high}}}{\ell_{D,\mathrm{low}}}.
Define
T=\frac{q\tau^p}{2}\bar h_{\mathrm{life}}, \qquad b=(1-\phi)B, \qquad y=W-(1+q)T-b. \tag{22}
The second bound in (19) ensures 
y>0
. The third ensures (F) for every owner type. It also verifies renters’ 
a'\ge0
: their analogous homogeneity bound is 
x_R\ge W/D
, and their old-age resources satisfy
qz_R\ge\beta(1+\omega_B)x_R>qT.
Thus none of the omitted financial constraints changes the constructed policies.
Finally choose the  actual , positive child costs
\kappa=\nu\bar n_{\mathrm{aux}}, \qquad \chi=g\kappa. \tag{23}
This preserves child space 
\kappa n
, child goods expenditure 
\chi n
, all other real allocations, and tenure probabilities. It subtracts the same constant 
\vartheta\log\kappa
 from both conditional young values. Actual fertility satisfies
\nu\bar n = \nu\frac{\bar n_{\mathrm{aux}}}{\kappa} =1.
Set
N=\frac{\bar H}{\bar h_{\mathrm{life}}}, \qquad Y=O=N,
and take the old distribution generated by these young decisions. Housing clears, the government budget holds by (22), and demography is stationary. The mapping from the two-dimensional 
(W,\ell)
 distribution to 
(y,b)
 is nonsingular, so this is genuine income–wealth heterogeneity. Estates do not enter the construction of subsequent entrants’ wealth. These are precisely the maintained equilibrium requirements. 
Pasted text
This proves joint compatibility analytically. The bounds remain uniformly finite as 
\beta
 ranges over any closed interval bounded away from zero and ending at one. No numerical reference equilibrium, vanishing neighborhood, limiting tenure share, or small-tax argument is involved.
5. A stronger falsification check: binding purchases can reflect future housing needs
The within-type result already identifies the obstacle. The following  zero-tax specialization  shows something stronger: current young-owner housing values can lie below  every  old household’s housing value despite strictly binding down payments.
Set 
\tau^p=0
, 
d=1-q
, and suppose
\beta\gamma>q(\alpha+\vartheta), \qquad d\omega_B>q\gamma, \qquad \phi>q^2. \tag{24}
Choose any positive 
g=\chi/\kappa
, and choose 
m
 in the nonempty interval
\frac{(1+q)d(\alpha+\vartheta)} {\alpha+\vartheta+\beta\gamma} <m<d. \tag{25}
Define
r(m)= \left[ \frac{\alpha}{m}+\frac{\vartheta}{g+m} \right]^{-1}.
At 
P=1
, choose a down-payment-limited house 
H
, and set
x=r(m)H,\qquad n=\frac{\vartheta r(m)H}{\kappa(g+m)},
c^2=\frac{\beta x}{q},\qquad h^2=H,\qquad e=\frac{\beta\omega_Bx}{q^2}.
Support these choices with
b=(1-\phi)H,
y+b=(1+q)dH+[1+\beta(1+\omega_B)]x+\chi n. \tag{26}
The young consumption and fertility first-order conditions give 
MV^Y=m
. Moreover,
r(m)>\frac{m}{\alpha+\vartheta},
so (25) implies
m+\beta\gamma r(m)>(1+q)d. \tag{27}
This is exactly the strictly positive purchase-ceiling derivative when future retention binds and the financial estate is positive.
It also implies 
\beta\gamma r(m)>qd
. Together with (24), this verifies both strict retention binding and estate-floor slackness. Gross bonds are
qa'+\phi H = \beta(1+\omega_B)x+(\phi-q^2)H>0.
Income is positive because 
(1+q)d=1-q^2>1-\phi
. Thus these are actual optimal constrained-owner choices, not an assumed pattern.
Choose nondegenerate intervals of 
m
 and 
H
 inside these explicit bounds. A rental cap satisfying
0<h_R^{\max}< \min\left\{ H_{\mathrm{low}}, \frac{\alpha W_{\mathrm{low}}}{dD} \right\}
makes young renters physically capped. Complete stationary demography and housing clearing as in the preceding construction.
Then
MV^Y_O=m<d,
whereas every old renter has 
MV^O_R\ge d
, and these old owners have 
MV^O_O>d
. Young renters cannot receive more housing because they are already at their physical cap. Concavity therefore rules out the proposed compensated transfer  from old housing to young housing  in this family.
This does not establish global efficiency; other reallocations may exist. It establishes that  financing constraints and tenure segmentation alone do not determine the desired direction . The cross-group separation in the main proposition is doing genuine economic work.
6. What Coven et al. do—and do not—supply
Your reading of  January 31, 2025  is correct. Sections 2.2–2.4, pp. 9–11, give welfare-derivative decompositions, not the requested primitive-parameter compensated Pareto theorem. In §2.1, pp. 7–9, old-age utility is a function of bequeathed total wealth, not your separate old consumption–retention–estate problem. Their equation (2) restricts initial income, assets and housing; footnote 2 explicitly identifies it as an initial-state restriction. Equation (3) imposes a collateral borrowing bound. That differs from your separate, pre-income 
b
-only closing constraint . 
The  August 1, 2026 text also contains analytical results : §3.1, pp. 15–18, has Lemma 1 and Proposition 2. Its baseline utility is 
\log C+\beta\log C'
; housing is an asset, the old sell it, and borrowing satisfies 
B\le\lambda PH
. Proposition 2 studies an uncompensated tax change and explicitly has the old lose. It is not a fixed-estate Pareto result. 
Their quantitative model allows repeated housing adjustment and tenure choice, with origination LTV and payment-to-income limits, an  owner minimum size , and uncapped renter quantities—rather than your owner/renter upper caps and prohibition on old-age purchases. Bequest incentives also interact with sale taxation and step-up basis. 
Thus your simplification introduced a substantive difficulty: current purchases limit later housing access even when gross financial saving is positive. Their tax-envelope results do not resolve that current-flow valuation problem. The main proposition resolves it without changing your model, by comparing different income–wealth groups.
7. Dated use along a transition
At date 
t
, the stationary donor-source calculation must be replaced by a check on  actual inherited households .
Define
u_t=(1+q\tau^p)P_t-qP_{t+1}, \qquad v_t= \max\left\{ u_t,\frac{\gamma}{\gamma+\omega_B}(1+q\tau^p)P_t \right\}.
For any compensation cutoff 
C_t\ge v_t
, an actual old owner satisfies 
MV^O\le C_t
 exactly when
\boxed{ \frac{a+P_tH+T_t}{H} \le \max\left\{ u_t+\frac{1+\omega_B}{\gamma}C_t,\, (1+q\tau^p)P_t+\frac{C_t}{\gamma} \right\}. } \tag{28}
This follows from the same old solution and includes binding-retention donors.
For current young recipients, put
H_i^{\mathrm{dp}}=\frac{b_i}{(1-\phi_t)P_t}, \quad W_{it}=y_i+b_i+T_t+qT_{t+1}, \quad D_t=1+\alpha+\vartheta_t+\beta(1+\gamma+\omega_B).
The corresponding sufficient checks are
H_i^{\mathrm{dp}}<h_O^{\max}, \qquad \frac{\alpha(W_{it}-u_tH_i^{\mathrm{dp}})} {(D_t-\alpha)H_i^{\mathrm{dp}}}>C_t, \tag{29}
and
qT_{t+1} +(qP_{t+1}-\phi_tP_t)_+H_i^{\mathrm{dp}} < \frac{\beta(1+\omega_B)}{D_t}W_{it}. \tag{30}
The proof above then establishes the dated improvement.
The remaining inputs are actual inherited 
a,H
, current and expected prices, and rebates—not a newly solved stationary old cohort. Initial donor eligibility must likewise be checked directly, or established from the decision date that generated the initial old cohort. Unexpected reforms preserve inherited contracts, as required by the packet. 
Pasted text
The completed result is an allocation theorem at high 
\beta
, with explicit income–wealth separation and financial-feasibility checks. It does not establish that a tax or credit policy implements the transfer, or that the eventual population level rises.
