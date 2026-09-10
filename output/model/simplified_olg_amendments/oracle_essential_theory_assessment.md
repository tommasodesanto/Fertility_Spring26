# Assessment of the September 9 Pro response

The first response contains a coherent conditional allocation and fertility argument. It is not yet a satisfactory replacement for the paper's theory: its mortgage restriction excludes the earlier long-period financing benchmark. One focused follow-up is running on that issue. No new model specification has been adopted.

[Full Pro response](oracle_essential_theory_response.md) · [Pro conversation](https://chatgpt.com/c/6aa1e33d-d328-83ea-971e-0541c9e3e359) · [Exact follow-up](../../../docs/prompts/oracle_essential_theory_finance_followup.md)

## What the proposed model explains

Two types have the same present-value resources. One receives \(W\) while young and a known, nonpledgeable payment \(Y\) when old. The other receives \(W+qY\) while young. Here \(q\) is the price of one unit of next-period goods. The first type therefore has fewer resources available for the down payment, without being poorer over its lifetime.

Both types choose housing, consumption, fertility and saving. A common divisible housing stock can be occupied by owners or renters; renters cannot choose more than \(r\) units. Old households can sell, resize and borrow against housing. The mortgage share is \(\phi\).

The proposed equilibrium has constrained young households renting at \(r\), liquid young households owning larger homes, and unconstrained old housing choices. This assignment is derived under sufficient primitive inequalities. It is not assumed to be every equilibrium of the model.

## The result, separated from its assumptions

Let \(x_i=c_i^y-\chi n_i\) denote consumption beyond children's goods needs, and let \(s_i=h_i^y-\kappa n_i\) denote housing beyond children's space needs. Equal current utility weights and the same housing coefficient \(\alpha\) at both ages imply that the dated planner pools adult goods and adult space:

\[
X=\frac{\bar x+\bar c_o}{2},\qquad
S=\frac{\bar s+\bar h_o}{2}.
\]

The planner chooses all current consumption and housing, initially fixes individual fertility, and preserves young continuation wealth and old net estates. Hence

\[
c_i^{y,F}=X+\chi n_i,\quad h_i^{y,F}=S+\kappa n_i,
\qquad c_j^{o,F}=X,\quad h_j^{o,F}=S.
\]

The equilibrium conditions establish

\[
\bar c_o\geq\bar x,\qquad
\bar h_o>\bar s.
\]

Consequently,

\[
\bar h_y^F-\bar h_y^E
=\frac{\bar h_o-\bar s}{2}>0.
\]

Every constrained young family receives more consumption and housing. The allocation may reduce liquid young families' utility and old households' utility; this is an equally weighted utilitarian comparison.

There is also a separate efficiency argument that does not rely on equal welfare weights. In this equilibrium, constrained renters value housing relative to goods more highly than unconstrained old households:

\[
\frac{\alpha x_C}{s_C}>p
=\frac{\alpha c_o}{h_o}.
\]

A small housing transfer can therefore compensate an old donor in goods while improving the constrained family's utility, keeping continuation wealth and estates fixed. The authority must relax the relevant private financial constraints. This is not a transfers-only policy theorem.

## Why the fertility result is useful

For the proposed Stone–Geary preferences, both additional consumption and additional housing raise the marginal utility of fertility. Each constrained family would therefore choose more children if offered its fixed-fertility planner bundle.

The joint planner result is stronger. Private fertility satisfies

\[
n_i=\vartheta
\left(\frac{\chi}{x_i}+\frac{\alpha\kappa}{s_i}\right)^{-1}.
\]

The reciprocal on the right is concave in \((x_i,s_i)\). Averaging private choices and using the fact that the planner gives the young more average adult goods and space yields

\[
\frac{\vartheta}{\bar n}
\geq \frac{\chi}{\bar x}+\frac{\alpha\kappa}{\bar s}
>
\frac{\chi}{X(\bar n)}+\frac{\alpha\kappa}{S(\bar n)}.
\]

The joint planner's fertility first-order condition then implies \(n^F>\bar n\). This proof does not require the child-goods cost to be negligible. The current parents' utility is the objective; unborn households receive no welfare weight.

## What the primitive conditions do

The source's inequalities have three separate roles:

1. The young with deferred income can afford the largest rental home but cannot finance buying that much housing. The unrestricted rental demand exceeds the rental cap.
2. Later resources and liquid households' resources are large enough to finance unconstrained old choices, make constrained young saving zero, and ensure that old average adult consumption is at least young average adult consumption.
3. Replacement fertility occurs inside the price interval where those household choices hold. A strictly decreasing fertility function determines rent; housing clearing determines the stationary population level.

These are genuine sufficient restrictions on primitives. They remain substantial. In particular, the price-window condition and liquid-young financing condition imply

\[
\frac{q-\phi}{1+\tau-q}(\alpha+\vartheta)>1,
\qquad
\beta(1+\alpha+\omega_B)>1.
\]

The response does not establish empirical mildness. It also does not establish that every old household occupies more total housing than every young household. The theorem compares old housing with young housing net of children's needs. Within the liquid type, old housing is smaller than its own earlier young home when \(\beta<q\).

The separate bound \(\beta<q\) is less fundamental. For the stationary allocation and dated fertility results, the proof extends beyond it if constrained young nonsaving is ensured directly by \(qY/(\beta K)\geq b\), retaining the other finance, resource and price-window conditions. Here \(K=1+\alpha+\omega_B\) and \(b=(q-\phi)W/(1+\tau-\phi)\). This extension has been checked for the core result only. It does not remove the \(\phi<q\) issue.

Uniqueness must refer to prices, real allocations and population within the specified regime. It need not include old tenure labels or portfolios: if an old household's optimal home is below the rental cap, renting and owning can tie.

## The decisive financing problem

The owner budget in Pro's proposed model is

\[
c+(1+\tau)Ph+a=w,\qquad a\geq-\phi Ph,
\qquad z=Y_{\rm next}+a/q+Ph.
\]

Let \(\lambda=u_c>0\), and let \(\mu\geq0\) be the multiplier on its collateral constraint. The housing first-order condition gives the exact identity

\[
\frac{u_h}{u_c}
=p+(q-\phi)P\frac{\mu}{\lambda},
\qquad p=(1+\tau-q)P.
\]

Thus the sign of the financial distortion depends on \(q-\phi\). If \(\phi>q\), a binding collateral constraint makes housing a way to obtain current liquidity and puts the marginal housing valuation below the user cost. This is a feature of these repayment and consumption dates, not a universal claim about mortgages.

The renter argument fails for the same reason. When \(\phi\geq q\), any renter plan can be reproduced by owning with \(a_O=a_R-qPh\), leaving consumption and continuation wealth unchanged and satisfying the mortgage limit. A household with a strictly binding rental cap then benefits from choosing a slightly larger owned home.

Consequently, \(\phi<q\) is not just an inconvenient proof bound for Pro's argument. It is necessary for the desired positive housing distortion in this model. At \(\phi=0.8\), it requires \(q>0.8\), and excludes the prior \(q=0.5\) benchmark. The single follow-up asks whether a conventional repayment structure or a minimum owner size can remove this problem while retaining the intended mechanism and a simple proof.

## What the transition section adds

Within the same strict regime, Pro derives a policy that raises the maximum rental size. An additional goods-versus-space restriction,

\[
(2\alpha+\vartheta)\chi\leq\kappa p_\ell,
\]

signs the stationary response: a higher rental limit raises constrained families' fertility and the terminal population. A lower fertility taste lowers the terminal population. Both destinations have mean replacement fertility.

The local stability calculation and the initial capital-gain boundary pass independent checking. The announcement comparison retains inherited bonds and houses, rather than fixing the initial old's market-valued wealth. For a taste change, their inherited saving choice must also be retained for the first settlement.

The signed path implication is higher cumulative reproduction between the same initial population and the higher terminal population. Neither a higher fertility rate at every date nor an unconditional impact increase is proved. The policy increases access to rental space; it does not implement the result through property taxes or increase constrained young homeownership.

## Model changes that remain author decisions

Pro replaces the morning/afternoon idea with exogenous later resources, permits old borrowing, removes the owner upper bound, makes housing tastes equal across ages, and uses property taxes for additively valued public services rather than household rebates. Its common stock can move freely between tenure categories. These changes must remain visible; none has been integrated into the main deck or manuscript.

The comparison with Coven et al. is broadly correct: their analytical model focuses on housing as an asset and tax capitalization, with old liquidation; Proposition 2 gives cohort-specific welfare effects. Their quantitative model contains richer lifecycle income, an owner-size minimum, separate tenure markets and borrowing against old housing. Their quantitative fiscal system also includes endogenous household transfers; the proposed public-service-only closure should not be attributed to that model. [August 2026 paper, sections 3–4](https://abdouecon.github.io/research/papers/Property_Tax.pdf).

## Verification record

The complete response was captured from the visible conversation, including 142 inline/displayed mathematical objects. Whitespace-normalized text matches the browser extraction exactly (20,980 characters; FNV-1a 360915802). The lead derived the financial-distortion identity and checked the principal planner and fertility inequalities. Two separate Astra/max reviewers checked the stationary household/planner block and the policy/transition block. Their technical reviews are reproduced below as supporting checks.


<details>
<summary>Households, planner and fertility — independent review</summary>

# Independent review of the Pro liquidity proposal

Reviewed September 9, 2026. Scope: stationary households and tenure, primitive conditions (4)–(6), the full dated planner and settlement, and the two dated fertility claims. No transition, literature, calibration, or alternative-model audit. Source: `output/model/simplified_olg_amendments/oracle_essential_theory_response.md`, left unchanged. Mandatory project context was read. This is a review of an unadopted proposal.

## Verdict

| Claim | Verdict | Qualification or smallest repair |
|---|---|---|
| Stationary household quantities and positive replacement equilibrium | **PASS** | Conditional on every stated inequality and the external-finance closure. |
| Liquid young necessarily own a home larger than the rental ceiling | **PASS** | This follows from the price window and \(L\ge b\); the response omits the short proof. |
| Old choices are financeable | **PASS** | Old households may rent or own when their optimal house is at most \(r\). Positive financial estates are not implied. |
| Equilibrium is unique within the regime | **FAIL literally; small wording repair** | Prices, housing/consumption quantities, fertility, and cohort mass are unique. Old tenure labels and associated portfolios need not be unique. |
| Conditions are jointly nonvacuous | **PASS** | They can be strong. The price window forces \(\beta K>1\) and a sharper restriction than \(\phi<q\). For fixed demographic conversion, the replacement bracket also imposes an upper bound on deferred income. |
| Full dated utilitarian planner raises aggregate young housing and each C family's consumption/housing at fixed fertility | **PASS** | Both goods and housing are optimized. Individual continuation wealth and estates are fixed. |
| Financial settlement preserves resources and promised wealth | **PASS** | Intermediaries' matching bond changes must be included. The planner explicitly relaxes private financing limits. |
| C families' fertility rises when offered their fixed-fertility planner bundles | **PASS** | Their consumption and housing are held at the offered bundle and continuation opportunities remain fixed. |
| Joint dated, parents-only planner raises mean fertility | **PASS** | The harmonic-mean Jensen argument has the correct direction; the full planner can be shown to choose common fertility. |
| Ordinary LTV makes the desired financial housing wedge general | **FAIL as a broader interpretation** | This budget's positive wedge requires \(\phi<q\). When \(\phi>q\), constrained ownership has the opposite housing wedge. This is a limitation of the specified mortgage/timing structure, not a general impossibility theorem. |

## 1. Household and tenure checks

Write \(t=(q-\phi)/d>0\), \(A=1+\alpha+\vartheta\), and \(M=A+\beta K=G+\vartheta\). An owner satisfies the consolidated budget and mortgage condition

\[
c+ph+qz=w+qY_{\rm next},\qquad
q(z-Y_{\rm next})\ge(q-\phi)Ph.
\]

For the old, the unconstrained shares are exactly (9). Ownership can implement them iff \(\omega_B\ge\alpha t\). When \(h_o\le r\), renting implements exactly the same allocation with saving \(qe>0\). Ownership instead uses

\[
a_o=\left(\omega_B-\frac{q\alpha}{d}\right)c_o,
\]

which can be negative even when all proposition conditions hold. Thus the old are free to select tenure; their optimal quantities do not depend on that selection. Their value is \(K\log z\) plus a constant, as used in the young problem.

The liquid young's unrestricted shares are correct. The second finance inequality is sufficient because

\[
ph_L=L\left[\alpha+\frac{\vartheta}{1+\chi/(\kappa p)}\right]
<(\alpha+\vartheta)L.
\]

For a C owner, \(c+\delta Ph\le W\) implies \(h<W/(\delta P)<r\) whenever \(p>p_\ell\). Every owner plan is then replicated by a legal renter with financial saving

\[
a_R=q(z-Y)=a_O+qPh\ge(q-\phi)Ph>0.
\]

At the proposed cap choice \(x_C=W-pr-\chi n_C<b\). Condition (5), with \(\beta<q\), gives

\[
f(Y/K-b)\ge(1-f)(1-\beta/q)L>0.
\]

Thus \(Y/K>b>x_C\), and the derivative toward positive saving is strictly negative. The cap binds because \(p<p_u\). Strict concavity in consumption, adult space, fertility and continuation wealth establishes the global rental optimum; replication excludes all owner deviations.

The omitted proof that liquid young own larger homes is short. Since \(\widehat h_C(p_\ell)>r\),

\[
t(\alpha+\vartheta)>1,
\qquad
b=\frac{Wt}{1+t}>\frac W{1+\alpha+\vartheta}.
\]

Consequently \(L\ge b>W/A\), so

\[
h_L=\frac{AL}{W}\widehat h_C(p)>\widehat h_C(p)>r.
\]

The fertility root is unique: \(F_p<0\), both bracket endpoints are in the feasible domain, and (6) places the root strictly between them. Equation (7) then gives the unique positive stationary cohort mass. Equal young/old masses are justified by replacement. External bonds and a residual rental intermediary are maintained closures, so no additional domestic bond-clearing condition has been omitted.

**Tenure-uniqueness repair.** Say “the equilibrium prices, quantities, fertility and population are unique within this regime, with old tenure potentially indeterminate.” The hypotheses do not force old ownership. For example, the exact internal witness

\[
\alpha=\vartheta=W=r=\kappa=1,\quad\chi=1/10,
\quad q=9/10,\ \phi=4/5,\ \tau=1/100,
\quad\beta=1/2,\ \omega_B=3,\ f=9/10,\ Y=14/5
\]

at \(p=3/5\), \(\nu=1/F(p)\), satisfies all inequalities. Here \(p_\ell=11/21<p\), \(\widehat h_C(p)=65/63>r\), and the old housing choices are \(14/15\) and \(16/27\), both below \(r\). Their financial positions satisfy \(a_o=-(57/11)c_o<0\). This witness only checks nonvacuity and refutes strict tenure uniqueness; it is not a recommended parameter family.

## 2. Exact economic force and nonvacuity

Let \(z_\ell=\chi/(\kappa p_\ell)>0\). The first price-window inequality is exactly

\[
\frac{q-\phi}{d}>
\frac{1+(1+\vartheta)z_\ell}{\alpha+\vartheta+\alpha z_\ell}
>\frac1{\alpha+\vartheta}.
\]

It therefore requires

\[
q>\frac{1+\tau+\phi(\alpha+\vartheta)}{1+\alpha+\vartheta},
\qquad
\beta K\ge\frac{q-\phi}{d}(\alpha+\vartheta)>1.
\]

The second implication does not contradict \(\beta/q<1\): a sufficiently large housing or estate weight can make \(\beta K>1\). But the response should not suggest that the impatient region is unrestricted. For \(\phi=.8\), \(q>.8\) is necessary, and the displayed sharper bound may be considerably higher. This is consequential when an age lasts many years.

Conditions (5) impose real restrictions on type composition and deferred income. At fixed other primitives they admit a finite \(Y\) iff

\[
f>f_{\min}:=\frac{(q-\beta)K}{A+qK}<1.
\]

Conditional on that inequality, their exact lower bounds are

\[
Y\ge\frac{Mb-W}{q},\qquad
Y\ge
\frac{fbM+(1-f)(1-\beta/q)W}
{fM/K-(1-f)(q-\beta)}.
\]

Thus there is a nonempty lifecycle-resource region, although small \(f\) cannot always be offset by arbitrarily large \(Y\): liquid types also receive the present value of that larger deferred endowment.

For a **fixed** \(\nu\), replacement adds a finite interval. Define

\[
Y_{\rm rep}(p)=\frac1q\left\{
\frac{M(\chi+\kappa p)}{(1-f)\vartheta}
\left[\frac1\nu-fn_C(p)\right]-W\right\}.
\]

Condition (6) requires \(Y_{\rm rep}(p_\ell)<Y<Y_{\rm rep}(p_u)\). This interval must intersect the lifecycle lower bounds; “sufficiently large deferred income” alone does not prove replacement for a predetermined \(\nu\). Joint nonvacuity is nevertheless exact: after satisfying finance and lifecycle restrictions, the nonempty fertility bracket permits a \(\nu\); for fixed \(\nu\), joint scaling of \(\chi,\kappa\) rescales fertility without changing housing demands or the price window. None of this makes the region empirically mild.

**Structural sign limitation.** If \(\phi\ge q\), any rental plan is replicated by ownership using \(a_O=a_R-qPh\). At a binding rental cap, an owner can then expand housing, so the capped-renter branch is impossible. More generally, for a borrowing-constrained young owner let

\[
\zeta=\frac1x-\frac{\beta K}{qz}\ge0.
\]

The exact housing first-order condition is

\[
\frac{\alpha x}{h-\kappa n}
=p+(q-\phi)Px\zeta.
\]

Hence \(\phi>q\) gives a negative housing wedge: housing collateral helps households obtain current resources. Old households are unrestricted in this region because their positive net estate automatically satisfies the consolidated mortgage constraint. The source's positive young housing wedge cannot be extended to this region by merely admitting C ownership. An equally weighted redistributive planner might still favor the young, but that would be a different argument.

**Small optional generalization.** The core does not otherwise need \(\beta<q\). Replace its use in proving C nonsaving by the direct condition \(qY/(\beta K)\ge b\), and retain finance, both resource inequalities and the price bracket. All stationary allocation and dated fertility proofs continue to hold. This statement does not extend the policy or transition proof.

## 3. Full dated planner and financial settlement

At fixed fertility and individual promised wealth, equal current household weights give common adult consumption and adult space. Strict concavity yields exactly

\[
X=(\bar x+B_o)/2,\qquad S=(\bar s+\bar h_o)/2.
\]

Here the resource condition is strictly stronger than the needed comparison because \(B_o\ge fb+(1-f)L>\bar x\). Also

\[
s_C<\alpha x_C/p<\alpha L/p=s_L,\qquad
\bar s<\alpha\bar x/p<\bar h_o.
\]

Therefore aggregate young housing rises by \((\bar h_o-\bar s)/2\) per young household, and every C household receives \(X>x_C\) and \(S>s_C\), hence more total consumption and housing at its own fixed fertility. Total old housing need not exceed total young housing before the intervention; the relevant comparison is old housing versus young adult space.

The financial settlement is valid. For each household, keep promised wealth fixed and set

\[
\Delta a_i=-qP\Delta h_i^{\rm owned},\qquad
T_i=\Delta c_i+p\Delta h_i.
\]

Let intermediary rental inventory be \(I=\bar H-\sum_i h_i^{\rm owned}\). Its bond position can be written \(a_I=-qPI\). Then

\[
\Delta a_I=qP\sum_i\Delta h_i^{\rm owned}
=-\sum_i\Delta a_i.
\]

Aggregate financial claims and goods are unchanged, and \(\sum_iT_i=0\). The tax base is the same total stock, so public spending remains fixed. Assigning ownership to allocations above \(r\) is feasible because the stock is common and owner housing has no upper bound. The planner relaxes private finance; this is not a transfer-only implementation with the original limits intact.

The separate small compensated housing transfer also has the correct derivative. It establishes a dated fixed-fertility gain under the relaxed financing permissions, not that the full utilitarian allocation is Pareto superior.

## 4. Both dated fertility claims

At a fixed offered bundle, the fertility derivative is strictly decreasing in \(n\), and its derivatives with respect to both \(c\) and \(h\) are positive. Every C family's fixed-fertility planner bundle therefore leads to a strictly larger privately chosen \(n\), holding the bundle and continuation opportunity fixed.

For the joint planner, common fertility is a result, not an extra restriction. Conditional on mean fertility, resources depend only on that mean; equalizing fertility maximizes the sum of \(\vartheta\log n_i\). Consumption and adult space are also equalized. Thus (12) describes the full optimum among current parents.

The function

\[
H(x,s)=\left(\frac\chi x+\frac{\alpha\kappa}s\right)^{-1}
=\frac{xs}{\chi s+\alpha\kappa x}
\]

is concave on positive \((x,s)\). Its Hessian is

\[
\nabla^2H=-\frac{2\chi\alpha\kappa}
{(\alpha\kappa x+\chi s)^3}
\begin{pmatrix}s^2&-xs\\-xs&x^2\end{pmatrix}\preceq0.
\]

Private first-order conditions imply \(\bar n\le\vartheta H(\bar x,\bar s)\). Since \(X(\bar n)>\bar x\) and \(S(\bar n)>\bar s\), the planner's fertility derivative at \(\bar n\) is strictly positive. Its derivative with respect to \(n\) is

\[
-\frac{\vartheta}{n^2}
-\frac{\chi^2}{2X(n)^2}
-\frac{\alpha\kappa^2}{2S(n)^2}<0.
\]

It has a unique interior zero, establishing \(n^F>\bar n\). This part requires no small child-goods-cost restriction. It remains a dated parental welfare comparison, not a dynamic population optimum.

</details>


<details>
<summary>Policy and transition — independent review</summary>

# Independent review: stationary policy and local transition

Reviewed 2026-09-09. Source: `output/model/simplified_olg_amendments/oracle_essential_theory_response.md`, particularly lines 360–415 and Appendices B–C, lines 490–671. Scope is policy, demographic dynamics and the surprise-price boundary only. The allocation/planner theorem and literature were not independently adjudicated here. No numerical examples, sweeps or simulations were used as proof.

## Verdict

**PASS, conditional on the stated strict household regime.** The stationary rental-cap and fertility-taste signs, Appendix B's sufficient inequalities, the demographic Jacobian and Jury tests, and the announcement-price derivative are analytically correct. The latter derivative correctly corresponds to *fixed initial cohort sizes*; the first perturbed future rent satisfies a different initialization from a generic demographic-state perturbation.

**One explicit presentation repair is necessary for a usable transition statement:** a change in the fertility preference changes the liquid young's desired saving. Consequently date-zero old wealth differs from the new stationary old-wealth coefficient even with no capital gain. Equation (14), if used literally with the inherited bond and housing holdings, already includes this correction. It must not be replaced by new-parameter stationary old wealth plus the capital gain alone. The two-dimensional autonomous population map applies from date 1 after a permanent date-zero change; date 0 needs the separate initial balance-sheet equation below.

The evidence supports local feasible convergent paths and the signed cumulative-reproduction comparison. It does not give a signed impact-fertility effect, all-date ordering, global uniqueness or convergence, large reforms crossing regimes, or an endogenous inheritance/type-transmission model.

## 1. Definitions and the stationary rental-access derivative

Write \(\theta=\vartheta\), \(K=1+\alpha+\omega_B\), \(G=1+\alpha+\beta K\), and \(L=(W+qY)/(G+\theta)\). The mean consumption resources of the old are

\[
B_o=fY/K+(1-f)\beta L/q.
\]

Per-young-cohort and per-old-cohort housing demand are

\[
D_y=fr+(1-f)\left(\kappa n_L+\frac{\alpha L}{p}\right),\qquad
D_o=\frac{\alpha B_o}{p},\qquad D=D_y+D_o.
\]

The source's household implicit derivatives, equations (503)–(506), follow directly from differentiating its fertility first-order condition. In particular,

\[
n_{C,p}=-\frac{\chi r}{x_C^2Q},\qquad
n_{C,r}=\frac{\alpha\kappa/s_C^2-\chi p/x_C^2}{Q}.
\]

At fixed price, \(D_r=f\) and \(F_r=fn_{C,r}\). Thus, on the replacement root \(F(p^*,r)=1/\nu\),

\[
p_r^*=\frac{fn_{C,r}}{-F_p},\qquad
\frac{dD^*}{dr}=f+D_pp_r^*.
\]

Therefore \(N_r^*>0\) is equivalent to the source's inequality (15), \((-D_p)n_{C,r}>-F_p\).

### Exact verification of the dimensionless bounds

Use the source's \(z,c,s,m,u,a\), where \(s+m=1\), and define

\[
\widetilde Q=\frac{\theta}{m^2}+\frac{z^2}{c^2}+\frac{\alpha}{s^2},\qquad
B=\frac{\alpha}{s^2}-\frac z{c^2}.
\]

The cap and fertility conditions imply

\[
\frac\theta m=\frac zc+\frac\alpha s,\qquad
\alpha c>s,\qquad
m<\frac\theta{\alpha+\theta},\qquad
c>\frac1{\alpha+\theta}.
\]

The last two follow because the positive term \(z/c\) implies \(\theta s>\alpha m\). The resource condition implies \(a>\alpha(fc+2u)\) because \(b>x_C\) and \(f>0\); the weak bound used by the source is sufficient.

The difference in (15), multiplied by the positive factor \(\kappa p\widetilde Q/r\), is exactly

\[
J=aB-\frac{fz}{c^2}
-\frac{u\theta}{(1+z)^2}
 \left(\frac\theta{m^2}+\frac{z(1+z)}{c^2}\right).
\]

Condition (13) and \(p>p_\ell\) imply \(z<1/(2\alpha+\theta)\). In particular \(B>0\), because \(\alpha c>s\) gives \(\alpha c^2>s^2/\alpha>zs^2\).

After inserting the bound on \(a\), the coefficient of \(f\) is

\[
\alpha cB-z/c^2
>\frac{c(1-\alpha z)-z}{c^2}>0.
\]

For the final inequality, \(1-\alpha z>0\) and

\[
c(1-\alpha z)-z
>\frac{1-z(2\alpha+\theta)}{\alpha+\theta}\geq0.
\]

For the coefficient of \(u\), put \(v=\alpha c/s>1\). Multiplication by \(c^2\) gives

\[
2v^2-\frac{(z+v)^2}{(1+z)^2}
-2\alpha z-\frac{\theta z}{1+z}.
\]

The first two terms form an increasing function of \(v\geq1\), with value 1 at \(v=1\). Hence the displayed expression exceeds

\[
1-z(2\alpha+\theta)\geq0.
\]

Thus \(J>0\), \(p_r^*>0\) and \(N_r^*>0\). The equilibrium cash-poor fertility derivative is also correctly reported:

\[
\frac{dn_C^*}{dr}
=n_{C,r}\frac{(1-f)n_{L,p}}{F_p}>0.
\]

Old consumption resources are fixed in this comparison, so old housing falls with rent. The source's expression for \(pD_y\) increases in both \(p\) and \(r\), whereas \(pD_o=\alpha B_o\) is fixed, proving a larger young housing share.

## 2. Stationary taste sign

Direct differentiation gives

\[
n_{C,\theta}=\frac1{n_CQ}>0,\qquad
n_{L,\theta}=\frac{(W+qY)G}{(G+\theta)^2(\chi+\kappa p)}>0.
\]

Thus \(F_\theta>0\) and \(p_\theta^*>0\). To verify the important inequality \(-D_p>\kappa(-F_p)\), observe that it is equivalent to

\[
a>\frac{fz}{c^2\widetilde Q}.
\]

But

\[
\frac{z}{c^2\widetilde Q}
<\frac{zs^2}{\alpha c^2}<\alpha z<\alpha c,
\]

where the final inequality uses \(z<1/(2\alpha+\theta)<1/(\alpha+\theta)<c\). The lower bound on \(a\) proves the result.

Since \(L_\theta<0\) and \(B_{o,\theta}=(1-f)\beta L_\theta/q<0\),

\[
D_\theta<(1-f)\kappa n_{L,\theta}.
\]

Consequently

\[
\frac{dD^*}{d\theta}
<(1-f)\kappa n_{L,\theta}-\kappa F_\theta
=-\kappa f n_{C,\theta}<0,
\]

so \(N_\theta^*>0\). A permanent taste decline lowers the stationary population, conditional on remaining on this branch.

## 3. Demographic dynamics and stability

With constant post-announcement parameters and all later choices anticipated, the liquid young choose continuation resources

\[
z_{L,t+1}=\beta KL/q
\]

independently of current or future prices. Cash-poor young choose zero financial saving and have old resources \(Y\). Because entrant type shares are fixed independently of parental type, the two cohort sizes suffice after the special initial date.

Linearizing static clearing in log cohort sizes gives

\[
\widehat p_t=A_y\widehat N_t^y+A_o\widehat N_t^o.
\]

The source's demographic matrix is therefore exactly

\[
M=\begin{pmatrix}1-\varepsilon A_y&-\varepsilon A_o\\1&0\end{pmatrix}.
\]

For cash-poor fertility,

\[
\varepsilon_C=\frac{z}{mc^2\widetilde Q}
<\frac{z}{(1+z)c}<1.
\]

Indeed \(mc\widetilde Q=z+\alpha c/s+mz^2/c+\alpha mc/s^2>1+z\); the last bound follows from the inequalities for \(c,z\) above. Liquid fertility has elasticity \(1/(1+z)<1\), and aggregate elasticity is their positive fertility-weighted average.

Also \(-pD_p>D_o\), hence \(0<A_o<1\). Finally,

\[
D_y-D_o-\kappa F
=fs_C+\frac{\alpha((1-f)L-B_o)}p<0.
\]

Together with \(-D_p>\kappa(-F_p)\), this yields \(\varepsilon(A_y-A_o)<1\). The characteristic polynomial is

\[
\lambda^2-(1-\varepsilon A_y)\lambda+\varepsilon A_o.
\]

All three Jury expressions are strictly positive:

\[
\varepsilon(A_y+A_o)>0,\qquad
2-\varepsilon(A_y-A_o)>1,\qquad
1-\varepsilon A_o>0.
\]

Hence the nonlinear autonomous demographic map is locally asymptotically stable, with geometric convergence in a sufficiently small neighborhood.

## 4. The initial balance sheet, including the saving correction

Let \(\theta_-\) describe a pre-shock stationary economy and \(\theta_+\) the permanent new preference. Let \(L_-\) and \(L_+\) be the corresponding expenditure shares. Fix the inherited bonds \(a_{L,-1}\), house \(h_{L,-1}\), and both initial cohort sizes. Then

\[
z_{L,0}=\frac{a_{L,-1}}q+P_0h_{L,-1}
=\frac{\beta KL_-}q+h_{L,-1}(P_0-P_-).
\]

Thus, **relative to the new stationary old resources**, the full discrepancy is

\[
\delta z_{L,0}
=\frac{\beta K}{q}(L_--L_+)
 +h_{L,-1}(P_0-P_-).
\]

The first term is the required one-period saving-state correction. A taste decline raises \(L_+\), so inherited saving alone leaves the initial old liquid type below the new stationary old resources; capital gains are a separate term with their own sign.

The exact initial market equation is

\[
\bar H=N_0^yD_y(p_0;\theta_+,r_+)
+\frac{\alpha N_0^o}{Kp_0}
\left[fY+(1-f)\left(\frac{a_{L,-1}}q+P_0h_{L,-1}\right)\right].
\]

Equivalently, use new-parameter \(D_o\) and add \(\alpha N_0^o(1-f)\delta z_{L,0}/(Kp_0)\). In a derivative with respect to \(\theta_+\), the \(-\beta KL_{+,\theta}/q\) part of this correction cancels the apparent direct change in initial old resources introduced by \(D_o(p_0;\theta_+)\). Omitting it would give incorrect impact derivatives, even though the post-impact stability calculation is unaffected.

At date 1 the new liquid young's optimized continuation resources are already \(\beta KL_+/q\), so no persistent saving state remains for a permanent surprise. For a changing preference sequence, the ordinary old resources instead use \(L(\theta_{t-1})\), not \(L(\theta_t)\). At a later surprise along a baseline transition, replace \(P_-\) by the inherited contract's previously expected date-\(T\) price and retain the actual inherited portfolio. These are direct consequences of equation (14), not new economic assumptions.

## 5. Announcement-price derivative with fixed initial cohorts

For the derivative with respect to the candidate \(p_0\), hold all primitives and inherited cohort sizes fixed. Write \(x_t=dp_t/dp_0\), evaluated at a stationary reference. The initial cohorts do not change, so

\[
x_0=1,\qquad
\frac{d\log N_1^y}{d\log p_0}=-\varepsilon,\qquad
\frac{d\log N_1^o}{d\log p_0}=0,
\qquad x_1=-\varepsilon A_y.
\]

In particular, \(x_1\) is **not** \(1-\varepsilon A_y\). For later dates,

\[
x_t=-\varepsilon(A_y,A_o)M^{t-1}(1,0)'\quad(t\geq1).
\]

Using the bounded/no-bubble price solution and \(\rho=q/(1+\tau)\),

\[
\frac{\partial P_0}{\partial p_0}
=\frac1{1+\tau}\left[1-\rho\varepsilon(A_y,A_o)(I-\rho M)^{-1}(1,0)'\right]
=\frac{1-\rho}{(1+\tau)\Delta},
\]

\[
\Delta=1-\rho+\rho\varepsilon A_y+\rho^2\varepsilon A_o.
\]

This is exactly the source's expression and lies strictly between zero and \(1/(1+\tau)\).

At the same reference, initial static rent feedback from an exogenous \(P_0\) change is

\[
k_P=\frac{\alpha(1-f)h_L}{K(-pD_p)}.
\]

Since \(h_L=(L/p)[\alpha+\theta/(1+z)]\) and

\[
-pD_p>(1-f)\frac Lp\left[\alpha+\frac\theta{(1+z)^2}\right],
\]

\[
k_P<\frac\alpha K\frac{\alpha+\theta/(1+z)}{\alpha+\theta/(1+z)^2}
\leq\frac{\alpha(1+z)}K<1.
\]

Therefore \(1-k_P\partial P_0/\partial p_0>0\): the initial market equation and forward discounted-price equation have a locally unique solution. This is a local implicit-function argument with the inherited *portfolio* fixed; it does not freeze inherited market-valued wealth.

## 6. Conditions and the cumulative comparison

The local continuation requires strict actual mortgage feasibility for old owners and liquid young, strict cash-poor exclusion from owning at the rental cap, strict cap binding and zero saving for cash-poor renters, and positive consumption/estate margins. The weak stationary inequalities (4) alone do not ensure a neighborhood with unchanged old financing status; use the response's stated strict-financing qualification. Along nonstationary paths the relevant owner constraint is

\[
q(z_{t+1}-Y_{t+1})\geq(qP_{t+1}-\phi P_t)h_t,
\]

so it must be checked by continuity from the strict stationary regime rather than replacing \(P_{t+1}\) with \(P_t\). Bounded convergent prices exclude the explosive asset-price solution. Finite types make preservation of uniform strict margins immediate for sufficiently small paths; no global shock radius is supplied.

The cumulative-reproduction identity is exact when baseline and policy share \(N_T^y\): summing their log entry laws telescopes to

\[
\sum_{t=T}^{\infty}\log\frac{F_t^{\rm policy}}{F_t^{\rm baseline}}
=\log\frac{N^{*,\rm policy}}{N^{*,\rm baseline}}.
\]

Local geometric convergence makes this sum well defined. Its positivity follows from the already-verified stationary rental-access sign. This statement includes the common inherited initial old through the announcement boundary, but supplies no separate impact or all-date fertility sign.

</details>
