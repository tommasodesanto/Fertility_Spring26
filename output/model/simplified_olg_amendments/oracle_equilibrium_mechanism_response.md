Source: https://chatgpt.com/c/6aa17908-01cc-83ea-ac9b-6193f4658c66
Captured September 9, 2026 after the response finished (33m35s). Advisory output; verify before adoption.

## Verdict

**Yes. The existing model can deliver a useful primitive equilibrium proposition—but the mechanism is forced equity accumulation, not old households’ reluctance to move.**

There is a region with comparable young cash resources and old incomes in which **renters save, while a positive mass of owners borrow at the mortgage limit**. Rental limits bind for both ages. Owners enter old age with enough income and housing equity to choose larger homes, despite being free to resize. The specified dated planner then reallocates housing toward young owners.

Crucially, this can happen when \(\beta/q<1\): unrestricted, uncapped households would choose *smaller* old homes. Mortgage finance and rental segmentation therefore do substantive work, rather than merely accompanying an income-timing result.

The proposition below is sufficient, not sharp. Its important economic restriction is a **tight LTV ceiling relative to lifecycle resources**. Its technical cap bounds are conservative. I have not established that a particular calibration satisfies them, nor a general result for conventional mortgage LTVs.

# Proposed main proposition

Write \(\tau=\tau^p\), \(H_R=h_R^{\max}\), \(H_O=h_O^{\max}\), and specialize to equal adult housing weights, \(\alpha=\gamma\). Retain the original budgets, mortgage repayment convention, endogenous rebate, free old-owner resizing, and exogenous entry-wealth distribution. In particular, the estate is not the source of entering households’ wealth.

Define

$$
D=1-q+q\tau,\qquad p=DP,\qquad
\ell=\frac{1-\phi+q\tau}{D},
$$

and retain

$$
E=1+\alpha+\vartheta,\qquad K=1+\alpha+\omega_B.
$$

Let \(w_i^0=y_i^y+b_i\) and \(v_i^0=y_i^o\).

### Assumptions

**Comparable lifecycle resources, with genuine heterogeneity.** The endowment distribution has bounded support, \(w_i^0\ge\underline w>0\), and

$$
\rho_-w_i^0\le v_i^0\le\rho_+w_i^0,
\qquad 0<\rho_-<1<\rho_+.
\tag{1}
$$

Current income, entry wealth, and old income may all vary. No proportional-income restriction is required.

**Saving motives and sufficiently tight mortgage credit.** Assume \(0<\tau<2\) and

$$
\omega_BD>q\alpha,\qquad
\frac{\rho_+E}{K}<\frac{\beta}{q}<1,
\tag{2}
$$

together with the primitive LTV restriction

$$
\boxed{
0<\phi<
\frac{1+\rho_-(1+q\tau)-(K/\alpha)D}
{\rho_-+1/q}.
}
\tag{3}
$$

The right-hand side must therefore be positive. These restrictions imply \(\phi<q\).

**Small rentals and sufficient owner capacity.** The explicit, price-free bounds in Appendix A hold. They make rental units small enough that both renter caps bind, while ensuring that not every young owner is physically capped. They also establish existence of a positive stationary equilibrium—not merely the desired behavior conditional on a price.

### Conclusions

A positive stationary equilibrium exists. **Every** positive stationary equilibrium under these assumptions has the following properties:

* Every renter saves strictly, \(a_i'>0\), and occupies \(H_R\) in both ages, with strictly positive housing-cap multipliers.
* A positive mass of young owners are below \(H_O\) and borrow at the maximum:

  $$
  a_i'=-\frac{\phi}{q}Ph_i^y.
  $$

  Their financing multipliers are strictly positive.
* Owners satisfy \(h_i^o\ge h_i^y\), strictly for every owner whose young home is below \(H_O\). Consequently,

  $$
  \bar h^{o,eq}>\bar h^{y,eq}.
  \tag{4}
  $$

For the packet’s full dated, fixed-fertility planner,

$$
\boxed{
\bar h^{y,F}>\frac{\bar H}{2N}>\bar h^{y,eq}.
}
\tag{5}
$$

Thus aggregate young housing increases, a positive mass of young owners receive more housing, and the equally weighted remaining-utility sum strictly increases.

This is **not** a Pareto result, an implementation by the property tax, or a claim that every young household gains utility.

## Why the economic argument works

### 1. The LTV restriction produces the age ordering

At stationarity put \(w_i=w_i^0+T\), \(v_i=v_i^0+T\). Since the rebate is common and \(\rho_-<1<\rho_+\),

$$
\rho_-w_i\le v_i\le\rho_+w_i.
\tag{6}
$$

For an owner, define its actual resources entering old age by

$$
z_i=a_i'+Ph_i^y+v_i.
$$

The original financing constraint and current budget imply, **whether or not the household borrows maximally**,

$$
z_i\ge v_i+\frac{\ell-1}{q}ph_i^y,
\qquad
w_i\ge c_i^y+\ell ph_i^y.
$$

Therefore

$$
z_i\ge
\rho_-c_i^y+
\left[\rho_-\ell+\frac{\ell-1}{q}\right]ph_i^y.
\tag{7}
$$

Restriction (3) is exactly

$$
\alpha\left[\rho_-\ell+\frac{\ell-1}{q}\right]>K.
$$

It follows that \(z_i>Kph_i^y/\alpha\).

The estate condition in (2) makes the old financial-estate floor slack, including when the old physical cap binds. The exact old housing choice is consequently

$$
h_i^o=\min\left\{H_O,\frac{\alpha z_i}{Kp}\right\}.
$$

Hence \(h_i^o\ge h_i^y\), strictly when \(h_i^y<H_O\).

This uses the old household’s **optimal resizing decision**. It does not treat its previous home as a lower bound. Nor does it infer the ordering merely from positive equity: old income contributes through the explicit lower bound (1), addressing the precise obstruction identified in the checked assessment.

### 2. Strict mortgage borrowing follows; it is not assumed

For an uncapped young owner, suppose its financing multiplier were zero. Its first-order conditions would give

$$
c_i^o=\frac{\beta}{q}x_i,\qquad
s_i=\frac{\alpha x_i}{p},
$$

where \(x_i=c_i^y-\chi n_i\) and \(s_i=h_i^y-\kappa n_i\).

Old optimality, including a possibly binding physical cap, gives

$$
h_i^o\le\frac{\alpha c_i^o}{p}
=\frac{\beta}{q}s_i<s_i<h_i^y.
$$

That contradicts the ordering just proved. Therefore its financing multiplier is **strictly positive**, and its mortgage is maximal.

The capacity assumptions guarantee positive mass of such uncapped owners. No conclusion about the financing multiplier of every *capped* owner is needed.

### 3. Renting is genuinely different

Under (2), an uncapped renter would choose positive saving: its saving condition is

$$
\frac{v_i}{w_i}<\frac{\beta K}{qE},
$$

which follows from (6). Its unconstrained age profile would satisfy

$$
h_i^{o,R}=\frac{\beta}{q}s_i^R<h_i^{y,R}.
$$

The rental cap prevents that allocation. Under the stated size bounds, renters instead save while occupying \(H_R\) in both ages. This is not the \(\phi=q\) benchmark in which conditional real choices coincide.

Indeed, at \(\phi=q\), \(\ell=1\), and the key mortgage restriction would require \(\alpha\rho_->K\), which is impossible here. **The result cannot survive by silently reverting to that benchmark.**

These observations identify the role of each institution. They are not claims about the aggregate general-equilibrium comparative statics of removing either constraint.

## What the planner result means

The packet’s planner preserves fertility, tenure, future real opportunities and current old net estates, but reallocates **both current goods and housing**. Its solution is

$$
X=\frac{\bar x+\bar c^o}{2},\qquad
c_i^{y,F}=\chi n_i+X,\qquad c_i^{o,F}=X,
$$

$$
h_i^{y,F}=\min\{H_i,\kappa n_i+t\},\qquad
h_i^{o,F}=\min\{H_i,t\},
\quad t=\alpha/\lambda_H.
$$

These are the existing planner and welfare weights, not a replacement criterion.

Because children require space, each young planner household receives weakly more housing than its matched old counterpart. Positive unused physical capacity makes that aggregate inequality strict. Equation (4) puts competitive young housing below half the stock; the planner puts it above half. In fact,

$$
\bar h^{y,F}-\bar h^{y,eq}
>
\frac{\bar h^{o,eq}-\bar h^{y,eq}}2>0.
$$

**Individual housing gains are more selective.** Exactly,

$$
\Delta h_i^y
=\min\{H_i-h_i^{y,eq},\,t-s_i\}.
$$

Renters are already capped and cannot receive a housing increase at fixed tenure. A positive mass of young owners must therefore gain housing, but not necessarily all owners.

For a young household, the remaining-utility change is

$$
\Delta U_i^y
=\log\frac{X}{x_i}
+\alpha\log\frac{h_i^{y,F}-\kappa n_i}{s_i}.
$$

This is the individual welfare test; receiving more housing does not suffice when consumption falls.

The original settlement preserves the fixed commitments through owner financial-position adjustments and balanced current transfers, with the rental-intermediary account included. Both individual financial lower bounds may need relaxation, exactly as allowed in the packet.

Children make the planner’s young housing share exceed one half. The financing restrictions generate the opposite competitive age ordering in this region. Nevertheless, an equal-weight utilitarian gain can contain ordinary redistribution as well as the consequences of constraints; the theorem does not decompose those components into separate welfare amounts.

## Fertility: the two experiments remain distinct

**Parents rechoosing after fixed-\(n\) bundles.** Let

$$
s_i^F=h_i^{y,F}-\kappa n_i.
$$

Evaluating the private fertility first-order condition at the original \(n_i\) gives the exact sign test

$$
\operatorname{sgn}(n_i^S-n_i)
=
\operatorname{sgn}\left[
\chi\left(\frac1{x_i}-\frac1X\right)
+\alpha\kappa\left(\frac1{s_i}-\frac1{s_i^F}\right)
\right].
$$

Thus the housing theorem does not establish higher individual or average privately readjusted fertility. This is the same goods-versus-space tradeoff emphasized in the packet.

**The dated planner choosing fertility jointly.** That planner must reoptimize goods, space and fertility together. The proposition above does not prove its fertility direction. Indeed, the explicit bounds below imply \(H_R<\bar h^{y,eq}\), so the packet’s sufficient joint-planner proposition—which requires the reverse strict capacity comparison—does not apply. This is an unresolved implication of the new housing result, not evidence that joint fertility falls.

# Proof appendix

## A. Explicit primitive capacity and existence bounds

All quantities in this subsection are functions of primitives. None uses an equilibrium price, rebate, multiplier or tenure share.

Let

$$
W=\int w_i^0\,dF,\qquad V=\int v_i^0\,dF,\qquad
M=E+\beta K,\qquad x_*=\frac{\underline w}{M}.
$$

Define

$$
\pi_*=
\operatorname{logistic}
\left(\frac{\bar\xi-E\log\ell}{\sigma_\xi}\right)>0,
$$

$$
T_*=\frac{\tau(W+qV)}{(1-q)(2-\tau)},\qquad
p_*=\frac{\nu\vartheta(W+T_*)}{\kappa(\alpha+\vartheta)}.
$$

For the renter calculation, write

$$
B=1+\omega_B,\qquad b_\beta=\beta/q,\qquad A=1+\vartheta,
$$

and define

$$
\eta_R=
\min\left\{
\frac{b_\beta B-A\rho_+}{b_\beta B-A},
\frac{\alpha b_\beta(1+q\rho_-)}
{A+\beta B+\alpha b_\beta(1+q)}
\right\}.
$$

Assumption (2) implies \(B>A\) and \(b_\beta B>A\rho_+\); hence \(0<\eta_R<1\).

The promised capacity assumptions are

$$
\boxed{
\begin{aligned}
0&<\chi<\nu\pi_*\vartheta x_*,\\
H_O&>
\frac{\kappa}{\nu\pi_*}
+\frac{\alpha\kappa}
{\nu\pi_*\vartheta-\chi/x_*},\\
0&<H_R<
\min\left\{H_O,\frac{\eta_R\underline w}{p_*}\right\}.
\end{aligned}}
\tag{8}
$$

The first two inequalities ensure that an economy in which every owner is capped would produce more than replacement fertility. The last ensures that rental caps bind even at the upper bound on any stationary service price.

### A.1 A uniform consumption bound

Conditional on tenure, write the lifetime problem in variables \(x,s,n,c^o,h^o,e\). Its lifetime budget is

$$
x+ps+(\chi+\kappa p)n+qc^o+qph^o+q^2e=w+qv.
$$

Its cash constraint is

$$
x+L_ds+(\chi+\kappa L_d)n\le w,
\qquad L_R=p,\quad L_O=\ell p.
$$

Let \(\Lambda,\mu\) be the lifetime and cash multipliers. Multiplying the first-order conditions by the corresponding quantities and adding yields

$$
M=\Lambda(w+qv)+\mu w+\eta_yH_d+\eta_oH_d.
$$

The estate-floor term vanishes by complementary slackness because that restriction is homogeneous. Since \(1/x=\Lambda+\mu\),

$$
M\ge(\Lambda+\mu)w=\frac wx,
\qquad x\ge\frac wM\ge x_*.
\tag{9}
$$

This bound allows binding young and old caps.

### A.2 A uniform lower bound on ownership

Take the optimal renter allocation. Scale its current \(c,h,n\) by \(1/\ell\), and allocate the released lifetime resources to old age. The resulting owner cash expenditure satisfies

$$
\frac{c+\ell ph}{\ell}\le c+ph\le w.
$$

The original renter’s old bundle is feasible for an owner: \(H_O>H_R\), and

$$
h^{o,R}\le\frac{\alpha c^{o,R}}p,\qquad
e^R=\frac{\omega_Bc^{o,R}}q
>
Ph^{o,R}
$$

by \(\omega_BD>q\alpha\).

Consequently,

$$
W_i^O\ge W_i^R-E\log\ell,
\qquad \pi_i^O\ge\pi_*.
\tag{10}
$$

This is only a bound. Actual tenure probabilities remain the original logistic functions of the two optimized value functions.

## B. Exact renter regime

Fix \(p\le p_*\), \(T\ge0\), and a type with resources \(w,v\). Consider the candidate with both homes equal to \(H_R\) and an interior saving choice.

Set

$$
J=w+qv-(1+q)pH_R.
$$

The candidate satisfies

$$
(1+\beta B)x+\chi n=J,\qquad c^o=b_\beta x,
\tag{11}
$$

and

$$
\frac{\vartheta}{n}
=\frac{\chi}{x}+\frac{\alpha\kappa}{H_R-\kappa n}.
\tag{12}
$$

These equations have a unique solution with positive adult consumption and adult space. The cap bound ensures \(J>0\).

Equation (12) implies \(\chi n<\vartheta x\), so

$$
x>\frac{J}{A+\beta B}.
\tag{13}
$$

Old resources and saving are

$$
z=Bb_\beta x+pH_R,\qquad a'=z-v.
$$

Using (13),

$$
(A+\beta B)a'
>
b_\beta Bw-Av-(b_\beta B-A)pH_R.
$$

The first term defining \(\eta_R\), together with \(v\le\rho_+w\), makes this strictly positive.

The second term defining \(\eta_R\), together with \(v\ge\rho_-w\), gives

$$
pH_R<\alpha b_\beta x=\alpha c^o.
$$

Therefore the old cap multiplier is strictly positive. Since \(b_\beta<1\),

$$
pH_R<\alpha x,
$$

which also makes the young cap multiplier strictly positive:

$$
\frac{\alpha}{H_R-\kappa n}-\frac px>0.
$$

Thus the candidate satisfies the original budgets and all first-order and complementary-slackness conditions. Concavity proves optimality. This establishes **strict saving and both binding renter caps**, rather than inserting capped households into an uncapped formula.

## C. Endogenous taxes, prices and stationary existence

For any trial \(p>0,T\ge0\), solve both conditional household problems and integrate using their actual logistic tenure probabilities. Denote the resulting means by \(\bar n(p,T)\), \(\bar h^y(p,T)\), and \(\bar h^o(p,T)\).

The stationary equations are exactly

$$
\bar n(p,T)=1/\nu,\qquad
T=\frac{q\tau}{2D}p\bigl[\bar h^y(p,T)+\bar h^o(p,T)\bigr].
\tag{14}
$$

Once solved,

$$
P=p/D,\qquad
N=\frac{\bar H}{\bar h^y+\bar h^o}.
$$

The old distribution is generated by those same stationary household choices; it is not separately imposed. These are the original housing, tax and demographic clearing conditions.

### C.1 Rebate and price bounds

Adding the current and discounted old budgets gives

$$
p(\bar h^y+q\bar h^o)<W+qV+(1+q)T.
$$

Hence the rebate map \(G(p,T)\), defined by the right-hand side of the tax equation, satisfies

$$
0<G(p,T)<
\frac{\tau}{2D}\bigl[W+qV+(1+q)T\bigr].
$$

Because

$$
2D-\tau(1+q)=(1-q)(2-\tau)>0,
$$

this map sends \([0,T_*]\) into itself. Every stationary rebate also lies below \(T_*\).

The fertility first-order condition gives

$$
n_i<
\frac{\vartheta h_i^y}{\kappa(\alpha+\vartheta)}.
$$

Since \(ph_i^y<w_i\),

$$
\bar n(p,T)<
\frac{\vartheta(W+T)}
{\kappa(\alpha+\vartheta)p}.
\tag{15}
$$

At \(p=p_*\), fertility is therefore strictly below replacement for every \(T\in[0,T_*]\).

### C.2 Low prices and positive uncapped-owner mass

If a young household were uncapped, its housing first-order condition and (9) would imply

$$
s_i\ge\frac{\alpha x_i}{\ell p}
\ge\frac{\alpha x_*}{\ell p}.
$$

For sufficiently small positive \(p\), this exceeds \(H_O\). Thus all young households are physically capped at sufficiently low prices, uniformly over the rebate interval.

When an owner is capped, its fertility is bounded below by the root \(n_*\) of

$$
\frac{\vartheta}{n_*}
=\frac{\chi}{x_*}
+\frac{\alpha\kappa}{H_O-\kappa n_*}.
$$

The first two bounds in (8) imply

$$
n_*>\frac1{\nu\pi_*}.
$$

Using (10), low-price mean fertility exceeds \(1/\nu\).

Household choices and integrated tenure-weighted demands are continuous on a compact positive-price rectangle. Combine the rebate map with a projected price adjustment

$$
p\longmapsto
\operatorname{proj}\bigl[p+\varepsilon(\bar n(p,T)-1/\nu)\bigr].
$$

Brouwer’s theorem supplies a fixed point. The strict fertility signs exclude the price boundaries, so it solves both equations in (14).

The same argument proves that **not all owners can be young-capped at any stationary equilibrium**: otherwise their fertility contribution alone would exceed replacement. Therefore uncapped young owners have positive mass.

This proves existence, not uniqueness.

## D. Completing the housing and planner proof

The old optimum used above follows directly from

$$
c^o=\frac{z-ph^o}{B},\qquad
h^o=\min\left\{H_d,\frac{\alpha z}{Kp}\right\},
\qquad e=\frac{\omega_Bc^o}{q}.
$$

Indeed \(h^o\le\alpha c^o/p\), so for owners

$$
e\ge\frac{\omega_Bp}{q\alpha}h^o>Ph^o.
$$

Thus the owner estate floor is strictly slack, including at a binding physical cap.

The resource argument (7) now proves the owner age ordering. Appendix C supplies positive uncapped-owner mass; the multiplier contradiction in the main argument proves maximal borrowing for that mass. Appendix B gives equality of the two renter homes. Integrating proves (4).

Because some young owners are uncapped, competitive total housing is strictly below total physical capacity:

$$
\frac{\bar H}{N}<2\int H_i\,dQ.
$$

At the planner allocation, if every old household were capped, its housing formula would make every young household capped too, contradicting this strict inequality. Some old planner housing is therefore uncapped. Since every \(n_i>0\),

$$
\min\{H_i,\kappa n_i+t\}>\min\{H_i,t\}
$$

on a positive-mass set. This proves the first inequality in (5); (4) proves the second.

The equilibrium allocation is feasible for the specified planner, and its housing differs from the unique optimum. Strict concavity therefore gives a strict utility-sum improvement.

Finally, (15) before applying the cash bound implies

$$
\bar h^{y,eq}>
\frac{\kappa(\alpha+\vartheta)}{\nu\vartheta}.
$$

But (8), \(\eta_R<1\), and \(\underline w\le W+T_*\) imply

$$
H_R<
\frac{\kappa(\alpha+\vartheta)}{\nu\vartheta}.
$$

This verifies the stated failure of the packet’s sufficient joint-fertility capacity test.

## E. Analytical nonemptiness

No numerical reference equilibrium or continuity around one is needed.

Choose \(q\in(0,1)\), \(\tau\in(0,2)\), and \(\alpha>D\). Select

$$
0<\phi<
\frac{q}{1+q}\left(1-\frac D\alpha\right).
$$

Then the interval

$$
\frac{q\alpha}{D}
<
\omega_B
<
\alpha\left[\ell+\frac{\ell-1}{q}\right]-1-\alpha
$$

is nonempty.

Choose \(\omega_B\) in that interval, then

$$
0<\vartheta<\omega_B,\qquad
\frac EK<b_\beta<1,\qquad \beta=qb_\beta.
$$

Next choose

$$
\max\left\{0,
\frac{K/\alpha-(\ell-1)/q}{\ell}\right\}
<\rho_-<1,
\qquad
1<\rho_+<\frac{b_\beta K}{E}.
$$

These are nonempty intervals and establish (1)–(3).

Choose any bounded heterogeneous distribution satisfying the resource-ratio bounds, splitting \(w_i^0\) into positive current income and entry wealth as desired. For any finite logistic-taste parameters and any \(\nu>0\), choose positive \(\chi\), then sufficiently large finite \(H_O\), and sufficiently small positive \(H_R\), according to (8). Any \(\bar H>0\) is admissible.

This produces an analytical family with positive taxes, heterogeneous incomes and wealth, saving renters, constrained mortgage borrowers, active rental limits, and the proved planner housing gain.

---

**The author’s decision is whether to make tight-credit-induced equity accumulation the explicit illustrative mechanism.** The existing model supports that conditional story. It does not yet justify presenting larger old homes as a general consequence of ownership, accumulated equity, or children’s housing needs alone.
