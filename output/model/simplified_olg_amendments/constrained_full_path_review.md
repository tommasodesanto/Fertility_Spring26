# Constrained efficiency with complete initial claims and a full equilibrium path

September 5, 2026. Second bounded investigation. This report supersedes no historical result: the earlier four-group conditional review remains unchanged. No model, manuscript, calibration, or policy file was edited.

**Result.** A complete local constrained Pareto improvement can be proved for an explicitly completed version of the illustrative model **if the planner can commit to a sequence of separately balanced cash transfers and future lump-sum taxes**. Every future housing market clears, the mortgage share and physical segmentation remain unchanged, every household's lifetime utility is preserved, and the fully accounted initial outside claimant gains. Housing shifts toward young households. The proof has an exact infinite path and a contracting tail; it does not freeze future household states.

The extra fiscal permission matters. Around the same illustrative steady state, **one-time cash redistribution has no nearby Pareto-improving equilibrium with a converging continuation** on the verified branches. A nonzero continuation makes some future renter cohorts lose; an exactly unchanged continuation fails the complete resource test. This is a local no-go theorem, not a global statement about every parameterization or large reform.

| Claim | Status and scope |
|---|---|
| Full-path Pareto improvement with committed date-by-date transfers, including future taxes | PROVED under the explicit ownership/fiscal completion and local conditions below |
| Conditions hold in the saved illustrative steady state | VERIFIED by closed-form arithmetic, all original budgets, and the lead's independent original-choice optimization |
| One-time cash gifts at that steady state | Local no-go for nearby converging equilibria on the maintained regular branches; not a global theorem |
| Automatic constrained inefficiency of the incompletely specified source economy | NOT CLAIMED: initial ownership, fiscal jurisdiction, and transfer commitment must be selected |

## 1. What is completed and what remains private

At date 0, current old owners hold \(m h\) units of title in the stationary illustration, where \(m\) is their mass and \(h\) is the house each acquired when young. Assign the residual initial title

\[
R=\bar H-mh
\]

to an explicitly modeled passive claimant. The claimant values consumption or financial wealth, demands no domestic housing services, trades the external bond at price \(q\), and owns the intermediary net of its matched debt. It has enough cash or ordinary outside financing to settle a sufficiently small transfer. Initial promised estate payments are liabilities to be honored; they are not unassigned windfalls. If there are several residual claimants, the aggregate transfer account can be divided in fixed shares, leaving each at least as wealthy. Existing mortgage and bond creditors continue receiving their contractual returns.

This is an **ownership and fiscal completion**, not an assertion already established in the source. The source allows external goods/bond trade and competitive rental intermediation, but does not assign every residual initial title or specify taxation of its ultimate owner. The planner must be allowed to make transfers to and collect transfers from this claimant. Calling it an outside investor does not by itself give the domestic government that jurisdiction.

Two transfer classes are distinguished:

1. **One-time redistribution:** a balanced vector of cash transfers at date 0, before closing; no extraordinary transfers at any later date.
2. **Committed sequence:** at date 0 the planner announces cash transfers for every subsequent date and predetermined household identity. Each date balances separately, including the passive claimant. Later payments can be negative. Government borrowing is zero at every date.

Both classes preserve the original mortgage share \(\phi\), nonnegative saving, owner/renter size limits, inability of an old owner to buy new housing or let retained housing, and the estate-composition constraint. The property tax remains positive and is rebated under the original rule. Identity-based targeting requires observing the relevant baseline type/taste or an equivalent predetermined identifier; the source does not establish that informational permission. A young household's net pre-closing transfer \(\ell_t^Y\) enters its budget once and changes the closing constraint to

\[
(1-\phi)P_th_t\le b+\ell_t^Y.
\]

An old transfer \(\ell_t^O\) enters the original old budget. The ordinary rebate remains unavailable at closing. Transfer schedules use a household's predetermined baseline identity; they are **not conditional on subsequently choosing ownership**. The same scheduled transfers must be used when checking a tenure deviation.

**What the sequence adds.** Some future owner cohorts receive a grant when young and pay a lump-sum tax when old. The tax remains due if they choose renting instead. It is not a loan contract secured by the house, but it is compulsory future taxation. Economically, commitment to that tax gives the planner a way to shift resources across a household's life. This permission is absent from the one-time-gift class. It cannot be described as merely distributing an ordinary current rebate. The passive claimant also advances resources at some dates against later fiscal receipts, using private outside finance. Its intertemporal budget and solvency are checked below.

## 2. Baseline and the local conditions

Keep every household's own baseline fertility fixed, with \(n^O\) for the baseline owner group and \(n^R\) for the baseline renter group. Keep all cohort masses fixed. Let each cohort have owner mass \(m>0\) and renter mass \(m_R>0\). The baseline has identical income and liquid wealth within each group, as in the saved illustration; idiosyncratic ownership tastes retain mixed tenure.

The relevant baseline branches are:

- Young owners have a strictly binding down-payment limit, positive saving, and a slack physical owner-size limit.
- Their old-age retention and estate constraints are slack.
- Young and old renters are at the rental size cap, with strictly positive saving and consumption.
- The corresponding conditions have strict local slackness or strict multipliers.

Use the original preference weights and define temporary coefficients

\[
A=1+\gamma+\omega_B,\quad
\rho=1+\beta A,\quad B=1+\beta(1+\omega_B),\quad
\eta=\frac{1+\omega_B}{A}.
\]

Let \(h,x,s=h-\kappa n^O\) be baseline young-owner housing, adult consumption, and adult space. Let \(h^O,c^O\) denote baseline old-owner housing and consumption. Let \(u=q r\) be the present-value housing user cost. The original interior saving and old-demand conditions give

\[
c^O=\frac{\beta x}{q},\qquad
h^O=\frac{\gamma\beta x}{qu},\qquad
e^O=\frac{\omega_B\beta x}{q^2}.
\]

The two substantive sufficient inequalities are

\[
MV^Y\equiv\frac{\alpha x}{s}>u,
\qquad
\lambda\equiv\frac{\alpha h^O}{\rho s}<1. \tag{2}
\]

The first is strict financing rationing in this branch. The second will make the compensation-induced adjustment die out. It is verified, not assumed numerically, for the saved example below. Global validity across different active constraints is not claimed.

## 3. Exact prices and allocations for all dates

Let \(a_P=1+q\tau^p\), and write the stationary price relation as \(u=(a_P-q)P\). For a sufficiently small \(\epsilon>0\), set

\[
u_0=(1+\epsilon)u,\qquad u_t=u\quad(t\ge1),
\]

and

\[
P_0=P+\frac{u}{a_P}\epsilon,\qquad P_t=P\quad(t\ge1),
\qquad
T_t=\frac{q\tau^pP_t\bar H}{2(m+m_R)}. \tag{3}
\]

Thus the user-cost relation and the original property-tax budget hold exactly. Only the initial house price changes. The future **allocations**, however, change and must be constructed.

Compensate initial old owners exactly at the new user cost:

\[
c_0^O=c^O(1+\epsilon)^{\gamma/A},\quad
h_0^O=h^O(1+\epsilon)^{-\eta},\quad
e_0^O=e^O(1+\epsilon)^{\gamma/A}. \tag{4}
\]

Their original utility is unchanged. Housing clearing determines initial young-owner housing:

\[
h_0=h+h^O-h_0^O.
\]

Preserve its lifetime utility by choosing adult consumption

\[
x_0=x\left(\frac{h_0-\kappa n^O}{s}\right)^{-\alpha/\rho}. \tag{5}
\]

For every \(t\ge1\), set recursively

\[
\begin{aligned}
h_t&=h+h^O\left(1-\frac{x_{t-1}}x\right),\\
x_t&=x\left(\frac{h_t-\kappa n^O}{s}\right)^{-\alpha/\rho}.
\end{aligned} \tag{6}
\]

Each owner cohort's old choices are

\[
c_{t+1}^O=\frac{\beta x_t}{q},\qquad
h_{t+1}^O=h^O\frac{x_t}{x},\qquad
e_{t+1}^O=\frac{\omega_B\beta x_t}{q^2}. \tag{7}
\]

Renters retain all their original real allocations. Equations (4), (6), and (7) give

\[
m(h_t+h_t^O)+2m_Rh_R^{\max}=\bar H
\]

at **every** date. No future housing demand is discarded or frozen.

The owner cohort's original utility, after substituting its optimal old allocations and keeping future user cost \(u\) fixed, differs from a constant only by

\[
\rho\log x_t+\alpha\log(h_t-\kappa n^O).
\]

Equations (5)–(6) preserve that expression exactly. The flow utilities and lifetime value functions themselves remain the original ones.

**Infinite-tail proof.** Let \(X_t=x_t/x\) and \(D=h^O/s\). The recursion is

\[
X_t=[1+D(1-X_{t-1})]^{-\alpha/\rho}.
\]

For \(0<X\le1\), its derivative is positive and at most \(\lambda<1\). Since \(X_0<1\), the path has \(X_t\uparrow1\), with

\[
0<1-X_t\le\lambda^t(1-X_0).
\]

Consequently young housing is higher, old housing is lower, and all deviations converge geometrically to zero. The estate/title/slackness conditions survive for sufficiently small \(\epsilon\). This is an actual infinite path, not an imposed terminal state or an IFT that assumes an improving direction.

## 4. Cash transfers implementing those private choices

Write \(\delta=1-\phi\). For each owner cohort, pay the pre-closing transfer

\[
\ell_t^Y=\delta P_th_t-b. \tag{8}
\]

It makes the original down-payment constraint bind at the desired privately optimal housing choice. Its saving is determined by the original budget:

\[
a'_t=y+b+T_t+\ell_t^Y-\chi n^O-x_t-(\delta+q\tau^p)P_th_t.
\]

At date \(t+1\), levy or pay the precommitted old-age transfer

\[
\ell_{t+1}^O
=A\frac{\beta x_t}{q}
-\left[R_fa'_t-R_f\phi P_th_t+P_{t+1}h_t+T_{t+1}\right]. \tag{9}
\]

This is simply the amount required by the original old budget to support (7). It is a stated fiscal instrument, not an unrecorded extra mortgage. Positive saving supplies the original Euler condition; the strictly positive housing multiplier remains positive locally; all other original KKT conditions hold. Thus households choose these allocations at the prices in (3).

For the initial old-owner state \((a,H=h)\), use

\[
\ell_0^O=A c_0^O-a-P_0h-T_0.
\]

Both initial renter groups receive

\[
\ell_0^R=(u_0-u)h_R^{\max}-(T_0-T).
\]

All later renter transfers are zero. Young renters' saving and lifetime allocations, and initial old renters' utility, therefore stay unchanged.

Finally set the passive claimant's transfer \(L_t\) to minus the sum of domestic transfers at date \(t\). Government borrowing is exactly zero. In particular,

\[
\begin{aligned}
L_0&=-m(\ell_0^Y+\ell_0^O)-2m_R\ell_0^R,\\
L_t&=-m(\ell_t^Y+\ell_t^O),\qquad t\ge1.
\end{aligned} \tag{10}
\]

This schedule includes future taxes. For cohorts \(t\ge1\), its first-order old-age tax is negative:

\[
q\,d\ell_{t+1}^O=(u-\delta P-MV^Y)\,dh_t<0.
\]

They gain housing and reduce nondurable expenditure while preserving lifetime utility. It is precisely this future fiscal incidence that a one-time transfer cannot reproduce.

## 5. The full resource and claimant account

The claimant's date-0 wealth change, including every residual initial title, is

\[
\Delta W_I=R(P_0-P)+\sum_{t\ge0}q^tL_t. \tag{11}
\]

Sum the original dated budgets, including estates and mortgage repayments, use housing clearing and the original rebate rule, and telescope the housing title terms. This gives the exact alternative expression

\[
\boxed{\displaystyle
\Delta W_I=-m\left[(1+\omega_B)(c_0^O-c^O)
+B\sum_{t\ge0}q^t(x_t-x)\right].} \tag{12}
\]

The coefficient \(B=1+\beta(1+\omega_B)\) counts a young cohort's current adult consumption, discounted old consumption, and discounted estate expenditure. It excludes housing payments, which are transfers once every title owner is included. Child goods are fixed and cancel. Property-tax payments and rebates cancel. This identity explicitly includes the initial title revaluation that invalidated the previous four-group exercise.

The initial old require more goods and estate expenditure when compensated for downsizing. The current and future young owners use fewer goods. Differentiate the exact path at zero:

\[
\frac{dh_0}{d\epsilon}=h^O\eta,
\quad
\frac{dx_0}{d\epsilon}=-\frac{MV^Yh^O\eta}{\rho},
\quad
\frac{dx_t}{d\epsilon}=\lambda^t\frac{dx_0}{d\epsilon}.
\]

Using \(h^O=\gamma\beta x/(qu)\) in (12) yields

\[
\boxed{\displaystyle
\Delta W_I'(0)
=\frac{m h^O\eta\,(MV^Y-u)}{1-q\lambda}>0.} \tag{13}
\]

The contraction supplies a summable uniform bound, justifying differentiation of the discounted series. Hence sufficiently small finite reforms leave **every household exactly indifferent** and give the fully accounted claimant strictly higher wealth. It can purchase extra consumption while preserving its original consumption plan. This is a Pareto improvement for the completed economy.

**A strict gain for young households is also feasible.** Fix one such sufficiently small \(\epsilon>0\). Raise \(x_0\) slightly above the indifference value in (5), keeping \(x_0<x\) and \(h_0\) unchanged, and continue to use (6) for every later cohort. The initial young-owner group now gains strictly. All later owners and all initial old households remain indifferent. Because \(x_0<x\), every subsequent housing shift remains toward the young. The positive claimant surplus survives by continuity of the uniformly contracting discounted series. Equations (8)–(10) implement the adjusted path as before. Thus the strict Pareto beneficiary need not be only the outside claimant.

**Solvency.** The transfer deviations decay geometrically. After initial titles are valued, the claimant's remaining fiscal obligations have finite market value. Its incremental marketable wealth can finance the remaining transfer stream at the existing world bond price; the required residual position tends to zero geometrically. Thus there is no unpaid terminal liability and no rollover based on population growth. The usual discounted no-Ponzi limit is zero. An unconstrained outside claimant with a positive liquidity buffer can implement the small perturbation. If the author instead forbids outside bridge finance or future fiscal transfers to this claimant, this implementation is not available; datewise treasury balance alone would not resolve that restriction.

**Optional real cost.** If each date uses \(C(h_t-h)\) goods per owner, with \(C(0)=0\) and \(C'(0)=K\ge0\), charged to the residual account, (13) becomes

\[
\Delta W_I'(0)=\frac{m h^O\eta\,(MV^Y-u-K)}{1-q\lambda}.
\]

Here the cost base is the date-specific extra young housing relative to baseline. This is an explicit optional convention, not a positive fixed moving cost or automatically the cost of changing housing between successive reform dates. The zero-cost theorem does not require choosing that extension.

## 6. Why one-time gifts fail locally in this illustration

The result in this section uses the same full ownership completion but permits transfers only at date 0. It applies to nearby converging equilibrium paths whose original strict branches remain valid. It does not exclude large changes that leave this neighborhood.

From date 1 onward, a young owner receives no extraordinary transfer. Let \(b/\delta\) be its fixed purchase-cash numerator. Its housing and adult consumption satisfy

\[
h_t=\frac{b}{\delta P_t},\qquad
\rho x_t=y+b-\chi n^O+T_t+qT_{t+1}-u_th_t.
\]

Old-owner housing at date \(t\) is \(\gamma\beta x_{t-1}/(qu_t)\). With both renter groups capped, these equations and housing clearing form an autonomous two-dimensional map for \((P_{t-1},P_t)\), smooth around the baseline.

To show the relevant directions, write \(p_t=dP_t/P\), \(X_t=dx_t/x\), and \(d=u/P=a_P-q\). Define

\[
A_x=\frac{T-qPh}{\rho x},\qquad
B_x=\frac{q(T+Ph)}{\rho x}.
\]

Then

\[
X_t=A_xp_t+B_xp_{t+1},
\]

and the linearized equilibrium recurrence is

\[
\frac qd p_{t+1}
+\left[B_x-\frac{a_P}{d}-\frac h{h^O}\right]p_t
+A_xp_{t-1}=0. \tag{14}
\]

At the saved baseline, its roots are

\[
3.30801404998,\qquad -0.05518475204.
\]

The second is the unique stable root. A nonconstant nearby converging path therefore eventually alternates around the baseline price. This conclusion is not limited to a first-order candidate: the exact smooth map has a one-dimensional local stable manifold, and its restriction has a strictly negative derivative at the stationary point. Every nonzero sufficiently close converging orbit on that manifold changes sides at successive dates.

For a future renter cohort, both housing quantities are fixed at the rental cap. Its lifetime utility increases with its adult consumption, whose budget change is

\[
B\,dx_t^R=dT_t+q\,dT_{t+1}
-h_R^{\max}(du_t+q\,du_{t+1}).
\]

Along the stable direction with root \(z=-0.05518475204\), this becomes

\[
B\,dx_t^R
=-\underbrace{[Ph_R^{\max}(a_P-qz)-T](1+qz)}_{0.1521144848>0}\,p_t.
\]

The nonzero derivative means exact renter welfare changes sign with the stable-manifold coordinate locally. A nonzero converging tail therefore hurts infinitely many future renter cohorts. It cannot be a Pareto improvement.

**The exactly zero-tail exception also fails.** The local map is invertible because \(A_x\ne0\). If a path reaches the stationary price pair exactly, the no-transfer recurrence forces all preceding post-reform price pairs back to the baseline. Thus \(P_t=P\) for every \(t\ge1\). Housing clearing at date 1 then fixes the mean adult consumption of the date-0 owner cohort at \(x\), since its old housing is proportional to that consumption.

If every date-0 young owner weakly gains, strict concavity/Jensen applied to \(\rho\log x_i+\alpha\log(h_i-\kappa n^O)\), with mean \(x_i=x\), requires mean young housing to weakly increase. Capped renters cannot change housing. Current old owners must therefore weakly reduce housing. If \(u_0<u\), exact old-owner compensation requires every old owner to increase housing, a contradiction. If \(u_0>u\), every old owner preserving utility requires strictly greater consumption and estate spending. Mean young-owner nonhousing expenditure cannot fall, because its mean \(x_i\) is fixed, and Pareto-preserved capped renters cannot supply a reduction. The complete present-value resource account then makes an initial claimant lose. At \(u_0=u\), strict concavity and budget balance leave no strict gain. This rules out the remaining finite-pulse case.

The no-go relies on actual signs verified for this illustrative baseline, rather than a generic incomplete-markets theorem. A different active constraint, heterogeneous income distribution, nonconvergent distant path, or altered instrument class needs a new argument.

## 7. Numerical receipt and three checks

The numerical input is row 0 of `output/model/simplified_olg_amendments/transition_phi80.csv`. Keeping its fertility and group masses fixed, stationary budget and market identities give \(P=0.3426835629750788\), only \(8.41\times10^{-13}\) below the saved rounded price. This removes harmless source residuals before the algebra check.

| Object | Value |
|---|---:|
| \(q,\tau^p,\phi\) | 0.5, 0.05, 0.8 |
| \(\alpha,\beta,\gamma,\omega_B,\chi,\kappa\) | 0.4, 0.4, 0.3, 0.4, 0.15, 0.5 |
| \(y,b\), rental/owner size limits | 1, 0.06; 0.45, 2 |
| \(m,m_R\) | 1.0893612826199552, 0.3659026465034622 |
| \(n^O,n^R\) | 0.5494856788378992, 0.3526717883084860 |
| \(u,T\) | 0.1799088705619164, 0.0058869658643552 |
| \(h,h^O,x,x^R\) | 0.8754432147123938, 0.6581964003552694, 0.4933973791493143, 0.5733917903339527 |
| Initial old \(a^O,a^R\) | 0.3651334697787124, 0.7172708310625342 |
| Residual initial title \(R\) | 1.0463260567599697 |
| \(MV^Y,\lambda,\Delta W_I'(0)\) | 0.3285480745093845, 0.2608845193070863, 0.1009349328350027 |
| One-time-gift \(A_x,B_x\) | -0.1738589858052817, 0.1845121018206285 |

For \(\epsilon=0.001\), selected exact path entries are:

| Object | Date 0 | Date 1 | Date 2 |
|---|---:|---:|---:|
| Price | 0.34285908382441 | 0.34268356297508 | 0.34268356297508 |
| Young-owner housing | 0.87598476505483 | 0.87558441801789 | 0.87548004710933 |
| Young-owner adult consumption | 0.49329153028593 | 0.49336976883710 | 0.49339017631968 |
| Owner young grant | +0.00006786679817 | +0.00000967761037 | determined by (8) |
| Owner old transfer at that date | -0.00003830717044 | +0.00001263795265 | -0.00006131841536 |
| Passive claimant transfer | -0.00008924073830 | -0.00002430971035 | +0.00006404795514 |

Initial old housing falls to \(0.6576548500128326\). Each initial renter receives \(0.00007794371655\). The claimant's initial title gain is \(0.00018365203816\); after all dated fiscal transfers its total wealth gain is \(0.00010067043878\). Independently summing the real-consumption account (12) gives the same number within \(7.0\times10^{-17}\). At \(\epsilon=0.0001\), the gain is \(0.00001009084427\), and gain divided by \(\epsilon\) approaches the positive analytical derivative.

Three distinct verification passes were completed:

1. **Derivation:** exact allocation recursion, original household implementation, complete title/estate telescoping identity, positive resource derivative, contraction tail, and the separate one-time-gift obstruction.
2. **Original-equation audit:** the closed-form path was reconstructed from saving and both original budgets, not from a frozen continuation value. Across 129 dates, maximum budget and housing errors were \(2.22\times10^{-16}\) and \(4.44\times10^{-16}\); owner lifetime utility errors were at most \(4.44\times10^{-16}\). Cash, retention, estate, and saving restrictions remain feasible. Estate slack exceeds \(0.0902\) in the larger perturbation. The code initially had an end-of-array indexing error in the diagnostic's final estate check; extending that diagnostic array corrected it before any result was accepted. No model solver or long computation ran.
3. **Adversarial and independent checks:** the full claimant cash-flow account and the goods account were computed separately; the first-order limit matches. The lead independently implemented three finite reforms, verified the infinite-tail bound, and optimized the original dated household problems for both tenure choices and both baseline identities. Its twelve optimization checks matched the proposed policies within \(1.8\times10^{-7}\), with original budget/fiscal/market residuals below \(9\times10^{-16}\). The lead's receipts are `planner_benchmark_checks.json` and `committed_cash_path.csv` in this output folder.

**Tenure is checked.** When testing a deviation, the household keeps its own frozen \(n^O\) or \(n^R\), and receives its own announced grants and future taxes under either tenure. The original marginal taste type then has a strict frozen-fertility margin because the two baseline fertility optima differ. In the finite checks the minimum margin is \(0.0357278\), so all logistic taste types retain their baseline tenure for a sufficiently small reform. This is not an assumption that the household may ignore an advantageous alternative tenure.

The useful remaining author decisions are now specific: select the ownership/fiscal completion and the information available for predetermined identity-based transfers; decide whether commitment to future lump-sum taxes and transfers to the passive claimant is permitted; and choose whether to retain the zero-cost theorem or add an explicitly defined real adjustment cost. Under those stated permissions, the full local positive result is established. Under a one-time-gift interpretation, the saved active regime instead has the local no-go result above.
