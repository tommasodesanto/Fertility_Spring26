# Estate and entry closure: decision for September 26

**Status:** lead recommendation for the author's decision today. No adoption,
target replacement, production change or new numerical run is made by this
sheet. Today's deadline supersedes the earlier Sunday working deadline.
This uses completed evidence only. The mortality decision remains with the
main integration task and must supply the household exit probabilities below.

## Recommended baseline

Use **liquidated net estates, uncertain adult receipts, an empirically measured
entry distribution, and an explicit external financial account**. Preserve the
current bequest taste function while changing its wealth argument consistently
to liquidatable resources. Replace the comparison of all model estates with
the child-directed empirical target by an explicitly constructed, age-weighted
wealth moment. This would not be an independently observed estate-flow target.

| Object | Decision recommended today | Evidence and approximation |
|---|---|---|
| Estate resources and utility | Liquidate the house at death, subtract the sale wedge once, settle liabilities from the estate, and apply the current bequest taste function to the remaining positive estate. Retain the externally fixed bequest shift and estimate the strength parameter. | Net valuation alone had small effects in the completed three-case experiment. This is a liquidation assumption: an in-kind inherited house is not modeled. |
| Donor and recipient coverage | Include positive estates from parents and childless households. Allocate the entire private estate pool to modeled adult households; fix estate taxation and non-household diversion to zero in this baseline. | Childless people can leave wealth to other relatives or unrelated recipients. Full household allocation and zero diversion are explicit simplifying restrictions, not estimated shares. Do not claim own-parent/own-child matching. |
| Receipt risk | At each supported age, use a zero/positive receipt lottery with the published probability and conditional mean, scaled to the available pool. Draw independently across receiving periods. | The matched experiment and its grid refinement show substantially different behavior from guaranteed payments. Positive-amount dispersion, repeated family links and a lifetime receipt cap are omitted. |
| Age and income allocation | Pool published usual-income groups using their 50/40/10 population shares. Pay at model ages 26, 30, ..., 78; set receipts at 18, 22, 82 to zero. | The Fed profile covers ages 25–80. Its usual-income groups are not model labor-income states, especially after retirement. Age pooling and unsupported-age exclusion are adopted approximations if approved, not empirical findings. |
| Timing | Pay at the start of the receiving period, before interest and choices, from estates generated at the end of the preceding period. A receipt is wealth, not additional earnings. | The tested wealth-jump operator uses the same interpolation weights in expectations and population transport. A dated funding path still needs implementation. |
| Entry wealth | Use the existing ages 18–24 childless-renter PSID sample, nonhousing net worth divided by the survey-wave mean annual gross working earnings, and the empirical joint distribution with earnings ranks. Preserve signed wealth and tied earnings ranks. | This common scale avoids dividing by very low own earnings and aligns with model annual-earnings units. The inherited B15 marginal/rank coupling is not an estimated empirical joint law. |
| Entry funding | Treat the measured initial financial position as resources and liabilities accumulated before the modeled adult life. Enter it once through an explicit external childhood/entry account; do not add an estate grant at entry or subtract later inheritances from the measured stock. | The retained entrant mean is 0.187 annual-earnings units; it contains both assets and debts. The common-scale candidate mean is 0.186. Similar means do not establish similar distributions. |
| Residual debt at death | Creditors recover what liquidation permits. Unpaid residual liabilities are written off to the external creditor account; they are not deducted from unrelated positive estates and are not inherited as negative lottery prizes. | All three retained estate cases have zero negative net estates. Positive future losses would be an explicit external credit subsidy at the maintained exogenous interest rate; this baseline does not claim competitive, zero-profit lending. |

The entry recommendation is a **distributional definition**, not approval of
the current 3-by-5 compression. That saved compression explains only 20.9% of
the raw weighted wealth variance and creates 0.277% mass at additional joint
combinations. Preserve the empirical weighted distribution through direct
projection onto the wealth grid and probability-rank overlaps with the B15
income nodes; pool exact earnings ties before splitting their probability
intervals. This introduces no new persistent state or estimated parameter.
The first model household-age node remains 18. The split 16/20 birth-entry queue
approximates the delay to that node; it does not turn this into measured wealth
at literal ages 16 and 20.

## Resource account and dated law

The signed estate is remaining financial wealth plus the house's sale value
after the sale cost. For death probability $d_{j,t}$ and post-choice mass
$g_t(x)$, define
\[
e_t(x)=b'_t(x)+(1-\psi)q_t h_t(x),\qquad
D_t^+=\sum_x d_{j,t}g_t(x)\max(e_t(x),0),\qquad
L_t=\sum_x d_{j,t}g_t(x)\max(-e_t(x),0).
\]
$D_t^+$ is the private recipient pool; $L_t$ is the creditor loss. If $C_t$ is
the sale cost at death, the signed gross-estate identity is
\[
D_t^{\mathrm{gross,signed}}=C_t+D_t^+-L_t.
\]
The current taste function receives $\max(e_t,0)$. Keep its existing
child-count dependence and level convention explicit. In particular,
subtracting a zero-estate utility value is an economic change when that value
depends on the number of children; it must not be smuggled in as a numerical
normalization. Do not add parent gating or per-child division in this estate
closure decision. Those would change the fertility incentives separately.

Let $p_j$ and $\mu_j$ be the published-profile probability and relative expected
receipt for model age $j$, and let $g^-_{t+1}$ be the surviving population
before receipts. The receiving-date scale and individual lottery are
\[
\lambda_{t+1}=\frac{D_t^+}{\sum_x g^-_{t+1}(x)\mu_j},\qquad
X_{j,t+1}=\begin{cases}
\lambda_{t+1}\mu_j/p_j&\text{with probability }p_j,\\
0&\text{otherwise.}
\end{cases}
\]
Thus total expected payments at $t+1$ equal $D_t^+$. Zero-probability ages have
zero payment. A positive pool with no eligible recipient mass is an explicit
failure, not a permission to discard wealth. Under a stationary distribution,
this becomes the joint price/estate fixed point. During a transition, both the
backward continuation from $t$ and forward advance into $t+1$ must use
$\lambda_{t+1}$ and the receiving-date law. Current stationary hooks have not
verified that dated indexing.

For entrant mass $M_t$ and initial wealth law $F_0$, record the external signed
entry position $E_t=M_t\int b_0\,dF_0$, together with its separate positive and
negative components. In the retained control those period components are
$+0.032$ and $-0.020$, with signed position $+0.012$ in model units. All death
estates are allocated through $D_t^+$; entry resources and any creditor losses
remain separate external accounts. This specifies the intended household-boundary
accounts; it does not yet verify the aggregate resource constraint.
A policy readout must report changes in
$E_t$ and $L_t$ alongside its tax account and household outcomes, so changes in
external support are visible.

The funding counterpart is an **outside donor/guarantor with resources outside
the modeled household sector**. It assigns the entrant financial position and
books its opposite position externally. External creditors retain the claims
associated with entrant debt. A loss guarantee pays creditors $G_t=L_t$ from
the outside donor's endowment; neither domestic taxes nor unrelated estates
finance it. The external financial sector supplies funds at the maintained
rate, and its net claims must reconcile to the opposite of modeled household
net financial positions. The implementation must carry the external balance
sheet and reconcile its changes with entry assignments, interest, ordinary
asset settlements and the guarantee. An exogenous interest rate alone does
not establish any of these identities. Policy welfare under this closure is
conditional on outside funding and excludes the donor's welfare. Approval
therefore accepts an economically substantive external resource assumption.

## The bequest calibration target needs an explicit replacement decision

The author-adopted 0.729% is **not an all-estate statistic**. Its saved builder
multiplies positive SCF 2007 `NETWORTH - TRUSTS` by head mortality and child
shares of 25% for married heads and 75% otherwise, and divides by signed
aggregate net worth. It mechanically applies those shares to all sampled
families because current rostered children do not identify lifetime offspring.
It also attributes whole-family wealth at head death. Consequently, neither
paying all model estates nor restricting model donors to $n>0$ reproduces its
definition. Multiplying by a guessed common child share does not fix this.

**Recommendation:** retain 0.729% as the documented child-directed empirical
proxy, and replace its calibration row with a **mortality-weighted positive
wealth / aggregate wealth moment**. Its interpretation is wealth held at ages
with higher exit risk, not independently observed estate transmission. Use
the same SCF 2007 sample, five implicates, survey weights, trust exclusion and
signed wealth denominator. For fixed empirical weight $w_i$, wealth $W_i$,
trust wealth $T_i$, and externally specified four-year exit weight $\bar d(a_i)$,
the proposed data object is
\[
m^{\mathrm{data}}=\frac{\sum_i w_i\bar d(a_i)\max(W_i-T_i,0)}
{4\sum_i w_iW_i}.
\]
The entire data-side object, including age bins, exit weights, coverage and
oldest-cell treatment, must be frozen before estimation. Its model counterpart
uses the same fixed age weights on endogenous household mass and wealth,
$m^{\mathrm{model}}(\theta)=\sum_x g_\theta(x)\bar d(a_x)
\max(W_\theta(x),0)/(4\sum_x g_\theta(x)W_\theta(x))$.
Never use fitted model wealth, population shares or parameter-dependent death
weights to rebuild the data target. Sharing an externally fixed age weight
does not by itself make the moment circular; choosing its weights to obtain
the desired fitted parameter would. The model has no separate trust asset;
the data numerator's trust exclusion remains a coverage approximation.

This construction still requires an explicit household mapping. A living
couple's whole SCF primary-economic-unit wealth cannot be assigned to the
first spouse's death and called household exit. A last-survivor model needs
an externally specified representative-household exit schedule, with a stated
mapping of head/spouse ages and surviving-spouse households. Current marital
status is not a lifetime joint-survival history. If that mapping cannot be
justified with the available data, the object must be described as an age-
weighted wealth statistic rather than a realized household-exit flow. The
existing independent-spouse/equal-ownership robustness row does not resolve
this problem automatically. SCF public head ages extend to a top-coded 95;
the current model ends at 82. Pooling older SCF families into the terminal
model cell or restricting empirical coverage are different measurement
choices and must be explicit, including their denominator effects. Keep sale
costs in utility and distributable resources; the wealth observer is
pre-liquidation.

This is a proposed target/observer change requiring the author's decision and
a new target fingerprint. **No replacement scalar or independent estate-flow
match has been established.** Nor has identification of bequest strength
$\theta_0$ been demonstrated. This moment may inform $\theta_0$ through the age
profile of wealth, but can be redundant with other wealth moments. After the
contract is fixed, require non-negligible sensitivity and a full-rank scaled
moment Jacobian for the estimated parameter block. Target count alone is
insufficient. $\theta_1$ remains fixed at 1% of median annual gross working
earnings.

There are two honest alternatives if the author rejects this constructed
moment or its identification check fails:

- Fix $\theta_0$ externally and remove it from the estimated parameter vector.
  A carried-forward fitted value is a maintained numerical restriction, not
  an independently estimated bequest preference. No independently validated
  external value has been identified in this review.
- Replace the row with a directly measured late-life wealth-retention moment
  (for example, mean wealth at older ages relative to middle ages), with
  fixed sample/weights and uncertainty, then test its identifying information.
  Such a target has not been constructed here; it cannot be declared ready.

An independently observed all-estate flow would also be useful, but none is
established by the completed evidence. Do not drop the current row and retain
a free $\theta_0$ without an explicit replacement restriction. The other fit
targets remain unchanged by this sheet.

## What is decided versus what still needs execution

The completed evidence supports uncertain adult receipts over their certain
conditional mean. Comparing certain payments with the receipt lottery on the
finer grid, ownership changes by $-2.094$ versus
$+0.105$ percentage points, the parent-ownership gap by $-4.948$ versus
$-1.615$ points, and fertility by $+0.091$ versus $+0.040$. These are fixed-price,
fixed-scale, fixed-parameter effects. The full fit/parameter tables and all
limitations are linked in the [experiment record](../../output/model/estate_receiver_probe/README.md).

The remaining author decisions are precise: approve or override the complete
baseline above; explicitly authorize replacement of the adopted child-directed
calibration row; and select the household mortality/exit definition with its
existing owner. In particular, approval means accepting age-only anonymous
receipts, zero external estate diversion, externally financed entry positions,
and external absorption of residual credit losses. None is an estimated fact.

Implementation estimates below are engineering estimates, not completed work:

| Work after the decisions | Bounded effort and acceptance condition |
|---|---|
| Constructed SCF wealth moment and model observer | One focused builder/observer pass, roughly 1–3 analyst-hours after the fixed age weights and household mapping are specified. Preserve sample/weights; save the scalar, uncertainty where available, source receipt and target fingerprint. Identification requires a separate local sensitivity/rank check; its cost depends on the accepted solver. No value can be guessed from the old 0.729%. |
| Common-scale joint entry projection | Roughly 2–3 analyst-hours plus one bounded Torch feasibility/solution check. Verify probability mass, signed wealth, joint support, tails and tied ranks. The old 3-by-5 current-resource check under 17.9% tax does not certify this law under the chosen utility/tax contract. Infeasible empirical support requires an explicit credit/measurement decision; no silent censoring or wealth floor. |
| Funded stationary receipt loop | Roughly 1–2 analyst-hours using the tested operator, then a capped Torch fixed-point smoke. Require $\mathrm{paid}=D^+$, market clearing, unchanged income measurement and the existing numerical gates. |
| Dated receipt and external accounts | Roughly 2–4 analyst-hours plus timing/adjoint tests and a constant-path replay. A hand example must show that deaths at $t$ fund receipts at $t+1$ exactly once; heterogeneous paths must use the receiving date's scale. |
| Creditor/entry reporting | Roughly 1–2 analyst-hours integrated with the borrowing owner. Test both signs of estates, net sale costs, entry positions and the external-loss identity. Fix the separately identified incumbent-owner borrowing adapter discrepancy before price transitions. |

These tasks can share code and partly run in parallel. Existing results do not
certify the combined specification. Wealth-grid convergence, consistent
mortality, the selected utility/tax rule, and transaction-specific treatment
of incumbent debt are integration requirements; no numerical tolerance is
relaxed by this recommendation. Today's model decision can be closed by an
explicit accepted contract, while its implementation checks must pass before
next week's calibration begins.

## Existing evidence

- [Fed 2018 note and published age profiles](https://www.federalreserve.gov/econres/notes/feds-notes/how-does-intergenerational-wealth-transmission-affect-wealth-concentration-20180601.html): pooled SCF 1995–2016, three-year receipt probability/conditional amount. The four-year mapping assumes constant within-cell arrival intensity. Published amount units are used only relatively and cancel in the estate-pool scaling.
- [Entry rationale and common-scale construction](../../output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/entry_wealth_rationale_review/report.md): deliberate ages 18–24 childless-renter sample, denominator correction, joint-distribution and support limitations.
- [Adopted child-directed target receipt](../../output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/bequest_flow_2007/calculation/verified_receipt.json): exact estimator and adopted 0.729% definition.
- [Signed estate and entry audit](../../output/model/estate_receiver_probe/results/recipient_evidence/accounts_v1/resource_account.json): negative net estates zero in all three retained cases; signed liquidation-cost and donor-observer identities verified.
- [Borrowing review](borrowing_negative_equity_resolution.md): incumbent-underwater constraint discrepancy and limits of the one-signed-asset representation.

No frozen run, manuscript, slide, empirical target or canonical status was
changed to prepare this decision sheet.
