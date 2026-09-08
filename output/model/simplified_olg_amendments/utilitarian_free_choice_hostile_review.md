# Hostile audit of the full-choice transfer construction

September 8, 2026. Independent algebraic derivation, followed by a complete
read of `utilitarian_free_choice_transfer_review.md`. No builds, simulations,
or numerical search.

**Verdict: PASS for the stated special case \(\phi=q\), with identity- and
taste-informed transfers. No mathematical or financing blocker found.**
This is a new full-choice construction: both young groups optimize fertility
and tenure. It preserves aggregate fertility through offsetting responses;
it does not establish a policy-induced increase in population or a general
result for given empirical primitives.

## 1. Household choices, tax ratios, and welfare — PASS

At \(\phi=q\), a binding owner mortgage gives \(a'=-Ph\), hence old
resources \(z=v\). Conditional on binding finance, the young problem is
exactly current utility maximization subject to \(c+ph=w\) and the owner
cap. The uncapped rules \(n_A=a_nw_A\), \(h_A=a_hw_A\), \(x_A=w_A/E\)
are correct.

For capped B, implicit differentiation of its fertility condition gives the
displayed \(b_B=n_{B,w}>0\). Writing \(t=\kappa n_B/H\), the inverse is
\[
w_B=pH+(\chi H/\kappa)f(t),\qquad
f(t)=t+\frac{t(1-t)}{\vartheta-(\alpha+\vartheta)t}.
\]
This equals the note's partial-fraction form. Its derivatives, \(t_*\),
\(t_\dagger\), and \(q_*\) are correct. The inequality
\(t_\dagger>t_*\) is equivalent to \(\alpha r>1\). Thus the advertised
strict interval supplies both \(0<q_B<1\) and \(q_B+d-1>0\).

The exact inverse formula for \(Q(G)\) cancels births, with B remaining
strictly capped for the stated finite interval. Since \(f''>0\),
\(Q'(G)\) decreases toward \(q_*\), remaining positive and below one;
\(S'(G)=Q'(G)+d-1>0\). B's adult consumption decreases when taxed, so
its finance restriction becomes more strictly binding. A's separate cap
and finance bounds preserve its branch. The complete log problem has affine
constraints and a strictly concave objective, making these conditional
solutions global optima.

The welfare derivative is exactly
\[
E/w_A-q_B/x_B-Kd/v_A+m_{OB}(q_B+d-1).
\]
The first strict inequality in welfare terms and positive old-B grant give
a finite gain by continuity. This construction does not require
\(\beta\ge q\); its explicit sufficient wealth bound is doing that work.
Fertility utility losses of young B are included through its optimized
value derivative \(1/x_B\), not omitted.

## 2. Original accounting and the finite tail — PASS

For A and B, all actual young choices satisfy
\(qa'+qPh=0\). A's extra mortgage repayment is exactly offset by its
larger inherited title; B's housing and mortgage remain fixed. Their old
resources, old portfolios, and estates therefore remain at baseline optima.
Old A's housing reduction \(gD=a_hG\) clears the additional young demand;
old B stays capped under a positive grant. Government balances immediately:
\(G+S=Q+D\).

The external resource check is exact. Across equal selected masses,
\[
\Delta C_0=(1-pa_h)G-Q-D/K+S/(1+\omega_B),
\]
\[
q\Delta E_1=-\omega_BD/K+\omega_BS/(1+\omega_B).
\]
Their sum is \(-pa_hG+\gamma D/K=0\). Initial old estates, including
old B's increased estate, are therefore fully financed. The change in old
death sales is offset by A's extra inherited title sales next date. There is
no new government debt or terminal refinancing.

Birth changes cancel exactly, preserving next-date cohort masses. Given the
model's exogenous entrant type distribution and estate-independent entrant
wealth, date-one entrants choose their original bundles. The entire
payoff-relevant inherited state is baseline from date two. This conclusion
would fail in a model where parental fertility or estates determined entrant
types or wealth.

## 3. Positive-tax stationary parameterization — PASS

The primitive construction is self-consistent. The stated bound on \(v_B\)
simultaneously makes its old cap and young finance strict. Every upper bound
for \(w_A\) is positive, so choosing it sufficiently small gives A's
strict branch and the welfare inequality. Both types have finite tenure
values, hence positive masses of strict owners under the original logistic
taste distribution.

Holding \(p,w_i,v_i\) fixed and setting
\(P=p/(1-q+q\tau^p)\) preserves all real choice menus. The stronger
condition \(\omega_B(1-q)>q\gamma\) ensures owner estate floors remain
slack as \(\tau^p\) increases from zero; the reduced owner problem and
all rental service costs therefore remain unchanged. The implied rebate
satisfies
\[
T\le q\tau^p pH/(1-q),
\]
so the note's explicit small-tax bound makes every \(w_i-T\) and
\(v_i-T\) positive. Splitting \(w_i-T\) into positive income and liquid
wealth respects the original primitives. Defining \(N\) from mean housing
closes housing and fiscal accounting exactly.

Common scaling of \((\chi,\kappa)\) rescales every fertility choice
inversely, preserving goods and housing allocations in both tenures. Both
tenure values shift by the same \(-\vartheta\log s\), so ownership
probabilities, housing, and the rebate remain unchanged. The proposed scale
therefore sets replacement fertility without disturbing the construction.
This supplies a constructed family of primitive economies, not a theorem
covering every stationary equilibrium of an already-fixed distribution.

## 4. Information and counterfactual timing — required clarifications

**Timing clarification resolved in the assembled appendix.** The policy is unexpected
at date zero, after current young ownership tastes are realized but before
their fertility, tenure, housing, and financial choices. The inherited old
state is the original stationary state. An anticipated intervention at an
earlier date would require redoing that inherited-state argument.

**RETAIN the strong targeting restriction.** Transfers are fixed cash
amounts for identified households, paid or collected regardless of actual
tenure. Selecting young households with a uniform positive baseline
ownership-value margin ensures that they still voluntarily own for small
transfers. This is not established for anonymous endowment-only grants.

Under taste-informed eligibility, use the individual tenure inequality and
integrate over the original logistic density. Conditioning on eligibility
does not itself make the truncated taste distribution logistic. In this
particular policy all selected tenure choices remain unchanged, while
unselected households face unchanged menus, so aggregate ownership also
stays exactly at baseline. Do not apply the original endowment-only logit
formula to a fictitious representative transfer.

Future cohort utility is unchanged in aggregate because masses and type
distributions are unchanged. Parentage can differ; no claim that every
future person's identity is unchanged is needed or established.

## 5. Assembled appendix follow-up

Read the complete new appendix “Voluntary fertility with unchanged cohort
size” in `simplified_olg_utilitarian.tex`. Its formulas, sufficient
conditions, explicit cap interval, welfare calculation, and stationary
constructor pass this audit. Its announcement timing is now explicit and
correct. It preserves the special-case and information restrictions.

One wording correction was sent to the lead: replace “The usual endowment-only
logit formula must condition on transfer eligibility” by an instruction to
compute tenure probabilities from individual choice inequalities and the
original taste distribution. Eligibility conditioning does not restore the
ordinary logit formula. This concerns the probability-description sentence;
the preservation of actual tenure choices used by the proof is valid.
