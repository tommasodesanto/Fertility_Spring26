# A finite demographic transition with renters and owners

You are an independent mathematical economic theorist. Please help settle one
specific proof question in a simple two-period overlapping-generations model
of housing and fertility. Work on the mathematics; do not rewrite the paper,
invent a different household problem, or give a literature survey. A short,
transparent result is more useful than an elaborate theorem with no clear
economic content. Treat all supplied claims as claims to check.

## The question

Can the finite transition argument in `transition_extensions.md`, section 4,
be extended to an economy with a substantial share of renters and strictly
positive child goods costs? Ideally give simple sufficient conditions that
cover a class of economies. A proved, explicit finite neighborhood around
the genuine mixed-tenure example in `mixed_transition_proof.md` would also
be useful. An existence statement for an unspecified sufficiently small
shock around that example is already available; repeating it would not
resolve this question.

The economic experiment has two steps. Start at a positive stationary
equilibrium with fertility weight \(\vartheta_0\) and financed share
\(\phi_0\). At date zero an unexpected permanent decline to
\(\vartheta_1<\vartheta_0\) starts a demographic transition. At a later date
\(t_p\geq1\), compare continuation of that baseline with an unexpected
permanent credit relaxation \(\phi_1>\phi_0\). Both continuations start
with the same actual inherited cohorts, financial claims, housing titles,
and tenure distribution at the intervention date. Household expectations
are correct following each surprise; previously chosen financial claims
remain those chosen under the earlier forecast.

The desired conclusions are:

1. Both infinite equilibrium paths exist and converge to positive stationary
   endpoints, with the original household constraints satisfied throughout.
2. The original preference decline lowers initial fertility and the baseline
   terminal population relative to the original stationary equilibrium.
3. Credit relaxation during that baseline raises fertility at the intervention
   date and leads to a larger terminal population than continuation without it.

Try first for a finite parameter range and a result valid at any later baseline
date. If that is too much, state precisely which of these conclusions you can
prove and why the rest fails or remains open. The priority is a useful theorem,
not defending every desired sign. If the proposed generality is false, give a
counterexample or an explicit additional restriction. Separate failures of a
proof method from failures of the economic claim.

## The model must stay fixed

The attached TeX gives the full model, budgets, value functions and equilibrium.
The following conventions matter:

- Young utility is
  \(\log(c-\chi n)+\alpha\log(h-\kappa n)+\vartheta_t\log n\).
  Old utility is \(\log c^2+\gamma\log h^2+\omega_B\log e\), discounted
  by \(\beta\) when young. Preserve the author's \(W^m,V^m\) functions.
- Completed fertility, housing, goods, saving and tenure are joint household
  choices. A single additive ownership-taste draw is observed before those
  choices. Both parameters of its logistic distribution stay fixed within
  each economy across both shocks. Do not reset them to keep a tenure share
  constant after a reform.
- The entrant distribution of income and liquid wealth is fixed. Homogeneous
  entrants are acceptable for this proof. Estates are warm-glow expenditure
  and are not automatically inherited as entrants' liquid wealth.
- Only pre-existing liquid wealth finances the down payment:
  \((1-\phi_t)P_t h\leq b\). Current income and rebates arrive after purchase.
  Rental and owner size limits are physical constraints. Old owners may
  downsize but cannot buy additional housing or rent out retained housing.
- The consumption good is the numeraire; goods and bonds trade with the rest
  of the world at a fixed bond price \(q\). Do not introduce domestic goods
  or bond clearing. Housing clears against a fixed \(\bar H\).
- The date-\(t\) housing service cost is
  \(u_t=q r_t=(1+q\tau^p)P_t-qP_{t+1}>0\). Property-tax revenue is rebated
  equally to young and old, with
  \(T_t=q\tau^pP_t\bar H/(Y_t+O_t)\). Zero tax is acceptable as a clearly
  stated intermediate theorem; positive child goods costs and material
  renting are the main extension being sought.
- Inherited old financial assets include repayment of the mortgage originally
  contracted. Do not recalculate that mortgage at the new financed share.
  Existing titles revalue at the new price. Consequently the initial old-owner
  housing coefficient \(M_0\) in the six-variable system is not fixed when
  the surprise changes \(P_0\). Retain the actual initial-old boundary.
- Demography is \(Y_{t+1}=\nu\bar n_tY_t\), \(O_{t+1}=Y_t\). Every positive
  stationary endpoint has \(\bar n^*=1/\nu\). The long-run comparison is
  population, not permanently above-replacement fertility. With a common
  inherited young cohort, the limiting young-cohort ratio equals the limit
  of the products of the relative fertility rates.

It is acceptable to retain the strict household branches used in the proofs:
young owners have restrictive down payments, positive saving and a slack
physical cap; young and old renters are at their rental caps; old owners have
slack retention and estate bounds. Do not assume these inequalities without
checking that your proposed path region implies them. Proving transitions
across changes in binding constraints is not required for a useful answer.

## What has already been obtained

Read the attachments in this order:

1. `simplified_olg_amendment_proposal.tex`: the original household and
   equilibrium specification and the intended two-stage comparison. The
   allocation proof is background; planner design is outside this question.
2. `mixed_transition_proof.md`: exact six-variable equilibrium map, initial-old
   boundary, a one-sided sequence argument, and a genuine mixed economy with
   owner share \(11/21\), \(\chi=3/20\), and \(\tau^p=467/9250\). The
   forward derivative has four stable and two unstable roots, including one
   zero stable root. Exact algebra and interval calculations certify the local
   credit-impact and terminal-population signs for a family of taste scales.
3. `transition_extensions.md`: the newest additions. Section 2 proves wider
   mixed-tenure stationary signs when \(\chi=\tau^p=0\). Section 3 supplies
   a simpler sufficient condition for convergence in the all-owner limit.
   Section 4 gives an explicit finite preference decline, later credit reform,
   and an infinite-sequence contraction with inherited claims handled exactly.
4. `local_transition_proof.md`: earlier all-owner derivations, the exact
   parameter example used by section 4 above, and limitations/counterexamples.
   Its older statements about unspecified neighborhoods should be read together
   with the newer finite construction, not as negating that construction.
5. `verify_simplified_olg_transition_extensions.py`: the latest symbolic and
   rational-interval verification calculations. It imports earlier repository
   helpers which are not included. You are not being given a self-contained
   runnable repository; this file supplies inspectable formulas and exact
   bounds, not independent authority for any claim.

Two limitations are already established. Positive child goods costs can make
stationary credit/population signs fail outside suitable restrictions. Also,
mixed equilibrium paths can oscillate around their endpoints. We therefore
do not ask for fertility to be higher at every future date or for monotone
convergence. A welfare improvement does not imply higher fertility, and the
competitive credit reform is not asserted to be a Pareto improvement.

## What would count as progress

Please begin by identifying any substantive error in the supplied equations
or arguments that changes the target. Otherwise spend most of the answer
trying to prove the extension. A better choice of variables or an
infinite-sequence boundary-value operator may help. The six-dimensional
forward map has unstable roots, so simply assuming its whole derivative is
a contraction cannot justify the desired equilibrium transition.

Give a precise proposition, define its assumptions, and provide the proof
with the estimates on which it depends. If using validated numerical bounds,
state the finite box, the uniform residual and derivative bounds, how the
infinite tail is controlled, and how the original household inequalities and
initial-old boundary are verified. A finite-horizon simulation with an imposed
terminal steady state, pointwise eigenvalues, or a plotted path is not an
infinite-horizon existence or convergence proof. Distinguish analytical
conditions, proved interval bounds, and numerical evidence.

Look for conditions that can ultimately be explained in a few lines to an
economist. Do not hide the result in an assumption that essentially restates
the desired signs, and do not insert an unexplained borrowing multiplier into
a purported primitive condition. A modest transparent theorem is preferable
to an unsupported global claim. Aim for a compact mathematical answer, with
one clear recommendation about what belongs in the illustrative model and
what is better left as an extension.
