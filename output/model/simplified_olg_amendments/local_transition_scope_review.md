## Verdict

The proposed stationary result is correct under its stated limiting branch. It follows from the original owner and old-owner problems, not from the transition proof. See the dated budgets in [simplified_olg_amendment_proposal.tex](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex:110) and the exact heterogeneous aggregation identities in [local_transition_proof.md](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/simplified_olg_amendments/local_transition_proof.md:384).

Let \(\ell=1-q\), \(\bar d=E[b_i]/(1-\phi)\), and
\[
K=\frac{\nu\vartheta}{\kappa(\alpha+\vartheta)},\qquad
\rho=1+\beta(1+\gamma+\omega_B).
\]
With \(\chi=\tau^p=0\), fixed \(F\), common \(\phi\), stationary prices, and the uniformly maintained owner branch,
\[
P^*=K\bar d,\qquad h_i^*=\frac{d_i}{K\bar d},\qquad
n_i^*=\frac{\vartheta h_i^*}{\kappa(\alpha+\vartheta)}
=\frac{b_i}{\nu E[b_i]},
\]
\[
x_i^*=\frac{y_i+b_i-\ell d_i}{\rho}.
\]
Thus a change in \(\phi\) leaves each young owner’s \(h_i^*\) and \(n_i^*\) unchanged, while
\[
\frac{\partial x_i^*}{\partial\phi}
=-\frac{\ell d_i}{\rho(1-\phi)},\qquad
\frac{\partial\log P^*}{\partial\phi}=\frac1{1-\phi}.
\]

After substituting the interior old solution into the original conditional owner value, the terms that vary are
\[
W_i^O=\rho\log x_i^*-\beta\gamma\log P^*+\text{terms constant in }\phi.
\]
Therefore
\[
\boxed{\quad
\frac{dW_i^O}{d\phi}
=-\frac{(1-q)d_i/x_i^*+\beta\gamma}{1-\phi}<0.
\quad}
\]
No correction is needed. If realised utility is written \(W_i^O+\xi_i\), the derivative is unchanged when the same household’s \(\xi_i\) is held fixed.

Old owner housing also falls:
\[
h_{i,\mathrm{old}}^*=\frac{\beta\gamma x_i^*}{q(1-q)P^*},
\qquad
\frac{d\log h_{i,\mathrm{old}}^*}{d\phi}
=-\frac{1+(1-q)d_i/(\rho x_i^*)}{1-\phi}<0.
\]
Finally, since mean young housing is fixed at \(1/K\), while the old/young housing ratio \(C\) falls with \(d=\bar d\),
\[
\frac{\partial\log N_{\rm hh}^*}{\partial\log d}
=\frac{Dr}{q(1-q)(1+C)}>0,
\]
under the maintained positive stationary branch. This is the stationary result recorded in [simplified_olg_simple_assessment.tex](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_suggestions/simplified_olg_simple_assessment.tex:510).

Exact conditions: \(q\) is fixed; \(F\) and \(b_i\) are unchanged across stationary equilibria; the purchase constraint binds strictly; saving is positive; owner physical caps, old retention, and estate constraints are slack; old choices are interior; all log arguments remain positive; and the \(\phi\) reform is sufficiently small to retain these branches. It is not a result about transitional cohorts, a variable-population social welfare criterion, or a compensated implementation. Nor does it take a finite limit of divergent ownership-taste utility levels: it concerns the conditional owner problem, with a fixed finite taste added afterward.

## Assumption assessment

Local smoothness assumptions are: uniform strictness of the conditional branches over the support of \(F\), a small permanent reform, and small positive renter mass, \(\chi\), and \(\tau^p\) when extending away from the limiting economy. They let the implicit-function/stable-path argument carry signs locally; they do not quantify an empirically relevant neighborhood. Sections 9–10 are appropriately explicit on this point.

The economically substantive restrictions are stronger:

- Near-all ownership removes tenure reallocation as the source of the transition. It is a clean demand limit, but not a claim about an economy with large renter responses.
- Slack old retention and estate constraints make old housing adjustable at the ordinary user cost. This is central: if old owners cannot release housing or must retain it for estates, the price and welfare formulas change.
- Fixed entrant endowments independent of estates shut down dynastic wealth feedback. That is important for interpreting the population recurrence, not merely for differentiability.
- A world bond price holds \(q\) fixed. The result is therefore a housing-credit/capitalization mechanism, not a closed-economy interest-rate or aggregate-saving result.
- Small \(\chi\) and \(\tau^p\) give a transparent zero-cost limit. With \(\chi>0\), stationary population can have either sign, as the proof itself shows in section 8.

The three restrictions that matter most economically are the binding young purchase constraint paired with slack old release constraints, estate-independent entrant wealth, and the world-bond closure. Near-all ownership is the principal scope limitation.

## Figure 2 and next discussion

Recommendation: retain Figure 2 as a local transition figure, placed after the conditional fertility discussion (main text only if the transition is part of the paper’s core; otherwise appendix). It should display the derived small-reform paths for fertility, young/old or total household population, price, and young versus old housing. Label it as an analytical response per unit rise in \(d\), in the \(\chi=\tau^p=0\) limit. Do not use it to depict welfare or a compensated allocation.

The next useful author discussion is whether the paper wants the stationary welfare derivative only as a disciplined caution: competitive credit relaxation raises price, lowers old housing, raises stationary household population, yet lowers each same-type owner’s conditional stationary lifetime value. That clarification cleanly separates four objects: conditional fertility, the full local transition, stationary population, and welfare. It does not require choosing public lending powers or claiming that the credit reform implements the fixed-fertility compensated reallocation.