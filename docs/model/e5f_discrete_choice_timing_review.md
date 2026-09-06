---
title: "When do households observe the choice shocks?"
subtitle: "First quantitative discussion before September 14"
date: "6 September 2026"
---

# The decision to settle first

**When deciding whether to attempt a birth, does the household already know its current housing-specific tastes or opportunities?** The current calculation averages over housing shocks before fertility choice. This information assumption is now the first quantitative priority and remains open for discussion.

The maintained code is consistent with the following sequence within each four-year period:

1. The household knows wealth, its inherited house, income state, children and prices. It observes fertility-choice shocks and chooses **wait or attempt a birth**.
2. Conception succeeds or fails. The household learns its resulting family state.
3. It observes the housing-choice shocks and chooses **rent or an owner product**. Renter rooms, consumption and saving are optimized conditional on the realized family state and housing choice.
4. Current choice shocks are averaged out of continuation values; they are not persistent state variables.

The Bellman solver works backward through these stages. Computational order does not imply chronological housing choice before fertility.

# Three different specifications

| Specification | What changes economically? |
|---|---|
| **Later housing information: current calculation** | Fertility responds to expected housing opportunities; the current housing-shock realization cannot select households into birth attempts. |
| **Earlier housing information** | Households know housing-specific shocks when deciding on a birth. They can still choose their actual house after conception is resolved. |
| **Housing committed before fertility** | The household chooses a house first, restricting how it can adjust space after a birth. This changes the feasible decisions as well as their information. |

**Knowing housing shocks early is different from buying a house early.** We should decide the information assumption before choosing whether housing must be committed in advance.

**Illustration only, not calibrated:** let housing information make the attempt gain either 2 or -1 with equal probability, and set fertility dispersion to 1. Later information gives the logistic of the average gain: **62.25%**. Earlier information gives the average conditional probability: **57.49%**. Moving the expectation changes behavior.

\clearpage

# The equation behind the timing

Let $Q_h(x)$ be the value of housing alternative $h$, after optimizing consumption and saving, for the realized family state $x$. The expected value before housing shocks is

$$
\mathcal H(x)=\mathbb E_{\varepsilon^H}\max_h\{Q_h(x)+\varepsilon_h^H\}
=\kappa_H\log\sum_h\exp\{Q_h(x)/\kappa_H\}.
$$

The equality uses independent, mean-zero type-I extreme-value housing shocks, the normalization consistent with the coded log-sum. The current rent/owner-product kernel returns both this value and its choice probabilities. Setting $\kappa_H=0$ takes the deterministic maximum; it does not move housing before fertility.

Let $x^+$ add a successful birth, $\pi_a$ be conception success, and $F_1$ the first-birth cost. The gain from attempting is

$$
\Delta(x)=\pi_a\left[\mathcal H(x^+)-\mathcal H(x)-F_1\mathbf 1\{n=0\}\right],
\qquad p_A(x)=\frac{1}{1+\exp\{-\Delta(x)/\kappa_F(n)\}}.
$$

Here $\kappa_F(n)$ is the first-birth or continuation-birth dispersion. The taste shocks attach to **wait versus attempt**; conception uncertainty is averaged inside the attempt value. The cost is paid only on first-birth success, and the birth probability is $\pi_a p_A(x)$. With housing information known earlier, the attempt probability would instead average conditional logits over that information. In general, an average of logits differs from a logit of an average.

# What is established, and what is still a choice

- **Code:** the Bellman nesting and forward population calculation agree on fertility/conception before current housing realization. The housing shocks cover rent plus individual owner sizes, not just a binary rent/own label. Persistent housing tastes or correlated fertility/housing shocks would require a different specification. The generic location-choice layer lies between fertility and tenure; the retained calibration has one housing market.
- **History:** the existing paper-issues ledger, **HT-1**, explicitly left fertility-before-housing, housing-before-fertility and joint choice open. Its surrounding description predates sequential births. Adding sequential births across ages did not by itself settle information within a period. A bounded historical search found no explicit author resolution of this information question.
- **Presentation:** the September 14 Bellman slide has housing expectations inside fertility choice. Its assigned-parameters table still lists zero tenure dispersion and older estimates; the retained calibration explicitly overrides tenure dispersion to **0.005**. These belong to different recorded specifications and need reconciliation before presentation.

**Suggested discussion order:** first state what the housing shocks represent and when households plausibly learn them; then decide whether a house must be committed before fertility; finally choose shock dependence/persistence and scales consistently. We have not quantified the effect of changing this timing, and cannot attribute the housing-response miss to it yet. Any adopted change requires consistent equations, code, measurement and re-estimation.

The checked solver and housing-kernel files exactly match the retained calibration's source hashes. Relevant evidence: [Bellman nesting](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/intergen_eqscale_seq_optimized/solver.py:2622), [forward sequence](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/model/tools/run_dynamic_population_transition.py:443), and [open timing issue HT-1](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/intergenerational_housing_fertility_paper_issues.md:90). No model, target, slide or protected manuscript was changed in this review.
