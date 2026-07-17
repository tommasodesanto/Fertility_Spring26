# Economics Writing Style Guide

Tool-agnostic version of the `writing-econ-papers` Claude skill. Paste or attach this
file in any LLM prompt (ChatGPT, Codex, Claude) that drafts or revises paper-facing
text: model sections, abstracts, theory notes, slide decks.


Audience for this guide: an LLM (or human) drafting or revising paper-facing text in
quantitative macro, housing, spatial, urban, or household finance. Distilled from real
revision rounds on theory drafts and decks; every rule below exists because its violation
was produced by a strong model and rejected by the author.

## 1. Section architecture (the Menzio rule)

Model sections follow this order, and the boundary between "environment" and "problems"
is strict:

1. **Environment** — labeled primitive paragraphs, in bold lead style:
   - **Agents.** Who exists, masses, types and their distributions.
   - **Preferences.** Utility functions for every agent, stated as primitives.
     Interpret parameters in words. If a structural identity makes the economics vivid
     (e.g. substituting residuals into a budget to show "the full price of a child"),
     show it here.
   - **Endowments / Technology / Markets.** What is traded, at what prices, what clears
     against what. Choose the simplest market structure that carries the point — one
     clearing condition if possible. Richer structure (multiple stocks, construction,
     segment-specific supply) is a footnoted extension, not the baseline.
   - **Institutions.** Financing limits, tax rules, assessment rules — as rules of the
     game, not as constraints inside someone's problem yet.
   - **Entry / demography** if relevant.
   - NO maximization problems, Bellman equations, FOCs, or value functions here.
2. **One section per agent's problem.** "The young household's problem" contains the
   problem, taking prices as given, and nothing that belongs in the environment.
3. **Equilibrium** — a numbered Definition listing the objects first, then the
   conditions. Lead the section with one paragraph saying what the equilibrium object is
   and what it does NOT pin down (e.g. "clearing pins aggregates, not assignment").
4. **Planner / efficiency** — open with the question the planner answers.
5. **Results, then policy.** Comparative statics before policy formulas; honest scoping
   sentences on what is local/fixed-price/conditional.

If the paper has a toy model and a full model (Coven et al. structure), the toy comes
first under its own heading, then "Setup", then the framework.

## 2. Primitives before derived objects

A derived object (a multiplier, gap, wedge, sufficient statistic) may not do work in the
text before it has been expressed in primitives and tied to observables.

**Bad (rejected in a real round):**

> Define the goods-equivalent value of relaxing the rental size limit by
> $\zeta_i^R=\theta_i^R/\lambda_i^R\ge0$.

**Good (accepted):**

> Converting the cap multiplier into consumption goods, $\zeta_i^R=\theta_i^R/\lambda_i^R$,
> the renter's marginal value of housing is $\alpha c_i^R/s_i^R = q+\zeta_i^R$.
> A renter against the size cap values one more unit of space above the rent it would
> command, and the gap is not an abstract multiplier: it is computable from the observed
> bundle as $\zeta_i^R=\alpha c_i^R/s_i^R-q$. An unconstrained renter has $\zeta_i^R=0$.

Corollaries of the rule:
- Convert every multiplier to goods units before interpreting it.
- Decompose composite gaps by the primitive constraint that generates each term
  ("the collateral margin... the debt-service margin... the menu").
- Name objects by their economics: "incumbent tax discount", "subsidy to staying put",
  "effective family-housing cap" — never decorative or purely technical labels.

## 3. Words around math

- Every display gets a one-line lead-in ending in a colon or verb, and a one-sentence
  interpretation immediately after. No orphan displays.
- Mechanisms are stated in numbered plain prose before any notation:
  "The mechanism has two parts. First, ... Second, ..."
- Facts first. Open sections and abstracts with the economic fact or question, not with
  notation or with the literature.
- Honest scoping is part of the result: say "holding prices, tenure, and active sets
  fixed", say what reactivates the omitted margin, say where the general case is handled.
- Do not overclaim: a conditional marginal comparison is not a welfare theorem; a local
  formula is not a GE counterfactual.

## 4. Remarks, footnotes, and apparatus

- Remark environments are for standalone formal caveats only. If it reads like
  commentary, fold it into prose; if it is bookkeeping, make it a footnote.
- Status notes, normalization conventions, and "we will extend this" content belong in
  footnotes, written as economics ("Two supply extensions are deliberately left out of
  the baseline and developed later: ...").
- Keep equation labels stable across revisions. Avoid notation collisions (x vs X, y vs
  Y); if a local construction needs symbols that collide, scope them explicitly in a
  footnote.
- Assumptions stated where used; iff results get both directions discussed; proofs
  include the nontrivial computation (no "it is easy to see").

## 5. Citation discipline

- Never describe another paper's model from memory. Open the PDF and read the model
  section, the equilibrium definition, and the institutional details before writing
  "as in X (year)".
- "As in X" must name exactly which feature is borrowed ("all floorspace clears in a
  single market, as in X"), and a footnote states where you differ ("X impose a minimum
  owned size; we cap rental size instead").
- Cite the version you actually read; update years and working-paper numbers.
- If a mechanism you use is one that a cited paper explicitly abstracts from, say so —
  that is positioning gold, not a problem.

## 6. Slides

- Conventional noun-phrase frame titles ("Young households", "Competitive equilibrium",
  "Marginal values and access gaps"). Not assertion sentences, not section numbers with
  no content.
- One idea per frame: at most one framing line, then the display(s), then at most one
  interpretation line. No symbol glossaries; define symbols in the line that uses them.
- Figures only where they carry a result (target ~2 per deck: the mechanism figure and
  the result figure). Diagrams of the economy, pipelines, and bar-chart metaphors are
  usually too much.
- Proofs and derivations go to appendix frames.
- Metaphors only if load-bearing and used once ("subsidy to staying put" earns its place;
  "misallocation thermometer" does not).
- Calibrate to "casual but technical reader": more normal than you think. When in doubt,
  copy the structure of a good job-market-paper deck, not a popular talk.

## 7. Workflow checks (LaTeX)

- Compile twice; require zero errors, zero undefined references, zero overfull boxes
  above ~2pt. Render pages to images and look at them before declaring victory —
  overlapping TikZ labels and clipped boxes do not show up in logs.
- After structural edits, grep for orphaned `\ref`/`\eqref` of deleted labels and for
  leftover symbols from the old structure.
- Commit a checkpoint before and after any restructure.

## 8. Red flags → corrections

| Thought | Reality |
|---|---|
| "I remember what that paper does" | Read the model section first; quote the clearing condition. |
| "This object is standard, no need to unpack" | Express it in primitives and observables first. |
| "A remark environment is safer" | Fold into prose or a footnote. |
| "The title can carry the claim" | Conventional title; the claim goes in the text. |
| "More structure shows rigor" | Simplest market structure that carries the point; extensions in footnotes. |
| "The math speaks for itself" | One sentence of economics before and after every display. |
| "Creative presentation helps the reader" | More normal than you think. One idea per frame. |
| "The preferences belong with the problem" | Preferences are primitives; they go in the environment. |
