# ChatGPT Pro Prompt: V4 Theory And Figure Review

You are advising on an economics theory/quantitative macro draft. Please take a
deep, skeptical pass. Do not summarize casually. I want an external referee/editor
style assessment of the latest iteration, with special attention to whether the
two new figures actually clarify the economics.

Attached/provided context:

- `latex/intergenerational_housing_fertility_v4.tex`: current full draft.
- `latex/intergenerational_housing_fertility_v4.pdf`: compiled current draft,
  if visible.
- `latex/figures/fig6_ce_planner_wedge.tex`: source for the first candidate
  figure.
- `latex/figures/fig6_ce_planner_wedge.png` / `.pdf`: rendered first candidate
  figure, if visible.
- `latex/figures/fig7_entry_fertility_decomposition.tex`: source for the second
  candidate figure.
- `latex/figures/fig7_entry_fertility_decomposition.png` / `.pdf`: rendered
  second candidate figure, if visible.

Project context:

The project studies a toy non-spatial overlapping-generations model linking
housing costs, old-household retention, young household housing constraints,
fertility, bequests, and entry through an outside option. The current v4 draft
tries to clean up the theory after several iterations. It should read like a
precise economics note: primitives, household problems, competitive equilibrium,
planner problem, and one or two clear efficiency results. Exposition should be
clear, disciplined, and close in spirit to JPE-style macro theory and Guido
Menzio's style: minimal naming, no tangents, clean timing, precise propositions,
transparent proof logic.

Important constraints:

- Do not restart from scratch.
- Do not propose a broad new model unless absolutely necessary.
- Do not invent weird terminology.
- Do not smooth over weak math. If a proposition is not right, say exactly why.
- Warm-glow bequests should remain in the model as true preferences.
- The outside option / endogenous entry should remain in the model, but it
  should be integrated cleanly rather than appended.
- The old-age retention wedge `L_j^p(q,x)` is not a true preference term by
  itself. The planner should not count it as utility or a resource saving unless
  its real/fiscal counterpart is modeled.
- The young down-payment or housing-finance friction creates a private shadow
  wedge `zeta_i`; be very careful about whether a constrained planner can remove
  it, or whether the result is only a local marginal reallocation with
  compensating transfers/instruments.

Tasks:

## 1. Figure review and redesign

Inspect the two candidate figures carefully. For each figure, answer:

1. What is the intended economic message?
2. Does the current visual successfully communicate that message?
3. Is anything conceptually wrong or potentially misleading?
4. Is the notation consistent with the model?
5. Is the figure too busy, too informal, or too graphically imprecise for a
   serious paper?
6. What exact changes would you make?

Then propose a better version of each figure. Be concrete. I want:

- a recommended title;
- a recommended caption;
- what should be on each axis;
- what curves/areas/arrows should appear;
- what labels should be removed or renamed;
- whether the figure should be one panel or two panels;
- whether it belongs in the main text or appendix;
- the exact conceptual claim the figure is allowed to support.

If useful, write replacement TikZ design instructions that a local coding agent
could implement.

## 2. Review the current theory iteration

Evaluate whether the compact theory block in v4 passes a serious review
standard. Focus on the current model, not the full quantitative codebase.

Check:

1. Are the primitives and timing stated cleanly enough?
2. Are young and old household problems well defined?
3. Is the competitive equilibrium definition correct and complete?
4. Is the planner problem the right object for the efficiency result?
5. Are the FOCs and wedge identities mathematically valid?
6. Is the no-wedge benchmark stated too strongly, too weakly, or correctly?
7. Is the young-old housing wedge proposition correct as written?
8. Does the outside-option entry margin belong exactly where it currently is?
9. Are warm-glow bequests integrated correctly?
10. Would a serious macro/urban/referee find the result clean, or would they
    object that the planner/instrument assumptions are unclear?

Please be especially hard on the constrained-efficiency logic. Distinguish:

- true preference MRS;
- private user-cost wedges;
- borrowing/down-payment constraints;
- fiscal or lender counterparts;
- local compensated reallocations;
- implementable policy changes;
- aggregate entry effects.

## 3. Referee-style bottom line

Give a concise but serious bottom line:

- Does this last iteration pass as a clean analytical note?
- What would the top two referee objections be?
- What is the smallest revision that would make the theory credible?
- Which figure, if any, should be kept?
- What should be cut or moved to the appendix?
- What is the one strongest result the note should emphasize?

## 4. Suggested edits

Provide a ranked list of actionable edits:

- P0: must fix before circulating.
- P1: important for clarity/polish.
- P2: optional.

Where possible, propose exact replacement text for:

- the planner paragraph;
- the no-wedge proposition;
- the young-old wedge proposition;
- the captions for both figures.

Do not provide a generic rewrite of the whole paper. Give precise diagnosis and
targeted repair instructions.
