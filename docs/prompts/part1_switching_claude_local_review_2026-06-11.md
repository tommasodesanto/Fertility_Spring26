You are reviewing a local LaTeX theory draft in:

`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

First follow the repository startup protocol:

1. Read `memory/AGENT_MEMORY.md`
2. Read the latest `memory/daily/YYYY-MM-DD.md`
3. Read `CALIBRATION_STATUS.md`
4. Read `SESSION_DIARY.md` only if historical detail is needed

Do not run broad directory scans. Do not edit files. This is review only.

Primary files to inspect:

- `latex/intergenerational_housing_fertility_part1.tex`
- `latex/intergenerational_housing_fertility_part1.pdf`

Optional context if needed:

- `docs/prompts/part1_switching_chatgpt_pro_review_2026-06-11.md`
- `latex/intergenerational_housing_fertility_part1_slides.tex`
- `latex/intergenerational_housing_fertility_part1_slides.pdf`

Task:

Act as a senior academic economist reviewing Part I of a quantitative macro / urban / household-finance paper. Treat the LaTeX source as the maintained draft. Do not rewrite the paper. Audit it.

Project context:

- The paper studies how housing costs and family-sized housing access affect fertility, tenure choice, and intergenerational housing allocation.
- The current file is a standalone Part I compact analytical model.
- The intended style is a clear macro theory paper: economy first, then households, choices, constraints, equilibrium, planner, marginal comparison, fertility, and policy implications.
- Avoid wedge-first exposition, decorative terminology, and unexplained reduced-form objects.
- Keep the economics fixed unless you identify a genuine error or a necessary qualification.

Current model architecture:

- Young households decide whether to enter the modeled housing market.
- Conditional on entry, they choose tenure: rent or own.
- Rental and owner segments have different family-size menus, with owner housing allowing larger family-capable units.
- Owners face down-payment and payment-to-income implementation constraints.
- Children require both nonhousing consumption and housing services via a Stone-Geary child-input structure.
- Old households are incumbent owners. The old-side mechanism, when active, is an incumbent property-tax discount: an old owner pays less tax if she stays in her current house than if that same house were taxed at current market value.
- Competitive equilibrium clears rental and owner segments in aggregate, with segment construction.
- The planner comparison distinguishes real menu constraints from private financial implementation constraints.
- Section 9 contains the fertility implications and a local policy decomposition.

Recent change to audit especially carefully:

- A new tail subsection, "Tenure switching at the boundary," has been added after the local policy decomposition.
- It studies deterministic tenure switching at the boundary \(W_i^O=W_i^R\).
- It claims that when \(q^O>q^R\), a household indifferent between renting and owning has a binding rental cap, chooses more occupied housing and more adult housing slack in the owner branch, and, under \(\chi\le \kappa q^R/\alpha\), has higher fertility in the owner branch: \(n_O>n_R\).
- It then adds a boundary term for the switching component of a fixed-price policy effect and claims owner-access policies make the fixed-tenure formula a lower bound under the same condition, while rental-menu policies can have an offsetting switching term.

Please produce a rigorous review with this structure:

1. Executive verdict: Is the updated Part I coherent as an economics model section? Is the new switching subsection worth keeping, moving to an appendix, or rewriting?
2. Mathematical audit of the new switching subsection:
   - Verify the \(X=\chi/c\), \(Y=\alpha\kappa/s\) transformation.
   - Verify the indifference-curve argument and the claim that fertility is single-peaked along a \(G(X,Y)\) level set.
   - Verify the proof that \(h_O>h_R\), \(s_O>s_R\), \(c_O<c_R\), and \(n_O>n_R\) under the displayed assumptions.
   - Check the boundary-integral switching term and the sign claims for owner-access and rental-menu policies.
   - Identify any hidden regularity conditions, missing cases, or sign mistakes.
3. Economic audit:
   - Is \(q^O>q^R\) the right assumption for the intended "ownership is costlier but offers larger family-space access" margin?
   - Does the result align with the core mechanism, or does it add too much functional-form-specific machinery?
   - Does this treatment of switching interact correctly with entry, the local policy decomposition, and the planner comparison?
4. Exposition/style audit:
   - Is the tail placement sensible?
   - Which pieces should be main text, appendix, or omitted?
   - Suggest minimal wording changes to make the economics more transparent and less notationally heavy.
5. Highest-priority patch list:
   - Give concrete minimal edits, preferably in LaTeX snippets, only where they materially improve correctness or clarity.

Be direct. Separate genuine errors from optional refinements. Do not praise the draft generically; focus on what is correct, what is risky, and what should change before this becomes paper-facing.
