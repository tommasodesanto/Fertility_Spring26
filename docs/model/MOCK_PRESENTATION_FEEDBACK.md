# Mock presentation feedback

Running record of the September presentation mock: slide fixes, substantive questions, and decisions.

First-pass review task: **Mock presentation: concise conceptual review** (`01a09818-86af-75a1-9ff4-3a396e394c4d`), max reasoning. Replies are advisory; only individually authorized fixes may be implemented.

## Open issues

### M02 — Adult aging and death

Do people die randomly and age randomly at any age? Clarify the actual timing.

**Status:** open; concise first-pass review requested. No changes authorized.

### M03 — Child versus adult aging

Children mature stochastically; do adults age deterministically?

**Status:** resolved for author-facing wording. The concern is preserved for review: the stochastic event is the exit from dependency, while adult household age advances deterministically. No implementation verification is claimed.

### M04 — Number of children: m versus n

Do we still need the distinction between children at home (m) and children ever born (n), now that bequests do not depend on n?

**Status:** open; concise first-pass review requested. No changes authorized.

### M05 — Utility specification and literature

The utility function is broadly unjustified and hard to parse. Is it common in the literature and quantitative work? Explain the role of every component: why divide by e(m) outside; why retain the Stone–Geary housing floor; why combine Cobb–Douglas and CRRA; why is the preference for children linear and outside the consumption/housing aggregator? Saverio suggests showing only a generic u and its shape properties in the model section, then specifying functional forms in Quantification.

**Status:** open; concise first-pass review requested. No changes authorized.

### M06 — Budget-constraint slide

The slide is messy; the constraints need to be made clear.

**Status:** open; concise first-pass review requested. No changes authorized.

### M07 — Fertility and child-aging slide

The slide is messy and confusing. Some material concerns calibration and should not be in the model exposition.

**Status:** open; concise first-pass review requested. No changes authorized.

### M08 — Age earnings profile

Is e_a exogenously fixed, and is it really needed?

**Status:** open; concise first-pass review requested. No changes authorized.

### M09 — Housing-supply timing

Is housing supply instantaneous? Confirm.

**Status:** open; concise first-pass review requested. No changes authorized.

### M10 — User-cost equation

Clarify the user-cost equation very deeply, for Tommaso as well as the audience.

**Status:** open; concise first-pass review requested. No changes authorized.

### M11 — Gumbel shocks, expectations, and constraints

If a household is constrained, it still receives the taste shock: is that appropriate? Clarify the Gumbel shocks and expectations.

**Status:** open; concise first-pass review requested. No changes authorized.

### M12 — Household-problem slides

The slides are messy, the notation is not transparent, and it is unclear whether m or n is needed. The preceding introductory slide may be unnecessary.

**Status:** open; concise first-pass review requested. No changes authorized.

### M13 — Parental death and children

Very important: how do survival s_a and the number of children interact? If parents die, do their children mature? If so, is market clearing/population accounting consistent?

**Status:** open; concise first-pass review requested. No changes authorized.

### M14 — Equilibrium and population accounting

The within-period equilibrium detail is unclear and probably unnecessary. Give a general equilibrium definition, presumably including the population law of motion. The preceding population slide may be unnecessary. Very important: clarify for Tommaso and the audience the distinction between households, agents, children, and the other population objects.

**Status:** open; concise first-pass review requested. No changes authorized.

### M15 — Speed of convergence

Why does convergence to a new equilibrium take so long? Why not one full generation, or two?

**Status:** open; concise first-pass review requested. No changes authorized.

### M16 — Equilibrium exposition as a whole

The sequence of slides around the equilibrium concept is confusing and needs substantial clarification.

**Status:** open; concise first-pass review requested. No changes authorized.

### M17 — Interest rate and calibration

Is i too low? More generally, the whole calibration needs a much clearer explanation and must be clearer in Tommaso’s own understanding.

**Status:** open; concise first-pass review requested. No changes authorized.

### M18 — Fertility measurement

Are the slides showing TFR or completed fertility?

**Status:** open; concise first-pass review requested. No changes authorized.


Each issue will have a stable number and contain:

- **Feedback / question:** preserve the concern as Tommaso states it.
- **Slide or topic:** identify the affected material.
- **Next step and owner:** Tommaso's hand edit, slides edit, or a specific question for another task.
- **Status:** open, investigating, awaiting an answer, or resolved.
- **Answer / decision:** record the resolution and supporting source beside the original issue.

**Current authorization:** record the issues and request one separate task, at max reasoning, to give concise first-pass answers to all 18 points. That task may flag problems and uncertainties only; no code inspection, implementation, slide edits, or model changes. Further investigation and fixes require individual authorization; M01 has been authorized and completed. Keep each eventual answer beside its original issue, without treating an unanswered concern as resolved.

## Resolved issues

### M01 — Earnings terminology

Replace “heterogeneous income” with “idiosyncratic earnings risk.”

**Status:** resolved; Tommaso authorized this edit with “do it.”

**Resolution:** On “This Paper,” replaced “Income heterogeneity, earnings risk, and bequests” with “Idiosyncratic earnings risk and bequests.” The model section still identifies permanent income groups separately.

## Related material

- [Working presentation](../../output/pdf/september_14_presentation.pdf)
- [Slides collaboration handoff](../prompts/HANDOFF_september14_slides.md)
- [Existing empirical follow-up ledger](POST_PRESENTATION_ISSUES.md) — earlier empirical issues, not yet classified as mock feedback.
