# Handoff: Fable audit of the simplified OLG theory note

Act as an independent senior macroeconomic theorist and a demanding journal
referee. Audit the simplified theory note described below. Do not assume that
the current draft is correct, but prefer the smallest defensible repair over a
larger model.

## Objective

The note is a deliberately minimal, self-contained theoretical exercise. It is
not the quantitative model and may not enter the final paper. It should be no
longer than five pages, should contain no introduction or workflow narration,
and should accomplish two things:

1. Show a precise sense in which housing can be misallocated across living
   generations when young buyers face a financing constraint and old
   incumbents do not.
2. Define equilibrium as a perfect-foresight transition between two positive
   steady states after one permanent change in fertility preferences.

The welfare result should use the weakest meaningful concept. In particular,
it should not require choosing welfare weights over unborn households. The
transition should not introduce a second shock merely to generate capital
gains. Capital gains must be zero at a constant real house price and may arise
only on appreciation dates along the transition.

## Files to read

Work from the repository root:

`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

Follow the mandatory startup in `CLAUDE.md`, but do not turn this into a review
of the quantitative model. Then read:

1. `docs/style/econ_writing_style_guide.md`
2. `latex/simplified_olg_transition_note.tex`
3. `latex/simplified_olg_transition_note.pdf` for the rendered presentation

Do not use older internal drafts as a source of authority. The object being
audited is the current note.

## Intended economic closure

Preserve these choices unless they are internally inconsistent:

- Agents live for two periods, but calendar time is infinite.
- Young households choose fertility, housing, saving, and either renting or
  owning. There is no hybrid tenure contract.
- A young owner uses a one-period mortgage. The financed share is `phi`, so the
  down payment is `(1-phi) P_t h`.
- An old owner either sells or retains the house. Housing retained until death
  receives stepped-up basis in the estate.
- The housing stock is fixed and clears every period.
- Property-tax and realized capital-gains-tax revenue is rebated through a
  balanced government budget every period.
- The draft conditions on two parameterizations that admit positive steady
  states. It need not establish global existence or uniqueness unless those
  claims are logically required for the transition definition.
- One permanent change from `vartheta^0` to `vartheta^1` moves the economy from
  the initial steady state to the terminal steady state.
- The misallocation exercise compares a financing-constrained young owner with
  an old incumbent at an interior retention choice.

## Audit tasks

### 1. Timing and accounting

Verify equations (3)--(10) line by line. In particular, check:

- the user-cost equation and discounting by `q`;
- the mortgage advance, down payment, repayment, resale proceeds, and old-age
  balance sheet;
- the timing and signs of `ell_t` and `q ell_{t+1}`;
- the stepped-up-basis treatment of retained housing;
- whether the renter and owner problems generate the stated housing first-order
  conditions.

Distinguish an actual accounting error from a normalization or an omitted
extension.

### 2. Steady states and transition equilibrium

Verify equations (14)--(19) and Definitions 1--2. Decide whether the equilibrium
sequence lists every state and condition needed to make the transition
well-defined, including the inherited old distribution, cohort masses, the
lagged house price used for tax basis, market clearing, government balance, and
the terminal condition. Check that capital gains vanish at both steady states.

Assess whether conditioning on the existence of `E^0` and `E^1` is sufficient
for this minimal note. Do not demand a global existence theorem simply for
completeness. If some additional local determinacy or terminal condition is
essential, state it exactly.

### 3. Conditional welfare and misallocation

Audit Definition 3, Proposition 1, its proof, and Corollary 1. The central
questions are:

- Is "conditional housing efficiency" a valid local Pareto/compensated
  reallocation concept among living households, or is it only a surplus test?
- Does the compensation argument correctly handle tax revenue and rebates at
  both `t` and `t+1`?
- Is the claimed surplus

  `zeta_it^{O,F} + q ell_{t+1} + ell_t - K_ijt`

  correct under the stated timing and fixed-price comparison?
- Does the steady-state corollary correctly reduce the gap to
  `zeta_i^{O,F} - K_ij`?
- Is any statement stronger than the proof permits once fertility, tenure,
  saving, cohort masses, and prices are held fixed?

The note should not become a social-welfare theorem over present and future
cohorts. If the current terminology overclaims, provide the narrowest accurate
replacement.

### 4. Fertility comparative static

Independently rederive equations (23)--(25). Check the lifetime resource
identity, the derivative of fertility with respect to constrained housing
access, and the stated sufficient condition. Verify that the text keeps this
comparative static distinct from the conditional reallocation result.

### 5. Figures

The author is unsure whether the note still has the right plots. Audit the two
current figures rather than automatically adding more:

1. The illustrative paths for `P_t` and `Y_t` between `E^0` and `E^1`.
2. The steady-state marginal-value wedge between a constrained young owner and
   an old incumbent.

Decide whether the minimum note should contain zero, one, or two figures. For
each retained figure, state the exact economic proposition it establishes.
Flag any figure that merely draws an unproved path. If the transition figure is
not defensible, choose between deleting it, replacing it with a timing/state
diagram, or deriving a genuinely solved illustrative example. Do not recommend
a numerical transition solely to decorate the note. Do not revive a
steady-state existence plot unless it is needed for the argument.

## Deliverable

Return a concise audit, not a rewritten paper. Use this structure:

1. **Bottom-line verdict**: whether the two intended theory pieces work.
2. **Claim-blocking issues**: errors that prevent circulation.
3. **Necessary corrections**: the minimum changes required before discussion.
4. **Figure decision**: exact number and content of figures.
5. **What already works**: objects that should be frozen rather than reopened.
6. **Minimal revised architecture**: section and result order, within five
   pages.

For every criticism, cite the relevant equation, definition, proposition, or
figure. Classify each recommendation as **claim-blocking**, **necessary**, or
**optional**. Supply replacement equations or short replacement language only
where needed to make a correction unambiguous.

Do not edit any repository files, run the quantitative model, launch searches,
or expand the model. Stop after the audit and wait for the author to decide
which corrections to implement.
