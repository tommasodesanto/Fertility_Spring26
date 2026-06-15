# ChatGPT Pro Review Prompt - Part 1 and Short Note (2026-06-15)

Attach, in this order:

1. `latex/intergenerational_housing_fertility_part1.tex`
2. `latex/intergenerational_housing_fertility_part1.pdf`
3. `latex/intergen_housing_fertility_short_note.tex`
4. `latex/intergen_housing_fertility_short_note.pdf`
5. `memory/project_part1_policy_stat_audit.md`
6. `latex/figures/example_misallocation.pdf`
7. `latex/figures/example_misallocation_only.pdf`
8. `latex/figures/example_fertility_cap.pdf`
9. `code/model/tools/make_part1_example_figure.py`

Do not attach the older advisor note. The current compressed mirror is
`intergen_housing_fertility_short_note`.

---

You are refereeing the analytical theory front end of a quantitative macro /
urban economics paper on housing access, intergenerational misallocation, and
fertility. The long document, `intergenerational_housing_fertility_part1`, is
the canonical Part 1. The short note, `intergen_housing_fertility_short_note`,
is a compressed mirror of the same core mechanism after a recent line-by-line
revision. Review both with top-five-referee rigor. Derive the equations and
proof steps; do not pattern-match.

The goal is not a prose rewrite. The goal is to find mathematical errors,
proof gaps, notation collisions, overclaims, internal inconsistencies, and
places where the short note and Part 1 now diverge in substance.

## What to verify hardest

1. Environment, stationarity, transition kernels, intergenerational books, and
   bequest/goods accounting in Part 1. Check that old housing wealth,
   inherited assets, bequests, gains-tax payments, and goods feasibility are
   not double counted.
2. Young household problem. Verify the combined choice over consumption,
   housing, fertility, savings, tenure, and access constraints. Confirm that
   owner financing uses the intended cap and that the notation treats
   `phi` as the financed share.
3. Owner effective price and access gaps. Check the derivation of
   `q^O = q + \bar\ell/(1+r)`, `q^R=q`, and the gap
   `\zeta_i^m = \alpha(c_i^m-\chi n_i^m)/(h_i^m-\kappa n_i^m)-q^m`.
   For `\zeta_i^{O,F}`, `F` means financing.
4. Savings margin. Verify the Euler inequality and the closed-form savings
   rule with `k=\beta(1+\gamma+b)` under the stated hypotheses. Check that
   all downstream `(1+k)` deformations reduce to the no-savings formulas at
   `k=0`.
5. Child-price proposition. Re-derive the fertility FOC displayed as
   willingness to pay for a marginal child equals the full child price, and
   check the claim that a binding housing constraint acts as a tax
   `\kappa\zeta_i^m` per child.
6. Efficiency proposition. Verify the local planner-equilibrium comparison:
   the variation resizes two owner-occupied houses and converts no tenures.
   Check the surplus display term by term: buyer effective price plus access
   gap, incumbent net-of-tax sale price, market price cancellation, tax terms
   as recoverable transfers, and `r_i^F` as the real implementation cost in
   Part 1.
7. Worked example and figures. Recompute all numbers in Part 1. In particular,
   verify the compensated-trade reading of the `0.88` surplus residual: pay
   the incumbent at `q-\bar\ell`, charge the buyer at `q^O+\zeta^{O,F}`, hold
   both fixed in utility, and compute the goods residual left over. The long
   document intentionally keeps the combined figure; the split panels in the
   short note are not a mandate to split Part 1.
8. Fertility response to relaxing the family-housing cap `H`. Verify the
   derivative, the two-force interpretation, and the sign threshold
   `\chi \le (1+k)\kappa q^m/\alpha`: extra space relaxes the child-space
   floor, while paying for that space tightens the goods floor, with savings
   weakening the second force.
9. Local policy decomposition. Check the lever-by-parameter framing:
   down-payment grants raise `a_i`, financial deregulation raises `\phi`, and
   family-sized rentals raise `\bar h^R`. Verify the claim that, with fixed
   aggregate stock, these are access reallocations rather than construction or
   aggregate supply shocks.
10. Tenure-switching proposition in Part 1. Re-derive the boundary formula,
    dominance conditions, fertility-jump decomposition, proof, and worked
    boundary numbers. The short note deliberately compresses this material;
    do not infer that Part 1 should drop it.

## Settled design decisions - do not relitigate

- Part 1 keeps the full `r_i^F` display, trichotomy, and proof.
- Part 1 keeps the full switching proposition, proof, and worked boundary
  numbers.
- Part 1 keeps the efficiency proof.
- Part 1 keeps the combined worked-example figure unless the author separately
  decides otherwise.
- Continuous completed fertility is the analytical object here. Do not convert
  the note into a childlessness or sequential parity-progression model.
- Tenure choice is a deterministic argmax in the theory. Logit entry terms,
  where relevant, enter through common factors that cancel in displayed
  ratios.
- There is no balanced growth path. Inflation/repricing language is a
  stationarity device for the gains-tax wedge, not a BGP.
- Notation is settled: `\beta` is the discount factor, `\vartheta` is the
  fertility weight, `k=\beta(1+\gamma+b)`, `\phi` is the financed share,
  `q^m` is the mode-specific effective price, and the primitives are written
  through `(c-\chi n)` and `(h-\kappa n)`.
- The planner comparison is intentionally local and conditional. Entry,
  tenure assignments, savings, and transition kernels are fixed unless a
  displayed statement explicitly varies them.

## Specific risk areas

- Prior June 12 prompts are stale: they refer to the advisor note and older
  repricing notation. Treat this prompt and the attached current files as the
  live review target.
- The short note may state the efficiency condition more broadly than Part 1.
  Flag any substantive mismatch.
- Verify that the worked-example prose clearly distinguishes both-tenures-price
  space at `q` from the owner lifecycle effective price `q^O`.
- Check whether Part 1's `\Theta_i^m` notation is used consistently when it
  moves between `\lambda_i^m\zeta_i^m` and normalized marginal-utility gaps.
- Scrutinize the inflation/stationary-gains device for conceptual coherence,
  but do not propose a BGP rewrite.
- Scrutinize the phrase "access, not construction" around `\bar h^R` with fixed
  stock: if it needs a feasibility qualifier, state the minimal qualifier.
- The measurement bridge for `\zeta^{O,F}` from household bundles and
  `\bar\ell` from deeds/assessments is only sketched. Flag what must be added
  for an empirical section, but do not rewrite it.

## Deliverables

Return numbered findings only. For each finding give:

- location: document, section, display, or nearby text;
- severity: must-fix, should-fix, optional, or reviewer-uncertain;
- diagnosis: the math/proof/notation/claim problem;
- minimal fix: the smallest textual or mathematical correction.

Separate the findings into:

1. mathematical or proof errors;
2. notation or definition problems;
3. claims stronger than proved;
4. inconsistencies between Part 1 and the short note;
5. figure/example numerical issues;
6. prose that obstructs the economics.

Do not rewrite sections wholesale. Do not restate the whole model back to me.
Do not recommend deleting settled Part 1 material merely because the short note
compresses it.
