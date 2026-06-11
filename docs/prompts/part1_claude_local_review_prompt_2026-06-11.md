# Claude Local Review Prompt: Part I Analytical Model

We are in:

`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

First follow the repository startup protocol:

1. Read `memory/AGENT_MEMORY.md`
2. Read the latest `memory/daily/YYYY-MM-DD.md`
3. Read `CALIBRATION_STATUS.md`
4. Read `SESSION_DIARY.md` only if historical context is needed
5. Run `git status -sb`

Treat existing dirty/untracked files as user or prior-agent work. Do not clean, revert, stage, or commit anything unless explicitly asked.

## Task

Audit the current standalone Part I analytical model as economics, not as LaTeX surface prose. This is an audit-only task: do not edit files yet. Provide a rigorous written assessment and compact LaTeX-ready patches only where useful.

## Files To Read

Primary maintained draft:

- `latex/intergenerational_housing_fertility_part1.tex`
- `latex/intergenerational_housing_fertility_part1.pdf`

Working audit/context:

- `docs/model/part1_deep_audit_policy_sufficient_stat_2026-06-10.md`

Figure concepts:

- `latex/part1_policy_sufficient_stat_figures_20260610.tex`
- `latex/part1_policy_sufficient_stat_figures_20260610.pdf`
- optional previews:
  - `output/part1_policy_figures_20260610/part1_policy_figures_page_1.png`
  - `output/part1_policy_figures_20260610/part1_policy_figures_page_2.png`

For comparison only if needed, not as source of truth:

- `latex/intergenerational_housing_fertility_v4.tex`
- `latex/intergenerational_housing_fertility_v4.pdf`

## Project Context

The paper studies intergenerational housing allocation, tenure choice, rental-owner segmentation in family-sized housing, and fertility. The model should read like a clean macro theory paper:

environment, households, choices, constraints, aggregate equilibrium, planner, CE/planner comparison, then fertility and policy implications.

Important preferences:

- Start from the economy, not from wedges or mysterious technical objects.
- Avoid decorative terminology and avoid wedge-first exposition.
- Keep the economics fixed unless you identify a genuine inconsistency.
- Use Coven et al.-style aggregate housing equilibrium as the reference point: tenure choices, segmented markets, construction/supply, market clearing, then planner comparison.
- The child-input primitive is Stone-Geary: children require both consumption goods and housing services.
- The core fertility mechanism should not depend on old-side property-tax lock-in. Property tax can be an old-side retention force, but young constrained access should stand on its own.

## Current Draft Features To Audit

The standalone Part I currently has:

1. Young households choosing tenure through rental and owner branch problems.
2. Rental and owner housing segments with different size menus, so large family-sized units are more available in the owner segment.
3. Segment construction/supply and aggregate market clearing.
4. A constrained planner with real menus and private financial implementation constraints distinguished.
5. A Section 9 fertility result based on the effective family-housing cap \(H_i^m\).
6. A local policy sufficient statistic:
   \[
   \frac{dN}{dz}
   =
   \bar M
   \int_{\mathcal B(z)}
   \left[
   \pi_i^E\Psi_i^m
   +
   n_i^m
   \frac{\pi_i^E(1-\pi_i^E)}{\kappa_E}
   \lambda_i^m\zeta_i^m
   \right]
   P_i^m(z)
   \,dG_Y(i).
   \]

The intended interpretation is: policy effectiveness is exposure to constrained households, times pass-through to the effective family-space cap, times fertility response, plus entry. If a policy does not relax the binding component of \(H_i^m\), then \(P_i^m(z)=0\) for this mechanism.

## What To Produce

Please respond in six sections:

1. **Executive Diagnosis**
   - Is this now a coherent Part I for a macro theory/quantitative paper?
   - What are the 5-8 highest-priority remaining issues?

2. **Mathematical Audit**
   - Verify the Stone-Geary child input structure.
   - Check the Section 9 derivative.
   - Check the equal-scaling condition \(\chi_i=\kappa q^m/\alpha\).
   - Check the sufficient-statistic formulas.
   - Correct any sign, notation, or active-set mistakes.

3. **Economic Audit**
   - Are tenure choice, rental-owner segmentation, construction, old retention, and planner feasibility coherently integrated?
   - Are real menu constraints and private financial constraints treated correctly?
   - Is the CE/planner result properly aggregate and conditional, or still too wedge/local-perturbation-like?

4. **Policy Statistic**
   - Is the sufficient statistic useful and correctly qualified?
   - Can it be made simpler, sharper, or closer to what a referee would accept?
   - What is true only locally/fixed-price/fixed-active-set versus in full equilibrium?

5. **Figures**
   - Are the two figure concepts the right figures?
   - Suggest better titles, labels, or alternative figures if needed.

6. **Revision Patches**
   - Provide compact LaTeX-ready replacement text only where useful.
   - Do not rewrite the whole document.
   - Do not edit files directly unless explicitly asked after the audit.

Be explicit about uncertainty. If the current draft overclaims, say exactly where and how to qualify it.
