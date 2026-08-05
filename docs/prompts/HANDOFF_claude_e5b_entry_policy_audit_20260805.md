# Claude Code prompt: audit the restored E5b calibration and policy mechanism

Use this entire document as the task. Do not ask me to restate the context.

Repository:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

## Role and objective

Act as an independent senior reviewer of a quantitative lifecycle model of
housing and fertility intended for a top economics journal. Audit the restored
pre-experiment E5b model, its calibration, and especially its policy mechanism.

The central concern from the author is:

> The weak fertility response to housing policy may not demonstrate that
> housing is unimportant. The model's entry mechanisms and equilibrium closure
> may make the policy comparison unbalanced or mechanically biased against a
> fertility response.

Treat this as a falsification exercise. Determine whether that concern is
correct, which entry mechanism is responsible, and what the smallest credible
repair or counterfactual design would be. Do not assume the existing agents'
interpretations are correct.

## Boundaries

- This is a read-only audit. Do not change model code, calibration targets,
  status files, or paper files; do not commit or push.
- You may write audit outputs only under
  `output/model/e5b_entry_policy_audit_claude_20260805/`.
- Small, targeted local solves are allowed. Do not launch a broad calibration
  or cluster search.
- Preserve the dirty worktree. Start with `git status -sb`.
- Report economically, not as a software review. Explain each object before
  using project shorthand.
- If you recommend a new economic mechanism, distinguish established evidence
  from conjecture and provide direct references. If claiming precedent, give
  priority to top-five economics papers from the last ten years and say clearly
  when a reference is instead a field journal, working paper, or demographic
  source.

## Mandatory orientation

Read in this order:

1. `AGENTS.md`
2. `memory/AGENT_MEMORY.md`
3. the latest `memory/daily/YYYY-MM-DD.md`
4. `CALIBRATION_STATUS.md`
5. `code/model/README.md`

Then inspect:

- `code/model/intergen_eqscale_seq_optimized/README.md`
- `code/model/intergen_eqscale_seq_optimized/IMPLEMENTATION_STATUS.md`
- `code/model/intergen_eqscale_seq_optimized/parameters.py`
- `code/model/intergen_eqscale_seq_optimized/production_profile.py`
- `code/model/intergen_eqscale_seq_optimized/calibration.py`
- `code/model/intergen_eqscale_seq_optimized/e5_profile.py`
- `code/model/intergen_eqscale_seq_optimized/run_e1_chain.py`
- `code/model/intergen_eqscale_seq_optimized/solver.py`
- `code/model/intergen_eqscale_seq_optimized/kernels.py`
- `code/model/intergen_eqscale_seq_optimized/build_e2_packet.py`
- relevant tests under
  `code/model/intergen_eqscale_seq_optimized/tests/`

Interpretive evidence, to verify rather than trust:

- `docs/model/e5_target_review_20260724.md`
- `docs/model/e5b_calibration_autopsy_20260726.md`
- `docs/model/entry_postponement_mechanism_memo_20260725.md`
- `docs/model/e6_decision_package_20260727.md`

The top of `CALIBRATION_STATUS.md` is authoritative: the author rejected E6a,
E6b, E6c, and all combined refits. They are experiments only. Do not promote
the late-age fecundity extension, permanent income classes, or readiness state
back into the maintained model merely because they lowered a calibration loss.

## Maintained model

The model to audit is the certified E5b specification:

- one housing market;
- four-year periods beginning at age 18;
- sequential births with literal parity states (0,1,2,3+);
- an imposed equivalence scale for children;
- a five-state persistent earnings process but no permanent earnings classes;
- deterministic tenure choice;
- renter housing continuous up to six rooms and owner housing on a discrete
  ladder;
- a financed housing share of (0.80), so the down payment is
  ((1-0.80)pH);
- normalized stationary population in the maintained calibration and policy
  packet;
- ten estimated parameters and twelve target moments.

The first-birth decision is an attempt decision. In schematic form, verify that
the code implements

\[
P_{\mathrm{try}}(s,a)
=\frac{\exp(V^{\mathrm{try}}(s,a)/\kappa_E)}
{\exp(V^{\mathrm{wait}}(s,a)/\kappa_E)
 +\exp(V^{\mathrm{try}}(s,a)/\kappa_E)},
\]

where

\[
V^{\mathrm{try}}(s,a)
=\pi(a)V(s',n=1)+(1-\pi(a))V(s',n=0),
\qquad
V^{\mathrm{wait}}(s,a)=V(s',n=0).
\]

Derive the exact ordering of fertility, tenure, location, consumption, housing,
and saving choices from the code. Establish which prices and housing
constraints enter (V^{\mathrm{try}}-V^{\mathrm{wait}}), and at what point.

The maintained biological schedule is

\[
\pi(a)=1-0.02\exp\{0.134(a-18)\}
\]

during the active fertility window, followed by the hard terminal closure.
Do not substitute the rejected E6a tail.

## Certified E5b calibration

Canonical artifacts:

- `output/model/eqscale_seq_e5b_recalibration_20260725/report/results.json`
- `output/model/eqscale_seq_e5b_recalibration_20260725/report/target_fit_full.csv`
- `output/model/eqscale_seq_e5b_recalibration_20260725/report/parameter_bounds.csv`
- `output/model/eqscale_seq_e5b_policy_packet_20260725/`

The selected result has reported loss `385.142880159`, strict market residual
`1.1091e-6`, and exact repeated evaluation. Recompute the objective directly.
The loss is a calibration criterion, not a formal J statistic: some rows use
measured standard errors and others use a synthetic five-percent-of-target
rule.

Complete target fit:

| Moment | Target | Model | Weight | Loss contribution |
|---|---:|---:|---:|---:|
| Completed fertility | 1.918000 | 1.903573 | 1425.739 | 0.297 |
| Childlessness | 0.188000 | 0.080188 | 17180.744 | 199.700 |
| Mean age at first birth | 25.310561 | 25.637433 | 16.000 | 1.710 |
| First births at age 30+ | 0.270062 | 0.332334 | 15625.000 | 60.590 |
| Rooms response to first child | 0.664435 | 0.657205 | 906.056 | 0.047 |
| Large-family rooms contrast | 0.367700 | 0.372240 | 2958.515 | 0.061 |
| Family ownership gap | 0.167662 | 0.165878 | 14229.591 | 0.045 |
| Ownership rate | 0.575472 | 0.563758 | 1207.846 | 0.166 |
| Mean occupied rooms | 5.779970 | 6.075505 | 11.973 | 1.046 |
| Wealth / annual gross earnings | 6.873100 | 5.454754 | 6.288 | 12.649 |
| Bequest flow / wealth | 0.008800 | 0.009078 | 5165289.256 | 0.399 |
| Old wealth p90/p50 | 3.448111 | 2.068370 | 56.960 | 108.433 |

Eight rows fit closely. Childlessness, first-birth timing, aggregate wealth,
and old-age wealth dispersion account for essentially all remaining loss.

Estimated parameters and search bounds:

| Parameter | Estimate | Bounds |
|---|---:|---:|
| Annual discount factor | 0.991097 | [0.94, 0.9995] |
| Per-child consumption/housing share shift | 0.012163 | [0, 0.25] |
| First-child share jump | 0.023896 | [0, 0.25] |
| Utility flow from a child | 0.585018 | [-3, 3] |
| First-birth logit scale | 3.552006 | [0.02, 50] |
| Later-birth logit scale | 1.254006 | [0.02, 50] |
| Owner housing-service premium | 1.038033 | [0.10, 5] |
| Housing supply intercept | 8.162391 | [0.20, 80] |
| Bequest strength | 0.168435 | [0, 8] |
| Bequest curvature shifter | 0.050046 | [0.02, 16] |

## Policy experiments to audit

The policy packet re-solves equilibrium prices at the fixed E5b parameters.
It contains two interventions:

1. A grant of `0.4` model wealth units associated with renter-to-owner entry
   into owner rungs corresponding to homes of at least six rooms while in an
   eligible child state.
2. The same grant plus a property-tax increase from one to two percent per
   year (`tau_H` rises from `0.04` to `0.08` in four-year units).

Reported outcomes:

| Case | Completed fertility | Childlessness | Ownership | House price |
|---|---:|---:|---:|---:|
| Baseline | 1.903573 | 0.080188 | 0.563758 | 0.815625 |
| Grant | 1.905159 | 0.079590 | 0.635743 | 0.819360 |
| Tax plus grant | 1.906023 | 0.079389 | 0.619109 | 0.661504 |

Thus the grant raises ownership by 7.2 percentage points but completed
fertility by only 0.0016. The combined policy lowers the house price by about
19 percent but raises completed fertility by only 0.0025.

These are currently mechanism diagnostics, not paper-ready policy estimates.
The grant is not financed through a government budget, no welfare calculation
is supplied, and there is no clean E5b tax-only case.

## Central audit: is the policy comparison structurally unfair?

Separate four different objects that may all be called “entry.”

### A. Entry into parenthood

- Verify the exact first-birth value difference and logit probability.
- Determine whether observed postponement and childlessness are generated
  mainly by the large logit scale rather than economic state dependence.
- Quantify the effect of house prices, rents, wealth, tenure, and the
  down-payment constraint on (V^{\mathrm{try}}-V^{\mathrm{wait}}).
- Compare those changes with the child utility flow and the logit scale.
- Test the existing claim that a 10 percent rent change moves the first-child
  surplus by only about `0.0033` against `psi_child≈0.585` and
  `kappa_fert≈3.55`.
- Explain whether the logit specification mechanically attenuates policy
  elasticities because it is simultaneously being used to create delay and
  non-entry.

### B. Entry into homeownership

- Trace the grant through the Bellman problem and forward distribution.
- Establish whether it is paid only on the birth transition, or whenever an
  eligible parent who remains in the child-at-home state later buys a home.
- Establish whether a household can receive it more than once.
- Verify that the grant is added before the down-payment and borrowing-floor
  tests in both the Bellman and forward equations.
- Determine whether the 7.2-point ownership response is a genuine behavioral
  response or a mechanical relaxation of an eligibility threshold.
- Check whether deterministic tenure choice creates discontinuous treatment
  responses that would differ materially under a small, empirically disciplined
  tenure-choice smoothing parameter.

### C. Entry of new households and population scale

- Verify whether E5b uses fixed entrant mass `1/J` and then normalizes the
  stationary population to one.
- Determine whether `entrant_conversion_factor=0.5` affects the maintained
  normalized closure at all, or merely records matured-child flows that do not
  feed back into new household mass.
- Explain whether a fertility policy can affect only family composition but
  not future cohort scale in the current packet.
- Determine how normalization affects housing demand, supply, prices, total
  births, and comparisons across policies.
- State whether this closure biases the fertility response, the price response,
  the ownership response, or only the interpretation of aggregate population
  effects.

### D. Fiscal balance and welfare

- Verify that the grant is an unfunded resource injection and property-tax
  revenue is neither rebated nor used to finance the grant.
- Compute the government-budget residual implied by each policy if possible
  from existing objects.
- Explain whether comparing a free grant with a tax-plus-grant case is
  economically meaningful without equalizing fiscal resources.
- Identify the correct welfare object and whose welfare should be compared:
  entrants, current households, newborn cohorts, or the stationary population.
- Explain how a revenue-neutral comparison should treat population scale and
  entry.

## Required targeted checks

Do not run a new calibration. Use existing caches where reliable and fresh
strict solves only where needed. At minimum:

1. Reproduce the E5b baseline and both policy cases, including all twelve
   moments and the equilibrium price.
2. Add a read-only tax-only evaluation at the same two-percent annual property
   tax if this can be done with the existing driver without modifying source.
3. Calculate treatment eligibility and take-up: mass eligible for the grant,
   renter-to-owner transitions receiving it, grant expenditure per household,
   and expenditure per additional birth.
4. Decompose the first-birth response by age, wealth, income state, renter versus
   owner, and whether the down-payment constraint binds.
5. Report the distribution of the first-birth value gap
   (V^{\mathrm{try}}-V^{\mathrm{wait}}), the implied attempt probability, and
   its derivative with respect to the house price at representative states.
6. Compare the current normalized-population closure with the logic of an
   endogenous entry/scale closure. A full new solver is not required: state the
   accounting equations and, if existing code supports it safely, run the
   smallest illustrative comparison.
7. Check whether the policy result is robust to evaluating household choices at
   fixed prices. This separates direct household incentives from the general-
   equilibrium house-price response.

For every new check, record the exact overrides, evaluator, convergence status,
and whether it is a diagnostic or a valid counterfactual.

## Economic questions

Answer directly:

1. Does the weak fertility response establish that housing policy has little
   effect in the economic model, or only that the current parenthood-entry
   specification insulates fertility from prices?
2. Is the large ownership response evidence for the housing mechanism, or an
   artifact of the grant's threshold implementation and deterministic tenure?
3. Is the stationary normalization innocuous for household-level fertility,
   or does it make the general-equilibrium comparison incomplete?
4. Which conclusion survives all reasonable interpretations of entry and
   fiscal closure?
5. What is the smallest model or policy-design change required before these
   experiments could support a paper claim?
6. Which additional moments would identify that change without weakening the
   existing twelve-target system?

Do not recommend permanent household types as a generic cure. The author has
rejected that device. If you believe heterogeneity is indispensable, explain
why, name the exact empirical variation that identifies it, and compare it with
alternatives such as externally measured child earnings costs, partnership or
employment transitions, and a housing-readiness constraint.

## Deliverables

Write:

`output/model/e5b_entry_policy_audit_claude_20260805/AUDIT.md`

with this structure:

1. **Executive verdict** — no more than one page.
2. **Exact model and policy equations** — code-to-math mapping with file and
   line references.
3. **Calibration audit** — full target table, all parameters and bounds,
   identification, and weight interpretation.
4. **Four entry mechanisms** — parenthood, homeownership, population, and
   fiscal entry/closure, clearly separated.
5. **Policy decomposition** — direct, price-mediated, eligibility, fiscal, and
   scale effects.
6. **Numerical checks** — exact commands, results, convergence, and limitations.
7. **Literature assessment** — only claims supported by verified references.
8. **Recommendation** — the smallest decisive next experiment, with equations,
   identifying moments, and outcomes that would change the verdict.

Also write compact machine-readable tables for:

- reproduced target fit;
- parameter estimates and bounds;
- baseline/grant/tax-only/tax-plus-grant outcomes;
- grant eligibility, take-up, and fiscal accounting;
- first-birth response decompositions.

End with one of these classifications, justified in plain language:

- **Policy result valid as currently interpreted**;
- **Valid only as a narrow mechanism diagnostic**;
- **Directionally informative but quantitatively biased by entry/closure**;
- **Not interpretable until the entry or fiscal mechanism is repaired**.

Do not implement the repair. Specify it precisely enough that the author can
decide whether it belongs in the model.
