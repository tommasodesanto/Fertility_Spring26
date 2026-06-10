# Scoring Rubric

Use this rubric after Claude returns its pass. Score each dimension from 0 to 4.

## 1. Theory Precision

- 0: summarizes prose without checking equations.
- 1: spots superficial notation/style issues.
- 2: checks some FOCs or definitions but leaves key assumptions vague.
- 3: gives a coherent CE/planner comparison and identifies proof gaps.
- 4: derives the result cleanly, distinguishes constraints from externalities,
  and proposes journal-quality statement/proof fixes.

## 2. Economic Judgment

- 0: treats every model ingredient as equally important.
- 1: gives generic macro comments.
- 2: recognizes housing/fertility mechanisms but misses identification or
  welfare subtleties.
- 3: diagnoses the role of outside option, bequests, tenure, liquidity, taxes,
  and market clearing.
- 4: gives a clear ranking of which assumptions are defensible, which are weak,
  and which need a different planner/instrument formulation.

## 3. Code Comprehension

- 0: does not inspect code.
- 1: lists files without explaining data flow.
- 2: maps major modules but cannot link equations to implementation.
- 3: identifies active solver/objective/calibration paths and concrete risks.
- 4: gives file/line-grounded findings, cheap verification commands, and useful
  test proposals.

## 4. Calibration Discipline

- 0: reports losses or old runs without context.
- 1: repeats `CALIBRATION_STATUS.md`.
- 2: checks targets but does not verify comparability.
- 3: distinguishes live/stale outputs, target systems, weights, closure rules,
  and feasible next experiments.
- 4: produces a rigorous target-versus-model table or explains exactly why it
  cannot be produced cheaply, then designs a recoverable calibration plan.

## 5. Evidence And Reproducibility

- 0: no file paths, no commands.
- 1: vague references to files.
- 2: some paths and command summaries.
- 3: file/line references, command log, and clear unresolved checks.
- 4: a reproducible output folder with report, logs, and proposed scratch
  artifacts.

## 6. Usefulness For The Next Work Session

- 0: interesting but not actionable.
- 1: too broad or too generic.
- 2: contains useful ideas but no prioritization.
- 3: gives a credible next-day work plan.
- 4: identifies the single highest-return next move and separates it from
  lower-priority cleanup.
