# Codex worker task

## Goal
Independently derive and review the natural-solvency feasibility logic for the exact block0506 household, supporting a lead-owned isolated implementation. No calibration judgment.

## Scope
Load mandatory current project context with bounded reads; prior reviewed context is in output/model/fixed_reference_transition_20260928/preparation_v1/transition_readiness.md and runtime_review.md. Inspect only solver.py (_savings_stage, _tenure_location_stage, solve_bellman_full_markov_income), matching compiled saving/choice kernels, child-aging and fertility mappings, and saved parameter manifest relevant fields. Model path code/model/intergen_eqscale_seq_optimized/. The frozen source will be authenticated on Torch before numerics. No unrelated scans.

## Context
Author-defined experiment removes artificial purchaser/down-payment, renter and incumbent-owner borrowing constraints while retaining repayment and no-default solvency, existing information order, income/survival/child transitions, death net liquidation, and every preference including psi. Current native_solvency_credit is an approximate prototype using sentinel values and first feasible grid node; do not certify it. We need an explicit backward Boolean feasibility/support recursion independent of utility magnitude, with continuous natural thresholds separate from numerical grid support. All household information already known at saving time vs later risks must be respected. Positive mortality makes net-estate solvency potentially binding even with future earnings. No subsidy, default or mortality insurance may be invented.

## Do not touch
Read-only. No model imports, numerical tests, hashing, rendering, SSH, jobs, active code edits or Git. No extra agents. Lead owns implementation, specification, tests and source pins.

## Required output
At most 900 words: exact feasibility recursion and timing; where constraints must change; pitfalls in interpolation, support holes and fertility probabilities; up to 12 meaningful boundary fixtures; concrete ambiguous economic choices if any. Review any credit_v1/natural_credit.py or build_overlay.py that exists near completion, without waiting beyond the budget.

## Verification
Trace actual call sites; distinguish native equations from proposed changes. Give file/line references. Lead will compare your recursion to the implementation line by line.

## Stop and report if
20-minute default profile cap; return best evidence and remaining uncertainty, no restart or scope expansion.
