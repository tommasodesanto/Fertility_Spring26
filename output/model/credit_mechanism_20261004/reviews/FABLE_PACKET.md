# Fable 5.1: adversarial review of credit, family space and fertility

## Goal

Tommaso asks for focused overnight research, not reassurance or a large report. A weak mortgage-credit fertility effect makes him doubt the paper. You are Fable5.1 independent reviewer, to debate an Astra max agent after first pass. Determine whether the unchanged model offers a credible sharp contribution; if not identify minimal evidence-grounded modification. Do not retune merely to produce larger effects.

## Scope

Read mandatory startup: memory/AGENT_MEMORY.md, latest memory/daily, CALIBRATION_STATUS.md, code/model/README.md. User claims: output/model/credit_mechanism_20261004/evidence/USER_CLAIMS.md. Work plan: output/model/credit_mechanism_20261004/README.md.

Read active defaults code/model/parameters/best_params.py, production/inputs.py, production/equilibrium.py, production/engine/household.py and shared.py/kernels.py only relevant birth/tenure/budget sections. Birth logit, common-state gain from credit, soft collateral/purchase rule, rental alternatives, equal saving/borrowing rate are key. Review exact theory public loan statement via latex/README.md and model manuscript sources READ ONLY. Supporting fit source: output/model/fixed_reference_economics_20260928/soft_timing_calibration_20261002_v1/collection/RESULTS.md (all14targets31parameters), production/default cached solution and output/model/experiments/credit_distribution_report_v1/README.md + manifest.json. Another Sol agent is preserving scratch data; do not duplicate whole archive/session searches. The claimed scratch numerical effects are unverified pending recovery, not grounds for adopting policy conclusions.

## Questions

1. Is incorrect income sorting necessary/sufficient to explain weak credit effect, or is small child-relative financing benefit structural? Identify evidence that distinguishes these. Check financial/state timing and distribution weights before economic explanations.
2. What in the existing model could deliver meaningful interaction of fertility and constraints even if unconditional mortgage response is modest? What result would refute this paper narrative? Separate positive mechanisms from normative theory/public-loan result and aggregation/GE closure.
3. Which empirical income gradient is relevant (household vs women, permanent vs current income, age/wealth/family sorting)? A means-tested child utility term changes preferences directly; fitting gradient by it is not a clean causal experiment. Do not assume time costs can never reverse gradient without proof for this actual utility.
4. Rank at most THREE decisive checks, computationally tiny where possible. Specify predicted outcomes under rival hypotheses and what decision each result settles. Minimal changes only if evidence demands them; name replacement identifying moments or external restrictions.
5. Property-tax gross-vs-net supplier accounting is separate from fertility-credit channel. Root closure does not force same price across policies and stationary endpoints do not prove feasible transition or welfare improvement.

## Do not touch

Read-only tools only. No source/manuscript edits, numerical jobs, tests, calibration/target changes, external messaging, Git, or subagents. Do not endorse a run or source identity not inspected. Do not silently adopt Estate-A, three-birth, new wealth or CES experimental inputs; reference remains working soft/post-interest old-target chain13 without these additions unless clearly comparison.

## Required output

At most1500words in your response: strongest diagnosis with source lines; what remains unknown;3specific falsifiable tests; unchanged-model paper contribution and limits; necessary-change decision conditional on evidence. Keep conjecture versus verified facts explicit. No full code or raw file dumps. Debate follow-ups will carry Astra critique and verified Sol results.

## Verification and stop condition

45-minute limit; stop after bounded evidence supports the diagnosis or required source is unavailable. Model must be claude-fable-5-1; never silently fallback. Use authenticated first-party CLI. No API billing fallback. The lead will independently check your math, decisive source passages, and eventual numerical results.
