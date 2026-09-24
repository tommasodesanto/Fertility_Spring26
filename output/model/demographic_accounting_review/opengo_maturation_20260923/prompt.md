# Task: population accounting for tomorrow's author decisions

## Goal
Perform a bounded read-only source review of births, dependent-child maturation, new household entry, and parental death. Produce a short, concrete decision aid for Tommaso tomorrow. You are the requested OpenCode Go worker; do the substantive review yourself. Maximum 30 minutes. Do not manufacture work to fill time.

## Scope and context
Full project startup required: memory/AGENT_MEMORY.md, latest memory/daily/2026-09-23.md, CALIBRATION_STATUS.md; SESSION_DIARY.md only for a specific historical question. Then code/model/README.md and the actual active source named by canonical status. Read narrowly, do not dump whole files or scan archives. Reuse existing demographic review in /Users/tommasodesanto/.codex/worktrees/d606/Fertility_Spring26/output/model/demographic_accounting_review/ if accessible, especially population_conversion_literature.md; do not redo its literature audit. Distinguish live code, frozen September14, proposed changes, and source defaults versus serialized run values. State exact source versions reviewed.

Author accepts conventional household abstraction: two adults (one man, one woman), one decision-maker, integer children individually counted, distinct first-child cost. No marriage/bargaining/single-household model requested. 2.1 mean completed fertility is settled as reference normalization. Existing entry factor1/2.1 and literal two-person conversion1/2 must be distinguished; no change adopted. User leans toward retaining 18-year average child departure, with an age-tail cutoff discussed but not adopted. Tomorrow starts with maturation, then parental death. Formation-factor robustness is low priority, not tonight's experiment.

## Questions
1. What exactly do implemented m (children at home), n (children ever born), maturation/departure probabilities, adult entrants, and household death count? Trace relevant code paths with file/line citations.
2. Trace one birth cohort through aging, maturation and household formation, while parents remain; separately trace dependent children when a parental household dies. Give small hand-worked count examples. Identify whether each child is counted once, omitted, duplicated, or represented only in a separate demographic process. Do not equate a household departure with biological maturation without evidence.
3. Explain precisely the implemented roles of1/2.1,1/2, survival and population normalization; distinguish stationary distributions from transition counts. Does normalization mask nonconservation? Only assert bugs when source proves them.
4. Supply a compact decision table: current rule, economic interpretation, established issue or uncertainty, minimal alternative, consequences for steady state and transition, what would need recalculation AFTER author choice. No recommended default based solely on numerical convenience. Separate accepted abstraction from outstanding accounting.

## Do not touch
No writes to repository, no model/household/GE solves, no code execution or shell commands, no tests, no cluster calls, no calibration, no new source estimates, no retries, no external delegation, no Git operations, no author draft/mock/slides/checklist changes. Use read/glob/grep tools only. No credentials. Existing dirty work is unrelated and must remain untouched. Do not implement any alternative.

## Required output
Return the complete review in your final message (supervisor will save it). Lead with at most six concrete findings, then the worked examples and one decision table. Cite exact paths and lines. Finish with missing evidence, files inspected, and a verification checklist distinguishing source-established facts from hypotheses. Keep under2200words. Avoid jargon and 'parity' in prose. No broad new literature review.

## Verification and stop
Cross-check household and aggregate renewal equations against actual source definitions, using hand arithmetic only. If source identity is ambiguous, report exact discrepancy rather than choosing silently. Stop at30minutes or once deliverable complete. If access/auth/tools block you, return exact blocker; no fallback model or repeated calls.
