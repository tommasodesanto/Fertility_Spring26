# Historical decisions checked before the policy comparison

The read-only `explorer_fast` worker completed its bounded evidence pass. The
lead checked the passages below in context. These are project conversation
records, not independent empirical verification. No new interpretation is
silently attributed to the author.

| Object | What the old record establishes | What it does not establish | Decision for this comparison |
|---|---|---|---|
| Dated 0.63 versus old-equilibrium 1.75 supply elasticity | The September 4 audit already reports both stages explicitly. An August 19 agent explanation describes re-anchoring a different supply elasticity at a fixed inherited price and quantity. | No recovered direct author statement selects the exact 0.63/1.75 pair. The old-equilibrium interpretation remains outstanding. | Preserve both values and date-0 normalization; do not label the distinction a newly discovered coding error. |
| Closed-national birth queue versus coherent-person demographic branch | August records distinguish the inherited birth queue from the later coherent-person construction. The author supports studying housing policies along a fertility-preference transition without exact house-price matching. | That broad author statement does not endorse every later migration, entry or terminal closure. It does not turn the queue into a validated person-population model. | Reproduce the existing finite-horizon mechanism comparison; label long-run population conclusions conditional and do not substitute another branch. |
| Tax experiment and denominator | The September 4 reconciliation explicitly distinguishes births per model household from total births. Its +0.498% rebated-tax result compares rebated 2% with rebated 1%. | It is not the effect of rebated 2% relative to the unrebated status quo, nor a welfare result. | Preserve each denominator and comparison baseline; show births per household and total births separately. |
| Owner premium and historical metadata | The September 4 parameter table calls chi the owner housing-service premium and theta0/theta1 bequest parameters. Current summaries retain old source metadata alongside explicit dated overrides. | Historical source defaults do not describe the final calibrated parameter vector. | Read actual selected theta, dated closure contract and prepared parameters, not isolated old metadata fields. |

## Checked passages

Paths are relative to the project root
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.

- `memory/transcripts/2026-09-04/combined_user_assistant.md:3857` and surrounding
  parameter table explicitly retain dated 0.63 and old 1.75. This is an agent
  audit statement, not a direct author choice of the pair.
- `memory/transcripts/2026-08-19/combined_user_assistant.md:48029`: the agent
  proposes a separate supply-elasticity experiment while retaining the inherited
  price and quantity. This is a precedent for the anchoring convention only;
  its alternative elasticity 0.232 is not the current 0.63 target.
- `memory/transcripts/2026-08-19/combined_user_assistant.md:52015` onward: the
  author describes starting from a steady state, changing preferences and
  studying housing policy during the transition. At line 52033 the author says
  imperfect price matching is acceptable because “it's an experiment.” The
  agent's surrounding discussion uses an open transition; it is not direct
  author approval of the later precise closed-national queue.
- `memory/transcripts/2026-08-19/combined_user_assistant.md:47685` records the
  agent's common fitted distribution/birth-queue plan with M=0, rho=1.
- `memory/daily/2026-08-26.md:87` distinguishes the observed household-head
  age profile, inherited births/2.1 queue and coherent-person-law development.
  Later integration supersedes its then-current implementation status; it does
  not retroactively change the frozen September 4 mechanism benchmark.
- `memory/transcripts/2026-09-04/combined_user_assistant.md:4830` onward gives
  separate 2063 births-per-household and total-birth effects. At line 4842:
  “the +0.498% rebate result compares rebated 2% taxation with rebated 1% taxation.”
  The same passage notes that rebates affect household resources; housing
  services are not the only fertility channel.

The worker did not find a September 5 transcript directory and did not search
unrelated private archives. Its direct-author evidence is limited as stated.
The lead keeps the remaining interpretation questions visible rather than
using the history pass as permission to change the economic specification.

## Reporting-driver difference reconciled

The old policy summary pins driver `628edeada...`; the new summary pins
`345f2b54...`. These match Git commits `7c32661` and `e72c238`, respectively.
The exact saved `policy_driver_comparison.diff` adds only a read-only population
accounting function and three output columns after state advancement. It changes
no policy equation, parameter, target, market solve or advancement. The five
policy definitions, fiscal treatment, population closure, lag, dates and scope
are identical. Numerical-gate summary values differ because they are realized
residuals and solve counts, not changed tolerances. The scientific source bundle
also contains the separately verified later-birth policy-reuse correction.
Neither historical nor current five-policy outputs use a tight-market fallback;
no affected retained-production result has been demonstrated. This comparison
is therefore described as the candidate on the verified repaired source, with
these source differences disclosed, rather than asserting byte-identical code.
