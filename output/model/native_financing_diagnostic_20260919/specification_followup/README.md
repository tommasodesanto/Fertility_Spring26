# Specification follow-up: September 20–21

## Authority and deliverable

On September 20, Tommaso explicitly authorized an active goal, check-ins,
continued experimentation, and longer cluster jobs while away until tomorrow.
The objective is a reviewed decision packet: a proposed provisional model,
evidence for the few choices that matter, unresolved questions, and the
implementation/measurement/calibration sequence after an author decision.
Experimental results do not adopt a new baseline, target system or policy.

Aim to deliver by **09:00 America/New_York, September 21**. This is a working
delivery target, not a deadline explicitly supplied by the author.

## Stages and initial resource envelope

1. Collect and review the existing Claude Max/Fable assessment. Session
   `fbfab0f8-fd81-4d5e-9453-23b6912efaf9` uses `claude-fable-5-1`, maximum
   effort, read-only tools, 60 turns and a 20-minute wall-clock cap. It began
   at 15:15:34 UTC September 20. Raw receipts are in
   `tmp/fable_specification_review_20260920/` at the repository root.
   A failure does not authorize an automatic restart.
2. Diagnose the credit/lifetime-cohort response using saved outputs first.
   Distinguish changed choices, evolving cohort composition, borrowing-support
   limits and measurement. The identical high-dose outcomes need a separate
   explanation. Do not presume either a bug or an economic sign theorem.
3. Establish which existing housing, ownership, wealth and fertility age
   profiles can be compared with consistent empirical definitions. Historical
   plotting shortcuts are not new calibration evidence. For example,
   `house_size_age_model_vs_data.py` uses a fixed renter-room proxy and cannot
   provide an exact realized-room comparison without correction.
4. Use these findings and Fable's assessment to specify at most two additional
   model variants worth testing. Changes must isolate a stated economic
   question and have an unchanged-source control. Limited conditional fitting
   can test feasibility; it is not serious calibration or adoption.

Initial envelope: **at most 16 new fixed-price lifecycle/cohort evaluations
and 32 new full stationary evaluations**, including controls, loop smokes and
selected-point repetitions. Replaying saved arrays does not count as a new
solve. At most eight single-threaded workers concurrently. Each batch must
state its actual dimensions, solve count and timing estimate from the latest
relevant measurement; these maxima are not an instruction to fill the budget.
No new batch after 02:00 September 21; planned numerical work must finish by
07:30 to leave time for collection and review. Each submitted allocation and
controller must enforce its own smaller applicable cap. Do not expand this
envelope silently.

## Required launch contract

Before submission, record a hypothesis, exact source/input/target hashes,
parameter overrides, population/entry contract, held-fixed objects, outputs,
per-case and total time budgets, case limits, and stop conditions in a launch
manifest in this directory. Estimate time from an observed comparable solve.
Use a fresh immutable remote snapshot and the exact-loop smoke; production
depends on successful smoke through Slurm, without laptop-dependent chaining.
No production job is authorized merely by appearing in a draft plan.

Keep the original objective, target definitions, weights, bounds and numerical
gates pinned for any comparable fit. A model variant may require a different
measurement contract; if so label it explicitly and do not compare losses as
though identical. Never add a free parameter without identifying variation or
an external restriction. All core numerical changes require lead review
against the model equations before they are used.

Progress/checkpoint output is required every case or five minutes, with a
latest-completed summary and best-so-far summary when selection is involved.
Investigate after 30 minutes without a heartbeat/checkpoint. On failure,
preserve evidence, stop the affected stage and cancel only never-started
dependents. A second attempt requires a documented changed hypothesis or
method; do not relax gates, silently change populations, or repeat a failed
search with a new label. Independent healthy branches may continue.

Each solved case must retain the standard 17 diagnostic plots and numerical
checks. A selected fitted candidate needs the complete target/parameter tables,
actual search bounds, source/objective fingerprints and two exact repetitions.
Keep snapshot total/first birth flows separate from explicit lifetime cohort
births, lifetime first-birth probabilities and the model's stationary
normalization. Across-checkpoint comparisons do not hold prices, preferences
or entrant wealth fixed unless explicitly verified.

## Check-ins and present state

### Reviewed Fable assessment and current work

The initial Max/Fable review completed in 517.5 seconds. A separately bounded,
focused clarification completed in 208.1 seconds in the same session. The
[assessment and clarification](../../../../docs/model/structural_model_specification_fable_20260920.md)
are retained verbatim, with the [review receipt](fable_review_receipt.json).
The clarification withdraws several initial overclaims. Counted moments do not
establish identification; bound hits do not prove target incompatibility; the
no-type candidate cannot be evidence that its nonexistent types caused the
wealth dispersion; a finite rental wedge only approaches a hard cap in a limit.
The plateau doses are four and twenty annual earnings amounts. Estate
liquidation and deterministic tenure remain substantive choices.

The lead does **not** accept importing BGM's transitory-variance correction into
the project's different income concept, or arbitrary economic success cutoffs
from the clarification. These are recorded suggestions, not experimental facts
or adopted restrictions. The next candidate intervention is the rental-size
price slope, with the owner premium retained and intercept zero; a cap-10,
zero-slope control must separate added rental support from price effects.
Implementation and launch remain conditional on source/accounting review.

Three bounded deliverables are being implemented by Luna with separate file
ownership: saved-array credit/debt/support diagnosis; exact-moment and simulated
annual-to-four-year earnings aggregation diagnosis; and early-ACS housing age
profiles. They do not themselves solve a new household model or change targets.
The earnings exercise measures the existing aggregation approximation and does
not silently install a new process. The empirical exercise uses the actual
active **2005–2006, 42-metro, cap-at-9** source, and separately labels national
diagnostic figures. The older 2012–2023 cache is metro-restricted and is not an
appropriate national or active-period replacement.

Current target provenance is the early-housing builder and candidate packet in
`output/model/e5f_matched_pf_20260909a/design_research/housing/`, with the exact
working target contract in `initial_calibration_contract/working_weights.csv`
under that experiment. The historical August 17 room receipt does not contain
the current capped-room values and must not be substituted.

The existing `earnings-refit-follow-up` heartbeat has been repurposed for this
goal. It checks every 30 minutes, gives a concise author-facing check-in about
every three hours between 08:00 and 22:00 Eastern, and reports substantive
findings, failures or needed decisions promptly. Otherwise unchanged polling
stays quiet. Last author-facing start update: approximately 11:22 Eastern,
September 20. End the follow-up after the final reviewed packet is delivered
and every submitted job is terminal or explicitly blocked.

At 11:24 Eastern, Torch authentication succeeded and the user queue was empty.
Three bounded Luna tasks cover credit-code grounding, empirical-profile
inventory, and deterministic collection of Fable's response. **No new model
job has yet been submitted.** The completed September 19 overnight jobs are
historical evidence; their failed/cancelled jobs must not be revived.

Keep task-owned changes separate from the repository's substantial unrelated
dirty work. The active decision ledger and author manuscript are not part of
this editing scope. Update this file and `CALIBRATION_STATUS.md` with reviewed
launches/results; commit and push only task-owned changes.
