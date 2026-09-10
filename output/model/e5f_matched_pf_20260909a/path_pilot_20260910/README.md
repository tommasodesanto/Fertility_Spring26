# Three-hour preference-path pilot

Authorized by Tommaso on 10 September 2026; review window ends at
**2026-09-11 00:30 UTC (20:30 EDT on September 10)**.

The baseline keeps the fertility-preference intercept constant after 2023.
Households know the complete path from 2007. Initial stationary fertility retains
the agreed 2.1 normalization. This pilot does not adopt a new production target
system or claim a completed initial-era calibration.

## Concrete experiment

Compare three fixed historical timing profiles, using the same inherited
structural parameters, initial state, terminal state and housing-supply anchor:

\[
\psi_t=\psi_0+\Delta\psi\{x_t+a x_t(1-x_t)\},\qquad
x_t=\min\{1,\max\{0,(t-2007)/16\}\},\quad
 a\in\{-0.5,0,0.5\}.
\]

The zero case reproduces the existing linear preference path; positive a moves
its decline earlier, negative a later. All three have identical endpoints and
remain constant after2023. They are fixed diagnostic cases, not new SMM estimates.
Changing the amplitude would require a newly solved terminal equilibrium and
is outside this first comparison.

The primary supplemental readout is aggregate births in2008–11,2012–15,
2016–19 and2020–23, compared with decisions2007,2011,2015 and2019. The
2023 decision supplies2024–27 births and is reported as continuation only.
The empirical source/aggregation packet is [fertility_data](fertility_data/).
This avoids falsely labeling household-based model rates as female TFR. Raw
birth levels, the top-bin addition and changes in each candidate's first block
remain visible beside normalized birth-count indices. The reviewed twelve-row
objective and every inherited parameter are retained as separate full tables.

## Budgets and gates

- Local preparation: three bounded25-minute supporting passes plus lead review.
- Exact candidate-loop smoke: all three timing cases on six dates, with the
  unchanged accounting, feasibility, source and target-fingerprint gates.
  This smoke is not a short-horizon equilibrium claim.
- Main comparison, conditional on a passed smoke and cluster access: at most
  three100-date candidates, one CPU/16GB each, up to three complete path
  evaluations per case including final reproduction if a root is attempted.
  Observed complete100-date mappings take about42–56minutes. Nine mappings
  imply6.3–8.4core-hours, roughly2.1–2.8hours with three concurrent cases,
  excluding queue/setup; therefore complete new roots before review are not
  guaranteed. A completed conditional comparison remains explicitly conditional.
- No launch exceeding the remaining review window; retain checkpoints for any
  incomplete case. No automatic expansion, amplitude search, new calibration,
  policy run or horizon extension in this pilot.
- Each case writes a heartbeat, latest completed evaluation and best-so-far
  record. Investigate30minutes without progress. Preserve the2e-4 market gate,
  all accounting/feasibility tolerances and exact reproduction checks.
- A converged parent's prices are only a numerical starting guess for a changed
  preference path; its market certification is not inherited.

## Initial access check (superseded by submission handoff below)

At21:31UTC the Torch connection reached the host but SSH authentication was
rejected. Tommaso has been asked to refresh access; local preparation continues.
No new cluster job has been submitted. The numerical source is isolated in
`tmp/e5f_matched_pf`, starting from96a41873; other project work is untouched.
The existing100-date equilibrium also has an outstanding terminal-rent and
horizon check, so this pilot cannot be promoted as a certified historical fit.

In parallel, [initial_fit](initial_fit/) maps every proposed initial restriction
to its model observer and prepares a readout from the existing normalized-old
checkpoint, without another equilibrium solve or invented weights.


## Submission handoff: connection interruption

Source47435e59 is committed and pushed on the isolated matched-PF branch.
The lead reviewed the optional preference hook and runner against the formula;
60focused tests passed locally, including independent birth-observer checks.
Source snapshot: `/scratch/td2248/projects/Fertility_Spring26_preference_shape_20260910a`.

Smoke array **17340342**, three cases, was confirmed RUNNING oncs618. A detached
submission coordinator was created so laptop closure would not interrupt setup.
The original foreground submission completed concurrently, creating a duplicate
submission risk. The coordinator's duplicate smoke was cancelled by reading its
saved job ID and preserving17340342. A main array **17340491** was already saved
by that coordinator. **Its dependency may still refer to the cancelled duplicate.**
An attempted repair/requeue was interrupted by loss of the SSH connection; its
completion is NOT verified. Do not say the main comparison is running or has
passed. At22:58UTC a fresh SSH attempt reached the host but authentication was
again rejected. No current smoke-completion receipt was collected.

On restored access, first inspect `smoke_job_id.txt`, `main_job_id.txt`,
`replaced_main_job_id.txt` if present, `logs/submission.log`, squeue/sacct and the
three output/smoke/case_N summary/contract files. Verify all three exact-loop
smokes and their artifact hashes. If17340491 is pending, repair its dependency
toafterok:17340342; if cancelled, submit exactly one replacement only after
checking that no replacement already exists. Record the actual ID. Preserve
source and all gates. Main jobs themselves verify all three smoke receipts
before computing, and the case builder refuses duplicate outputs/contracts.
The absolute review-window deadline is2026-09-11 00:30UTC; if it is already past,
collect the completed smoke and prepare the next decision rather than launch
another unapproved round. Even before the deadline, disclose incomplete roots.

The existing thread follow-up was reactivated at23:02UTC to inspect this handoff every15minutes, notify only meaningful changes, and stop at the00:30UTC review deadline. The first scheduling attempt timed out; the retry returned an explicit ACTIVE confirmation.
