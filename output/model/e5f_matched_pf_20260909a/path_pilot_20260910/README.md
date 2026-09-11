# Three-hour preference-path pilot

**Completed:** see [final results and full tables](RESULTS.md). All numerical paths
finished; the original baseline file-check failure is diagnosed and preserved.
The monitor is paused. The dated progress notes below are historical.

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


## Current handoff: cluster execution verified at 23:27 UTC

Corrected source **e399c90e** is committed and pushed on
`codex/matched-perfect-foresight`. The immutable remote snapshot is
`/scratch/td2248/projects/Fertility_Spring26_preference_shape_20260910b`.
The source and numerical formula retain the original 60-test verification.

All three six-date smoke cases, array **17347123**, completed with exit code zero
in 3m12–3m25. Their `afterok` dependency released main array **17347124** at
23:25:35 UTC. All three main cases are independently confirmed RUNNING, with
fresh `historical_backward_and_forward` heartbeats. Each main case verifies all
three smoke summaries, gates and artifact hashes before computing. At 23:58 UTC the lead independently verified all three local smoke receipts:
508 source pins per case, all 15 saved-artifact hashes, all 36 numerical gates,
18 dated household budgets, every fit-row calculation and unchanged parameters.
These are six-date plumbing tests; they do not pass an equilibrium market test.

The remaining main budget is now **one 100-date mapping per case**: the linear
case exactly replays the baseline, while the earlier/later cases hold the
inherited parent price path fixed. `E5F_PILOT_CONDITIONAL_ONLY=1` makes that
restriction explicit. This supersedes the optional three-mapping root budget
above. At the previously observed 42–56 minutes per mapping, the three concurrent
cases require about 2.1–2.8 core-hours and 42–56 minutes of running wall time,
plus setup. This is a conditional diagnostic, not a newly solved equilibrium or
calibration. Do not inherit the parent's market certificate for changed paths.

The fixed review deadline remains **2026-09-11 00:30 UTC**. The case builder
sets its internal budget from that absolute deadline, and the cluster allocation
has a separate 70-minute limit. Preserve unfinished checkpoints at the review
cutoff; do not automatically start another round. The jobs continue independently
of laptop sleep. The app's 15-minute follow-up now names the corrected jobs and
may resume when the app is available; remote progress does not rely on it.

## Earlier launch failure and correction

Original source 47435e59 was staged under the separate `...20260910a` snapshot.
Original smoke array 17340342 timed out after roughly 20 minutes per case.
The submission inherited `NUMBA_DISABLE_JIT=1` from the local-test environment,
so the model ran without compilation. This was a launch mistake, not evidence
that the model's numerical method needs to change. The corrected shell explicitly
exports `NUMBA_DISABLE_JIT=0`, checks Numba's runtime setting and prints the result
before the model starts; that check passed in the new cluster logs. A compiled
Numba probe, Python compilation and shell syntax checks also passed.

During the original connection interruption, a detached submission raced with
the foreground submission, creating duplicate smoke 17340490 and main 17340491.
Those jobs are cancelled. **Do not repair or relaunch them.** The corrected
snapshot has one smoke array and one main array, submitted once and verified.
Old outputs remain intact for diagnosis. No targets, weights, scientific gates,
model equations or production benchmark were changed by the launch repair.

## Local verification and historical comparison

The main paths were still progressing at 00:00 UTC: cases 0 and 1 had completed
57 of 100 forward dates, and case 2 had completed 68. No completed main output
is collected yet. All three runtime compilation checks passed in their logs.

The existing, exactly reproduced baseline path has now been matched to the
four empirical birth-count blocks. Normalizing 2008–2011 births to 100 gives:

| Birth years | Data | Existing baseline |
|---|---:|---:|
| 2008–2011 | 100.00 | 100.00 |
| 2012–2015 | 97.06 | 88.79 |
| 2016–2019 | 93.93 | 82.75 |
| 2020–2023 | 89.04 | 79.07 |

The existing baseline therefore has a 20.93% decline versus 10.96% in the
observed aggregate birth count. This is a new diagnostic readout of an existing
solution, not a result of the running alternatives or a female-TFR comparison.
Own-first-block normalization discards the initial level; historical household
counts/ages are externally conditioned, and national birth data versus the
model's housing geography remains an approximation. Additional births in the
3+ parity bin are imputed at entry into that bin. All raw levels and denominators
are retained in `computation/baseline_birth_comparison.json`.

[Full inherited fit table](../design_research/computation/final_replay/evaluation_003/target_fit.csv)
and [all parameter estimates, bounds and current fixed roles](computation/all_inherited_parameters.csv)
remain available. No target or parameter was changed by this diagnostic.

After copying completed remote `output/main/` cases into `computation/main/`,
regenerate verification and birth comparisons from the project root with:

```sh
python3 tmp/e5f_matched_pf/code/model/tools/collect_e5f_matched_pf_preference_pilot.py \
  --pilot-root output/model/e5f_matched_pf_20260909a/path_pilot_20260910 \
  --source-root tmp/e5f_matched_pf \
  --parent-evaluation output/model/e5f_matched_pf_20260909a/design_research/computation/final_replay/evaluation_003
```

The collector performs no model solve. It verifies each completed case before
writing `computation/verified_receipts.json`, full main fit and birth-comparison
CSVs, the common inherited-parameter table, and each candidate's detailed birth
comparison. Missing completed summaries remain explicitly pending. The unchanged
main case reproduces fit, parameter and measurement tables byte for byte. The
collector separately diagnoses exactly four relocated ACS path labels in the
transition table and requires every other field to remain identical; the original
failed file check remains in the record.
