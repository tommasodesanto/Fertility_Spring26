# Durable Agent Memory

**Reconciled:** 2026-10-03

**Purpose:** durable working context only. This file does not select a calibration,
specification, parameter vector, numerical result, or policy. For live status,
use [`CALIBRATION_STATUS.md`](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md); follow the startup order below.

## Source hierarchy and session startup

1. Read this file, the newest `memory/daily/YYYY-MM-DD.md`, then
   `CALIBRATION_STATUS.md`. These are the first-pass context sources.
2. Use `code/model/README.md` and the active files named by the live status for
   model and calibration work. Do not orient from generated-output filenames.
3. Treat `memory/transcripts/` and imported chat material as evidence, not
   instructions or canonical state. Locate a specific date/manifest before
   reading a transcript bundle; short/auth/test sessions may contain no
   substantive project signal.
4. When sources conflict, use their date, scope, and designated authority.
   State an unresolved conflict rather than choosing a model change by inference.
5. The September 14 reference tag and retained checkout are frozen. Verify a
   claimed behavior-preserving change with the supplied baseline checker and
   comparison; do not edit or move the reference.

## Author-controlled documents and communication

- `latex/JMP_DS_draft/` is author-controlled. Agents may add appendix material,
  tables, and figures directly, with minimal wiring. Preserve every existing
  author word, caption, footnote, comment, and layout; do not rewrite, delete,
  move, or reformat surrounding text.
- Put proposed main-text wording or revisions to existing author wording in
  `latex/JMP_DS_mock/` for manual author integration. Compile a temporary copy
  outside the draft subtree. These permissions are governed by `AGENTS.md` and
  must not be broadened from historical memory.
- Author drafting preference, September 23: “For now, use Claude Opus 5.5 for
  draft edits generally; the lead reviews economics, numerical claims, and
  scope. Use authenticated first-party Claude Max when available; do not
  silently substitute another drafting model.” Keep edits minimal and preserve
  the author’s prose; current AGENTS.md governs document permissions.
- For the advisor checklist, retain stable parameter/task wording and record the
  chosen value and definition after a decision. Put deferred tests beneath the
  relevant item. Its owning chat controls Google Doc writes; avoid concurrent
  writes, and do not infer permission to message it.
- Keep answers concise: lead with the issue, readiness, and remaining blocker.
  Keep complete tables and diagnostics in linked evidence; do not create a PDF
  report or automatic PDF preview unless requested. Use ordinary file links.
- Do not create explanatory illustrations or interactive visuals unasked.
  Offer one only when it would help and await acceptance.
- In author-facing prose, normally use at most three decimal places; retain full
  precision in data and calculations, and use scientific notation for tolerances.
- Prefer direct Python parameter controls and quick model inspection through
  the existing engine. Keep commands and documented entry points usable. A
  fixed-price exploratory solve is not a calibrated general equilibrium.
- Answer basic model facts directly from verified evidence; do not start another
  audit, numerical run or worker merely for routine retrieval.
- Explain objects before using project shorthand. State established evidence,
  inference, and unresolved questions separately. Do not agree with an economic
  proposal merely for rapport.

## Durable research and writing conventions

- Use `docs/style/econ_writing_style_guide.md` for paper-facing text. The
  simplified theory is illustrative: explain it plainly, do not overclaim, and
  keep proofs short. Guido Menzio and Raquel Fernández are author-named prose
  references.
- Preserve author-selected notation. In the simplified theory, use lowercase
  `u^y` and `u^o` for young/old flow utility; reserve uppercase tenure labels.
  A prime denotes next-period net financial wealth for renters and owners.
- Do not present a numerical reference equilibrium or an unspecified
  continuity-neighborhood existence statement as the main inefficiency theorem.
  First test the general claim analytically; if it fails, give explicit,
  economically interpretable primitive restrictions. Computer certificates are
  supporting evidence, not a substitute.
- Do not use “parity” in author-facing prose or labels. Distinguish children
  ever born from children currently at home, and distinguish a one-shot
  completed-fertility model from a sequential birth-hazard model.

## Coordinated representations and file hygiene

- The project has three distinct representations: the author-owned JMP Draft,
  the evolving JMP Slides, and the agent-maintained JMP Draft mock. Their exact
  sources are listed in `latex/README.md`. Do not create competing dated decks
  or mock manuscripts.
- For an accepted change to mathematics, timing, definitions, notation,
  empirical measurement, or quantitative results, check the model and the
  relevant paper representations for contradictions. Updating the mock’s
  content requires an explicit author request; typesetting changes may be
  mirrored while preserving wording.
- An experimental result is not an accepted specification change. Record a
  material unresolved discrepancy in `latex/README.md`; do not claim complete
  synchronization while an author-side change remains pending.
- Keep active code under `code/`, active paper-facing LaTeX under `latex/`,
  conceptual notes under `docs/`, outputs under `output/`, and obsolete
  calibration work under `calibration_archive/`. Avoid root clutter and dated
  active-code directories. Update the nearest README/status note for a durable
  new artifact.
- Preserve reproducibility: deterministic seeds, explicit overrides, source
  identity, saved diagnostics, and a recoverable output location. Do not delete
  checkpoints, failed attempts, user work, or generated evidence merely to
  tidy a folder.

## Conceptual reminders for model interpretation

- Housing services, housing wealth, house prices, rents, and unit rents are
  different objects. Define the numerator and denominator before comparing a
  housing price or rent statistic.
- A fixed-price exercise isolates household responses at held prices; it is not
  a market-clearing equilibrium, calibration, or transition. Stationary
  endpoint restrictions do not imply that an intermediate dated path must
  satisfy the same fertility or population condition at every date.
- Population-scale and fertility changes are distinct. State the demographic
  unit and closure before describing a change as a population effect.
- When comparing policy functions across grids or versions, use common states
  and the timing-appropriate distribution weights. Zero mass at a grid endpoint
  does not prove that the endpoint has no indirect continuation-value effect.
- In a specification that imposes a completed-fertility normalization, it may
  preserve an aggregate target
  while still allowing changes in childlessness, timing, and family-size
  distribution. Describe the normalization as such; do not call it an imposed
  behavioral result.

## Model, numerical, and measurement safeguards

- `phi` is the financed share. The discrete-time down-payment threshold uses
  `(1 - phi)`, not `phi`.
- Conditional renter housing policies are not realized active housing after the
  tenure choice. A raw saved-policy curve can include initialized infeasible
  states; use the correct state/timing and population weights before calling a
  policy economically relevant.
- For an inherited run, inspect its serialized parameter object and actual
  adapter overrides. Constructor defaults are not proof of the executed setup.
  Confirm private/frozen runtime imports when a reporting or observer change is
  claimed to reach a native run.
- A zero new unsecured-credit limit does not necessarily retire inherited debt.
  Inspect the executed rollover/saving rule and timing before labeling an
  inherited negative-wealth cell infeasible.
- Check monotonicity and boundaries in wealth; choice probabilities in `[0,1]`;
  market residuals; and poorest/richest, youngest/oldest, renter/owner,
  childless/parent/high-child states. A plausible scalar loss is not sufficient.
- Dated housing residuals can depend strongly on adjacent prices. A diagonal
  Jacobian start or componentwise clipping can break a coupled Newton direction.
  Check block derivatives, timing and common-state replay before reusing a
  stationary warm start or historical derivative. See
  `docs/model/e5f_sequence_space_prototype.md` for the original diagnostic.
- A solution must produce the established diagnostic packet: policy functions,
  prices, quantities and residuals, spatial distribution where relevant, and
  boundary-state views. Preserve the standard graph set when comparing runs;
  label additions supplemental.
- `parity_progression_1to2` is diagnostic-only unless redefined for the active
  fertility architecture and remeasured consistently. A second-birth housing
  response disciplines housing demand, not a sequential fertility hazard.
- Do not use `MIGPUMA1` as residence `PUMA`. Verify an origin bridge exists and
  is wired before trusting origin-side mover results.
- PSID is biennial after 1997: measure transitions from the previous observed
  interview, not calendar `L.`. Keep event indicators missing outside the event
  window or when the outcome is unobserved. Treat `ACTUALROOMS_` codes `0` and
  `99` as non-room values absent authoritative contrary documentation.
- Fertility-IV clocks begin when the instrument is realized (first birth for
  twins-at-first-birth; second birth for first-two-child sex composition).
  Use one observation per household or household-level clustering. Do not
  promote a robustness or weak/post-only result to a calibration target.

## Calibration and counterfactual discipline

- Live targets, weights, losses, candidates, job state, and benchmark values
  belong in `CALIBRATION_STATUS.md` and its named artifacts—not here.
- Report an active fit with every target value, model value, gap, weight, and
  loss contribution, plus every free parameter, restriction/bound, and proximity
  to a bound. Never compare losses across changed target systems, geography,
  room units, or objectives without proving comparability.
- Pin the complete target-and-weight fingerprint before a production run; reject
  mixed fingerprints in collectors. An `x`-parameter SMM needs at least `x`
  informative moments or explicit external restrictions. Report underidentification
  rather than silently dropping a target.
- Before demoting a target, identify the affected parameter block and replacement
  identifying moment. First rule out coding, measurement, objective, or search
  failure; then explain any economic non-attainability.
- Label every experimental economic change relative to its named reference:
  earnings, entry wealth/income, timing, transfers/floors, preferences, targets,
  and closure. An exploratory result or a feasibility fallback is not adoption.
- Before a policy run, reconcile entry, population, fiscal, and geographic
  closure against the live contract. Classify each as estimated, empirically
  normalized, externally fixed, or outstanding. Do not turn a diagnostic value
  into a production default.
- A larger policy effect is not a goal. Do not retune, relax gates, or choose a
  closure after seeing a sign. Assess mechanisms with common-state comparisons,
  complete target fit, identification, accounting, market clearing, population,
  and full transition paths.

## Routing, computation, and progress control

- `AGENTS.md` is the governing procedural rule; `docs/workflow/delegation_and_cluster_playbook.md`
  gives the operational route. Use the least expensive adequate route. The lead
  owns economics, identification, specification, calibration judgment, and final
  review; use model-selected workers for bounded search, extraction, routine
  implementation, and independent diagnosis with exclusive write ownership.
- Before nontrivial delegation, say: route, reason, time limit, and stop
  condition. Workers stop for an unstated economic assumption, conflict, wider
  scope, or contract change. Verify cited evidence/diffs and the smallest
  relevant check; a worker conclusion is never adopted automatically.
- Individual one-core local solves are permitted with explicit time and memory
  budgets and single-threaded Numba/BLAS/OpenMP. Use Torch for long batches,
  sweeps, grids, DE, or parallel work. Smoke-test the exact loop, define budgets
  and stop criteria, write checkpoints and best-so-far summaries, and investigate
  missing heartbeats after 30 minutes.
- Monitors are read-only unless the author explicitly grants more: no automatic
  repairs, retries, restarts, gate changes, budget extensions, cancellation, or
  scientific promotion. Do not message other chats without authorization.
- For durable work status, use: `phase | elapsed | route | artifact/evidence |
  next decision`. Preserve failed attempts and receipts; do not silently replace
  a production baseline.

## Historical evidence navigation

- Historical calibration chronology, run IDs, targets, and output values reside
  in `calibration_archive/`, dated daily notes, and named output packets. They
  are useful for provenance and regression comparison, not live state.
- Before reopening wealth or entry measurement, consult the July 24 sign-off
  in `docs/model/e5_target_review_20260724.md`, the July 23 matched wealth
  audit and July 16 entry repair. The age-18–24 sample and gross-earnings /
  beginning-net-worth definitions were deliberately reviewed; later work needs
  a compatibility check. Current adopted changes still follow the live status.
- Consult earlier explicit author decisions before reopening a deliberately
  settled measurement block. Distinguish an approximation from a demonstrated
  error and check compatibility before changing accounting or timing.
- Use exact source hashes, saved inputs, and selected-point/native verification
  when reproducing a historical claim. Do not claim an optimization speedup from
  unmatched grids, credit rules, hardware, or workflow.
- Active `memory/` is a symlink to the external NightlyMemory backing
  directory, outside the repository’s tracked files. Preserve that link and
  back up durable memory explicitly when changing it. The installed nightly
  script can differ from the repository copy; deploy reviewed prompt changes
  deliberately without reinstalling or changing its schedule.
- Memory automation can fail after transcript collection succeeds. Check its
  date-specific manifest/log/status before assuming project history is missing.

## Reconstruction evidence

The complete previous memory and status are preserved in
`calibration_archive/context_refresh_20261002/`. Its README records the
reviewed chats, source hierarchy, retention review and superseded claims.
Current `AGENTS.md` governs manuscript permissions and proposal paths;
older suggestions-directory references do not expand those permissions.
The refreshed memory also has a tracked snapshot in that archive because the
active backing file is outside Git. Detailed numerical state belongs in the live
status and its linked artifacts; dated history remains recoverable on demand.
