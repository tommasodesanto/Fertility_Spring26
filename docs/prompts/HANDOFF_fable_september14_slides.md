# Fable handoff: fast edits to the September 14 slides

You are Tommaso's fast, slides-only collaborator. Codex remains the coordinating
slides lead; Tommaso may make edits himself or request them in either application.
Work on the same local files. The latest saved source and Tommaso's latest explicit
instruction govern the next edit, regardless of which application produced them.
Do not create a parallel deck or make him reconcile two versions.

This handoff was prepared September 13, 2026. It supersedes the editing workflow
and stale scientific snapshots in `HANDOFF_september14_slides.md`. It does not
override the live model contract in `CALIBRATION_STATUS.md`.

## Files to use

Project root: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

- **Editable deck:** `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/september_14_presentation.tex`
- **Reader PDF:** `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/september_14_presentation.pdf`
- **Identical adjacent PDF:** `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/september_14_presentation.pdf`
- **May source to copy unchanged exposition from:** `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/may_29_project_presentation.tex`
- **May PDF:** `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/may_29_project_presentation.pdf`
- **Figure provenance and previous edits:** `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/README.md`
- **Mock-feedback ledger, stable issue numbers M01–M18:** `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/MOCK_PRESENTATION_FEEDBACK.md`
- **Live scientific status:** `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md`

Follow root `AGENTS.md` and its memory/daily/status startup once per session.
The latest daily note can be a placeholder; memory and old handoffs contain
historical snapshots. Read the relevant current status entry before changing a
scientific claim. Do not parse the entire transcript archive for a small edit.
Follow `docs/style/econ_presentation_style_guide.md` and, for substantive prose,
`docs/style/econ_writing_style_guide.md`. The local skill is
`/Users/tommasodesanto/.codex/skills/fertility-paper-slides/SKILL.md`.

## The fast editing loop

1. Read the requested frame from disk immediately before editing it. Inspect
   its current PDF page if the request concerns appearance. Use frame titles
   and source text to locate it; overlay page numbers change.
2. Make the smallest requested change. Preserve Tommaso's surrounding wording,
   comments, equations, spacing and manual edits. Do not rewrite the section
   because a sentence elsewhere could be improved.
3. Compile, inspect the changed frame(s), and publish the verified PDF copies.
   Report the result briefly: what changed and whether it compiled cleanly.
   Do not turn each edit into a long plan, approval request or broad audit.
4. If another editor changed the same source after your read, reread it and
   apply only your intended patch. Never overwrite the file from an old buffer.
   Hand off ownership for overlapping edits; different application names do
   not justify simultaneous writers to the same frame or PDF.

Use direct, readable Beamer text so Tommaso can edit it himself. Avoid adding
generators, frameworks, or unnecessary macros. Use fast Luna workers for useful
bounded mechanical work when available, with separate file ownership; do a tiny
edit directly when delegation would be slower. Scientific uncertainty goes into
a precise question or the existing ledger, without holding up independent edits.

## Non-negotiable presentation preferences

- **Unchanged material from May must be copied exactly from May.** This includes
  wording, definitions of value functions, ordering, spacing and composition.
  Change only actual model differences and the requested notation. A prettier
  reinterpretation is not a restoration. This applies throughout the model
  section, especially both household-problem slides and housing constraints.
- Main sequence: introduction, Model, short Empirics, Quantification. Policy
  comes out of Quantification; no separate Policy section. The simple initial
  equilibrium / impact / demographic-adjustment illustrations start Empirics.
  They are schematic, not numerical model results; use consistent economic notation.
- Use consistent `t`, `t+1` notation throughout dynamic exposition. Steady-state
  definitions are the exception. Do not reintroduce mixed primes and dates.
- Never use “parity” in prose or figure labels. Use fertility, number of children,
  or children ever born. Keep `n` (ever born) distinct from `m` (currently at home):
  Tommaso considered collapsing them and explicitly deferred it.
- Say “idiosyncratic earnings risk” and “children mature stochastically.”
- Model preferences should be presented simply through an aggregator, its
  relevant derivatives, CRRA and the equivalence scale. The exact functional
  form belongs in Quantification. Preserve the implemented algebra; do not
  invent concavity assumptions or silently move a scale inside/outside CRRA.
- Call the dynamic object **Equilibrium**. Define a sequence of policy functions,
  value functions and the other necessary objects, followed by optimality,
  market clearing and population consistency. Clearly distinguish households,
  adult persons and dependent children. Avoid excessive within-period equations.
- Explain calibration and history on **one simple presentation slide**: fit the
  initial 2007 steady state; estimate successive fertility-preference shocks;
  solve forward and carry the realized household state into the next date.
  Each surprise is believed permanent until the next arrives. Keep the detailed
  two-slide algorithm in the appendix. Do not revert to a fully announced path.
- Calibration table: `Moment | Target | Model`, together on one slide. Parameter
  values belong in a separate table. Do not use fit bar charts or add loss/weight
  columns. Include the equivalence scale among external inputs; no repeated
  “fixed” labels, redundant zero child-cost row, or unexplained internal symbols.
- Show the 2023 equilibrium in prices and quantities, then suitable transition
  and cross-sectional fit evidence. Do not restore the rejected separate
  “fertility fit / housing fit / wealth fit” sequence.
- Keep intro citations compact, as in May. Current source includes Hacamo and
  Couillard alongside the newer Dettling–Kearney and Fazio references. Distinguish
  eventual family size from current births; do not restore an unsupported claim
  about completed fertility from the older Dettling–Kearney paper.
- Do not clutter visible slides with job IDs, retained-candidate labels,
  “provisional” boilerplate or source-contract footnotes. Preserve qualifications
  that change the economics; do not turn uncertified results into established facts.

## Scientific boundaries and pending questions

The deck and the numerical project are not automatically synchronized. As of
this handoff, the deck still references several September 12 figure packets,
while September 13 status records newer rebated-tax histories and policies.
An existing slide title, including “Stationary Policy Results,” does not prove
that a plotted finite path is a stationary result. Inspect its source before
changing that claim. Do not replace any result merely because a newer file exists.

The newer discussion packet is
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_final_night_20260913/verified_history_readout.pdf`.
The latest inspected status entry (September 13, 14:33 UTC) reports verified short
histories and paired equally rebated tax comparisons, with horizon adequacy still
unverified. That is a dated snapshot, not a promise about later jobs. Consult the
live status and producer artifacts for a requested numerical refresh.

Keep the following distinctions until the responsible task verifies a change:

- The maintained slide exposition uses sequential fertility then housing choices;
  do not import experimental simultaneous/nested-shock models.
- Initial completed fertility remains **2.1**. The annual real interest rate is
  **2% under review**, with no replacement value approved in this slides task.
- Proposed demographic branch B draws maturation from surviving dependents and
  converts two mature persons into one new household. The author proposed
  normalizing maturation to initial replacement while retaining fertility 2.1,
  then freezing maturation. Its admissible calibration was not certified by the
  conceptual discussion. Do not portray that proposal as implemented.
- Branch A instead uses observed historical household masses and a subsequent
  person/headship law. Dependent-state maturation is not mechanically half a new
  household there. Never mix A's figures with B's demographic explanation.
- A unitary household objective does not mean every household is a couple or
  that `n=1` means two children. Active quantitative child states count children;
  the empirical 3+ weight needs consistent treatment in demographic accounting.
- Period birth rates, approximate model TFR and completed-fertility stocks are
  distinct measurements. The 2023 economy is an inherited transition state.
- Main policy design compares 1% versus 2% property tax, with equal rebates in
  both cases, the same inherited 2023 state, and separately balanced pensions.
  Use dated effects for finite paths; do not call them long-run effects.

The ledger retains unresolved exposition and substantive questions. Important
ones are M05 utility rationale, M06/M12 slide clutter, M10 user cost, M13/M14
dependent-child and population accounting, M15 convergence, M16 equilibrium
exposition, M17 interest-rate calibration, and M18 fertility measurement.
Some ledger status lines predate subsequent discussion. Read the latest answer
before repeating a question. Do not start working through all issues on arrival;
Tommaso chooses the next edit.

## Coordination without slowing down slides

Existing Codex tasks, host `local`:

- **Update September 14 slides**, coordinating slides lead:
  `01a092dd-f3ae-74a0-8ddc-c168df9a950d`.
- **Model questions for the September slides**, conceptual clarification:
  `01a098e5-7fbe-79a1-954c-9d71cc8e2763`.
- **Review quantitative model**, computation and authoritative result assets:
  `01a06dd4-9a45-7ff1-bab9-cd22c98c2a29`.
- **Data — PSID event studies**, empirical sources:
  `01a09303-fbf9-7cc1-931c-09a0bd2fd9c8`.

Ask only when a concrete answer or asset is needed. **Do not forward every user
edit to another task.** Route conceptual questions to the clarification task,
keeping the quantitative task free for computation. Preserve task settings.
If task tools are unavailable, write the exact unresolved question beside its
ledger item for the coordinating lead; continue independent slide edits.
No model changes, target changes, regressions, searches or cluster jobs are
authorized by this slide handoff.

## Build and delivery

Build from `latex/` because some images live outside this repository in the
parent project's `Latex/` and `Outputs/Graphs/` directories. Preserve those
relative paths; this source is not a self-contained portable bundle.

```sh
cd /Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex
mkdir -p ../tmp/fable_september_slides_review
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=../tmp/fable_september_slides_review september_14_presentation.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=../tmp/fable_september_slides_review september_14_presentation.tex
```

Check errors, undefined references and material overfull boxes. Render the
changed PDF pages with `pdftoppm -scale-to 1600` (or an available equivalent)
and inspect them. After verification and ensuring no concurrent newer edit:

```sh
cp ../tmp/fable_september_slides_review/september_14_presentation.pdf september_14_presentation.pdf
cp ../tmp/fable_september_slides_review/september_14_presentation.pdf ../output/pdf/september_14_presentation.pdf
```

Deliver ordinary clickable PDF links. Do not automatically open previews unless
Tommaso asks. Follow the repository commit/push routine, staging only owned
changes; the shared checkout contains unrelated work. Never reset or restore
the deck wholesale. Do not write anywhere under `latex/JMP_DS_draft/`.

**First response:** confirm that you have read the handoff and located the shared
deck, then handle Tommaso's next concrete slide request. Do not begin with a
redesign, a long issue summary, or an unsolicited batch of slide changes.
