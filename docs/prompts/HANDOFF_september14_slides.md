# September 14 presentation: slides-only collaboration

Work with Tommaso only on the presentation: its argument, sequencing, wording,
equations, figures, tables, speaker notes, and the compiled PDF. Tommaso wants a
separate place to discuss slides without managing the technical project again.
The presentation is on September 14, 2026. Be concise, concrete, and economical
with tokens. Do not begin with a broad audit or redesign the entire deck.

## Workspace and startup

Use the existing project folder:
`/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`

Follow its `AGENTS.md`, including the mandatory memory/daily/status startup and
`git status -sb`. Memory contains older snapshots; `CALIBRATION_STATUS.md` owns
the live quantitative status. Read the `fertility-paper-slides` skill at
`/Users/tommasodesanto/.codex/skills/fertility-paper-slides/SKILL.md`, then both
project writing guides before drafting presentation prose:

- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/style/econ_presentation_style_guide.md`
- `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/style/econ_writing_style_guide.md`

## The single working deck

- Source: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/september_14_presentation.tex`
- Reader PDF: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/pdf/september_14_presentation.pdf`
- Adjacent build copy: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/september_14_presentation.pdf`
- Scope, figure provenance, archive and build instructions: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/README.md`
- Visual reference: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/may_29_project_presentation.tex`
- Existing reading notes: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/docs/model/model_reading_notes.txt`

The September 10 refocus retained the quantitative-model/empirical/policy
structure and parked the separate simplified-model exposition. Preserve that
decision unless Tommaso asks to reopen it. Do not create a competing working
deck. The prior task **Update September presentation slides**
(`01a08d25-2b8a-7823-81d9-adfeeedbb56f`, host `local`) can be read for background;
do not restart it as a second simultaneous deck editor.

## Technical coordination: do not make Tommaso relay everything

Tommaso authorizes you to bring technical questions to these existing tasks or
to bounded specialist subagents. You own presentation decisions; the technical
tasks own scientific verification and computation.

1. **Review quantitative model** — this quantitative/calibration task.
   ID `01a06dd4-9a45-7ff1-bab9-cd22c98c2a29`, host `local`.
   Ask it about implemented utility/budgets/timing, calibration and target
   definitions, dated fertility measurement, pensions, numerical validity,
   policy results, and the exact files/numbers approved for a slide.
2. **Review two-period OLG model** — analytical/theory task.
   ID `01a06dc4-5f03-7351-b7df-2b0e80a427e1`, host `local`.
   Ask it about propositions, welfare criteria, assumptions, proofs and the
   relationship between the analytical argument and the full lifecycle model.
   Its remit has recently included the full-model welfare argument despite its
   older title; do not assume its task title determines the current theory.

Where available, use task read/message tools directly. Send a bounded request:
the frame and proposed claim/equation, the precise uncertainty, and a request for
a concise verdict with an authoritative file/line or result artifact. Read
existing context first to avoid repeatedly asking an answered question. Preserve
the destination task's model/settings. Do not spawn a duplicate quantitative
lead. If cross-task tools are unavailable, put the ready-to-send question in
your handoff notes rather than guessing or turning the slides conversation into
a long technical detour.

Use low-cost subagents for genuinely independent reference/figure checks or
bounded technical verification. Keep ownership separate and review their work.
Changes to the numerical model, empirical targets, equilibrium methods, or
cluster searches go through the quantitative task; a slide edit must not silently
change the science. Continue independent presentation work while awaiting a
technical answer.

## Immediate corrections and evidence boundaries

This handoff was prepared on September 11 after the commute. Refresh the live
status before using numerical results. These are internal editing instructions,
not material to paste onto audience-facing slides.

- **Family Space currently has an obsolete equation** at approximately line138:
  `h_1 1{m>0} + h_m m`. The implemented specification is a single parenthood
  requirement, `h_P 1{m>0}`, applying only while dependent children are present.
  The equivalence scale remains; do not remove it or invent a nonhousing floor.
  Update the explanatory sentence and the housing-parameter row in **Empirical
  Discipline** consistently. Verify any further affected notation.
- The general CRRA expression on **Preferences** is algebraically compatible
  with the approved scale. With sigma2, flow utility is
  `-e(m)/(c^alpha [s-h_P 1{m>0}]^(1-alpha)) + psi_t m`, with
  `e(m)=((2+0.7m)/2)^0.7`. Do not rewrite surrounding prose unnecessarily.
- **Earnings and Bequests** mentions retirement income but does not explain
  the repaired Social Security budget. Obtain the precise concise equation from
  the quantitative task: fixed payroll tax, benefits balancing the actual dated
  worker/retiree budget, and those same benefits anticipated by households.
- Ask the quantitative task to verify the financial-return/transaction timing
  on **Budget Constraints**. The theory task has flagged a possible mismatch;
  this handoff does not certify a replacement equation.
- The maintained quantitative experiment starts from an approximate pre-2007
  stationary economy with author-selected model fertility normalization2.1.
  Households learn the preference path in2007; its intercept is flat after2023.
  The2023 economy is a transition state. This is conditional scenario analysis,
  not an identified explanation of the U.S. fertility decline.
- A revised initial calibration and balanced short transition trials exist.
  A fitted historical preference path, horizon validation, and revised policy
  results are not yet established. Do not transplant older policy percentages
  or calibration tables into this specification. The active deck already omits
  obsolete numerical tables; preserve that restraint until verified replacements
  arrive. Keep missing-content decisions in notes, not fabricated slide results.
- Three short trial paths now have reproduced age-specific **household-based**
  fertility diagnostics. These are not automatically empirical female TFR.
  Birth counts, period rates and completed fertility are different objects.
  A2019 model decision covers births in2020–2023; a2023 decision covers2024–2027.
- The intended main policy compares1% and2% property tax with equal rebates in
  both paths, a common inherited2023 state, and a separate balanced payroll
  pension budget. A description of this experiment is not a computed result.
- The latest theory task reports a conditional rental-space welfare result;
  the intended owner-housing result remains unproved. Get its exact statement
  and permission to use it rather than reviving a broad older misallocation claim.
- Existing August transition diagrams are schematic, not simulated paths. Keep
  their interpretation distinct from numerical results. The reviewed PSID
  Sun–Abraham first-birth evidence is an empirical source; do not replace its
  estimator, sample or uncertainty through a cosmetic plotting edit.

Quantitative starting points:

- Live status: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/CALIBRATION_STATUS.md`
- Approved design: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/OVERNIGHT_PLAN.md`
- Initial fit, every target and parameter: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/initial_calibration_contract/extended_refinement/collected_17378993/READOUT.md`
- Return-home comparison, including fixed-beta profiles and completed short-path diagnostics: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_matched_pf_20260909a/current_candidate_transition/return_home_20260911/READOUT.md`

## Presentation workflow

Keep a conventional economics seminar style, close to the May deck. One idea per
frame, short bullets, selective equations, consistent notation, and economic
interpretation. Default calibration table: **Moment | Target | Model**. Keep
weights, losses, job IDs and debugging out of the main deck. Do not turn an
unverified result into a factual slide by hiding its qualification in a footer.
Do not make new illustrations unless Tommaso explicitly requests them.

Begin by inspecting the actual deck and PDF. Give Tommaso a short map of its
sections and the three most urgent slide issues. Then work through the slides
with him; do not spend the first turn re-litigating the technical research plan.

Use LaTeX/Beamer and deliver a PDF, not PowerPoint. Compile twice to the external
build directory named in `latex/README.md`; inspect warnings, links and every
changed frame. Reuse verified data/figures, and record exact evidence sources in
internal notes. Commit only your coherent owned edits, preserving unrelated work.

The entire `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/latex/JMP_DS_draft/`
subtree is strictly read-only, including build products. Do not edit it as part
of preparing slides.
