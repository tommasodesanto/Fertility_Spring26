# LaTeX Workspace

Active documents:

- Slides-only collaboration handoff: `../docs/prompts/HANDOFF_september14_slides.md`.
  It identifies the single working deck, technical task contacts, and the
  September 11 utility/pension updates that still need presentation edits.

- `september_14_presentation.tex`: the single working September 14 seminar
  presentation. The reader PDF is `../output/pdf/september_14_presentation.pdf`;
  `september_14_presentation.pdf` is an identical build copy. The September 10
  refocus follows the May deck's appearance and model/empirical/policy structure.
  It presents the quantitative household environment, sequential choices,
  population accounting, equilibrium, and a policy comparison introduced along
  the same inherited transition. The three original August
  `housing_fertility_stage_{initial,impact,adjustment}.pdf` figures are reused
  unchanged around the equilibrium exposition. They are explicitly schematic:
  their fixed-composition housing curves and replacement-one household units
  are not the quantitative demographic law or a computed transition. Exact
  attribution to the Raquel meeting remains unverified. The separate
  `demographic_adjustment.pdf` is redundant, and `figures/example_misallocation.pdf`
  is excluded because it carries the parked allocation argument.

  The complete pre-refocus source, including all simplified-model main and
  appendix frames and the earlier numerical readout, is preserved verbatim in
  `archive/september_14_before_slide_refocus_20260910.tex` (original source
  SHA256 `90e94b4bd20dbf780a4eacbde49eb2099313a29926e409e8cf068a77faf9e1ea`).
  It is historical source, not a second working deck. All separate theory
  notes and the September 10 planner work remain untouched and parked.

  **Scope of the presentation.** The AHS 2023 tenure-by-bedroom figure is
  reused from May. The PSID rooms profile comes from the corrected
  household-aligned Sun--Abraham event table, via the saved September 5
  measurement-review figure; the pre-birth pattern and pointwise intervals
  remain visible. Only its internal figure title is trimmed in LaTeX. No
  regression, calibration, model solve, or new mechanism figure was generated.
  The old calibration tables and computed transition claims are preserved in
  the archive but omitted from the active deck because they do not establish
  a calibrated history under the latest intended specification.

  **Unresolved choices (kept outside audience-facing slides).** Sequential
  choice is retained and explicitly named for exposition; the comparison with
  simultaneous fertility nests does not itself promote either arm. Initial
  fertility 2.1 is author-selected. Initial-economy-first versus joint
  transition estimation, the final measurement/calendar mappings, and the
  initial/dynamic supply, demographic, and fiscal closures remain governed
  by `../CALIBRATION_STATUS.md`; this edit does not adopt a new target contract.
  The rebated-tax frame defines a conditional policy experiment with one common
  inherited state and supply schedule, not a computed or authorized production
  policy result. Initial-equilibrium conditioning must still be reconciled
  with the historical demographic inputs before presenting a numerical fit.

  Build twice from `latex/`, writing auxiliary files outside the active folder:
  `pdflatex -interaction=nonstopmode -halt-on-error -output-directory=../tmp/september_slides_review september_14_presentation.tex`.
  Copy the verified PDF to the reader path and the adjacent build copy.
  The refocused deck has 29 numbered main frames and six appendix frames
  (42 PDF pages including overlays and dividers). Two final compilation passes
  have no warnings or overfull boxes; all frames were visually inspected and
  all figure paths and the remaining appendix link resolve.

- `JMP_DS_draft/`: author-controlled source for the new job-market-paper draft.
  Its main file is `JMP_DS_draft/JMP_DS_draft.tex`, with separate section and
  appendix files. The entire subtree is read-only for agents after its initial
  creation; proposed text belongs in `JMP_DS_suggestions/` for manual
  copy-and-paste by the author.
- `JMP_DS_suggestions/`: agent-writable staging area for passages, equations,
  or revisions proposed for the author-controlled draft. Nothing placed here
  is part of the manuscript unless the author copies it into the draft.

- `simplified_olg_paper_theory_package.tex` /
  `simplified_olg_paper_theory_package.pdf`: definitive paper-facing theory
  package. Its first five pages are the compact analytical section; the
  appendix gives complete household derivations, positive-steady-state and
  conditional transition results, the intergenerational reallocation proof,
  the exact fertility derivative, and a high-accuracy terminally closed
  mixed-tenure transition approximation. The same section and appendix are
  input into `intergenerational_housing_fertility_paper_draft.tex`. The claim
  ledger is `../docs/model/simplified_olg_theory_claim_ledger.md`.
- `full_quantitative_model_analytical_note.tex`: compact analytical map for the
  quantitative model. It derives the stationary person law, the renter,
  tenure-product, and sequential-fertility margins, the long-run and finite-
  transition implicit systems, and the appropriate money-metric reallocation
  test without introducing a continuous owner-housing choice. The verified
  reader PDF is written to
  `output/pdf/full_quantitative_model_analytical_note.pdf`.
- `archive/simplified_olg_owner_only_development_20260830/`: superseded
  owner-only development note, technical appendix, PDFs, and figure. They are
  retained for provenance; the mixed-tenure package above is the only current
  paper-facing simplified theory.
- `transition_closure_update_presentation.tex` /
  `transition_closure_update_presentation.pdf`: 20-slide advisor update. It
  defines the timing and choices in the two-generation population--housing
  model, presents the initial equilibrium, impact response, and demographic
  adjustment in three separate stages, maps the closure into the existing
  lifecycle model, explains the 2007--2023 dated calibration, and reports the
  no-policy continuations and the absence of a positive closed stationary root
  under the current estimate.
- `aug_07_model_closure_presentation.tex`: reference deck for the production
  paper architecture following the August 15 author decision. Its main model
  is the one-market sequential-fertility `E5F` floor arm with independent child
  maturation. The parameters and policy rows shown there remain provisional;
  architecture promotion is not calibration promotion.
- `model_writeup.tex` / `model_writeup.pdf`: current model writeup.
- `main_note.tex` / `main_note.pdf`: larger paper-style writeup.
- `april_20_project_presentation.tex` / `april_20_project_presentation.pdf`:
  recent slide deck.
- `may_29_project_presentation.tex` / `may_29_project_presentation.pdf`:
  reference presentation deck.
- `intergenerational_housing_fertility_note_slides.tex` /
  `intergenerational_housing_fertility_note_slides.pdf`: expanded academic
  presentation of the full circulated July 2026 paper draft, with separate
  expositions of the simplified analytical model and full quantitative
  lifecycle model, followed by calibration, results, and policy.
- `distributional_empirics_report.tex` /
  `distributional_empirics_report.pdf`: data-vs-model distributional discipline
  report.
- `fertility_slice_diagnosis_report.tex` /
  `fertility_slice_diagnosis_report.pdf`: fertility diagnostics by age,
  location, tenure, income, and wealth.
- `intergenerational_housing_fertility_v4.tex` /
  `intergenerational_housing_fertility_v4.pdf`: revised full non-spatial
  intergenerational housing-fertility draft with a cleaned compact analytical
  model and the quantitative blueprint retained from v3.
- `intergenerational_housing_fertility_paper_draft.tex` /
  `intergenerational_housing_fertility_paper_draft.pdf`: full paper draft with
  the mixed-tenure analytical section and proof appendix integrated; the
  quantitative lifecycle material remains the existing draft block.
- `housing_efficiency_proof_audit.tex`: standalone reconstruction and proof
  audit of the simplified model's local first-best, constrained-efficiency,
  tenure, and planner results. The verified reader PDF is written to
  `output/pdf/housing_efficiency_proof_audit.pdf`.
- `fertility_population_housing_transition_note.tex`: unified, self-contained
  report on reproductive equilibrium, demographic momentum, and
  intergenerational housing. It begins with the transparent representative-
  family model, adds young and old cohorts to derive the exact conditions under
  which population and prices can rise during a low-fertility transition, and
  then maps both layers into the full lifecycle model. It also states the
  interpretation of the fertility path, the appropriate initial distribution,
  the identification requirement, policy incidence, a minimal quantitative
  roadmap, and the current limitation of the maintained calibration. Figures
  are regenerated by
  `code/model/tools/build_fertility_population_housing_transition_figures.py`
  and
  `code/model/tools/build_demographic_momentum_intergenerational_housing_illustration.py`;
  the verified reader PDF is written to
  `output/pdf/fertility_population_housing_transition_note.pdf`.
- `dynamic_intergenerational_housing_fertility_model.tex`: self-contained
  research update on endogenous population scale in the paper's sequential-
  fertility model. It presents the transparent two-generation theory,
  demographic momentum, the open and closed stationary conditions, the exact
  fertility--price schedule, policy incidence with endogenous population, and
  a 2007-normalized temporary-equilibrium fertility transition with fixed-stock
  and static-elastic housing-supply cases, together with its historical
  population, rent, house-price, and price-to-rent comparison. The Stone--
  Geary one-shot implementation remains a preserved fallback rather than the
  paper-facing quantitative model. It also reports, without promoting, an
  inherited-cohort existence check showing what extra demographic initialization
  would be needed to match the interim population and housing-cost signs.
  Figures and numerical checks are regenerated by
  `code/model/tools/audit_closed_reproductive_closure.py` and
  `code/model/tools/run_e5f_open_population_transition.py`; the verified reader
  PDF is written to
  `output/pdf/dynamic_intergenerational_housing_fertility_model.pdf`.

The three immediate development notes subsumed by this report are retained
together under `archive/reproductive_equilibrium_development_20260813/`; they
are working history, not alternative active drafts.

Build/support files:

- `.latexmkrc`
- `main_bib.bib`
- `lit_review_extra.bib`
- `new_proposal_refs.bib`

Candidate figure drafts:

- `figures/fig6_ce_planner_wedge.*`: standalone draft figure for the
  competitive-equilibrium versus planner housing wedge.
- `figures/fig7_entry_fertility_decomposition.*`: standalone draft figure for
  the outside-option entry margin and aggregate fertility decomposition.

Archived material:

- `archive/legacy_2026-05-07/`: old drafts, old slides, build artifacts,
  diagnostic images, and theory experiments moved out of the active folder.
- `archive/intergenerational_housing_fertility_v3_20260609/`: archived v3
  source, PDF, and bundle readme used to initialize the v4 draft.
