# LaTeX Workspace

Active documents:

- Slides-only collaboration handoff: `../docs/prompts/HANDOFF_september14_slides.md`.
  It identifies the single working deck, technical task contacts, and the
  September 11 utility/pension updates and their presentation status below.

- `september_14_presentation.tex`: the single working September 14 seminar
  presentation. The reader PDF is `../output/pdf/september_14_presentation.pdf`;
  `september_14_presentation.pdf` is an identical build copy. The September 10
  refocus follows the May deck's appearance. The September 12 structure follows
  May: Quantitative Model, a short Empirics section, then Quantification.
  Empirics contains the data overview and the AHS and PSID evidence frames.
  Quantification contains calibration strategy, empirical moment construction,
  complete initial target and parameter tables, fit plots, lifecycle equilibrium
  profiles, the historical transition design and illustrations, and the policy
  comparisons; policy has no separate section divider.
  It presents the quantitative household environment, sequential choices,
  population accounting, equilibrium, and a policy comparison introduced along
  the same inherited transition. The three original August
  `housing_fertility_stage_{initial,impact,adjustment}.pdf` figures are reused
  unchanged within Quantification, following the dated transition strategy. They are explicitly schematic:
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

  **September 12 return to May's model layout.** Unchanged primitives use
  May's organization: Environment with a named choice list; one Preferences
  frame with flow utility, family needs and housing services; a separate Bequests
  frame; Earnings and Pensions; and one Housing, Tenure, and Budget Constraints
  frame with saving written explicitly for renters and owners. The transaction
  helper was removed from the main slides. Current return timing is preserved,
  and the actual collateral/unsecured-debt floors are defined in a linked
  Borrowing Limits appendix. May's unconditional renter nonnegativity rule is
  not valid for the maintained debt-carryover specification. The quantitative
  task checked the expanded equations and the exact debt-floor function against
  its pinned `intergen_eqscale_seq_optimized/parameters.py:622` and
  `solver.py:112--163`. The household-problem overlays retain sequential timing
  with shorter state and Bellman notation.

  **Choice specification.** The current calibration and transition calculations
  enforce `joint_nested_choice=False`: see
  `../tmp/e5f_matched_pf/code/model/tools/run_e5f_initial_revision_probe.py:77--81`
  and `../tmp/e5f_matched_pf/code/model/tools/e5f_approved_initial_state.py:32`.
  Fertility choice integrates over conception and subsequent housing decisions;
  the housing stage has its separate extreme-value taste shock. The simultaneous
  nested alternative was tested, but no permanent author rejection of that
  alternative was verified. The deck describes the active calculation rather
  than claiming that every nested alternative has been definitively abandoned.

  **September 12 policy overview and transition figures.** The two policy
  frames now describe the common inherited 2023 economy, permanent annual
  property tax of 1% versus 2%, equal per-head distribution of each path's own
  revenue, common preference/supply/demographic environment, and separate
  balanced pension budget. They describe an experiment, not computed effects.
  The quantitative task reconfirmed this contract directly. An accepted
  historical 2023 state and the entry/demographic closure remain unresolved.
  At the author's request, **Review quantitative model** supplied numerical
  transition panels and is investigating a bounded fertility-preference IRF.
  The separate review figures are
  [fertility and demography](../output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/seminar_transition_panels/fertility_demography.pdf)
  and [housing](../output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/seminar_transition_panels/housing.pdf).
  Their directory contains the plotted CSV, source hashes, and the IRF design.
  Saved short forecasts are responses to a permanent preference step, with
  finite market/fiscal closure but uncertified terminal distance; they cannot
  be labeled a fitted history or a temporary, baseline-subtracted impulse
  response. Their initial calibration also differs from the latest provisional
  initial-fit tables. No such plot has yet been inserted into the deck.
  The first no-shock terminal probe failed at its initial demographic-renewal
  guess before completing a root evaluation. This is not a nonexistence result;
  a bounded diagnosis under the unchanged closure is requested. No matched
  baseline/pulse IRF has been computed.

  **Two-slide algorithm walkthrough.** The first two numbered frames in
  Quantification explain the actual nested loops in literal steps. For each
  structural candidate, the initial pension is derived from the stationary
  age--earnings marginal (not separately iterated); a bracketed preference
  normalization wraps complete housing-equilibrium solves. The outer search
  uses finite-difference moment slopes and bounded ridge-regularized joint
  proposals. The historical slide distinguishes the conditional price/pension
  path iterations, the scalar preference search, and advancement to the next
  calendar window only after acceptance. The quantitative task verified this
  nesting and the stated tolerances directly on September 12.
  `../output/pdf/september_14_calibration_walkthrough.pdf` is a two-page extract
  of these frames for sharing; it has no independent slide source. Regenerate
  it with pypdf by locating the two pages titled `1. Calibrating the Initial
  Economy` and `2. Fitting the Historical Transition` in the verified full deck.
  These replace the former generic Calibration Strategy and From the Initial
  Economy to 2023 frames, leaving the full deck length unchanged.

  **September 12 model and quantification revision.** The Environment and
  household-problem overlays return to May's age/choice/heterogeneity exposition,
  retaining the active one-market sequential model. Household state is
  `(b,h,z,n,m)` at date `t`, age `a`: liquid wealth, inherited owner housing,
  persistent Markov earnings, lifetime parity, and dependent-child count.
  The quantitative task confirmed that permanent-income extensions are off and
  the active child state is a count, not a vector of child ages. Values now show
  their state arguments, fertility outcomes, earnings and maturation expectations,
  and the associated attempt probability. The owner service multiplier is
  positive; values above one imply a premium, consistent with its search range.
  The classroom equilibrium definition and explicit population accounting remain
  at the end of Model. The stationary definition is in the appendix.

  Initial tables and the two supplementary PGFPlots fit graphics use the same
  September 12 provisional candidate, verified by **Review quantitative model**:
  `../output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/verified_final/`.
  [Complete target-fit CSV](../output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/verified_final/selected_target_fit.csv)
  preserves every target, model value, gap, weight, loss contribution, provenance,
  sample and measurement caveat; the [complete parameter CSV](../output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/verified_final/selected_parameters.csv)
  preserves estimates, bounds and external restrictions. The source
  [readout](../output/model/e5f_matched_pf_20260909a/current_candidate_transition/overnight_20260912/verified_final/README.md)
  contains both complete human-readable tables. Fingerprint:
  `c0e266d3a0d430343c469d780d1aedb45fa87f8763c9c938889e0c37daa31de2`.
  All 13 targets appear in the main tables; all 12 scored targets appear in the
  fit plots as model/target ratios. The unscored 2.1 normalization is separate.
  The main parameter tables contain all 17 free/fixed/normalized/derived rows;
  the appendix gives every free coordinate's bounds and position. Near-bound
  means within 1% of the allowed range, not proof of a binding restriction.
  The existing `fertility_by_age.png` and `ownership_by_age.png` in
  `selected_standard_diagnostics/` are included unchanged as model equilibrium
  profiles, not empirical-fit overlays. The stable 17-graph packet is unchanged.

  The AHS 2023 tenure-by-bedroom figure remains the May asset. The existing
  September 5 PSID profile remains visible but is explicitly provisional.
  **Data — PSID event studies** confirmed no approved replacement figure or target:
  the retained 0.720246 room contrast is the old computational -1/+3 contract,
  with timing and reference support under review. It is not the newly requested
  -2 comparison. No diagnostic .403018 or .238478 estimate was promoted.
  The housing-target row and fit plot flag this limitation. No regression,
  calibration, model solve, or new mechanism figure was generated in this task.

  **Historical transition and evidence boundary.** Quantification explains
  successive unanticipated preference innovations at 2007/2011/2015/2019,
  matching the subsequent four-year fertility windows from
  `../output/model/e5f_matched_pf_20260909a/current_candidate_transition/inputs/empirical_blocks.csv`.
  Each level is believed
  permanent; its conditional forward-looking equilibrium is solved, only the
  first four-year period is realized, and the next innovation inherits household
  states and birth queues. Preferences remain at the final level after 2023.
  This replaces the old single announced-path description. The retained
  age-specific births per adult-household diagnostic is an approximate analogue
  of female TFR; exposure mapping remains unresolved. The pre-2023 demographic
  bridge conditions on external information and is not demographic validation.
  A complete historical fit, horizon validation and policy results remain absent.
  Existing historical workers use an earlier pinned initial point: the displayed
  latest initial diagnostics are not claimed to generate their trial paths.
  The policy frames define comparisons from a common inherited 2023 state.
  Live scientific closure and numerical status remain in `../CALIBRATION_STATUS.md`.

  Build twice from `latex/`, writing auxiliary files outside the active folder:
  `pdflatex -interaction=nonstopmode -halt-on-error -output-directory=../tmp/september_slides_review september_14_presentation.tex`.
  Copy the verified PDF to the reader path and the adjacent build copy.
  The deck has 36 numbered main frames and nine appendix frames
  (52 PDF pages including overlays and dividers). Two final compilation passes
  have no warnings or overfull boxes; all frames were visually inspected and
  all figure paths and all three appendix links resolve.

  **September 12 equilibrium wording.** The main definition is titled
  `Equilibrium` and names sequences of policy functions, value functions,
  distributions, populations, prices and fiscal objects in prose. Its classroom
  reference is `/Users/tommasodesanto/Documents/5_Stationary_Equilibria.pdf`,
  section 1.2, pp. 5--6: name the equilibrium objects, then require household
  optimality, market clearing and distribution consistency. The stationary
  invariance condition is adapted to this model's demographic evolution.
  Population accounting is explicit: fertility policies produce births, the
  annual person law includes surviving cohorts/newborns and net migration, and
  fixed headship rates determine head counts. Household mass reflects choices
  plus demographic entry/exit and equals heads over represented adult ages.
  The named quantitative task verified these statements against
  `person_cohort_law.py`, `four_year_bridge.py`, and `household_head_bridge.py`
  in its active `demographic_transition` package. The annual law uses `y`;
  household conditions use model dates `t`. Physical-room clearing is written
  simply as demand equals supply; the supply expansion remains on its own slide.
  The quantitative task also confirmed the neutral expectations wording.
  The successive-surprise information assumption was subsequently integrated
  into Quantification in the revision described above. Its implementation is
  `e5f_successive_surprises.py`, `evaluate_forecast`, in the quantitative task's
  active isolated model tools; the technical task verified the contract directly.

  **September 11 slides-only pass (retained in the consolidated layout).**
  The former `Family Space` and `Empirical Discipline`
  now use the parenthood-only requirement `h_P 1{m>0}`, while dependents are
  present. The equivalence scale and surrounding preferences are unchanged.
  Sources: the handoff's immediate corrections and the September 11 approved
  utility entry in `../memory/daily/2026-09-11.md`, confirmed by the named
  quantitative task against `../tmp/e5f_matched_pf/code/model/tools/e5f_parenthood_utility.py:67`
  and `../tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/solver.py:2262`.

  `Earnings and Pensions` explains the fixed payroll tax and benefits that balance the actual
  dated household budget, with the same benefits anticipated by households.
  Pension notation is `varpi_t`, to avoid conflict with choice probabilities
  `pi_t` and liquid wealth `b`. Disposable income excludes property-tax rebates.
  The combined budget frame applies the financial return to liquid wealth after
  the housing transaction. Equilibrium notation and the appendix solution steps
  include pensions consistently. The named quantitative task supplied these
  equations, and the slide editor checked the cited source excerpts:
  `../tmp/e5f_matched_pf/code/model/tools/e5f_social_security.py:29,70,99`,
  `../tmp/e5f_matched_pf/code/model/tools/e5f_balanced_history.py:198`,
  `../tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/solver.py:228,2522`,
  and `../tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/kernels.py:709,722,737`.
  These are implementation-consistent equations, not a historical-fit or policy
  certificate. All changed frames were visually checked after two clean builds.

  The named theory task confirmed a conditional equal-weight lifetime-welfare
  improvement from a small rental-space transfer to a parent: both renters must
  be interior, with positive consumption and `c_y < sqrt(e(m_y)) c_o`.
  The direct allocator preserves aggregate resources and future states, but
  need not preserve individual rental expenditure entitlements. This is not an
  owner-housing or mortgage-constraint theorem and has not been added to the deck.
  Proof source: the September 11 completed reply in
  `https://chatgpt.com/c/6aa4182b-2214-83ea-8184-efef74d558bb`, assistant message
  `92769ac8-590a-48df-a646-98f6941c04b8`, equations (4)--(6), (10), and the
  full-path feasibility argument, independently checked by the theory task.

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
