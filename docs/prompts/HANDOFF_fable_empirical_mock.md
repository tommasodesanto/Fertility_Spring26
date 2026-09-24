# Fable: write the empirical evidence section of the JMP mock

Tommaso wants you to turn the empirical work you have been doing together into
one coherent paper section. This is a writing assignment using existing results,
not another regression campaign. Use Claude Opus 5.5 through the authenticated
Claude Max route. Do not silently substitute another drafting model.

Project root: `/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26`.
All paths below are relative to this root. Read the actual saved sources; this
handoff is a map, not a replacement for estimator metadata or author decisions.
Snapshot: September 24, 2026. Later explicit author decisions take precedence.

## Deliverable and ownership

Write into the EXISTING `latex/JMP_DS_mock/sections/02_empirical_evidence.tex`.
The file already has uncommitted prose: read it immediately before editing,
preserve anything useful, and replace obsolete empirical claims deliberately.
Do not create another suggestions folder or manuscript. Update the EXISTING `latex/JMP_DS_mock/data_appendix.tex` where needed
and wire it into `JMP_DS_mock.tex`; it currently exists but is not included.
Do not create a competing empirical appendix. Use existing figures/tables. Keep main-file changes limited to required wiring.
Use the existing bibliography; add only verified missing entries if necessary.

Do not edit `latex/JMP_DS_draft/`, the slides, the model, data, regression scripts,
or `sections/04_quantification.tex`. Read the author's current model section for
voice and notation, but do not propagate its content into other sections.
No new estimation, cluster jobs, calibrations, or policy exercises. Do not change
an empirical target merely by putting a revised estimate into prose. Preserve
unrelated dirty files and any concurrent edits.

Complete one grounded writing pass, one source check, and one typesetting pass.
Stop after delivering the section, necessary appendix, compiled PDF, and a compact
receipt. If a source or decision is missing, finish the independent text and
list the exact unresolved choice outside the paper; do not launch new analysis.

## Read first

Follow the repository startup, then read:

- `docs/style/econ_writing_style_guide.md` in full.
- `latex/README.md` for the three document roles and current discrepancies.
- `latex/JMP_DS_draft/sections/03_model.tex` (read only).
- `latex/JMP_DS_mock/sections/02_empirical_evidence.tex` and
  `latex/JMP_DS_mock/sections/04_quantification.tex` (the latter read only).
- September 24 entries in `CALIBRATION_STATUS.md` and `memory/daily/2026-09-24.md`.
- `docs/model/ACTIVE_DECISION_LEDGER.md` for author decisions, not a license to
  settle open choices yourself.

The mock quantification was updated on September 24 to reflect September 23
choices. It is not an authoritative source for today's empirical estimates:
its first-birth row still says 0.600, and its bequest row still says 0.88%.
Today's results/decisions have moved on. Flag cross-section discrepancies in the
receipt, without editing quantification in this assignment. Its blank model-fit
cells must remain blank. The mock README's original setup inventory is historical.

## Paper question and empirical logic

The paper asks how housing costs and access to family-sized housing interact
with fertility, tenure, and housing allocation over the lifecycle. It is not
being redirected into a collateral-only paper or a narrative about failed tests.
The evidence should make a small number of economic facts understandable:

1. How housing size is distributed across tenure, and what that implies about
   the housing bundles families actually occupy.
2. How rooms, ownership, and mobility change around family formation, including
   adjustment before childbirth and differences by initial household position.
3. What complementary cross-sectional and instrument-based exercises add, and
   what variation they actually identify.

These are proposed organizing questions, not predetermined empirical conclusions.
Let verified evidence determine the claims. Explain how facts motivate model
ingredients without claiming that they uniquely identify a mechanism. Neither
large-home ownership shares nor housing adjustment around births proves that
housing constraints cause low fertility. Model policy counterfactuals are distinct
from empirical event studies and fertility instruments; keep them out of this
section. Do not resurrect the old spatial model or center–periphery evidence as
a required paper component simply because old packets contain it.

## Evidence map: start with the latest work you did with Tommaso

### PSID: main longitudinal evidence

Read the complete README, `window_path.csv`, `summary.csv`, and verification/fit metadata in:

- `code/data/psid_followup_mar2026/output/sa_rooms_first_birth_v2/`
- `code/data/psid_followup_mar2026/output/sa_first_birth_outcomes_v3/`

Use the actual saved coefficients, covariance matrices, support and run receipts
for each selected curve. The older chronology and measurement reconciliation are
in `docs/model/first_birth_rooms_timing_memo_20260923.md` and
`code/data/psid_followup_mar2026/README_psid_followup_mar2026.md`.

The author's chosen household design uses PSID individual weights, biological
birth histories, entry at first adult observation, and a −3/−2 reference window.
It selects a woman who is reference person or spouse within the household-year.
The later status-at-baseline comparison instead fixes household position before
treatment; do not call these the same sample. Read the ending updates in the
rooms README: its early interpretation of anticipation is superseded where the
later diagnosis attributes the household-design pre-rise to row selection.
The latest status says the final headline sample, baseline, and model horizon
still need reconciliation. Do not choose them by significance or model fit.
Present well-defined comparisons and place the exact author decision outside
paper prose rather than silently choosing a new headline.

Explain the biological event clock, biennial interview windows, reference group,
fixed effects, weights, clustering, and support from the actual selected estimator.
Do not describe a cohort-interacted estimator as generic two-way fixed effects.
The current estimator uses person and survey-year fixed effects, age and education
covariates, longitudinal IW weights and person-clustered uncertainty. Verify
these in each selected fit, including outcome-specific complete-case samples.
Missing observations outside a window are not zeros. The latest official rooms
reconstruction retains valid zero room counts and excludes vintage-specific
non-room codes; do not apply an older blanket rule that every zero is missing.
Rooms and moving variables
were rebuilt/checked against official yearly PSID items; ownership has a different
dating convention. Use the audited crosswalk, not a plausible variable label.

Distinguish: moving at all; moving for space as a share of all observations; and
moving for space conditional on moving. The unconditional response can be small
while the conditional share rises because overall mobility falls. Neighborhood
moves are a separate outcome. Explain this economic distinction plainly. Conditioning on movers selects on an
outcome: the conditional share is not a substitute for the unconditional effect
or causal mechanism evidence. Separate outcome paths do not prove that the same
household buys a larger home in a single observed joint event.

The old 0.72 sentence in the mock and historical 0.60/0.73/0.93 summaries are
not interchangeable with the final window/sample estimates. Do not import an old
second-birth result without checking its current measurement and validation.

### Housing stock and tenure: descriptive motivation

- `docs/model/evidence_tenure_segmentation_20260918.md`
- `code/data/ahs_supply_snapshot/README.md`
- `code/data/ahs_supply_snapshot/output_ahs_family_unit_menu_national/AHS_2023_FAMILY_UNIT_MENU.md`
- The tables and figure packet in that same national output directory.

Verify denominators: share of owners in large homes differs from share of large
homes rented. This describes occupied stock, not a causal supply elasticity or
proof that a particular household cannot rent a large home. Keep rooms and
bedrooms distinct. Do not silently mix metro, national, or period definitions.

### ACS matched pseudo-panel: complementary descriptive evidence

- `code/empirical/acs/kleven_pseudo/README.md`
- `first_birth_housing_run_report_20260920_job18079576.md`
- `first_birth_sensitivity_bundle_receipt_18080591.md`
- `second_birth_housing_run_report_20260920_job18080900.md`
- `second_birth_proxy_design.md` and the saved support/fit receipts those notes link.

The report filenames above are under `kleven_pseudo/`. National continuation
results also exist at `output/national_acs_comparison/national_continuation_20260921b/`:
read `continuation_receipt.json` and the selected `*_receipt.json`, contrast and
covariance files. Do not repeat the older README's claim that national matching
is still pending; independently distinguish completed fits from unavailable ones. Check exact geography
and coverage; regional descriptive curves are not national estimates. Constructed
prebirth donors are not repeated observations of the same women. Roster-based
birth ordering is not a complete biological history. Explain the match and
common-cohort support briefly. Use only exercises that add something to the PSID
story; do not inventory every exploratory fit in the main text.

### Fertility instruments: a separate identification exercise

- `output/acs_fertility_iv/data_appendix/SOURCE_AUDIT.md`
- `output/acs_fertility_iv/national_128g_results/review_summary.json`
- `output/acs_fertility_iv/national_128g_results/national_18row_table.csv`
- `output/acs_fertility_iv/samesex_diagnosis/report.md`
- `latex/JMP_DS_draft/sections/appendix_acs_fertility_iv.tex` and its bibliography
  (read only; reuse reviewed substance by adapting it into the mock if useful).

For the PSID instrument exercises, read
`code/data/psid_followup_mar2026/output/iv_housing_reaudit_20260809/README.md`:
it supersedes the March first-pass interpretations. These are exploratory
triangulation, with limited twin-like support and weak same-sex first stages.
For ACS read `code/empirical/acs/kleven_pseudo/acs_twins_samesex_housing_contract.md`.
Roster proxies do not establish biological motherhood or exact twin births.
Do not restrict the sex-composition sample to exactly two children, which
conditions on a post-instrument outcome. The same-sex housing pattern raises
exclusion concerns; it does not formally reject exclusion or identify its cause.

The national IV results and same-sex checks ARE completed; some older inventories
say otherwise. Use final receipts. Distinguish twins from first-two-child sex
composition, their realization dates, samples, first stages, reduced forms, IV
estimates and weak-instrument uncertainty. Housing can respond directly to sibling
sex composition; take the completed exclusion diagnostics seriously. Do not label
a descriptive pseudo-panel or a weak/exclusion-sensitive IV estimate a clean
fertility experiment. Put detailed robustness and identification checks in an
appendix where useful, without framing the paper around negative diagnostics.

## Writing and presentation

Write paper prose in the first person, matching the author's draft. Start from
the economic fact, explain the comparison, then give the evidence. Use ordinary
language and the restrained exposition of Menzio/Fernández and the author's
Boar–Gorea–Midrigan reference. Consult actual source passages before claiming to
follow a particular paper; do not copy wording. Do not change the paper's voice
into a project report. Preserve author-chosen notation and explain new objects.

Use compact subsections and the mock's existing run-in subheaders. No workflow
footnotes, job IDs, internal specification nicknames, calibration receipts, or
repeated caveat lists in the paper. State economically material limitations once,
where they affect interpretation. Keep pending author choices in the handoff
receipt. Do not use “parity.” Display at most three decimal places; use meaningful
units, confidence intervals and sample sizes. Cite estimators and data sources
accurately from verified references, not memory.

Choose a small number of existing figures/tables that support the narrative.
Retain their economic definitions; do not generate new scientific results.
A table should identify the sample, window and outcome rather than conceal
incomparable estimates in a single “effect” column. Do not promise an appendix
that does not exist: resolve the current missing `app:empirical-data` reference
by updating and wiring the existing `data_appendix.tex`, after checking its label.

## Completion check

Check each quantitative claim against its saved source and record a compact
claim-to-source map outside the manuscript. Compile the existing mock twice,
keeping build products outside its source directory; visually inspect changed
pages and deliver `output/pdf/JMP_DS_mock.pdf`. Do not open PDFs automatically.
Report the files changed, the narrative selected, unresolved decisions, and any
unverified claims or compilation issue. Stop there and return control to Tommaso.
