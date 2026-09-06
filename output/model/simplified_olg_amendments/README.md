# Simplified OLG amendment checks

## September 6: transition with substantial renting and positive child costs

Start with the revised **four-page** `output/pdf/simplified_olg_paper_core.pdf`.
Its allocation and conditional-fertility pages are unchanged. Figure 2 now
uses the analytical mixed-tenure response. The **eight-page** supporting proof
is `output/pdf/simplified_olg_mixed_transition_proof.pdf`; both LaTeX sources
are under `latex/JMP_DS_suggestions/`. The earlier limiting figure, ten-page
assessment and all-owner proof are preserved below.

The new exact six-variable representation retains original dated budgets,
positive child costs, a positive rebated property tax, finite logistic tastes,
endogenous tenure and the actual initial old's assets and housing titles.
An explicit economy has `11/21` owners, child goods cost `3/20`, and property
tax `467/9250`. The same stationary bundles support a family of taste scales.
Within each credit reform both taste parameters remain fixed.

- For every positive taste scale: four stable/two unstable roots, a positive
  stationary population derivative, and negative stationary conditional
  lifetime-value derivatives in both tenures. Exact polynomial coefficient
  signs prove these statements throughout the positive half-line.
- For every scale in `[1,4]`: the actual initial-old boundary is transverse
  and initial fertility rises. Rational intervals cover the entire interval,
  with 913 adjacent cells; this is not just a list of parameter-point solves.
  The sequence implicit-function argument gives locally unique nonlinear
  converging paths after small permanent credit relaxations.
- Near scale `1`: a nonzero dominant complex mode and a uniform nonlinear
  remainder bound prove that small finite reforms eventually approach
  replacement fertility, and final population, from both sides. The initial
  fertility increase and final population gain do not imply monotonicity.

Small primitive and compact-support entrant-distribution changes preserve the
strict conditions. The reform radius and a multidimensional parameter box are
not quantified. This remains an illustrative local result, with fixed entrant
endowments and world bond price. Population growth is distinct from the
compensated allocation welfare gain and from a social population ranking.
No planner permission, author convention, or quantitative-model change is made.

Evidence:

- `mixed_transition_proof.md`: full derivation, general local theorem,
  positive-tax inversion lemma, exact construction, finite-reform oscillation
  proof, and the wider family result.
- `mixed_transition_map_review.md`, `mixed_transition_certificate_review.md`,
  `mixed_transition_family_review.md`: three sequential bounded read-only
  reviews of distinct new steps. The second requested an explicit uniform
  tail lemma; the third confirms that the added lemma closes that issue.
  The last review audited source and receipt but could not replay them in its
  model environment because `sympy` was absent. The lead executed the complete
  checks in the working Python environment.
- `mixed_transition_certificate.json`: exact rational matrices, root bounds,
  boundary determinant, response signs and complex-mode signal; also records
  the analytical Figure 2 values.
- `mixed_transition_family_certificate.json`: all-positive-scale Routh,
  discriminant and stationary-value proofs; complete `[1,4]` interval coverage;
  original-equation comparisons at four declared scales. Completed in about
  91 seconds with a 600-second/12,000-subdivision fail limit.
- `mixed_transition_smoke.json`, `mixed_transition_checks.json`: exact original
  branch checks, four short original-equation paths, 24 household optimizations,
  central derivatives and 24/40-date horizon comparison. Maximum equilibrium
  residual is below `2.9e-15`, budget error below `1.6e-15`, and initial/final
  derivative discrepancies below `3e-9`. Finite paths check arithmetic; they
  do not supply the convergence proof or the plotted curves.

One new driver reproduces the work; the old limiting driver is unchanged:

```sh
python3 code/model/tools/verify_simplified_olg_mixed_transition.py --smoke
python3 code/model/tools/verify_simplified_olg_mixed_transition.py
python3 code/model/tools/verify_simplified_olg_mixed_transition.py --certificate-only --figure
python3 code/model/tools/verify_simplified_olg_mixed_transition.py --family-only
```

The verified interpreter is `/Library/Developer/CommandLineTools/usr/bin/python3`,
with NumPy 2.0.2, SciPy 1.13.1, SymPy 1.14.0 and Matplotlib 3.9.4. Use that
interpreter for these commands if the active model virtual environment lacks
SymPy. No package or environment changes were made for this proof.

All four receipts match the final new source hash. The supplemental figure
files are `mixed_transition_analytical_figure.pdf/png`; the compact note uses
this response as Figure 2, with the allocation figure unchanged.

## September 5 late-evening scope and prose pass

Start with `output/pdf/simplified_olg_paper_core.pdf`: **four pages**, comprising
three pages of proposed results after the existing household setup and one page
of economic assessment. Both existing figures are retained. Source:
`latex/JMP_DS_suggestions/simplified_olg_paper_core.tex`. The ten-page assessment
below remains the fuller reading note; neither document is an adopted paper edit.

The new section 11 of `local_transition_proof.md` explains the limiting
stationary credit result for heterogeneous entrants. Price rises in proportion
to common borrowing capacity, each young type's housing and fertility remain
unchanged, and every old type uses less housing. Stationary population rises,
but conditional owner lifetime utility falls for every same-type entrant.
This holds at the all-owner, zero-child-cost/tax limit with the original
uniformly strict branches and fixed income, wealth distribution and world bond
price. It is separate from transitional-cohort welfare, variable-population
social welfare, and the compensated allocation result.

- `local_transition_scope_review.md`: one bounded read-only second opinion,
  independently deriving the general welfare formula and assessing the economic
  assumptions. The lead qualifies its estate wording: the original estate
  inequality is not a minimum-retention rule and does not prohibit a sale.
- `local_transition_welfare_checks.json`: exact symbolic utility and housing
  identities, original budgets and constraints for five fixed types at three
  financed shares, ten original household optimizations, and utility derivatives
  checked to `8.3e-9`. Regenerate with the existing driver:

```sh
python3 code/model/tools/verify_simplified_olg_local_transition.py --welfare-only
```

Recommendation, pending author discussion: retain the allocation and conditional
fertility arguments as the core; treat the second figure as a separate local
demographic extension, with the price and welfare interpretation explicit.
Near-all ownership, small child goods cost, fixed entrant endowments and
unchanged household constraints limit its economic coverage. The simple
stability inequality does not make all these restrictions innocuous. No new
planner power, ownership preference, population convention, or model is chosen.

## September 5 continuation: analytical transition and primitive conditions

The longer reading note is **ten pages**, at
`output/pdf/simplified_olg_simple_assessment.pdf`; its LaTeX source remains
`latex/JMP_DS_suggestions/simplified_olg_simple_assessment.tex`. Pages 1–4 keep
the allocation and conditional fertility arguments. Page 5 replaces the
prescribed demographic curve with the analytical equilibrium response.
Appendix B gives primitive stability and initial-fertility conditions;
Appendix C gives a stationary population condition with positive child costs.
No preference or planner-instrument decision is adopted.

The all-owner limit with zero child goods cost and property tax has an exact
aggregate recurrence, including with heterogeneous entrants on uniformly
strict household branches. The condition `ell + (3+q) C > 4D`, with each
coefficient defined from primitives in the proof, gives exactly two stable
roots and one unstable root. The predetermined initial cohort and actual old
assets select a local converging equilibrium. The exact initial-fertility test
uses a primitive polynomial evaluation; a simpler sufficient inequality is
also available. These extend locally to positive renter mass, child costs and
property tax through the full dated equilibrium operator.

The original model's initial-old accounting, endogenous tenure, extra price
leads, normalized type distribution and all private constraints are retained.
The positive-parameter neighborhood and reform must be small; their numerical
size and a global transition theorem remain unproved. Initial fertility and
terminal population rise under the stated conditions. All-date finite-reform
fertility and monotone population are proved only for the plotted limiting
example. An admissible second economy instead has lower initial fertility and
a larger final population. The uncompensated credit reform lowers stationary
household welfare in the plotted limiting example; it is not the compensated
Pareto allocation proof.

Supporting files:

- `local_transition_proof.md`: full derivation, general primitive conditions,
  exact aggregation, uniform infinite-tail proof, welfare distinction, and
  counterexamples. Sections 10.4 and 10.2 contain the broad stability condition
  and exact initial-fertility test.
- `local_transition_independent_review.md`: first review of the explicit
  recurrence, household choices and original initial-old boundary; identifies
  the subsequently closed mixed-extension and uniform-tail gaps.
- `local_transition_closure_review.md`: independent check of those two newly
  completed proof steps and the negative welfare derivative.
- `local_transition_heterogeneity_review.md`: independent aggregation and
  positive-child-cost stationary-condition review; required probability,
  inherited-state and uniform-operator qualifications are incorporated.
- `local_transition_general_review.md`: independent general convergence review
  under `w>d`, the exact/sufficient fertility tests and negative initial-sign
  example. This report predates the sharper condition in section 10.4.
- `local_transition_sharp_stability_review.md`: final narrow independent check
  of that sharper root condition, including zero/repeated stable roots and
  carryover of the original initial-old boundary and primitive fertility tests.
  Its requested standalone operator/decay-weight qualification is incorporated.
- `local_transition_checks.json`: eight original-equation cases, original
  household optimizations, a central difference versus the analytical path,
  and 24/40-date horizon comparison. Maximum equilibrium residual is below
  `2.3e-14`; derivative error is `1.48e-9`. Finite paths are arithmetic support,
  not the proof of infinite convergence or an admissible-neighborhood bound.
- `local_transition_stationary_checks.json`: exact symbolic derivative and
  three original stationary cases with both population signs at positive cost.
- `local_transition_heterogeneity_checks.json`: zero Hessians of eight owner
  quantity/asset maps, four-corner constraints and 24 original household
  optimizations, plus exact aggregate dated-demand and inherited-old checks.
  Its dated prices are diagnostic inputs, not equilibrium paths.
- `local_transition_general_checks.json`: six symbolic polynomial/transform
  identities, three initial-boundary identities, 270 declared coefficient
  cases (including unstable cases), and an original-equation counterexample
  with initial fertility derivative `-0.0597632` and final population derivative
  `0.0952381` per proportional credit-cap increase. This grid is an arithmetic
  check, not a calibration or an economically admissible-parameter map.
- `local_transition_smoke.json`: the exact zero/tiny-reform loop smoke check.
- `local_transition_analytical_figure.pdf` / `.png`: Figure 2, obtained directly
  from the analytical formula, never from the finite-path solver.

All verification modes take seconds locally and record the current driver hash.
From the repository root:

```sh
python3 code/model/tools/verify_simplified_olg_local_transition.py --smoke
python3 code/model/tools/verify_simplified_olg_local_transition.py
python3 code/model/tools/verify_simplified_olg_local_transition.py --stationary-only
python3 code/model/tools/verify_simplified_olg_local_transition.py --heterogeneity-only
python3 code/model/tools/verify_simplified_olg_local_transition.py --general-only
python3 code/model/tools/verify_simplified_olg_local_transition.py --figure
```

The first unscaled optimizer check at the deliberately tiny rental cap was
ill-conditioned and did not report success. Scaling the optimizer's choice
units resolved it; the original objective and constraints were retained. All
final modes pass with independent household optimization and strict private
conditions. No failed check was accepted as proof.

## September 5: earlier simple allocation assessment

The initial reading note was the eight-page
`output/pdf/simplified_olg_simple_assessment.pdf`, with source
`latex/JMP_DS_suggestions/simplified_olg_simple_assessment.tex`. It restores
the author's priority: first establish housing misallocation simply; then
assess fertility, population and constrained-finance extensions. It does not
adopt new preferences or public collection powers.

- `simple_assessment_checks.json`: four explicit mixed-tenure equilibria,
  nine independent original household optimizations, 36 finite compensated
  reallocations, owner-to-renter replication, conditional fertility signs,
  demographic accounting and stationary stock scaling. Driver/specification
  hashes are recorded. The family varies income and ownership-taste location
  with the tax rate; it is not a tax comparative static.
- `simple_assessment_direct_review.md`: independent direct-proof and
  constraint-role review. Two transcription corrections are disclosed at the
  top. This review preceded the new equilibrium constructions.
- `simple_assessment_new_claims_review.md`: separate independent derivation
  of replication, the equilibrium family, the slack-rental equilibrium, and
  stationary scaling. Its scope and incorporated clarifications are explicit.
- `simple_allocation_figure.pdf` / `.png`: exact compensated reallocation for
  the eligible pair in the constructed equilibrium; the crossing is the best
  allocation on this comparison, not a full economy-wide first best.
- `simple_population_figure.pdf` / `.png`: a prescribed fertility sequence
  and the resulting demographic identity. This is not a solved equilibrium
  transition or a prediction of adjustment speed. This earlier figure is retained
  but the latest note uses the analytical Figure 2 described above.

Regenerate this bounded packet in a few seconds, without calibration or
equilibrium search:

```sh
python3 code/model/tools/verify_simplified_olg_simple_assessment.py --plot
```

The full derivations and assessment are also recorded near the top of
`docs/model/simplified_olg_overnight_work.md`. A complete equilibrium example
now establishes that a positive group of eligible owner pairs can exist.
A second example establishes that the binding rental cap is not necessary
with the existing ownership taste. The continuation above supplies local
analytical credit-reform results; global transitions and anonymous public-loan
implementation remain unproved.

## September 5: external discussion assessment

`external_discussion_review.md` preserves the independent agent's response
returned by the author. It recommends public down-payment finance and revising
the interpretation of the existing conditional proof. It is a discussion input,
not acceptance of planner powers or verification of its proposed anonymous
pure-loan implementation. The lead's qualifications are in the live decision
ledger; the existing mathematical sources and receipts remain unchanged.

## September 5: focused constrained-efficiency result (latest)

The author asked to settle constrained inefficiency before moving to other
issues. The two existing Astra reviewers and the lead completed a new pass.
The current four-page synthesis is
`output/pdf/simplified_olg_constrained_efficiency.pdf`, from
`latex/JMP_DS_suggestions/simplified_olg_constrained_efficiency.tex`.
This result supersedes the earlier statement below that the complete-path
constrained comparison is entirely unproved.

- `constrained_full_path_review.md` proves a complete infinite-path Pareto
  improvement with committed household-specific young and old transfers,
  enforceable future taxes, and an explicit passive initial title owner with
  outside bond access and participation in the fiscal account. Young housing
  rises at every date. Initial young owners can gain strictly. Every market
  clears and all initial claims are included. It also proves a local obstruction
  to one-time gifts around the verified stationary regime, not a global no-go.
- `constrained_instruments_review.md` gives a separate finite-support fallback
  using a current-income advance and two targeted housing incentives. It is an
  expanded-instrument result, not a claim of globally minimal instruments.
- `planner_benchmark_checks.json` now also records three finite committed
  reforms, a fourth case with a strict initial-young gain, original-budget and
  full-resource identities, the infinite contraction bound, 12 independent
  household optimizations, full tenure-deviation checks at fixed individual
  fertility, and the one-time-gift price roots.
- `committed_cash_path.csv` contains dated receipts for the analytical path.
  `committed_cash_diagnostics.png` is a supplemental six-panel internal check;
  the agreed theory and earlier numerical figures are unchanged.

Regenerate the complete receipts and supplemental diagnostic with SciPy and
Matplotlib installed (about five seconds locally):

```sh
python3 code/model/tools/verify_simplified_olg_planner_benchmarks.py --original-optimizers --plot
```

The positive theorem is conditional on the stated transfer and ownership
permissions, stationary group structure, strict private branches and tenure
margin, and a contraction condition. Author acceptance of those permissions
remains OPEN in the decision ledger. In particular, public commitment provides
financing across ages; an unchanged mortgage share does not remove this extra
economic power. No new fertility externality is used. The author manuscript
and quantitative model were not edited.

## September 5 planner benchmarks

The author requested two Astra agents at maximum reasoning, one for direct
allocation and one for cash transfers followed by markets. The lead reviewed
both arguments and performed independent algebra and finite-allocation checks.
The five-page synthesis is `output/pdf/simplified_olg_planner_benchmarks.pdf`;
its source is `latex/JMP_DS_suggestions/simplified_olg_planner_benchmarks.tex`.

- `first_best_review.md`: exact current-goods compensation, original title and
  estate accounting, conditional marginal-value diagram, and global choices.
- `constrained_efficiency_review.md`: exact cash-transfer envelopes, a local
  obstruction and four-household-group improvement. The latter moves housing
  toward the old and omits future market clearing and intermediary-owner welfare.
  Neither example establishes constrained inefficiency of the complete OLG model.
- `planner_benchmark_checks.json`: independent receipts, including source hashes,
  27 direct perturbations, three exact conditional cash cases, and the obstruction.

Reproduce these small checks without equilibrium or calibration solves:

```sh
python3 code/model/tools/verify_simplified_olg_planner_benchmarks.py
```

The live decisions and remaining ownership/transfer/estate questions are in
`docs/model/ACTIVE_DECISION_LEDGER.md`. The earlier amendment and its evidence
remain below as a separate completed pass.

## Earlier amendment checks

This is the supporting output folder for the September 4–5 theory development.
The live phase/checkpoint record is `docs/model/simplified_olg_overnight_work.md`.
The proposal is `latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex`;
its readable version is `output/pdf/simplified_olg_amendment_proposal.pdf`.

Regenerate the bounded analytical checks from the repository root:

```sh
python3 code/model/tools/verify_simplified_olg_amendments.py
```

- `verification_summary.json`: check counts, errors, source/code hashes,
  constructed original-owner bundles, and common-cap correction.
- `fertility_checks.csv`: direct scalar household solutions and finite differences.
- `welfare_checks.csv`: exact compensation, uniform premium, and ordinary
  market-price comparisons; dated financing and utility changes.
- `analytical_checks.png` / `.pdf`: the stable two-panel diagnostic for these
  analytical examples. No equilibrium or calibrated result is shown.
- `independent_review.md`: completed bounded read-only reviewer report. The
  lead adopted its tax-reserve and seller-receipt clarity corrections. Its
  original line references refer to the reviewed checkpoint, not later edits.
- `independent_review.log`: transient wrapper log; not a canonical claim source.

The primitive sufficient condition and allocation proof retain their explicit
stationary/fixed-price/branch qualifications. These checks do not establish a
GE policy sign or broad transition theorem. No production calibration is run.

Regenerate the small illustrative equilibrium checks and the full six-panel
packet with:

```sh
python3 code/model/tools/verify_simplified_olg_amendments.py --transitions
```

The predeclared experiment is four financed shares (80%, 81%, 85%, 90%), with
identical initial cohorts, fixed fertility preferences, zero gains tax, and
property tax 5%. A horizon-40 comparison, second solution method, original
household optimization, and endpoint derivative checks are included. Expected
local runtime is about three minutes, with a ten-minute stop. This is a theory
example, not a calibration or production quantitative policy run.

To rebuild the compact report figure and check population accounting from
existing paths, without solving another equilibrium:

```sh
python3 code/model/tools/verify_simplified_olg_amendments.py --saved-transitions
```

- `transition_verification.json`: all policy summaries, parameters, original
  optimizer discrepancy, longer-horizon and second-solver differences.
- `transition_phi80.csv`, `transition_phi81.csv`, `transition_phi85.csv`,
  `transition_phi90.csv`: full path quantities and original-equation checks.
- `transition_phi85_h40.csv`: longer-horizon comparison.
- `transition_stability_checks.json`: initial and final endpoint linearizations,
  root counts, inherited-state projection, and derivative step comparisons.
- `original_household_optimizer_checks.csv`: 24 direct original-budget
  constrained optimizations, separately checking the reduced choices.
- `population_and_composition_checks.json`: cohort-product identity and exact
  initial-tenure-share decomposition of the fertility change.
- `credit_policy_transitions.pdf` / `.png`: full six-panel diagnostic packet.
- `credit_policy_summary.pdf` / `.png`: compact fertility/population figure for
  the consolidated proposal, with the full packet retained above.
- `transition_independent_review.md`: bounded independent mathematical review
  of the continuation theorem and finite-history argument.
- `final_independent_review.md`: final independent synthesis review; the lead's
  resolution is in the work record. Line numbers refer to the reviewed version.
- `delivery_manifest.json`: source/output hashes and the final verification
  record. Numerical solver-function hashes are distinct from presentation edits.

The 90% reform has a slightly negative initial fertility response despite a
larger final household population. The owner borrowing limit is slack in that
case. These findings are retained; the example does not establish a positive
fertility response at every date or verify interval-wide continuation hypotheses.

Build the note from the repository root, directing all intermediate products
outside the protected author manuscript. Run this twice:

```sh
/Library/TeX/texbin/pdflatex -interaction=nonstopmode -halt-on-error -output-directory=tmp/pdfs/simplified_olg_amendments latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex
```

The final PDF is copied to `output/pdf/simplified_olg_amendment_proposal.pdf`
after rendering and visual inspection. Wrapper logs, progress JSON, runtime
records, and temporary build products are not canonical evidence.
