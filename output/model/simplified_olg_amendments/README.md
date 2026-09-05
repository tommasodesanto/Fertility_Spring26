# Simplified OLG amendment checks

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
