# Simplified OLG amendment checks

This is the supporting output folder for the September 4–5 theory development.
The live phase/checkpoint record is `docs/model/simplified_olg_overnight_work.md`.
The proposal is `latex/JMP_DS_suggestions/simplified_olg_amendment_proposal.tex`;
its readable checkpoint is `output/pdf/simplified_olg_amendment_proposal.pdf`.

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
