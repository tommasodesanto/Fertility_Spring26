# Verified first-shock candidate: diagnostic packet

Open `recovery_review.pdf`. Five overview figures and the unchanged 17 saved
standard household diagnostics are included, with all 13 initial empirical
review rows and all 17 parameter/restriction rows. `shock_fit.csv` shows all
four historical targets; `forecast.csv` records all 104 candidate forecast dates.

This is the September 15 recovery's first candidate, not a completed historical
fit. Its preference is 0.12891531457859182, held permanently from 2007. Market
and fiscal convergence and exact replay pass; first-period fertility misses
its target, and later surprises remain unestimated. The separate terminal
steady state is verified, but terminal distance and horizon independence are
not certified. No model solves or new cluster jobs were run for this packet.

The grey comparator is the actual curve linked in the saved September slides.
It comes from four announced shocks, despite the old permanent-shock caption.
The recovered candidate has one smaller shock, so differences combine shock
size, expectations and numerical convergence; this is not a pure solver test.

The full initial fit table is the current provenance review of the retained
parameter vector, not a new SMM estimate. Its 2.1 normalization is unscored.
Blank weights are not zeros; all metadata and warnings remain in target_fit.csv.

Saved standard household diagnostics refer to 2007 after the first shock.
The source filename lifecycle_2023.csv is historical and does not date that
snapshot. Native 2023 fertility-by-age observations are available; 2023 housing,
wealth, ownership and intergenerational allocation by age were not saved for
this candidate. No old 2023 profile is substituted. Legacy diagnostic fertility
choice statistics and the summary field named tfr are not native period TFR;
the latter is a completed-fertility stock. Main fertility figures use the
verified native birth-flow diagnostics. Reconstructed completed fertility refers
to model age 42–45, not an empirical ages-40–44 comparison.

`review_input.json` freezes all plotted series, initial and terminal references,
full tables, and SHA-256 source hashes. `source/recovered_candidate/` retains the
unaltered downloaded receipts and standard graphs. `verification.json` records
native-flow, cohort-stock and plotted-array checks. The lead additionally
matched native price/fiscal coordinates exactly to the accepted root receipt.

Rebuild the packet without model imports or solves, from the repository root:

```sh
/Users/tommasodesanto/miniconda3/bin/python code/model/tools/build_e5f_recovery_review.py --mode all
```

The plotting runtime uses matplotlib/numpy; PDF assembly uses the bundled
Codex Python with reportlab/pypdf. Individual `--mode plots` and `--mode report`
commands remain available. Existing presentation files are unchanged.
