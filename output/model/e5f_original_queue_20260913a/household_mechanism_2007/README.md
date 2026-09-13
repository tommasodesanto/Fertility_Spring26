# 2007 household mechanism

Supplemental saved-policy figure; no model solve or calibration change. The CSV retains the full wealth grid; the figure displays valid nonnegative renter wealth through at least the slice's 99th percentile, extending to the six-room down payment if needed. Quantile and source details are in verification.json.

Birth probability is fecundity times the saved first-birth attempt probability. Housing averages the actual tenure choice conditional on a realized birth or no birth from the same initial childless state. Current income varies within the same permanent-income group. Dotted lines indicate down-payment eligibility, not proof of a binding desired borrowing constraint. This is a state-conditioned policy comparison, not an income-fertility cross-sectional regression. Curves join the existing grid points without smoothing.

Re-render with:
```
python code/model/tools/build_e5f_ss2007_household_mechanism.py --render-only --outdir output/model/e5f_original_queue_20260913a/household_mechanism_2007
```
