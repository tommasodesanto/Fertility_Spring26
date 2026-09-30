# Entrant wealth distribution

Reference: **2007 stationary reference — block0506, September 28 verified export**.

This is a small extraction from authenticated saved entrant arrays, not a new
household solve or equilibrium. No checkpoint download or model import occurred.
See `entry_distribution_summary.json` for paths, array IDs, source hashes and
checks; `entry_distribution_summary.csv` supplies compact tabular results.

The entrant joint probability is the saved conditional wealth matrix multiplied
by the saved income-state weights. Column sums and total probability equal one;
the conditional matrix hash and wealth grid match the frozen manifest. Wealth
units are equal-working-age mean annual gross earnings. Current income covers
four years. These are new age-18 renters, not the full stationary population.

Negative / zero / positive financial wealth shares are 26.2509% / 28.3397% /
45.4094%; mean wealth is 0.186520. The two occupied states with nonpositive cash
have wealth index 44 and income indices 0 and 1. Their joint probability is
0.00008021960141084162, or 0.00802196% of entrants. The lead independently
recomputed these results from the small saved arrays. Zero-percent quantiles
were removed because the original extraction selected zero-mass grid points;
the reported interior quantiles are unaffected.

The five empirical wealth-to-income bin means are input approximations, not
individual microdata or empirical quantiles. Their negative-bin weight is not
the actual on-grid negative-wealth share. The wealth–income coupling is retained
from the reference and is not an estimated empirical joint distribution (see
`code/model/tools/e5f_earnings_wealth_contract.py`). The tiny failing probability
therefore depends on that approximation. It does not measure the total behavioral
impact of a zero borrowing limit.

A proposed 0.139535 top-up for the two failing states is one existing grid step,
not a minimum continuous transfer. No top-up, entry revision or new economic
case has been implemented or adopted by this extraction.
