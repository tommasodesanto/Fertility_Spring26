# Matched saved-policy purchase response

Frozen winner31; q=.719168368828958. Compare saved baseline buyer/stayer financed shares.8/.8 with purchase-only1/.8, retaining all scalar preferences, room menu[2,4,6,8,10], parent floor2.3, targets, entry, units and incumbent-owner financing. Source folders: [baseline](../../purchase_ltv_v1/local_run/retry5/results/baseline_80_80/closure.json) and [purchase-only100](../../purchase_ltv_v1/local_run/retry8/results/purchase_100_stayer_80/closure.json). Both impact packets use exactly the same inherited PRE-fertility distribution SHA459ab9229e7ce58376bda8c933b583ec059104ede31c57627ce7368c41e5a51b.

This bounded diagnostic uses saved fertility/location/tenure policies and native calendar algebra only. It performs zero Bellman, lifecycle, GE or optimizer calls. The whole-PRE replay reproduces saved ownership and births within2e-10; every group conserves its inherited mass with zero feasibility projection. Policies were optimized separately under each financing rule, so continuation values also respond. This is not an isolated same-value purchase affordability check.

Groups are inherited PRE-fertility renters, before this period's fertility and location/tenure decisions. Young means model age-cell starts18,22,26,30,34,38,42; it is supplemental, not an annual-age or target-clock change. n=0 means no children ever born. Current parent means m>0 children at home; this differs from all ever-parents. Empty former-parent states account for the remaining young renters.

| Initial group | PRE mass | Baseline ownership | Purchase-only ownership | Change, pp | Birth-child flow change per initial household |
|---|---:|---:|---:|---:|---:|
| Young renter, n=0 |.186376|.167809|.202672|+3.48627|−.00019142|
| Young renter, current parent |.054793|.325702|.329579|+.38776|+.00000579|
| All young renters |.265355|.219089|.245100|+2.60104|−.00013236|
| Entire inherited population |1.000000|.659374|.668716|+.93424|−.00003242|

Ownership is realized current ownership after fertility and location/tenure transactions, on each group's unchanged inherited mass. n=0 birth-child flow is first-birth probability for this period; current-parent flow includes continuation births. The reported differences are net expected outcomes, not the fraction of individuals whose choice switches under a coupled taste-shock draw.

Young renters with n=0 account for about69.55% of the aggregate net ownership increase, whereas young current-parent renters account for about2.27%. This saved-policy comparison shows a substantially stronger financing response among prospective parents than among current parents. It does not support a claim that the original down payment prevented first births: first births actually decline slightly in the n=0 group, while continuation birth-child flow for current parents changes negligibly upward.

The saved arrays lack pre-mask buyer-versus-renter option values, desired house size absent financing restrictions, purchase constraint shadow prices and a coupled individual taste-shock distribution. Consequently, neither a zero ownership probability nor a small ownership response establishes infeasibility or an unfulfilled purchase desire. More leverage also changes saving/financing preferences and continuation values. Earlier bunching at the six-room rental cap under another calibration is not a purchase-affordability measure for this winner.

[Paired group outcomes](paired_group_outcomes.csv), [source/identity and zero-solve receipt](summary.json), and [saved-data extractor](extract_saved.py) retain definitions, source receipts, weights and limitations. No new solve or model change.
