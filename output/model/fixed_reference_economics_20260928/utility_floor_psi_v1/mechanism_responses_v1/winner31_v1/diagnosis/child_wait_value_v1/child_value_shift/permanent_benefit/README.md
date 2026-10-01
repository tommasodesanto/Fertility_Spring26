# Permanent child-benefit sensitivity at fixed prices

Frozen winner31, q=.719168368828958, original120×9 grid and identical inherited PRE population. Four one-core native backward-policy calls permanently change psi_child to.10 or.25, under baseline financed share.8 or both buyer/stayer share1. continuation_V is OMITTED: all future household policies and values are fully reoptimized under the respective permanent benefit/financing treatment. The original psi=.171561928 comparison reuses certified saved policies; no fifth call.

Four calls completed in14.03 numerical seconds (16.50 seconds including authentication/extraction), within120 seconds. No KFE, new stationary population, GE, search or calibration. Current and future benefit is psi*m^(1−curvature), with unchanged curvature; no benefit normalization. Original financial, possible-death estate and housing-floor constraints remain active.

Matched inherited never-parent renters at model age-cell starts18–42, total mass.1863762:

| Permanent child benefit | Baseline first-birth probability | Expanded-credit probability | Credit effect, pp |
|---|---:|---:|---:|
|.10|11.1013%|11.0817%|−.01957|
|.171561928, saved winner|22.2374%|22.1857%|−.05168|
|.25|34.6777%|34.6131%|−.06464|

Higher permanent benefit raises the birth level materially, but financing relief remains mildly negative at every tested benefit. Raising this parameter therefore does not reverse the response over these three points at the frozen price/common PRE. This does not establish global impossibility, settle GE, or cover other parameter combinations. Mean credit-induced attempt gaps become less negative as psi rises (−.003888,−.003155,−.001958), while birth-sensitive weighted changes become more negative (−.000382,−.000689,−.000790); the nonlinear probability response depends on which states are near the birth-choice margin, not only the mean value gap.

For psi=.10/.25, positive first-birth-change contributions are+.00002287/+.00002612 of the full population and negative contributions−.00005934/−.00014658. Net changes are−.00003647/−.00012046. Common interior value coverage remains98.9161%; excluded states contribute exactly zero actual birth change. Original PRE weights are retained for birth probabilities; value means use common interior states without invented logit levels.

Every call authenticates31 effective scalar parameters. Baseline variants change psi_child and its derived child_benefit_CRRA_coefficient=(1−curvature)*psi, leaving29 unchanged. Credit variants additionally change financed_share, leaving28 unchanged. All other preferences, first-birth cost, taste-shock scales, rooms/owner premium, earnings, entry, fiscal quantities, targets/weights and prices remain fixed. Source identities match the reviewed original household and estate-patched credit module. No new forward lifecycle/estate ledger or17 diagnostic packet is claimed for this policy-only sensitivity.

[Four new policy-call rows with wait/child values and attempt gaps](four_policy_calls.csv), [three-point credit comparison and positive/negative contributions](comparison_with_saved_reference.csv), [full31-parameter/source/runtime receipt](summary.json), [driver](run_sensitivity.py). Existing certified14-target/31-parameter/17-plot packets: [baseline](../../../../purchase_ltv_v1/local_run/retry5/results/baseline_80_80/), [expanded credit](../../../../purchase_ltv_v1/local_run/retry9/results/both_100_100/). No further run authorized by this packet.
