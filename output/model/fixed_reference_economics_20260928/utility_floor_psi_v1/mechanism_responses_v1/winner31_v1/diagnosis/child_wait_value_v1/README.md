# First-child versus postponement: saved-value diagnostic

Reference: frozen winner31, q=.719168368828958, common inherited PRE population. Saved arms: baseline buyer/stayer financed shares.8/.8, purchase-only1/.8, and both1/1. Housing floor2.3, menu[2,4,6,8,10], all preferences/entry/units/targets/grids remain fixed; only both changes reported global financed share. No CES, calibration, GE or new lifecycle solve. Existing complete14/31/17 packets: [baseline](../../purchase_ltv_v1/local_run/retry5/results/baseline_80_80/), [purchase-only](../../purchase_ltv_v1/local_run/retry8/results/purchase_100_stayer_80/), [both](../../purchase_ltv_v1/local_run/retry9/results/both_100_100/).

## Exact recovered values

Let W denote the value of postponing a first birth now, including optimal later fertility, housing and saving. Let C denote the successful-child destination value net of the unchanged first-birth fixed cost. Let A denote the value of attempting now, pi its conception success probability, kappa the first-birth taste-shock scale, and F the saved fertility inclusive value. The executed source uses

\[
 A=\pi C+(1-\pi)W,\qquad F=\kappa\log(e^{W/\kappa}+e^{A/\kappa}).
\]

On feasible interior two-choice states, saved probabilities recover both branch levels:

\[
 W=F+\kappa\log p_W,\quad A=F+\kappa\log p_A,
 \quad G=A-W=\kappa\log(p_A/p_W),\quad C=W+G/\pi.
\]

This is attempt versus postpone, not child versus never-child. There is no Euler/menu normalization constant in the executed logsumexp. kappa=.13056430591307258; first-birth fixed cost=.4009349064519595. pi declines with age, from.98 at18 to.501436 at42, and is unchanged across arms. Subsequent births use another scale and are not inverted here. The actual newborn-exempt VI_ex branch is inactive, so saved child-success housing policies correspond to the successful first-child value branch.

Source identity: original household SHA2a34f5f28c0c63ca3d24d9e759cd1b5b92aadaddb5fe04795e65abbc0b448082 matches both saved source receipts; buyer/death-floor override SHAa13e9b3944bc69e0a5229191a8f8133e8542d60f4464d4eaed51a0f9220abdb6 matches relaxed receipts. Their first-birth blocks are exactly identical. Each arm's31 effective parameters matches the authenticated candidate, allowing only financed_share=1 in both1/1. Prices and common PRE hashes agree.

## Matched results

Population: inherited renters with no children ever born, model age-cell starts18–42. Total mass.1863762. Interior value coverage is98.9161%; excluded mass.00202013 has exact saved attempt probability zero in all arms. Its actual first-birth change is zero. No probabilities are clipped and the original population is not renormalized for birth accounting. Value means condition on the common interior; actual birth effects retain all original group weights.

| Relaxation relative to baseline | Purchase-only1/.8 | Both1/1 |
|---|---:|---:|
| Mean postponement value gain |+.00247302|+.00745334|
| Mean successful-child net-cost value gain |+.00122300|+.00396818|
| Mean attempt-minus-postpone gap change |−.00111169|−.00315503|
| Birth-sensitive weighted gap change |−.00025414|−.00068907|
| Conditional postpone ownership change, pp |+3.52212|+9.56526|
| Conditional successful-child ownership change, pp |+.61909|+1.27062|
| First-birth probability change, pp |−.01914|−.05168|

Thus both branches gain, but postponement gains more on average. Financing relief raises expected ownership much more when postponing than after a successful first child. This is consistent with the relative value of waiting rising; it is not a decomposition into immediate affordability versus future continuation channels. Both are already inside each reoptimized branch value, and the means alone do not explain aggregate changes.

Exact first-birth accounting is Delta B=sum g*pi*(p_A,new−p_A,base). Purchase-only has positive contributions+1.54430e−5 and negative contributions−5.11191e−5, net−3.56761e−5 of the whole population; both has+2.99237e−5 and−1.26244e−4, net−9.63206e−5. Birth-sensitive weights g*pi*p_A*(1−p_A)/kappa produce linear gap approximations−3.58413e−5 and−9.71784e−5, close to the exact net changes. Purchase-only declines most at age-cell starts22 and18; both declines most at22,18 and26. The baseline attempt-probability5–25% group accounts for the largest negative contribution, rather than the states with effectively zero birth probability dominating mean gaps.

The conditional baseline WAIT renter policy reaches the six-room cap for10.42% of interior weight; the successful-child conditional renter policy reaches it for35.68%. These are conditional renter policies even when ownership might be chosen. The WAIT-cap subset has a small positive birth response (+9.6865e−7 purchase-only; +3.4337e−6 both). The negative aggregate response is concentrated below that cap. This is rental-cap pressure, not proof of a down-payment-constrained desired purchase. Purchase affordability and pre-mask desired owner values are not measured here.

## Artifacts and limits

[Summary/source hashes/invariant checks](summary.json), [weighted value and probability distributions](weighted_distributions.csv), [age contributions](by_age.csv), [conditional-cap/probability/ownership-response splits](by_baseline_pressure.csv), [extractor](extract_value.py). The positive expected-ownership-change-weighted gap changes are−.00201099/−.00614142; these weight expected marginal changes, not individually blocked purchaser counts. Conditional ownership is taken from saved tenure probabilities on the same renter-origin wait/success states; the single-location probability is checked for both branches. Earlier native saved-policy replay independently reproduces matched total ownership/birth outcomes. Supplemental plot uses these existing age contributions; no target clock is changed.

Pre-mask desired house values and constraint shadow prices remain unavailable. A total branch-value comparison cannot isolate immediate purchase constraints from later saving/ownership/fertility opportunities. It nevertheless completes the requested child-versus-postpone comparison: credit benefits both choices, benefits postponement more, and the exact birth response is modestly negative in this matched population. No new solve is needed or authorized by this packet.
