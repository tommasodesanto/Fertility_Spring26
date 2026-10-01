# Current-period versus permanent financing diagnostic

Frozen winner31, q=.719168368828958 and identical inherited PRE population. This experiment relaxes current buyer and incumbent-owner financing from.8 to1 for ONE FOUR-YEAR MODEL PERIOD, while anticipating the baseline's future financing rules and values. Permanent both1/1 uses the existing saved checkpoint. Preferences, children costs/floors, housing menu/units, entry, survival, income transitions, targets and grids are unchanged. It is a policy-horizon diagnostic, not a stationary calibration or market-clearing experiment.

## Controlled native evaluations

The existing native household hook accepts continuation_V. At every age j, the temporary call evaluates current choices against saved baseline V at j+1, using unchanged native income-transition expectations, survival, bequests and child aging. Terminal bequest handling is unchanged. The exact baseline continuation tensor—including200097 native DEAD nodes—is preserved byte-for-byte (SHA9277a9de998b3c57a9757b723549efc82050f0d262ec1680636fd5b311a30ddb). No future support is projected or relaxed, no new bailout/support operator is added, and the current possible-death buyer estate floor remains active through the previously reviewed estate-patched household module. Thus borrowing today faces continuation values under the future baseline limits.

Two one-core native POLICY evaluations completed in9.63 seconds total (3.64/3.59 seconds evaluation). The baseline-P call with saved baseline continuation reproduced the FULL saved V, fert_probs and fert_value arrays with maximum errors EXACTLY ZERO. The second uses current P.phi=1 and no stale purchase override, with the same baseline continuation tensor. No forward lifecycle/KFE, equilibrium root, optimizer or additional numerical run was performed. An independent occupied-state death-estate ledger was not regenerated for this policy-only packet; original native constraints and continuation/interpolation were retained. No full support certification beyond the retained native contract is claimed.

## Matched first-birth and housing results

Same inherited never-parent renter states at model age-cell starts18–42; total PRE mass.1863762. Common interior weight.1843560 (98.9161% coverage). Excluded mass has saved attempt probability zero and contributes exactly zero actual birth change; no population dropping or probability clipping. Values are reconstructed from the exact two-choice logit as documented in[the saved-value packet](../README.md): W=F+kappa log(p_wait), A=F+kappa log(p_try), C=W+(A−W)/pi. W postpones a first birth and preserves the option to have children later; C is successful-child value net of the fixed first-birth cost.

| Relative to baseline | Temporary current1/1, baseline future | Permanent1/1 |
|---|---:|---:|
| Mean postponement value gain |.00052521|.00745334|
| Mean successful-child net-cost gain |.00009776|.00396818|
| Mean attempt-minus-postpone gap change |−.00035852|−.00315503|
| Birth-sensitive weighted gap change |−.00009748|−.00068907|
| First-birth probability change, pp |−.00735|−.05168|
| Positive birth-change mass |+.00001433|+.00002992|
| Negative birth-change mass |−.00002803|−.00012624|

Current financing relaxation already benefits postponement more than having a first child now and yields a small negative birth response. The permanent policy yields a larger change. Permanent minus temporary includes future financing opportunities, expectations, induced future-state choices and duration interactions; it is not an additive causal decomposition into primitive channels. Holding future value functions fixed does not hold actual future housing or savings states fixed.

Conditional ownership probability changes by owned rooms, percentage points on common interior PRE:

| Branch/arm |2 rooms|4 rooms|6 rooms|8 rooms|10 rooms|
|---|---:|---:|---:|---:|---:|
| Postpone, temporary |+2.0063|+1.2566|+.2629|+.0896|+.0052|
| First child, temporary |0|+.0043|+.2515|+.3414|+.0244|
| Postpone, permanent |+7.3227|+1.6554|+.4814|+.0966|+.0092|
| First child, permanent |0|+.0094|+.8676|+.3635|+.0301|

Two-room purchases account for76.55% of the permanent increase in conditional postponement ownership. Two-room ownership is exactly zero in the first-child branch: the physical parent floor2.3 rejects H=2 before applying the owner service premium. This identifies a concrete asymmetry in current purchase opportunities, not the fraction of individual households who wanted to buy but were blocked.

## Zero-solve current two-room choice ablation

As a diagnostic only, remove the two-room OWNED option from the CURRENT decision under both saved baseline and permanent-credit branches, while preserving their ORIGINAL future menus/values. The tenure logit is unnormalized. If p2 is the conditional two-room probability, the exact theoretical identity gives

\[
 \Delta W=\kappa_{tenure}\log(1-p_2),\quad C'=C,\quad
 G'=G-\pi\Delta W.
\]

The updated try probability uses the corresponding finite odds multiplier; no implicit normalization of utility is introduced. Saved tenure probabilities are float32, so numerical results inherit that finite precision. Positive-weight p2 never equals1 (max.391456 baseline/.610166 permanent), child p2 is exactly0, and the largest odds multiplier is1.08137. Zero saved attempt probabilities remain zero; no branch levels are invented for underflow states. Unoccupied cells are masked only to avoid undefined0/0, retaining all positive PRE weight.

The permanent-minus-baseline first-birth change moves from−.05168pp to−.04637pp after this CURRENT-option ablation: about10.27% attenuation, still negative. Therefore two-room availability explains a meaningful current-choice contribution but does not account for most of this matched birth response, despite dominating the ownership increase. This is not a proposed deletion of two-room housing, a unit correction, full-menu counterfactual, or recalibration.

## Retained artifacts

[Driver](run_temporary.py), [control/source/continuation receipt](summary.json), [value/birth comparison](value_comparison.csv), [conditional ownership by size](ownership_by_size.csv), [age contributions](birth_by_age.csv), [zero-solve ablation source](ablate_current_two.py), [ablation receipt](two_room_ablation.json). Source hashes match the reviewed original and buyer/death-floor module. Existing full14-target/31-parameter/17-plot packets: [baseline](../../../purchase_ltv_v1/local_run/retry5/results/baseline_80_80/), [permanent both](../../../purchase_ltv_v1/local_run/retry9/results/both_100_100/). The present diagnostic changes policy duration only; no preferences, targets, model-unit compensation or weights changed. Two numerical policy calls exhausted this task's allowance; no further run authorized.
