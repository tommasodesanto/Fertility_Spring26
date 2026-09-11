# Early CPS/NCHS fertility contract: verified points, explicit remaining gates

The four empirical points reproduce. **The new early target-and-weight system is not activated.** Annual CPS uncertainty is recoverable from official generalized-variance formulas; the pooled covariance convention still needs a decision. Existing model code already computes dated first-birth flows, but the scored historical observer uses a different cohort statistic. This review preserves all original samples and sets every active weight to null.

Author authorization is the September 11, 00:11 EDT entry in the main `CALIBRATION_STATUS.md`. The parenthood-only utility and initial normalization 2.1 remain as chosen. No model solve, cluster launch, empirical-builder edit, core edit, commit or push was performed here.

The complete machine-readable assessment is [fertility_target_contract.json](fertility_target_contract.json); exact numbers and raw-partition checks are in [verification.json](verification.json). These main-checkout files supersede the recovered worker's assessment under `tmp/e5f_matched_pf/output/.../fertility_contract/`, which remains unchanged. In particular, the earlier blanket claim that no defensible CPS uncertainty can be recovered was too strong.

## Empirical definitions and verification

| Candidate | Reproduced value | Population and estimator |
|---|---:|---|
| Childless at 40–44 | 0.198278751006843 | CPS June 2004/2006 women, valid live births ever had, positive supplement weight; pooled weighted zero-birth indicator |
| Exactly one among mothers at 40–44 | 0.213655325220146 | Same records; weighted count with one birth divided by weighted count with at least one |
| Mean first-birth age, midpoint coded | 25.976263860992496 | NCHS 2003–2006 first-birth records; count-weighted four-year midpoint labels |
| First births at 30+ | 0.249278013041067 | Same NCHS records; maternal age at least 30 divided by all first births |

CPS uses the original IPUMS `FREVER` and `FRSUPPWT` definitions: women ages 40–44, June 2004/2006, `FREVER` 0–20, exclusion of 999, and positive weights divided by 10,000. There are 10,872 selected records and pooled population weight 22,769,521.99. We re-read and rehashed only the original indexed partitions, totaling 71,309,637 bytes, and verified every record's year/month and the schema hash. Recomputed shares differ from the old receipt only by floating-point summation below 2e-14. The source weights/populations also reproduce the official examples' rounded 2004/2006 populations and uncapped completed fertility (11.535/11.235 million; 1.895/1.862). This is a strong aggregate check, not a record-by-record official-file crosswalk.

NCHS uses the pinned age/order cache and authoritative `build_first_birth_timing.R`: 6,611,269 first births, `lbo_rec=1`; 2003 maternal age is decoded from `mager41`, later years use `mager`. The observed age support is nominally 12–49, with the 2003 under-15 group mapped to 14. There are no fixed effects, survey weights or regression clustering in these count ratios. The current raw-file manifest is retained; no multi-gigabyte natality files were rescanned or freshly hashed in this bounded review.

**Residence caution:** the builder reads only age and birth order and applies no residence filter. Its 2006 all-order count at ages 12–49 is 4,272,728, exceeding the published all-age US-resident total 4,265,555. This proves those totals are not identical; nonresident inclusion is a plausible explanation, not a newly certified decomposition. The existing exact sample is preserved and should be described as first births in the US-named raw files, not certified US-resident counts. [Published residence convention and counts, Table 1](https://www.cdc.gov/nchs/data/nvsr/nvsr58/nvsr58_24.pdf).

## CPS uncertainty that is actually available

Census supplies the approximate generalized-variance formula
\[
s_y(p)=\sqrt{b_y p(1-p)/X_y},\qquad b_{2004}=b_{2006}=2016,
\]
where \(p\) is expressed as a proportion and \(X_y\) is its weighted population denominator. Use the fertility characteristic, not the number-of-births parameter: both scored outcomes are indicators with at most one event per woman. The denominator for exactly-one conditional on motherhood is mothers. These are official **approximate sampling SEs**, not exact replicate estimates. [CPS June 2004, Attachment 17 Formula (2), Table 3](https://www2.census.gov/programs-surveys/cps/techdocs/cpsjun04.pdf), [CPS June 2006, Attachment 16 Formula (2), Table 4](https://www2.census.gov/programs-surveys/cps/techdocs/cpsjun06.pdf).

| Moment | 2004 approximate SE | 2006 approximate SE | Pooled, zero covariance candidate | Pooled, correlation-one upper candidate |
|---|---:|---:|---:|---:|
| Childless | .005215790 | .005396663 | .003751255 | .005305036 |
| Exactly one among mothers | .006046408 | .006137682 | .004307419 | .006091130 |

The pooled candidates treat observed annual base shares as fixed. If \(\alpha=X_{2004}/(X_{2004}+X_{2006})\), their variance is
\[
\alpha^2s_{2004}^2+(1-\alpha)^2s_{2006}^2
+2\alpha(1-\alpha)\rho s_{2004}s_{2006}.
\]
The table uses \(\rho=0\) and \(\rho=1\). Two years apart removes the standard rotating household-panel overlap, but does not establish independence of shared area sampling. Random mother-denominator weights also need treatment in a full pooled ratio variance. Therefore neither pooled candidate is called a certified design SE or silently adopted. The extract lacks replicate weights and sampling strata/PSUs; a person bootstrap would not recover them. A documented generalized-variance approximation can support a diagnostic scale now, with these qualifications.

## NCHS: distinguish process precision from discrepancy scales

Complete natality registrations have no probability-sampling error conditional on the recorded population. This does **not** imply that every useful statistical SE must be zero: NCHS explicitly discusses event-process random variation under stated assumptions. [Births: Final Data for 2006, printed pages 96–99](https://www.cdc.gov/nchs/data/nvsr/nvsr57/nvsr57_07.pdf).

| Scale candidate, not adopted | Midpoint mean age | Share 30+ | Interpretation |
|---|---:|---:|---|
| Conditional multinomial process SE | .002249892 | .000168244 | Independent event-process approximation conditional on total first births |
| Sample SD across 2003–2006 annual statistics | .084567369 | .008492262 | Actual temporal dispersion, without division by two or trend removal |
| Maximum leave-one-year-out shift | .039691174 | .003709611 | Deterministic window sensitivity |
| Inherited old-cohort scale | .15 | .01 | Old rounded cohort-window difference; diagnostic borrowing only |

For the process calculation, \(\operatorname{Var}(\hat\mu)=\operatorname{Var}(A)/N\), \(\operatorname{Var}(\hat p)=p(1-p)/N\), and their covariance is also saved. This mean-age extension is our calculation, not an NCHS published formula. It captures neither registration error nor the approximate stationary model's discrepancy. The annual SD is not a sampling SE. A finite SMM scale must state which of these concepts it represents and be pinned in a named contract; no automatic inverse-variance choice has been made.

There are 71,956 records with unknown birth order in the same four-year age sample. Allowing all of them to be first births at either boundary yields conservative mean bounds [25.911919,26.170319] and share bounds [.246594,.257361]. These are sensitivity bounds, not SEs. They are deliberately conservative and do not alter the authoritative known-order sample.

## Implementable model observers for a default-off diagnostic

**Period timing already has the necessary primitive.** In the isolated checkout, `run_e5f_transition_calibration.py:876` computes first-birth flows from the pre-choice distribution, first-birth attempt probability and conception success. `period_fertility_diagnostics` at line 918 independently reconciles them with all-order birth accounting and already returns the desired flow-weighted mean/share. Reuse these functions on the newly solved stationary evaluation. Do not use `HistoricalMomentObserver`, whose final scored timing comes from `cohort_timing_moments` and a synthetic cohort ending at age 42.

For flow \(F_j\) at age-cell start \(a_j\), return \(\sum_j(a_j+2)F_j/\sum_jF_j\) and \(\sum_{a_j\ge30}F_j/\sum_jF_j\). No female-exposure denominator is needed for timing conditional on first birth; female exposure is needed for a separate all-order fertility rate. A stationary initial analogue pools repeated identical flows, which cancels from the ratios. It still assumes that the model household's reproductive member represents the maternal population.

The empirical bin map is ages ≤21→20, 22–25→24, …, 42+→44. It moves **7.73136% of early first births below age 18** into the first model cell; the above-45 share is only .047585%. Preserve this data operator, but disclose that it represents those births by boundary decisions, not actual pre-18 motherhood. This matters for interpreting an age-18 entrant initialized childless.

**CPS needs a declared within-cell stock projection.** Integer ages 40–44 correspond to [40,45), overlapping cells [38,42) for two years and [42,46) for three. Two concrete diagnostic alternatives are:

1. **Constant cell stock.** Use parity mass from a specified phase of each cell, with weights .5 on start-age 38 and .75 on start-age 42. Then calculate the ratio from aggregated numerator/denominator masses. This is simple but treats the cell's stock as constant despite births inside the interval.
2. **Uniform birth-time parity interpolation (preferred diagnostic).** Let \(q_j^-\) and \(q_j^+\) be the normalized parity distributions before and after that cell's fertility step. Assume the one possible parity transition occurs uniformly within its four-year interval. At elapsed fraction \(u\), set \(q_j(u)=(1-u)q_j^-+u q_j^+\). The [40,42) subwindow uses mean \(u=.75\); [42,45) uses \(u=.375\). Combine these with age-cell population masses times .5 and .75, then compute \(p_0\) and \(p_1/(1-p_0)\). Validate equal pre/post age mass and one-step parity flows first. This preserves the birth midpoint convention under the stated interpolation; it is an observation approximation, not an additional within-period optimization.

To separate age composition from interpolation, a supplemental option replaces projected model age weights by the fixed CPS subgroup weights: ages 40–41 carry 0.403248458 of all women (calculated from the saved totals), and ages 42–44 the remainder. An even finer version uses each exact age's saved empirical weight and midpoint fraction within its model cell. This holds the population age mix at the sample's mix; it must not be relabeled as the endogenous model cross-section. The exact subgroup share is provided by the saved totals rather than the displayed rounding.

Pure tests should verify: constant parity across ages gives unchanged moments under every projection; zero birth flow makes pre/post projections identical; a unit first-birth transition gives the stated elapsed-fraction zero/one shares; continuation births affect exactly-one in the correct direction; aggregated one-child numerator uses mothers as denominator; age labels 20/24/…/44 reproduce the data's ≥30 classification. Nonpositive denominators must fail, not become artificial zeros. Compare both projections at each retained diagnostic candidate to display approximation sensitivity.

For dated work, the existing bridge allocates decision-date \(t\) births to \(t+1,\ldots,t+4\); 2020–2023 thus maps to decision 2019. Calendar allocation and maternal-age reconstruction must be checked jointly before calling this an exact annual observer. A label change is insufficient. The static interpolation above is a bounded initial diagnostic, not a certified history reconstruction.

## Normalization and activation boundary

The initial 2.1 is the author's stationary model completed-fertility normalization. The inherited literal states are 0,1,2,3+, with representative top-bin count 3.602359422009; it is not the actual female period TFR. The early CPS capped top-bin count is separately recovered in the receipt. The birth-to-entry conversion and queues require their own accounting check; this note does not modify them.

If 2.1, the two CPS shares and the inherited top-bin weight were applied to exactly the same parity distribution, their algebra would require a 3+ share .416791, versus the observed .286255. This is feasible but a substantial high-parity discrepancy. It is an illustrative conditional calculation; actual model age projections differ. Preserve the full parity distribution as a diagnostic and do not describe exactly-one conditional on motherhood as a second-birth hazard.

Activation still requires: a tested early observer; selected age/population approximation; a named pooled CPS variance convention; an explicit NCHS scale and sample/tail declaration; separate 2.1 population identities; and the complete target/weight fingerprint with identification checks. None of these findings authorizes dropping a moment. Documented provisional weights and both stock projections can support useful default-off diagnostic panels; they cannot become certified SMM or production policies automatically.

## Reproduce this bounded verification

From the main repository run:

```bash
python3 output/model/e5f_matched_pf_20260909a/design_research/fertility_contract/verify_fertility_contract.py
```

Observed runtime was about .52 seconds. This verifies selected CPS source bytes, all four points, annual/sample/tail arithmetic, and saved provenance without importing the model. JSON and all active-weight-null assertions were checked separately. Web text verified the official formula/table inputs; PDF screenshots were unavailable from the browser cache, so no visual PDF check is claimed.
