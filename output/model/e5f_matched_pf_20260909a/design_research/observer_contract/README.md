# Early housing and wealth observation contract

The approved early-target design has usable empirical values, but its model observation contract is **not yet certified for SMM**. This review preserves all identifying housing rows and gives the exact blockers rather than making the target system appear coherent through relabeling. The author approved the design on September 11, conditional on verification; the old late-target objective remains distinct.

The full machine-readable mapping is [observer_contract.json](observer_contract.json); [observer_contract.csv](observer_contract.csv) provides ten rows, including the two validation rows requested in the approved plan. JSON contains exact filters, formulas, candidate weights, source paths and SHA-256 pins, actual age masks, provenance limits, verification results and outstanding decisions. No model or target builder was edited, and no regression, model solve or cluster job ran.

## Verified empirical values

All ACS rows use survey waves 2005 and 2006, direct membership in the same 42 MET2013 city codes, and the common occupied-household-head sample. Their geography differs from the old admitted-PUMA footprint. Room outcomes are capped **before** aggregation. Ownership rows, including both parent/control groups, additionally require UNITSSTR in 3:10. The PSID wealth rows use survey waves 2003 and 2005.

| Observation | Target | Uncertainty | Role |
|---|---:|---:|---|
| Mean occupied rooms, heads 18–85, capped at nine | 5.5610973761 | 0.0883812008 | Initial restriction |
| Ownership, heads 30–55 | 0.6483340343 | 0.0206752729 | Initial restriction |
| Recent-parent minus no-resident-child ownership | 0.1628955092 | 0.0060795247 | Initial restriction; group blocker |
| Rooms, 3+ versus 1–2 resident own children | 0.3470669318 | 0.0597051546 | Initial restriction; group blocker |
| Wealth / annual gross labor earnings | 6.1458613945 | 0.3628551534 | Initial restriction |
| Old wealth/income p90 / median, ages 76–84 | 3.5159350865 | 0.3069107851 | Initial restriction; age/income blocker |
| Annual bequests / wealth | 0.0088 | 0.00044 **synthetic scale** | Maintained external restriction |
| Pooled first-birth rooms contrast, −1 to +3 | 0.7202462624 | 0.0852600513 | Preserved structural restriction |
| Ownership, heads 25–34 | 0.4311583769 | 0.0210233346 | Validation only |
| Old wealth/income median, ages 76–84 | 7.2857930019 | 0.5225311594 | Validation only |

Uncertainty is a standard error except the explicitly synthetic bequest scale. ACS errors come from 1,000 metro-resampling draws, not the official ACS survey variance. PSID aggregate wealth uses 999 reference-person bootstrap draws; old wealth uses 499 draws from the broader living reference-person population aged 65–84, followed by the exact outcome filters. These separate bootstraps do not provide joint ACS/PSID covariance. Candidate inverse marginal variance weights are saved, but there is no new complete target-and-weight fingerprint in this folder.

## What must change in the observation code

The isolated `HistoricalMomentObserver` explicitly accepts only the old twelve-row fingerprint and measures the cross section in 2023. A new initial observer must use the pre-announcement stationary solution. Reusing a 2023 statistic and replacing its empirical target value would fail the timing contract.

Let \(g(x)\) be living household mass after current housing choices and \(s(x)\) actual occupied rooms. The new ACS room statistic is \(\sum_x g(x)\min\{s(x),9\}/\sum_x g(x)\), with the relevant sample selection. For renters, use the policy at the **full income state** before income collapse; for owners, use the chosen housing product. Capping an income-averaged policy is wrong because \(\min\{E[s],9\}\ne E[\min\{s,9\}]\) in general. Existing uncapped mean/family-room statistics cannot simply be relabeled. All 17 four-year cells already span ages 18–85 under the code's interval convention, so aggregate rooms needs a geometry check, not arbitrary extra exclusions.

A common age observation rule remains necessary for narrower samples. The existing prime-age indices select labels 30, 34, …, 54, covering intervals 30–57; old wealth selects labels 78 and 82, covering 78–85; young ownership selects 26, 30 and 34, covering 26–37. These differ from data ages 30–55, 76–84 and 25–34. JSON gives explicit uniform-within-cell overlap weights as an implementable **approximation**, not an adopted solution. For example, the old sample would use half the 74–77 cell, all of 78–81, and three quarters of 82–85. Exact age-conditioned wealth/housing differences within a cell are absent from the state.

## Family groups are not exact model observables

ACS recent parent means `NCHILD>0` and `ELDCH<4`: **the oldest resident own child is under four**, so every resident own child is young. Its comparison group has `NCHILD=0`, which can include former parents. The model's independent-count `newparent_cs` instead includes every positive at-home child count, regardless of recency; its control group has lifetime parity zero. Child ages and time since first birth are not state variables. Using all dependent parents, or merely switching controls to zero dependents, does not reproduce this empirical contrast. Structure type also has no exact model analogue for the DUE restriction.

For the family-room contrast, ACS `NCHILD` counts all resident own children and `YNGCH<18` requires at least one minor. The model bins current at-home count at three; that implementation is verified, but independent binomial departure with mean 18 years does not encode actual child ages or an under-18 condition. The required resident-child mapping is therefore still an explicit approximation/observer decision. It is not solved by verifying the cutoff or by calling lifetime parity the resident count. Preserve this row: the parenthood-only utility deliberately leaves it as an informative test of larger-family housing behavior.

## Wealth and income timing

The aggregate wealth observer uses coherent beginning-period wealth \(W=b+pH\) for owners, or \(b\) for renters. The labor-earnings denominator divides the four-year after-payroll income flow by four and by \(1-\tau\); the 12 working cells cover ages 18–65. This matches the intended RP/spouse gross earnings concept, rather than total family income. Apply it to the initial balance-sheet distribution, not post-transaction housing paired with pre-transaction assets.

The old-age quantile routine correctly uses living households, not decedents; its weighted-quantile algorithm agrees with the empirical builder. However, the retirement-state denominator is model pension income, while PSID `INCFAMR` contains pensions, transfers, asset income and other family components. The model also does not impose the empirical $1,000 real-income cutoff or missing completed-child-history selection. These are additional sample/income reconciliation requirements. Wealth valuation dates versus earnings/income reference years within the 2003/2005 survey waves were not newly audited. No exact same-calendar stock/flow alignment is claimed.

The annual bequest ratio 0.0088 remains external with its inherited synthetic weight. The model uses annualized, nonnegative, post-saving estates, actual death probabilities and forced terminal death, divided by beginning-period aggregate wealth. No early PSID microdata sample or empirical standard error exists for this restriction.

## Preserved pooled PSID birth response

The primary receipt remains 0.7202462623815278 rooms with SE 0.0852600513385958. It is the Sun–Abraham \(+3-(-1)\) contrast, with omitted event time −2, person and year fixed effects, age and education controls, longitudinal IW weights, and clustering by the selected woman's stable ID. The sample is one current reference/spouse woman per single-family-unit household-year: 49,457 estimation rows and 4,112 women. Its rooms are aligned forward one observed interview and remain uncapped.

The correct stationary analogue can reuse the **dated pair** of `begin_dated_first_birth_housing_branch` and `finish_dated_first_birth_housing_branch`, supplying the stationary policy at both dates. It selects successful first births, advances equal treated/control copies one period, permits treated continuation births at the destination, and keeps the control childless; no Census age bridge enters. Do not substitute the sequential `first_birth_housing_response` legacy helper, whose treatment branch suppresses destination continuation births.

This pooled regression is not an early-only estimate. Applying it to the initial stationary economy requires stability and structural-mapping assumptions. Non-flat earlier leads, different empirical cohort weights and the documented reference-normalization sensitivity remain disclosed. The September 5 retained-target decision stands; no regression, normalization or target was replaced.

## Verification and limits

Independent checks reproduced all five early ACS points from saved metro totals and all 1,000 paired draws, their standard errors and covariance. Aggregate PSID wealth reproduced from yearly weighted totals; all early saved marginal bootstrap standard errors reproduced; old p90/median arithmetic passed. The old quantiles themselves were not reestimated from raw records. The PSID birth contrast reproduced from its two coefficients, and its SE from both marginal variances and covariance within the original table's display tolerance. Primary output hashes matched. Actual small age-index and child-bin functions were extracted and executed without importing or running the model.

All 34 source/evidence files are hashed. Large raw-source hashes are explicitly inherited from authoritative receipts; this review checked ACS size/modification time and PSID size, without rereading either multi-gigabyte source. The original recovered worker copy is preserved in the isolated checkout; this main-checkout version corrects its oldest-child wording, mistyped hash, omitted age mismatch, cap order and missing validation rows.

Certification still requires the initial observer and pure/compiled observation checks, a declared age and family/structure/income mapping, the full new weight fingerprint, and lead-owned local identification and complete target-fit/parameter tables. No target should be dropped to bypass a mapping failure, and no production policy should be launched from a partially certified contract.
