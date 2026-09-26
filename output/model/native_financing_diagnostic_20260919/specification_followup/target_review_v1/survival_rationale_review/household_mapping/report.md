# A 2007 age-specific household-exit proxy and fixed-distribution exposure

**Question.** What does a defensible age-specific adult/household exit schedule look like when pre-retirement mortality is allowed, and how many model households carrying a dependent-child state are exposed to an adult death under the selected saved distribution? This is a bounded demographic mapping and fixed-distribution accounting exercise. It does not change or solve the model.

## Finding

The official 2007 U.S. life table implies appreciable adult death well before retirement: under an equal male/female count illustration, the 18-to-66 death risk is 16.62%. An age-specific household-exit proxy that lets a household persist after its first adult death is much lower than the corresponding first-adult-death exposure. At ages 62–66, the SCF-weighted four-year rates are **8.05% for at least one adult death** and **2.02% for household exit**. Treating first death as dissolution would therefore overstate exits by about four times in that bin.

Applying the candidate exit schedule to the selected checkpoint's fixed age-18 entrant mass lowers total adult household stock by 1.62%; the working-to-retiree mass ratio moves from 2.85805 to 2.84817 (−0.35%). This is a small change in that particular ratio, but it does not establish that mortality is economically immaterial: the calculation freezes household choices, child-state composition, saving, fertility, births, and the entry law. It measures the mechanical age-mass response only.

The checkpoint contains positive mass in a model-dependent-child state from age 22 onward. Under the candidate schedule, the normalized expected mass of those households experiencing a **nonterminal** first adult death is 0.00971 per four-year model period; the mass with a nonterminal household exit is 0.00303. These are dependent-state household events, not child-unit events. For coupled households, the first adult death does not imply household exit when the second adult survives. The age-82 terminal row is outside these hazard totals: it has positive candidate (cs=1) household mass (0.002890), which exits the model age support through its terminal boundary. The reported mortality exposure therefore excludes that terminal exit and is not total lifecycle exit exposure.

## Source and household mapping

I combine the NCHS sex-specific 2007 exact-age one-year death probabilities (q_{x,s}) with the weighted 2007 Survey of Consumer Finances (SCF) cross-sectional distribution of reference-person ages, sex, partner ages/sex, and single/partnered status. The mortality CSV has 101 rows for ages 0–99 and 100+, retains source spreadsheet precision, and was visually checked against the printed NCHS tables. The official SCF public file contains 4,417 primary economic units (PEUs), represented in five implicates; I pool rows with the official simple-estimate weight (X42001/5). Head age, sex, and living arrangement are invariant across implicates; partner age (X19) can vary and is retained by implicate. “Partnered” includes married and cohabiting resident partners under SCF codes.

For each person of sex (s) and age (a), the four-year death probability is computed from the annual life-table risks as

\[
q^{(4)}_{a,s}=1-\prod_{k=0}^{3}(1-q_{a+k,s}).
\]

For a single-adult household, exit is the reference adult's (q^{(4)}). For a two-adult household with reference adult hazard (q_h) and partner hazard (q_p), first-adult-death risk is (1-(1-q_h)(1-q_p)), while household exit is (q_hq_p). The last-survivor formula assumes independent spouse death times conditional on ages and sexes and treats the household as continuing after the first death. It is a transparent approximation, not a source-prescribed joint-death law. I average these single and coupled probabilities by SCF weights within four-year reference-age bins. I apply the age-specific SCF partnered share equally to all model households at that head age, including (cs=1) dependent-child households; this assumes child households have the same couple mix as all households of that age. Parent/couple status is not cross-tabulated in the checkpoint and the resulting child-household exposure is therefore approximate. Partnered shares range from 28.3% at ages 18–22 to 61.3% at 62–66 and 34.0% at 78–82.

This is a one-state, age-mixed cross-sectional schedule: it does not update a couple into a single/widowed survivor after first death, remap the surviving adult's age/earnings, or account for remarriage. Independence can understate joint risk if spouses share fatal shocks; the candidate is also sensitive to the household-versus-person unit. Public SCF partner age is top-coded at 95. Mapping 95+ to (q_{95}), (q_{99}), or the terminal (q_{100+}) does not change any of the displayed four-year bins in this sample; the exact qx95/qx99/qx100+ schedule sensitivities are saved in `demographic_profiles.json`.

| Reference ages | Partnered share | First adult death, q4 | Household exit, q4 |
|---|---:|---:|---:|
| 18–22 | 28.27% | 0.519% | 0.286% |
| 22–26 | 49.15% | 0.572% | 0.188% |
| 26–30 | 61.91% | 0.643% | 0.142% |
| 38–42 | 63.76% | 1.299% | 0.282% |
| 42–46 | 61.15% | 1.779% | 0.445% |
| 50–54 | 64.11% | 3.364% | 0.762% |
| 58–62 | 62.93% | 5.905% | 1.411% |
| 62–66 | 61.26% | 8.045% | 2.022% |
| 66–70 | 58.27% | 10.350% | 3.123% |
| 70–74 | 61.61% | 14.877% | 4.161% |
| 74–78 | 48.83% | 19.848% | 8.094% |
| 78–82 | 34.00% | 27.496% | 15.739% |

All 16 four-year bins, single-only, partnered last-survivor, head-death and top-code schedules are in [candidate_schedule.csv](candidate_schedule.csv) and [demographic_profiles.json](demographic_profiles.json). The empirical inputs and provenance are in the upstream [mortality receipt](../../overnight/bequest_flow_2007/mortality_source_receipt.json) and [SCF calculation receipt](../../overnight/bequest_flow_2007/calculation/candidate_receipt.json). Primary sources: NCHS, *United States Life Tables, 2007* ([report](https://stacks.cdc.gov/view/cdc/221925), [official data tables](https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/59_09/Table02.xls)); Federal Reserve, [2007 SCF files and codebook](https://www.federalreserve.gov/econres/scf_2007.htm).

## Age-profile comparison and normalization

The selected serialized checkpoint uses four-year periods, age 18 entry, (J=17) age cells, and (J_R=12). Its survival switch is on; conditional survival is one for ages 18–62 and uses the saved combined-sex 2023 post-retirement probabilities for 66–70 through 78–82. The candidate replaces all 16 transition hazards with the 2007 cross-sectional household-exit schedule. I reweight the checkpoint's age marginals by the ratio of candidate to saved survival along each age path, holding age-18 entry mass and within-age state shares fixed. The reweighted marginals match the checkpoint's compact lifecycle output before reweighting to (1.91\times10^{-14}) maximum absolute error.

| Fixed age-18 entrant mass | Current saved schedule | Candidate 2007 household exits | Change |
|---|---:|---:|---:|
| Working household mass, ages 18–62 cells | 0.74080 | 0.72811 | −1.71% |
| Retiree household mass, ages 66–82 cells | 0.25920 | 0.25564 | −1.37% |
| Total age-cell household stock | 1.00000 | 0.98375 | −1.62% |
| Working / retiree mass | 2.85805 | 2.84817 | −0.35% |
| Mass at age-66 cell | 0.06173 | 0.05720 | −7.35% |
| Mass at age-82 terminal cell | 0.03912 | 0.04112 | +5.14% |

This fixed-entry result is not a closed stationary population: no births, entrant scale, fertility-weighted birth flow, death replacement, or age-state policies were recomputed. In the selected checkpoint, the saved birth-based adult entry is 0.061733465 and the normalized entrant flow is 0.061733456. The adopted birth-to-household conversion (1/2.1=0.4761905) is a one-time conversion; it cannot simply be reused as an all-age death-replacement rule. Once mortality is introduced, entry must be reconciled to births and deaths under a stated cohort/period timing and dependency-support law. Otherwise only a fixed-entry sensitivity such as this one is defined.

For context, the separate one-person upper-bound schedules in `demographic_profiles.json` show why household unit matters: applying reference-adult death to every household gives a higher W/R ratio (3.656), and treating first adult death as exit gives 4.138, because these alternatives remove more post-working-age mass. They are bounding mappings, not preferred household laws. Neither is a preferred alternative to the last-survivor candidate; the current model has a single household state and no spouse transition, so neither mapping is fully represented.

## Dependent-child-state exposure in the selected checkpoint

The checkpoint is immutable and identified by SHA256 `d5ef71bdaf9960273035c722a2428a55f14bab160e0596881c8a981e67b8ead1`; its frozen-source manifest SHA256 is `237904131d159f775c7ae89d1bbf1e8d1dd70c79ad012c80a36d948658d6f9c6`. The calculation read `stationary_g_pre` on the Torch cluster from the compact, matching estate-control packet; it did not download a giant checkpoint locally. It verified source-manifest identity, selected reference identity, serialized survival settings, and age marginals against the control lifecycle CSV. The saved tensor has separate parity and child-state axes (sizes 4 and 4 after the income axis); local/frozen parameter source defaults (n_{parity}=3), so the compact exposure receipt does not establish the serialized top-bin override or weights.

Under the current selected distribution, (cs=1) mass is zero at the age-18 entry cell and about 1.2–2.4% of normalized mass in ages 22–42 cells, then declines with age. The model keeps children-ever-born category (n) separate from child stage (cs). With one dependent stage ((K=1)), (cs=1) is the active dependent state and (n) supplies the model child-count category used in child costs and housing; the joint ((n,cs)) state therefore recovers model-coded (m), subject to the saved top-bin convention/weights. The exposure calculation sums over (n), so the results below count households, not (m) child units. It does not encode individual child ages or biological-minor status. Holding within-age state shares fixed and reweighting by candidate survival gives:

- Expected (cs=1) households experiencing a nonterminal first-adult death in the next four-year period: **0.009711 per normalized household-mass unit**.
- Expected (cs=1) households exiting through a nonterminal mortality transition: **0.003030**.
- Current-schedule (cs=1) household nonterminal mortality exits for comparison: **0.002701**.
- Candidate (cs=1) household mass at the age-82 terminal cell: **0.002890**; terminal age-support exit is separate from the mortality-event totals above.

The first-adult-death event is over three times the nonterminal household-exit event because a first death in a surviving couple is not household dissolution. The age-82 terminal exit is additional to the reported hazard exposure and includes positive dependent-state mass. The parity-by-stage state can recover model (m) in principle, but this saved age-only aggregation has already summed out (n); computing child-unit exposures would require the joint ((n,cs)) mass, the serialized top-bin mapping/weights, and the same event hazards. The compact receipt does not preserve those cross-tabs, and this task did not repeat checkpoint computation. Actual biological child ages, minor status, and caregiver/orphan outcomes remain unidentified even if model (m) is recovered. The SCF age data are cross-sectional and are not linked to the model households.

Age-specific dependent-state mass and exposure are in [checkpoint_age_exposure.csv](checkpoint_age_exposure.csv); full calculation and exact checkpoint/source/target identity are in [checkpoint_exposure.json](checkpoint_exposure.json). This is a fixed-distribution exposure diagnostic: no household or equilibrium solve, policy/fertility reoptimization, widow succession, or entry-law revision was run.

## Accounting implications and assessment

This mapping confirms that first adult death, household exit, and orphan event are three distinct objects. The model currently treats survival as household attrition; it has no spouse state, remarriage or widow's age/earnings succession, child support/caregiver state, or age-specific children-at-home roster. A reduced-form household hazard can be implemented without those optional detailed states, but its exit unit and the dependent-state/entry accounting must be internally consistent. For estates, a first spouse death would also require a rule for asset ownership and transfer; at last-survivor exit, housing, debt, and estate timing need explicit handling. For PAYGO, age-specific exits alter contributor and beneficiary stocks, and working-age exits also change the entry/birth denominator. The close W/R ratio here is only a frozen-entry accounting fact, not a validation of the fiscal closure.

Under the author’s stated acceptance of a reduced-form mapping, this candidate admits positive pre-retirement deaths, distinguishes first-adult death from last-survivor exit, and uses observed age-specific single/couple shares. It can be considered as a reduced-form household-exit law without adding a full widow or genealogical state; that choice does not make its cross-sectional marital-mix/independence assumptions empirical facts. Before adopting it, the author still needs to state the unit of exit, how household mortality interacts with dependent-state and birth-entry accounting, the terminal boundary treatment, estate/debt handling, and survivor-weighted PAYGO closure. Detailed spouse succession and child-level support are optional richer states if the research question requires them, not prerequisites to a reduced-form attrition law. The existing certain-survival-before-66 restriction remains a transparent model boundary, not a claim that deaths are empirically negligible.

## Reproduction and scope

Run `python3 build_household_exit_schedule.py` locally to recreate the SCF-weighted schedule and fixed-entry profile using only the CSV inputs. `measure_checkpoint_exposure.py` is run on Torch against the immutable compact control packet; its Slurm script is retained as `submit_exposure_job.sh`. Job 18566254 completed successfully in 10 seconds after a single concrete `PYTHONPATH` correction to include the frozen `source/code/model/tools` directory; the prior job failed before loading the checkpoint and produced no result. No automatic retry chain was used. No model/equilibrium solve, source change, calibration, checkpoint download to the Mac, or adoption/status edit was made.
