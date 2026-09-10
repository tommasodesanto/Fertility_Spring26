# Initial-period wealth target feasibility

Completed September 10, 2026. Diagnostic empirical construction only; no active target, sample builder, parameter or model code changed.

**The initial wealth block can be measured today.** The raw PSID is local. One selected-column pass now delivers the exact existing wealth definitions over 2005–2007 and 2003–2007, including fresh person-clustered uncertainty. The full-vintage estimates and uncertainty reproduce the authoritative builders, providing a direct definition check. A comparable old-age wealth target is therefore not a data-availability obstacle to an initial-period calibration.

**Announcement-timing refinement:** If the preference path is announced at the start of 2007, the 2007 survey wave should not be used as strictly pre-announcement evidence. A separately authorized second selected-column pass adds 2003/2005 and 2005-only windows below; it reproduces every previous result and bootstrap-draw row exactly. This does not select a target window.

## Constructed observations

Standard errors are person-cluster bootstrap standard errors; brackets show bootstrap percentile 95% intervals. Pooled observations retain original survey weights, and persons retain all of their included family-year observations when resampled. These are new initial-window empirical estimates, not approved calibration targets or survey-design standard errors accounting for all PSID sampling strata and PSUs.

| Window | Aggregate wealth / annual gross labor earnings | Old-age wealth/income p90 / median | Old-age wealth/income median |
|---|---:|---:|---:|
| **Pre-announcement: 2003 and 2005** | **6.145861**, SE **0.362855**, [5.469085, 6.908434] | **3.515935**, SE **0.306911**, [2.986597, 4.156280] | **7.285793**, SE **0.522531**, [6.418606, 9.035234] |
| **Pre-announcement: 2005 only** | **6.452448**, SE **0.394280**, [5.713560, 7.231378] | **3.380016**, SE **0.464589**, [2.625279, 4.583395] | **7.876782**, SE **0.827502**, [6.924940, 9.514655] |
| 2005 and 2007 | **6.926584**, SE **0.417310**, [6.143310, 7.774446] | **3.369136**, SE **0.395227**, [2.770989, 4.312729] | **7.991506**, SE **0.702557**, [6.890337, 9.514655] |
| 2003, 2005 and 2007 | **6.580858**, SE **0.361876**, [5.893710, 7.331343] | **3.571120**, SE **0.285920**, [3.081571, 4.191757] | **7.462511**, SE **0.535265**, [6.533799, 8.416353] |
| 2007 only | 7.364028, SE 0.513215, [6.424380, 8.461104] | 3.406323, SE 0.513788, [2.801671, 4.705160] | 8.047822, SE 0.826639, [6.386608, 9.599205] |
| Existing long pool | 6.873077, SE 0.398836 (2005–2019) | 3.448111, SE 0.132477 (1984–2019) | 6.501316, SE 0.231958 (1984–2019) |

The old-age median is a **supplemental empirical object**, already defined and constructed by the authoritative old-wealth builder; it is not currently one of the twelve hard targets.

| Window | Aggregate wealth: family-years / unique people | Gross earnings: family-years / unique people | Old wealth statistic: family-years / unique people |
|---|---:|---:|---:|
| **Pre-announcement: 2003–2005** | **10,828 / 6,018** | **9,265 / 5,215** | **644 / 416** |
| **Pre-announcement: 2005** | **5,565 / 5,565** | **4,769 / 4,769** | **350 / 350** |
| 2005–2007 | 11,324 / 6,349 | 9,750 / 5,521 | 670 / 425 |
| 2003–2007 | 16,587 / 6,714 | 14,246 / 5,885 | 964 / 490 |
| 2007 | 5,759 / 5,759 | 4,981 / 4,981 | 320 / 320 |

## Exact definitions and source provenance

Raw source: `/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta`. The diagnostic reads ten columns and 3,533,123 rows once. The pass and all bootstraps finished in approximately 20 seconds. Raw file size and modification time, exact selected columns and package versions are saved in `wealth/run_metadata.csv`; a full hash of the 5.9 GB raw file was not computed.

**Aggregate ratio:** weighted total family net worth `NETWORTHR` among reference persons (`RELTOHEAD_ == 10`) ages 18–85, divided by weighted total reference-person/spouse gross labor earnings `EARNINDRRC` for reference persons ages 18–65. Require finite positive longitudinal individual weights `IW`, finite net worth, and finite nonnegative gross labor earnings at ages 18–65. The denominator is not total family income. The ratio is the ratio of pooled weighted totals, not the mean of household ratios or the mean of yearly aggregate ratios. This precisely retains `code/data/psid_followup_mar2026/audit_aggregate_wealth_earnings_ratio.R` apart from the explicitly varied year window. Bootstrap: 999 draws, seed 20260723, resampling unique reference-person IDs.

The date labels identify PSID **survey waves**. Detailed within-wave wealth dates and the income variables' reference-year metadata are inherited from the existing shelf and were not newly audited here. A final calendar-specific measurement contract should retain that distinction.

**Old wealth distribution:** living reference persons ages 76–84 with finite positive `IW`, observed completed-child count `RELCHINUM`, finite `NETWORTHR / INCFAMR`, and annual family income `INCFAMR > 1000`. Living means `DEATHYEAR` missing or survey year no later than death year. Both wealth and family income retain the shelf's existing real-unit normalization. The inherited observed-children filter remains, even though this moment pools fertility groups. The weighted quantile is the smallest sorted observed ratio where cumulative weight reaches the requested fraction. This is a **living-household** wealth/income distribution, not realized estates or decedents' bequests. Exact source: `code/data/psid_followup_mar2026/audit_intergen_bequest_family_size_targets.R`, `make_sample(reference_pooled_num_7684)` and `bequest_target_values`.

The old-tail bootstrap uses the existing builder's broader sampling population of weighted living reference persons ages 65–84, before outcome/children filters; then applies the exact outcome sample to each draw. It has 499 draws, seed 20260715. Sampling-population counts are saved separately from final valid-outcome sample counts. The point estimate **and bootstrap SE** of the existing 1984–2019 p90/p50 reproduce exactly to their saved precision.

The aggregate builder's complete 2005–2019 yearly weighted totals and counts reproduce the saved authoritative `yearly_ratios.csv`; the reconstructed long-pool aggregate estimate and SE reproduce its report. `wealth/verification.json` records checks and SHA-256 source/output fingerprints.

## What this resolves, and what remains

1. **Use an early window if the author chooses an initial-economy calibration.** The strictly pre-2007 old-tail SE is about 8.7% of the 2003/2005 point estimate, compared with 13.7% for 2005 alone. The 2005/2007 old-tail SE is about 11.7% of the point estimate; extending through 2003 reduces this to about 8.0%. These are usable empirical restrictions with appreciable uncertainty, not absent data. The longer pool's apparent precision must not be carried over to the early target.
2. **Dates matter for wealth levels.** The old-age median is about 7.99 in 2005/2007, compared with 6.50 in the long pool. The aggregate ratio is 7.36 in 2007 alone versus 6.93 in 2005/2007. Choosing a window is an economic initial-state approximation, not merely a precision choice. The matched initial economy must reproduce that chosen empirical window's interpretation.
3. **Parameter identification still needs assessment.** Wealth/earnings, old wealth level and old wealth dispersion all jointly constrain patience and bequest preferences. Their availability does not establish an invertible three-parameter mapping. In particular, the inherited weak sensitivity to the bequest wealth-shift parameter cannot be fixed by sample size alone. No new model Jacobian was run.
4. **Annual bequests/wealth 0.0088 remains external.** The target-provenance ledger identifies it as a Gale–Scholz historical normalization transmitted through De Nardi–Yang, with no project microdata builder, initial-period sample or measured standard error. A strict initial-data design must label it a maintained external restriction. Alternatively, replacing its identifying role with the newly available old-age wealth level would be an explicit model-identification decision for the bequest block and would require checking the joint Jacobian; it cannot be adopted automatically.
5. **No new joint covariance matrix is claimed.** The saved marginal bootstraps reproduce the existing separate-builder conventions. If the new objective uses a covariance-based joint weighting matrix, a coordinated bootstrap across the aggregate and old-age samples is additional work.

## Reproduction

Run from the repository root:

```sh
Rscript output/model/e5f_matched_pf_20260909a/design_research/wealth/build_initial_wealth.R output/model/e5f_matched_pf_20260909a/design_research/wealth
```

Deliverables are `wealth/aggregate_wealth_results.csv`, `wealth/old_wealth_results.csv`, the two bootstrap-draw tables, `wealth/aggregate_yearly_components.csv`, `wealth/run_metadata.csv`, `wealth/run.log` and `wealth/verification.json`. The diagnostic script is self-contained and writes only this evidence folder.
