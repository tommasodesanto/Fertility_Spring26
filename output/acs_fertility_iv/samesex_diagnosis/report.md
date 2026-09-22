# Same-sex diagnosis (national ACS 2005–2019): collected results

Source: job 18284841 (COMPLETED 0:0, 19m05s), remote outdir
`/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_samesex_diagnosis_20260922a`.
Only small JSON/CSV/log files were copied (`remote_collection/`); the analytic frames and metadata RDS were not.
No IV or AR fits. The only derived statistic is the BB − GG contrast, computed from each joint fit's saved clustered covariance.
All specifications: RF/FS on `samesex` (or BB, GG) + i(mother age at event) + i(RACE) + i(survey year) + i(event age), PERWT, household clusters.
RF units: rooms, bedrooms; ownership in pp. FS in pp for D = 1{3+ coresident children}.

## Validation
- 31/31 case receipts. 27 regression cases, all `full_fit`, 0 warnings. The single log line matching "warn" is a unit-test PASS message.
- A (baseline reproduction): pass. Source path/size/md5, original production lib md5 and saved-receipt identity all match; rf/fs coef, SE and nobs match exactly.
- C (coverage): all 51 state files. 656,986 cached eligible = 656,986 metadata eligible = 656,986 one-to-one joined rows. 0 duplicate, missing or unmatched keys; 0 recompute mismatches of HH/Z/D/event age (`one_to_one_and_recomputed_equal`). VT smoke passed.
- Pool audit (1,132,295 women with 2+ linked children): 0 invalid sexes, 0 age-order violations, 20,719 with second = third child age, 228 NCHILD–link-count mismatches, 2,486 households with more than one pool mother.
- D sample N equals baseline N (656,985 rooms; 656,986 bedrooms/owner). Boy1/Boy2 enter as non-constant regressors. Mixed reference = 327,011 (BG 164,253 + GB 162,759 − 1 missing-rooms mother).
- E excludes 11,500 of 656,986 eligible mothers (second = third child age): N 645,485 rooms, 645,486 others.

## Findings (associations; none of this identifies a mechanism)
1. **Boy1/Boy2 controls change nothing.** Rooms RF −0.0376 → −0.0380 (SE 0.0051); bedrooms −0.0341 → −0.0342; FS 3.27 pp unchanged. The sex-of-child coefficients are small (rooms: first boy +0.005, second boy +0.009).
2. **Both same-sex compositions are negative relative to mixed.** Rooms: BB −0.0306 (0.0061), GG −0.0454 (0.0063). Bedrooms: BB −0.0322 (0.0033), GG −0.0363 (0.0035). The first stages are identical (BB 3.28, GG 3.27 pp). BB − GG, from the full V: rooms +0.0148 (SE 0.0072, 95% CI [0.0006, 0.0289]); bedrooms +0.0041 (0.0039, CI includes 0); ownership +0.09 pp (0.21).
3. **By second-child age.** FS is ≈ 0 at age 0 (0.09 pp, F 0.9) and rises to 1.20, 3.10, 4.67, 5.58, 6.04 pp at ages 1–5. At age 0 the RFs are also ≈ 0 (rooms −0.009, bedrooms −0.004, both CIs include 0). At age 1 the bedroom RF is already −0.028 (0.0066) while FS is only 1.2 pp. From age 1 to age 5, the bedroom RF grows about 1.6× while the FS grows about 5×. The RF is not proportional to the FS across ages. That pattern is consistent with, but does not establish, a composition effect that does not run through a third child. Ownership RFs are imprecise at every age (all CIs include 0).
4. **Second/third tie exclusion changes nothing.** Rooms −0.0378, bedrooms −0.0347, owner −0.23 pp; FS 3.27 pp.
5. **Raw weighted means by cell** (table below) show the same ordering: same-sex cells have fewer rooms and bedrooms than mixed cells and more third children (D ≈ 0.245 vs 0.213).
