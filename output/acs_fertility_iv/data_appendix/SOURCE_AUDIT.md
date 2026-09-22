# Source audit: ACS sibling-composition housing appendix (final integration)

Appendix source: `latex/JMP_DS_draft/sections/appendix_acs_fertility_iv.tex`
(+ `.bib` alongside). Standalone build: `acs_iv_data_appendix.tex` here
(pdflatex, bibtex, pdflatex x2), 12 pages. Tables and figures come from
`code/empirical/acs/kleven_pseudo/build_acs_iv_appendix_tables.py`, which the
appendix `\input`s through `\acsivtabdir` and `\acsivfigdir`. It is not wired into
the author's main draft. Status: prose, tables and figures accepted by the lead
(2026-09-22). Five precision edits were then applied: the tie paragraph (the event age
is unchanged), the age-0 wording, the AR large-sample qualification, the Table 6
weighted-count note, and the closing persistence sentence. Ready for lead delivery.

Author: claude-opus-5-5, 2026-09-22. Local only in the integration pass: no SSH,
no jobs, no re-estimation.

## 1. Numbers and where each comes from

Baseline (Tables 2–3, Figure 1): `output/acs_fertility_iv/national_128g_results/national_18row_table.csv`.
The builder asserts equal rf/fs/iv nobs, `full_fit`, and 0 AR errors. The companion
file `acs_iv_18row_estimates_companion.csv` here gives all 18 rows' conventional
RF/FS/IV 95% CIs (±1.96 SE) and the tested-grid AR summaries, in display units.

Same-sex checks (Tables 4–6, Figures 2–3): `output/acs_fertility_iv/samesex_diagnosis/samesex_diagnosis_all_cases.csv`
(from `collect_samesex_diagnosis.py`, commit f72fc94c, accepted by the lead) and
`samesex_diagnosis/remote_collection/outdir/D_sex_cell_counts.csv`. The source is
national job 18284841 (COMPLETED 0:0), remote outdir
`/scratch/td2248/projects/kleven_acs_pilot_20260917/output/acs_samesex_diagnosis_20260922a`.
The lead independently verified 27 cases, 54 stage coefficients, full V, SE, CI and N
(`samesex_diagnosis/lead_numeric_review.json`).

| Paper statement | Value | Source |
|---|---|---|
| Sex-of-child controls, rooms RF | −0.037987 vs baseline −0.037627 | D_additive receipt; national table |
| Tie exclusion | 11,500 excluded; rooms −0.037827; N 645,485/645,486 | E receipts |
| BB, GG rooms RF | −0.0306013 (0.00613), −0.0453563 (0.00634) | D_joint receipt |
| BB, GG first stage | 0.032794, 0.032690 | D_joint receipt |
| BB−GG rooms | 0.01475499, SE 0.00719948, CI [0.0006443, 0.0288657] | joint V (bb, gg, cov); matches lead review |
| BB−GG bedrooms / owner | 0.00406 (0.00395) / 0.00092 (0.00208) | lead review |
| Age 0 | FS 0.09 pp, F 0.9; rooms, bedrooms RF CIs include 0 | B receipts |
| Age 1 bedrooms | RF −0.0278 (0.0066), FS 1.20 pp | B receipts |
| Age 1→5 ratios | FS 6.04/1.20 = 5.0; bedroom RF 0.0450/0.0278 = 1.62 | B receipts |
| Mixed reference | 327,012 (BG 164,253 + GB 162,759); rooms fit 327,011 | sex cells; D_joint |
| Coverage | 51 states, 656,986 keys one-to-one, 0 recompute mismatches in HH/Z/D/event age; baseline reproduction exact | C_full_coverage_gate, A_reproduction_gate (receipts only, not in paper) |
| Sample counts (Table 1) | see first-phase ledger below | `sample_gate_receipt.json`, `counts.json`, `analytic_frames_identity.json` |
| 91% / 4% RELATE; 1,692 invalid links | 11,777,336 and 543,774 of 12,975,041 | `child_relate_audit.json`, `link_audit.json` |
| ROOMS code 28 | 107 records, 2009 | `code/empirical/acs/kleven_pseudo/source_audit_extract27_20260919.md` l.66 |

AR implementation, as described in the paper: `ar_confidence_set` in
`twins_samesex_iv_lib.R` runs three regressions (Y, D, Y+D). For each grid value it
forms π̂_Y − β0·π̂_D with variance V_Y − 2β0·C + β0²·V_D, where C = (V_{Y+D} − V_Y − V_D)/2,
and compares the Wald statistic with F(1, df2), using fixest's clustered df2. The
paper describes this statistic and makes no claim beyond the tested grid. The old
3.84 boundedness argument was removed.

## 2. Literature: versions actually read

| Work | Version read | Pinpoints | Supports |
|---|---|---|---|
| Angrist & Evans (1998) AER 88(3):450–477 | Published article, JSTOR-distributed PDF (dpipe.tsukuba.ac.jp mirror), journal pagination checked by running headers | 452–453; 454; 460–461; 461; 465 fn12; 466–467 | Sample design; twin-age overcount; sex controls; FS 6.2 pp; excluded covariates; invariance to sex controls |
| Angrist, Lavy & Schlosser (2010) JOLE 28(4):773–823 | Warwick-hosted Oct 2010 PDF with journal pagination (warwick.ac.uk/fac/soc/economics/staff/vlavy/angrist_lavy_schlosser_qq_jole_october_2010.pdf) | 775; 799 | Twin non-randomness; room/clothes-sharing concern, which they attribute to R&W 2000; no-first-stage subsample check |
| Jones (2015) NBER WP 21391 | NBER PDF (nber.org/system/files/working_papers/w21391/w21391.pdf) | pp.2–3 (no pinpoint in prose) | Direct effects on non-compliers bias the IV estimand |
| Andrews, Stock & Sun | **Nov 20, 2018 manuscript** (par.nsf.gov/servlets/purl/10142670). The bib cites this version, with a cross-reference to ARE 11 (2019):727–753. | p.3; pp.16–17 | AR recommendation with one instrument; robust F = KP with one endogenous regressor |
| Anderson & Rubin (1949) AMS 20(1):46–63 | Bibliographic record only (Project Euclid, DOI 10.1214/aoms/1177730090) | none | Citation of the test |
| IPUMS USA variable pages MOMLOC, RELATE, ROOMS, BEDROOMS, OWNERSHP, PERWT | usa.ipums.org/usa-action/variables/<VAR>, 2026-09-22 | — | Link semantics; code frames; 2007 topcodes (rooms 9, bedrooms 5; bedrooms = code − 1, confirmed by lead) |

**Secondary attribution only:** Rosenzweig & Wolpin (2000) JEL 38(4):827–874 was not
opened. The paper cites it only as the source to which ALS (p.775) attribute the
room-sharing concern. The bib entry was checked against the AEA article page.

**Not read, not cited:** Bhalotra & Clarke (REStat 2019); Rosenzweig & Wolpin (1980).

## 3. Corrections to earlier project reports (made by minimal edit, authorized by lead review)

1. `output/acs_fertility_iv/RESULTS.md` l.16 said the robust F is "**not** a
   Kleibergen–Paap statistic". It now says that with one instrument and one
   endogenous regressor this is the statistic to which KP reduces (ASS 2018
   pp.16–17), and that no separate KP routine was run. Nothing else changed.
2. `national_18row_table_with_pp.csv`, row Twin1_pooled0_5 / BEDROOMS_out: the
   `ar_honesty_note` called the set a single grid point. It now reads "no tested grid
   value accepted (empty on the tested grid); this does not establish the location
   or width of the continuous AR confidence set". This is a one-line diff; estimates
   and line endings are unchanged.
3. Earlier notes called the 4,024,238 base count "mothers". They are women aged
   21–35, including women with no linked child. Table 1 says so.
4. The 32,186 same-age ties are counted over all ≥2-child base women, not within the
   event window, and are not a sequential step. Table 1 is restructured to avoid
   implying attrition.
5. The ASS preprint's own reference list gives AR 1949 as 21:570–582. The bib uses
   20(1):46–63.

## 4. Deviations from Angrist–Evans (stated in the appendix)

No CEB-consistency filter (the ACS lacks CEB); no allocation-flag drop (the extract
lacks the flags); event age 0 included; race indicators but no separate Hispanic
control; baseline without Boy1st/Boy2nd, reported as a check; whole-year ages with
no quarter of birth.

## 5. Remaining limitations (not identified by these data)

- Same-age proxy: the twin share is unmeasured. IPUMS lists BIRTHQTR for ACS 2005+,
  but it is not in the extract.
- Clock: the event age is the age of the oldest or second-oldest linked resident
  child, not a verified birth date or birth order.
- Linked children are social/coresident links (step and adoptive children included).
  Non-resident children are unobserved. There are no imputation flags. The pool has
  228 resident-count–link mismatches and 2,486 households with more than one pool
  mother; SEs are household-clustered.
- The negative same-sex IV cannot be classified as an exclusion failure or a genuine
  effect. No sign restrictions were tested, and the room-sharing channel is not
  identified.
- The BB−GG rooms contrast is marginal, one of several, and unadjusted for
  multiplicity.
- Across-age non-proportionality is descriptive: the ages are separate
  cross-sections, and composition and heterogeneous effects may differ across them.
- AR: only the finite grid was evaluated. There are no continuous endpoints, and no
  claim about the Twin1 pooled bedrooms set beyond "no grid value accepted".
- There is no joint covariance for rooms vs bedrooms. Capped rooms minus capped
  bedrooms is not a decomposition.
