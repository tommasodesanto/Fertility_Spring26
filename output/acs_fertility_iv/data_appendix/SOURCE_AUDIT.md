# Source audit: ACS sibling-composition housing appendix

Appendix source: `latex/JMP_DS_draft/sections/appendix_acs_fertility_iv.tex`
(+ `.bib` alongside). Standalone build: `acs_iv_data_appendix.tex` here
(pdflatex, bibtex, pdflatex x2). Status: first-phase draft, settled material
only. Not wired into the author's main draft; not final.

Author: claude-opus-5-5, 2026-09-22. No estimation, SSH, or large RDS load was
done for this phase.

## 1. Numbers: where each comes from

All from `output/acs_fertility_iv/national_128g_results/`:

| Appendix object | Source file | Check |
|---|---|---|
| Tables 2–3 (all 18 rows) | `national_18row_table.csv` via `code/empirical/acs/kleven_pseudo/build_acs_iv_appendix_tables.py` | Builder asserts rf/fs/iv nobs equal the pre-fit usable N, `status==full_fit`, `ar_n_errors==0`. The table blocks typed into the appendix are byte-identical to the builder's fragments (checked by substring match). |
| Figure 1 | same CSV, same builder | RF coef and 95% CI from `rf_ci_lower/upper`. |
| Table 1 source/gate rows | `sample_gate_receipt.json` | 59,046,776 → 46,373,936 (12,672,840 out-of-year; 0 non-1-year). |
| Table 1 base, design counts | `counts.json`, `analytic_frames_identity.json` | 4,024,238 women / 3,783,862 HH; Twin1 857,607 / 325,171 / 16,496 / 388; SameSex 1,132,295 / 32,186 / 656,986 / 329,974 (BB 173,813, GG 156,161) / 148,775. |
| 91% / 4% RELATE shares; 1,692 invalid links | `child_relate_audit.json`, `link_audit.json` | 11,777,336 / 12,975,041 = 90.8%; 543,774 / 12,975,041 = 4.2%. Computed after the year gate, before any age restriction on mothers. |
| ROOMS code 28: 107 records, 2009 | `code/empirical/acs/kleven_pseudo/source_audit_extract27_20260919.md` l.66 | Cap/valid codes match `twins_samesex_iv_lib.R` l.32–35. |
| −0.0376/0.0327 = −1.15 | CSV: rf −0.037627, fs 0.032745, iv −1.14909 | 1/π = 30.5. |
| AR boundedness claim | CSV `ar_n_components==1`, `ar_fully_interior_bounded==True` in 17/18 cells | Just-identified AR set is bounded iff first-stage F > 3.84; min F in tables = 174. Twin1 pooled bedrooms: 0 components on grid (set lies strictly between 0.3 and 0.4; contains β̂ = 0.333 because AR stat is 0 at β̂). |

## 2. Literature actually read (full text, local copies in /tmp/lit)

| Work | Copy read | Pinpoints used | Supports |
|---|---|---|---|
| Angrist & Evans (1998) AER 88(3):450–477 | JSTOR-distributed published PDF (journal pagination), via dpipe.tsukuba.ac.jp mirror | p.452–453 (sample: women 21–35, oldest < 18, CEB-consistency and allocation drops); p.454 (age-only twins inflate twin rate 35%, twins restricted to 1980); pp.460–461 (Boy1st/Boy2nd controls, sex correlation, fn10); p.461 (FS 6.2 pp, 1980); p.465 fn12 (education, husband's earnings endogenous); pp.466–467 (invariance to sex controls) | Sample design, controls rationale, first-stage comparison |
| Angrist, Lavy & Schlosser (2010) JOLE 28(4):773–823 | Warwick-hosted Oct 2010 PDF, journal pagination | p.775 (twins vary with maternal age/race, spacing, health; sex-composition room/clothes sharing, attributed to R&W 2000); p.799 (no-first-stage subsamples, no reduced-form relation for two-boy instrument) | Twin non-randomness; room-sharing exclusion concern |
| Jones (2015) NBER WP 21391 | NBER PDF | pp.2–3 (direct effects on non-compliers bias the IV estimand; R&W hand-me-down example) | Direct-effect amplification argument |
| Andrews, Stock & Sun (2019) ARE 11:727–753 | **Nov 20, 2018 preprint** (NSF PAR purl 10142670), not the published version | Preprint p.3 (report AR with one instrument), pp.16–17 (robust F = KP with one endogenous regressor), p.8 (ratio) | F/KP equivalence; AR recommendation. Pinpoints removed from prose because they are preprint pages. |
| Anderson & Rubin (1949) AMS 20(1):46–63 | Bibliographic record only (Project Euclid, DOI 10.1214/aoms/1177730090) | none | Citation of the test only |
| IPUMS USA variable pages: MOMLOC, RELATE, ROOMS, BEDROOMS, OWNERSHP, PERWT | usa.ipums.org/usa-action/variables/<VAR>, fetched 2026-09-22 | — | Link semantics (step/adoptive), code frames, 2007 topcodes, universes |

**Cited but not opened:** Rosenzweig & Wolpin (2000) JEL 38(4):827–874. Bib entry
checked against the AEA page; the appendix attributes the room-sharing concern
to it only through ALS p.775 and Jones pp.2–3.

**Not read, not cited:** Bhalotra & Clarke (REStat 2019, twin births and maternal
condition; MIT Press 403); Rosenzweig & Wolpin (1980).

## 3. Discrepancies found (not in prose)

1. The 4,024,238 base count is women aged 21–35, including those with no linked
   child. Earlier reports call them "mothers".
2. `national_18row_table_with_pp.csv` `ar_honesty_note` says the Twin1 pooled
   bedrooms AR set is a single grid point. It is empty on the grid (0 components).
   The appendix states the correct reading.
3. `output/acs_fertility_iv/RESULTS.md` l.16 says the robust F is "**not** a
   Kleibergen–Paap statistic". With one instrument
   and one endogenous regressor the two coincide (ASS preprint pp.16–17).
4. The 32,186 same-age tie count for the same-sex design is taken over the whole
   ≥2-child pool, not after the event-age window. It is not a sequential step
   down to 656,986; Table 1's note says so.
5. The ASS preprint's reference list gives the wrong AR 1949 volume/pages (21:570–582).
   The bib uses 20(1):46–63 (Project Euclid).
6. IPUMS extract version/DOI unknown to this author; the bib cites variable pages
   with an access date only.

## 4. Deviations from Angrist–Evans, stated in the appendix

No CEB-consistency filter (ACS lacks CEB); no allocation-flag drop (extract lacks
flags); event age 0 included (AE drop second child < 1); race indicators but no
separate Hispanic control; no Boy1st/Boy2nd in baseline; no quarter of birth in the
extract, so the same-age proxy uses whole-year ages.

## 5. Evidence gaps

- **BIRTHQTR**: IPUMS lists it for ACS 2005+. It is not in the extract. Adding it
  would sharpen the same-age proxy (the AE p.454 35% overcount is the reason).
- No Boy1st/Boy2nd robustness row.
- No joint covariance for rooms vs bedrooms, or BB vs GG. The appendix makes no
  formal contrast.
- No imputation flags for child age/sex.
- Ages 3 and 5 are described as the reported windows, not as pre-registered,
  because no pre-registration record was checked.

## 6. Dependency: pending same-sex diagnostic subsection

The same-sex diagnostic results (matrix A–E) are not in the appendix. No
placeholder text was written. When the lead collects the national diagnosis
results, add a subsection after "Results" (label `app:acs_iv_diagnostics`) with
its own table, built by extending the builder. Until then, the room-sharing
reading in the prose stays "consistent with, not established".
