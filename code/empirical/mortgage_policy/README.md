# Mortgage Policy Empirical Pipeline

This folder implements the first-stage setup for the mortgage-policy natural
experiment note:

- source manifest for official HUD, FHFA, and CFPB/HMDA files;
- downloader for small loan-limit files and optional HMDA LAR files;
- builder for county and CBSA loan-limit exposure measures;
- optional HMDA aggregation script for the first-stage lending response.

## Recommended Run Order

From the repository root:

```bash
Rscript code/empirical/mortgage_policy/download_mortgage_policy_sources.R
Rscript code/empirical/mortgage_policy/build_loan_limit_exposure_panel.R
```

This downloads only the small HUD loan-limit files by default and creates:

- `output/hud_loan_limits_county_long_2008_2009.csv`
- `output/loan_limit_exposure_county_2008_2009.csv`
- `output/loan_limit_exposure_cbsa_2008_2009.csv`
- `output/loan_limit_exposure_summary.md`

## Optional HMDA First Stage

The HMDA historical LAR files are large: roughly 330-570 MB per year for the
filtered nationwide first-lien owner-occupied 1-4 family files. Download them
only when ready to run the first-stage aggregation:

```bash
DOWNLOAD_HMDA=1 Rscript code/empirical/mortgage_policy/download_mortgage_policy_sources.R
YEARS=2007,2008,2009 Rscript code/empirical/mortgage_policy/build_hmda_first_stage_county.R
Rscript code/empirical/mortgage_policy/summarize_hmda_first_stage.R
```

The first-stage script aggregates originated home-purchase loans by county and
CBSA, with counts for:

- FHA originations;
- conventional loans sold to Fannie Mae or Freddie Mac;
- loan amounts in the GSE 2009-vs-2007 newly eligible band;
- loan amounts in the FHA 2009-vs-2008 band.

The summary script writes:

- `output/hmda_first_stage_exposure_bins.csv`
- `output/hmda_first_stage_exposure_bin_deltas.csv`
- `output/hmda_first_stage_summary.md`

Initial run for `2007-2009`: above-median treated counties show a visible GSE
first stage. Conventional loans sold to Fannie/Freddie in the newly eligible
GSE band rise from `0.48%` of originations in `2007` to `2.73%` in `2008` and
`3.84%` in `2009`.

## Natality Merge Scaffold

After the first stage is visible, build the preliminary natality merge:

```bash
Rscript code/empirical/mortgage_policy/build_natality_exposure_panel.R
```

This writes:

- `output/natality_gse_exposure_panel_2005_2012.csv`
- `output/natality_gse_exposure_summary_2005_2012.csv`
- `output/natality_gse_exposure_summary.md`

The initial default uses
`/Users/tommasodesanto/Desktop/Projects/Datasets/BirthsHybridMSAHpiPop_collapsed.dta`.
Treat this as a merge scaffold only. The generated birth-rate summaries need a
separate panel audit before they can support a fertility-response claim.

## ACS Housing Outcomes

Build ACS young-household ownership and room outcomes:

```bash
Rscript code/empirical/mortgage_policy/build_acs_gse_exposure_panel.R
```

This writes:

- `output/acs_gse_exposure_panel_2005_2012.csv`
- `output/acs_gse_exposure_summary_2005_2012.csv`
- `output/acs_gse_event_study_diagnostics.csv`
- `output/acs_gse_exposure_summary.md`

Initial run for household heads ages `25_34`: the first-stage mortgage channel
does not translate into a clean positive ACS housing response. Ownership
interactions are small, while rooms and the share of households in `6+` room
units move negative after `2008` in the broad metro panel. Treat this as a
warning sign before running fertility regressions.

## Refined HMDA-Takeup Screens

After building the HMDA, ACS, and natality scaffolds, run:

```bash
Rscript code/empirical/mortgage_policy/run_refined_acs_takeup_design.R
Rscript code/empirical/mortgage_policy/run_refined_natality_takeup_design.R
```

The ACS script constructs a narrower treatment screen from the HMDA first
stage:

\[
\Delta Share^{GSEsold,newband}_{c,2008/09-2007},
\]

then marks high-takeup metros as positive-exposure CBSAs in the top quartile of
that takeup change. It adds local HPI relative to `2007` and MSA unemployment
controls.

It writes:

- `output/hmda_refined_takeup_design.csv`
- `output/acs_refined_takeup_panel_2005_2012.csv`
- `output/acs_refined_takeup_summary_bins.csv`
- `output/acs_refined_takeup_event_study.csv`
- `output/acs_refined_takeup_summary.md`

The natality script applies the same takeup screen to the preliminary natality
merge and writes:

- `output/natality_refined_takeup_panel_2005_2012.csv`
- `output/natality_refined_takeup_summary_bins.csv`
- `output/natality_refined_takeup_event_study.csv`
- `output/natality_refined_takeup_summary.md`

Current read: the refinement does not rescue the mechanism. The HMDA first
stage is real, but ages `25_34` ACS outcomes still do not show a positive
ownership or housing-space response. The natality output is diagnostic only:
birth-rate levels are unstable and placebo age groups do not behave cleanly.
Do not use the refined natality screen as a headline fertility estimate.

## First-Time Homebuyer Tax Credit Screen

As a separate crisis-era policy screen, run:

```bash
Rscript code/empirical/mortgage_policy/run_fthc_tax_credit_screen.R
```

This uses ACS owner-reported 2007 median home values for owner household heads
ages `25_34` to construct:

\[
CreditShare_c = \frac{8000}{MedianHomeValue_{c,2007}}.
\]

The continuous treatment is scaled per `10` percentage points of credit share;
the binary treatment marks top-quartile exposed CBSAs. The script writes:

- `output/fthc_exposure_2007.csv`
- `output/fthc_acs_screen_panel_2005_2012.csv`
- `output/fthc_acs_summary_bins.csv`
- `output/fthc_acs_event_study.csv`
- `output/fthc_acs_summary.md`

Current read: this is a better housing-mechanism screen than the FHA/GSE
loan-limit design. More exposed cheaper metros show positive ACS responses for
mean rooms and `6+` room units among ages `25_34`, but ownership effects are
small and ACS birth-in-past-year does not move positively. Treat the screen as
supportive evidence that down-payment relief can shift family-sized housing
access, not as a fertility reduced form.

## Identification Note

The cleanest public exposure measure in the current build is:

\[
\log L^{GSE}_{c,2009} - \log L^{GSE,proxy}_{c,2007},
\]

where \(L^{GSE,proxy}_{c,2007}=417000\) outside statutory high-cost areas and
\(625500\) in AK/HI/GU/VI. FHA 2009-vs-2008 is included, but it is a narrower
comparison because the 2008 Economic Stimulus Act had already expanded FHA
limits.

Do not move to ACS or natality outcomes unless HMDA shows a first-stage lending
response in the newly eligible loan-amount bands.
