# Kleven ACS Northeast employment and housing diagnostics

## Current validated package — 2026-09-20

The current package uses the verified Northeast ACS extract and the preserved author matching lineage. It contains descriptive matched pseudo-panel diagnostics for employment and housing; it does not establish a causal fertility shock or a national replication.

| Stage | Status and authoritative artifacts |
| --- | --- |
| Source audit | Torch job `18078493`; [audit report](source_audit_extract27_20260919.md) and [source manifest](output/source_audit_extract27_20260919/source_audit_manifest.json). Exact source-key uniqueness and zero sex/age/ownership mismatches hold on 2,190,987 shared keys, with shared source years 2005–2019. |
| Employment benchmark | The preserved 20.26% Northeast benchmark receipt is [`ne_benchmark_receipt.json`](/Users/tommasodesanto/.local/share/kleven_execution_review/results/ne_benchmark_receipt.json). Treat its rounded-label corroboration as descriptive, not causal; the older VT-only Torch receipt is `/scratch/td2248/projects/kleven_acs_pilot_20260917/vt_employment_glm_20260917/vt_employment_receipt.json` and must not be substituted for the Northeast artifact. |
| First-birth housing, full window | Torch job `18079576`; [run report](first_birth_housing_run_report_20260920_job18079576.md), with saved full-fit and support receipts under `output/first_birth_housing_20260920_job18079576/`. |
| First-birth housing, short window | Torch job `18080591`; [bundle receipt](first_birth_sensitivity_bundle_receipt_18080591.md), with the common-cohort support gate and saved curves under `output/first_birth_sensitivity_bundle_18080591/`. |
| Second-birth housing | Torch job `18080900`; [run report](second_birth_housing_run_report_20260920_job18080900.md), with corrected saved curves and fit receipts under `output/second_birth_housing_20260920_job18080900/`. |
| Second-birth readout | Torch job `18082204`; [readout outputs](output/second_birth_housing_readout_20260920_job18082204/readout_receipt.json) and nested specification, FERTYR, composition, and common-cohort receipts in the same directory. |
| PSID comparison | [Saved PSID reference](psid_comparison_reference.md), used as external descriptive context only. |

The second-birth readout used 5,047 primary gap≥2 anchors and 1,734 joint-negative anchors; all selected anchors have `FERTYR_status=yes`. The common implied-event-year 2007–2016 gate passed all 60 cohort-year by event-time cells for rooms, bedrooms, and ownership before fitting. The nested readout retained 29,546 observations in every specification. These results are conditional on the constructed coresident roster proxy, exact donor matches, person weights, and source-household clustering.

The verified raw input and output locations are `/scratch/td2248/projects/kleven_acs_pilot_20260917/vendor/` and `/scratch/td2248/projects/kleven_acs_pilot_20260917/output/kleven_acs_pilot/`. The readout job output is `/scratch/td2248/projects/kleven_acs_pilot_20260917/output/kleven_acs_pilot/second_birth_housing_readout_20260920_job716202e6/`. The reviewed source version is local `ae144970` (the readout driver and launcher were introduced in ancestor `716202e6`).

The completed package can be checked locally without compute:

```bash
Rscript code/empirical/acs/kleven_pseudo/test_run_second_birth_housing.R
Rscript code/empirical/acs/kleven_pseudo/test_second_birth_housing_readout.R
```

To reproduce a stage on Torch, copy the reviewed source files at commit
`ae144970` into `/scratch/td2248/projects/kleven_acs_pilot_20260917/code/`,
then run the launchers sequentially from
`/scratch/td2248/projects/kleven_acs_pilot_20260917/code/cluster`:
`run_kleven_acs_source_audit.sbatch`, `run_first_birth_housing.sbatch`,
`run_first_birth_sensitivity_bundle.sbatch`, `run_second_birth_housing.sbatch`,
and finally `run_second_birth_housing_readout.sbatch`. Set a fresh
`OUTDIR`, `SHORT_WINDOW_OUTDIR`, `SECOND_BIRTH_OUTDIR`,
`SECOND_BIRTH_HOUSING_OUTDIR`, or `SECOND_BIRTH_READOUT_OUTDIR` through
`sbatch --export=ALL,...` for any rerun; the saved output directories are
receipts for the completed jobs and should not be overwritten. The final two
stages consume the verified source and saved matching checkpoints; they do not
replace the source transfer or rematch contract.

The empirical limits are explicit. ACS repeated cross-sections provide a coresident two-child newborn proxy rather than observed birth histories. The strict second-birth population requires `NCHILD=2` and exactly two linked children, so observed third-child households are excluded; this changes the negative-composition population. The source does not establish nonresident children, biological parentage, exact twinning, or the realization clock for a twins or same-sex fertility instrument. The matched event curves therefore do not validate an IV, identify an exogenous fertility shock, or support a causal housing claim. PSID comparisons use different populations and clocks, so they provide directional context rather than empirical validation. `ROOMS=28` remains unknown; bedrooms use the reviewed (x-1) cap-5 rule; ownership is a 0/1 proportion and is displayed in percentage points only in figures.

## Historical preparation notes

The historical pilot began with the author's original ACS data and R matching/cleaning code. Its initial Vermont slice and CPS donor preparation are retained below as background; they do not define the current Northeast package.

All raw data loading, cleaning, matching, estimation and plotting run on Torch in `/scratch/td2248/projects/kleven_acs_pilot_20260917`. The laptop handles only source inspection, transfer and small checks. No model calibration or paper text changes.

Preparation job **17917426**: partition `cs`, account `torch_pr_570_general`, 8 allocated CPUs and 128 GiB (thread counts capped at one), 30-minute walltime. Two sequential original-file downloads, approximately 4.2 GB total; each has a 3 GB / 10-minute cap and atomic rename. Data are loaded sequentially on a compute node to inspect schemas and retain Vermont observations with labels. Expected preparation is minutes, but no observed throughput yet; hard caps govern the job. Job 17917334 was cancelled while pending because its default partition was unsuitable; no scientific computation ran there. A partition-specific test established the accepted CPU/memory allocation before replacement.

The author's full pipeline reports roughly 63 hours on a 1.5 TB RAM server. Do not launch `MASTER.R`, all figures, or national matching from this pilot. The next stages must clean, run the exact matching/regression path on a smaller supported demographic slice, and only then run the full Vermont sample. Freeze national resources from measured pilot time/memory before launching national work. No automatic production promotion or repeated retries.

## Files and evidence

- `source_contract.json`: immutable author-source hashes and intended estimator parameters.
- `prepare_inputs.py`: original author URLs, bounded acquisition and input hashes.
- `inspect_inputs.R`: compute-node schema/package inventory and pilot extraction.
- `prepare.sh`: compute-node preparation entry point.
- `../../../cluster/submit_kleven_acs_pilot.sh`: preparation submission wrapper (relative path from this directory).

Remote evidence goes to `output/kleven_acs_pilot/`: `input_receipts.json`, per-source schema/year/sample counts, package inventory, session information, completion status, Slurm logs, and download heartbeat. Inspect those files before running downstream work. Raw microdata stay on Torch and are never sent to the model service. Author sources are extracted unchanged into remote `vendor/`.

The OpenCode/Kimi source-map pass identified exact required fields, weights and early CPS donor behavior. Implementation is separately reviewed before use. Housing comes after employment checks; whether the author's raw extract contains housing fields remains an explicit schema question. Matching on housing outcomes is not permitted to force agreement.

The verified extract27 source audit completed as Torch job `18078493`. Compact
receipts are in `output/source_audit_extract27_20260919/`, and the full reduced
Northeast packet remains on Torch. The authoritative digest and overlap values
are in `output/source_audit_extract27_20260919/source_audit_manifest.json`;
this README does not duplicate the hash. The audit found exact-key uniqueness
and zero sex/age/ownership mismatches on 2,190,987 shared unique keys. It also
found literal `ROOMS=28` values, which are preserved as raw codes and remain
missing for outcome analysis under the reviewed coding contract.

The NE-only second-birth roster diagnostic completed as Torch job `18078758`.
It produced 6,872 strict event-0 anchors, 1,771 full-pre anchors, 6,325
reference-period anchors, 70,063 post rows, 211,056 one-child donors, and
21,919 negative-time donor targets. This paragraph records the historical
roster-only diagnostic; later transformed-support matching and first-birth
housing results are indexed below. Compact support receipts are in
`output/second_birth_proxy_diagnostic_20260920/`.

## Detailed receipts and historical context

The source audit is documented in `source_audit_extract27_20260919.md`. Its
authoritative SHA-256 and byte count are recorded in
`output/source_audit_extract27_20260919/source_audit_manifest.json`; the
verified extract has exact source-key uniqueness and zero sex, age, or
ownership mismatches on 2,190,987 shared keys. Raw `ROOMS=28` values remain
explicitly unknown for outcome analysis; they are not recoded to nine.

The completed first-birth housing fit is summarized in
`first_birth_housing_run_report_20260920_job18079576.md`. Torch job 18079576
used source-household clustered uncertainty with heteroskedastic sensitivity
and produced state-separated rooms, bedrooms, and ownership curves. Its
matched pseudo-panel estimates are descriptive conditional on the constructed
source-key matches.

The corrected short-window and transformed-support bundle is recorded in
`first_birth_sensitivity_bundle_receipt_18080591.md` and the compact readout
under `output/first_birth_sensitivity_bundle_18080591/`. Job 18080591 passed
the six-event support gate for all 360 requested outcome/state/gender/cohort
groups and saved event -2 weighted baselines, Kish effective sample sizes, and
the comparable +3 minus -1 contrast. The full-versus-short table compares
that common contrast only; post-minus-pre estimates use different windows and
are not treated as a replication check.

The transformed second-birth stage verified raw-field concordance and donor
support. It is still a coresident roster proxy with conditional matching
uncertainty, not an observed biological second-birth event. The estimator and
readout are complete; the saved anchor, gap, donor-reuse, and source-weight
receipts remain the authoritative inputs for review, and no further allocation
is active.

The saved PSID comparisons in `psid_comparison_reference.md` provide external
descriptive context with different longitudinal samples and event clocks.
They do not validate the ACS proxy or establish a causal effect.

The original motivation included twin and same-sex fertility shocks. This ACS
package does not estimate or validate those instruments. The verified extract
supports coresident sex-composition and same-age roster proxies, but it does
not establish historical first-two-birth links, exact twinning, nonresident
children, biological parentage, or the instrument-realization outcome clock.
A dedicated design must define and audit those objects before any IV claim.
