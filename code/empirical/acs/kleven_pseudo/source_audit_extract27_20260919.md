# ACS extract27 source receipt — 2026-09-19

Status: source verified and atomically promoted on Torch; source-audit driver is being repaired for the approved Northeast memory scope.

## Local source

- Path: `code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta`
- Size: 9,919,999,546 bytes (Stata header metadata)
- SHA-256: `edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e`
- Observations: 59,046,776
- Variables: 69
- Header timestamp: 2025-07-07 17:46
- Required housing/link fields present in the header: `YEAR SAMPLE SERIAL PERNUM HHWT PERWT STATEFIP PUMA GQ RELATE MOMLOC POPLOC SEX AGE MARST RACE EDUC NCHILD NCHLT5 ELDCH YNGCH FERTYR ROOMS BEDROOMS OWNERSHP OWNERSHPD`.
- `extractor27.do:30` and `:33` define rooms and bedrooms; `:113-126` label rooms, bedrooms, own children, and oldest/youngest child ages; `:1627-1630` define `FERTYR==2` as “Yes” (children born within the last year).
- Existing metadata receipts show one-year ACS sample codes `200701`, `201101`, `201501`, `201901`, and `202301`; the complete year/sample distribution remains unverified until the remote audit can read the staged file.
- Sample labels in `extractor27.do:201-236` distinguish ACS one-year codes ending `01`, ACS five-year codes ending `03`, PRCS one-year codes ending `02`, and PRCS five-year codes ending `04`. A matching key is valid only within the same product code; year overlap alone does not establish source identity.

## Remote staging state

The intended non-overwriting destination was:

`/scratch/td2248/projects/kleven_acs_pilot_20260917/inputs/ACS/local_extract27_20260919/extract27.dta.gz`

The transfer was compressed from the local file after a sample showed substantial compression. It terminated before completion with:

`Permission denied (gssapi-keyex,gssapi-with-mic,password,keyboard-interactive)`

The remote partial file was not checksum-verified and must not be used. No `sbatch` job was submitted. The planned audit script was staged under `forensics/source_audit_extract27_20260919/`, but it was not run. No source request, credential action, or retry was attempted after the new authentication rejection.

## Recovery recheck — 2026-09-20

The remote decompression produced `/scratch/td2248/projects/kleven_acs_pilot_20260917/inputs/ACS/local_extract27_20260919/extract27.dta.new` at the expected size of 9,919,999,546 bytes. The adjacent `sha256.new` records the expected SHA-256, `edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e`, for that `.dta.new` path; a direct remote `sha256sum` reproduced the same hash. This corrects an earlier 63-character transcription that omitted the final `e`.

The staged `watch_stage_submit.sh` process is no longer present. Its `watch.log` is empty, and the input directory has no `source_ready.txt` or promoted `extract27.dta`. The planned audit output directory is absent, and Slurm accounting shows no `acs27-source-audit` job or other matching source-audit job. Therefore no year/sample, key, or housing-field compatibility findings are available yet, and no source-audit result should be inferred from the completed `.dta.new` file.

## Current recovery state — 2026-09-20

The verified `.dta.new` was promoted to `/scratch/td2248/projects/kleven_acs_pilot_20260917/inputs/ACS/local_extract27_20260919/extract27.dta`; `source_ready.txt` remains absent because the original watcher exited before its promotion step. The durable 64-character source hash is `edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e`; an earlier receipt transcription omitted the final `e` and must not be reused.

Audit jobs `18075798`, `18077232`, `18077677`, and `18077695` failed before producing audit outputs. The first two failures were launcher/source-hash setup failures; the fourth exposed the 63-versus-64-character transcription in the launcher. The corrected durable launcher and driver use JSON receipt validation, exact 64-hex checking, compact national fingerprints, and Northeast-only source-key concordance against the prepared `ne_acs.RData` object. No source compatibility result is established until that corrected audit completes.

## Consequence for the housing adapter

The local file proves that rooms and bedrooms exist, but it does not yet prove exact overlap with the Torch Kleven ACS source. Before joining housing outcomes to the preserved first-birth employment panel, verify the full `(YEAR,SAMPLE,SERIAL,PERNUM)` key, `SAMPLE` product labels, year coverage, serial namespace, and row-level key uniqueness. If `SAMPLE` is absent from the panel or the product codes differ, stop instead of joining on `(YEAR,SERIAL,PERNUM)`. No matching rerun or housing estimate is valid before this receipt is complete.

## Completed source audit — 2026-09-20

The corrected continuation audit completed as Torch job `18078493` on `cs618` in
6:13 with exit code `0:0` and peak memory `39,626,376 KB`. The launcher passed
the exact 64-hex SHA receipt and byte-count checks, the case-validated header
smoke, and the selected-field read of 59,046,776 rows by 29 fields. The audit
saved the reduced Northeast packet immediately after the read:

`/scratch/td2248/projects/kleven_acs_pilot_20260917/output/kleven_acs_pilot/source_audit_extract27_20260919/ne_extract27_housing_key_packet.rds`

The packet has 2,784,052 Northeast rows and 38 columns, including raw housing
codes, code-validity flags, outcome-safe fields, source keys, weights, and the
matching/concordance fields. The exact source key `(YEAR,SAMPLE,SERIAL,PERNUM)`
is unique in both extract27 and the prepared Northeast panel. The unique-key
intersection is 2,190,987 rows, with zero sex, age, or ownership mismatches.
Year/sample support matches exactly for 2005–2019 ACS 1-year samples. The
prepared panel additionally contains 2000–2004 samples, while extract27
additionally contains 2020–2023 samples; these are explicit coverage differences,
not a key or value mismatch.

The audit found 107 unknown literal `ROOMS` codes, all `28` in the 2009 ACS
1-year sample (`200901`), and all within the Northeast subset. The official
IPUMS ROOMS page lists `N/A`, `1–27`, and `30`, with no category 28 or 29:
<https://usa.ipums.org/usa-action/variables/ROOMS>. The raw value is retained,
its validity flag is false, and its future outcome field is missing. No rows,
keys, weights, or samples were dropped or imputed. Housing outcome use therefore
requires lead review of this explicit 2009 code anomaly.

Compact durable receipts are under
`code/empirical/acs/kleven_pseudo/output/source_audit_extract27_20260919/`.

## NE proxy readiness diagnostic — 2026-09-20

The existing strict roster builder ran on the verified NE packet as Torch job
`18078758` and completed in 51 seconds with exit code `0:0` and peak memory
`5,630,184 KB`. It used explicit match covariates `SEX`, `EDUC`, `MARST`,
`RACE`, and `STATEFIP`, with `FERTYR` codes yes=`2`, no=`1`, and unknown=`0,8`.
No matching or housing estimation ran.

The diagnostic found 1,439,607 female rows with valid age, 722,337 valid
MOMLOC links, 6,872 strict event-0 anchors, 1,771 anchors supporting the full
−5 through −1 pre-window, 6,325 supporting the reference period −2, 70,063
strict post rows, 211,056 one-child donors, and 21,919 negative-time donor
targets. At event 0, `FERTYR` counts were 7,931 yes, 496 no, and 16 unknown;
the strict anchors all come from the observed-yes group after the other roster
and age rules. The support and link tables are saved under
`code/empirical/acs/kleven_pseudo/output/second_birth_proxy_diagnostic_20260920/`.

The historical sample label `200004` is recorded as `ACS 2000` from
`extractor27.do`, rather than inferred as a PRCS 5-year code. This metadata
correction changes labels only; it does not change row counts, keys, or
concordance.
