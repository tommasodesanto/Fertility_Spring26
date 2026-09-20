# Kleven ACS employment pilot

This task starts with the author's original ACS data and original R matching/cleaning code. It does not implement a new pseudo-panel estimator. The first empirical run is Vermont (state FIPS 50), not a national replication claim. CPS 1995–1999 childless observations provide the early donor years used by the original ACS branch.

## Scope and budgets

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
Northeast packet remains on Torch. The audit found exact-key uniqueness and zero
sex/age/ownership mismatches on 2,190,987 shared unique keys. It also found 107
literal `ROOMS=28` values, all in the 2009 ACS 1-year sample; these are preserved
as raw codes with missing outcome-safe values pending code review.

The NE-only second-birth roster diagnostic completed as Torch job `18078758`.
It produced 6,872 strict event-0 anchors, 1,771 full-pre anchors, 6,325
reference-period anchors, 70,063 post rows, 211,056 one-child donors, and
21,919 negative-time donor targets. Matching and housing estimation remain
unrun. Compact support receipts are in
`output/second_birth_proxy_diagnostic_20260920/`.
