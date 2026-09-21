# First-birth rooms and ownership: PSID-to-ACS mapping

Status: specification and source-contract map for lead review. This note does
not submit a Torch job. The latest user instruction authorizes a national ACS
rematch; the source and matching stages below therefore describe that path
explicitly. The prior Northeast housing-only/no-rematch restriction is no
longer the national contract.

## Authoritative PSID specifications

The corrected rooms target is defined in
[`sa_rooms_first_birth_household_aligned_v1.do`](../../../data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do).
The header fixes the unit and event baseline (lines 6--12); the source load and
one-wave `ACTUALROOMS_` alignment are lines 43--59; the first biological birth
is the earliest `TYPE==1` year across `RELCHI1`--`RELCHI20` (lines 61--72).
The sample keeps current women age 18+, sex 2, reference person or
spouse/partner, positive `IW`, and a single current family unit, then keeps one
woman per household-year (lines 74--115). Full relationship histories define
confirmed zero-child controls and exclude untimed or unknown histories (lines
117--134). Event time is `K = year - first_child_year`, with -2 omitted and -6
estimated (lines 135--147). The estimator and contrast are lines 164--209:

```text
rooms ~ event indicators [pw=IW], absorb(person ID, survey year),
        covariates(age, education), vce(cluster person ID),
        Sun--Abraham cohort(first_child_year),
        controls = confirmed zero-child women.
target = coefficient(+3) - coefficient(-1), using the full covariance.
```

The durable receipt is
[`target_receipt.csv`](../../../data/psid_followup_mar2026/output/sa_rooms_first_birth_household_aligned_v1/target_receipt.csv):
49,457 estimation woman-years, 4,112 women, 16,031 confirmed never-treated
woman-years, estimate `0.7202462623815278`, clustered SE
`0.0852600513385958`, and pre-event mean rooms `4.792846924194579`. The
outcome is the aligned interview-level `ACTUALROOMS_`; vintage-specific
non-room codes are removed before estimation. This is a descriptive/calibration
mapping and has documented pre-event movement; it is not a causal claim.

The saved primary ownership replication is
[`sa_replication_own_only.do`](../../../data/psid_followup_mar2026/sa_replication_own_only.do).
It loads `HOMEOWN` and recodes 2 to renter/0 and 3 to missing (lines 21--30),
uses only `RELCHI1BYEAR` as the first-birth date, drops missing first-birth
years and births before panel entry (lines 38--47), and uses the same event
window with -2 omitted (lines 49--64). Its model is lines 78--81:

```text
own ~ event indicators, absorb(survey year),
     covariates(age, education), vce(cluster person ID),
     cohort(first_child_year), control_cohort(last observed cohort).
```

Although the script loads `IW` and defines a `$weight` macro, the actual
`eventstudyinteract` command has no `[pw=IW]`; the saved ownership result is
therefore an unweighted year-FE replication. The log reports 252,343
observations and 23,761 person clusters. The saved estimates are
[`own_f_c_y_all_repl_estimates.dta`](../../../data/psid_followup_mar2026/output/sa_replication/own_f_c_y_all_repl_estimates.dta)
and the event-study figure is in the same output directory. At -1 the saved
coefficient is `-0.014130468480289`; at +3 it is `0.00709169870242476`, so the
point difference is `0.0212221671827137` (2.12 percentage points). The script
does not save the covariance between those coefficients, so an exact contrast
SE or CI must not be reported from this output. The saved -2 pre-birth mean is
`0.503281533718109`.

The ownership object is not the corrected rooms object: it has no person fixed
effect, no probability weight in the estimation command, no full-history
never-treated group, and a narrower first-child field. Any ACS comparison must
label these as distinct PSID specifications rather than claiming an identical
replication.

## ACS contract and feasibility

The frozen Kleven matching contract is in
[`source_contract.json`](source_contract.json), lines 11--32: first-birth age
25--45, event window -5 through +10, reference -2, matching on `doiy`, age,
gender, education, marital status, race, and state, with original cleaned
childless CPS donors from 1995--1999. The persisted v5 panel is the matching
object; the housing adapter must not rematch or use housing variables as
matching covariates. Its event-time support label is `doiy - binned t_es_lw`,
not a biological birth year.

The authoritative housing adapter is
[`estimate_first_birth_housing.R`](estimate_first_birth_housing.R), lines
321--399. It joins on `(YEAR,SAMPLE,SERIAL,PERNUM)`, retains `ROOMS`,
`BEDROOMS`, and `OWNERSHP` as separate interview-level outcomes, uses the
author matching weight, fits one state/event specification by gender,

```text
outcome ~ state + event:state | age_factor + doiy_factor,
weights = author matching weight,
primary variance = source-household cluster;
heteroskedastic variance is a sensitivity.
```

The +3 minus -1 contrast is computed from the fitted covariance (lines
277--319). This ACS estimator is a matched repeated-cross-section diagnostic;
it is not numerically the PSID Sun--Abraham estimator.

The full national housing pull must freeze these fields before launch:

```text
YEAR SAMPLE SERIAL CBSERIAL HHWT CLUSTER STATEFIP PUMA STRATA GQ
OWNERSHP OWNERSHPD ROOMS BEDROOMS PERNUM PERWT MOMLOC POPLOC SPLOC
NCHILD NCHLT5 ELDCH YNGCH RELATE SEX AGE MARST FERTYR RACE EDUC
```

These housing and matching fields are present in the local `extract27.dta`
definition at
[`extractor27.do`](../../../data/Spatial_aggregate_withmicrodata/raw_data/extractor27.do):
identifiers and weights are lines 8--18 and 35--43, housing fields are
lines 23--34, and roster/demographic fields are lines 44--58. The local
extract is the housing bridge, not the author matching input. The author raw
ACS/CPS inputs are staged separately at
`/scratch/td2248/projects/kleven_acs_pilot_20260917/inputs/ACS/raw_acs.RData`
and `.../inputs/CPS/raw_cps.RData`, with the unchanged cleaner/matcher files
under `/scratch/td2248/projects/kleven_acs_pilot_20260917/vendor/`.

The archived author cleaner proves that raw `hispan` is used to construct the
four-category `race`/`race.num` matching variable
(`clean_acs.R` lines 88--119 inside
[`replication_code.zip`](../../../literature/kleven_pseudo_event/replication_code.zip)).
`BIRTHQTR` is not referenced by the archived cleaner or matcher. The national
source inventory must nevertheless record both names explicitly: `hispan` is a
required raw input for the author race recode, while `BIRTHQTR` is an audit
field requested by the current full-pull contract and must either be present or
be recorded as absent. It must not be substituted for `hispan`, and the absence
of `BIRTHQTR` must not silently change the matching cells. The national rematch
must run the unchanged `clean_acs.R`, `clean_cps.R`, `setup.R`, and
`matching.R` lineage, then attach housing outcomes by a proven full key.

The bridge key is `(YEAR,SAMPLE,SERIAL,PERNUM)`. If the author raw object does
not retain `SAMPLE`, the driver must stop and produce a field-level receipt;
it may not infer sample from year or serial. A successful bridge must report
one-to-one key uniqueness, matched/unmatched rows by ACS/CPS provenance and
state, and zero changes to author matching weights, event time, or donor roles.

The verified extract27 source fingerprint is recorded in
[`source_audit_manifest.json`](output/source_audit_extract27_20260919/source_audit_manifest.json):
9,919,999,546 bytes and SHA-256
`edb1afe53d4b6e6c5c5b8075bb83b81e1569c3cd9b619fe030af2fba0d33324e`, with
unique source keys `(YEAR,SAMPLE,SERIAL,PERNUM)` for the audited extract. The
audit found 2,190,987 unique key intersections for the prepared New England
packet and exact source year/sample counts for 2005--2019. That is six-state
New England
(`CT, ME, MA, NH, RI, VT`), not the full Census Northeast and not national
panel concordance. The current adapter therefore proves the join only for the
verified packet; the digest does not establish national key uniqueness or
national panel coverage. The national rematch must produce its own national
key-coverage receipt and source-household support table before housing fitting.

Housing coding must retain raw codes and missingness. The current adapter uses
valid `ROOMS` 1--27 and 30 with primary cap 9 and unknown literal code 28;
valid `BEDROOMS` 1--22 transformed by `x-1` with primary cap 5 and unknown 23;
and `OWNERSHP` 1/2 mapped to 1/0 with 3/9 unknown. The national report must
show year-specific top-code support, because the meaning of a literal rooms 9
changes over the source vintages. No row, person, weight, or unknown housing
code may be dropped merely because its outcome is unavailable.

## Freezeable national plan (no submission yet)

The national ACS driver should be a separate wrapper around the unchanged author
cleaner/matcher and the existing housing estimator. It must not alter the
vendor matching logic or estimator interface. Its stages are:

1. **Author-source gate.** On Torch, inventory the two existing national raw
   objects and vendor code, assert that `hispan`, all frozen match fields,
   identifiers, weights, and the explicit `BIRTHQTR` status are known, and write
   row/year/sample/schema receipts. Do not reload or alter the source archive.
2. **Author matching gate.** Run the unchanged national `clean_acs.R`,
   `clean_cps.R`, `setup.R`, and `matching.R` path on a bounded state smoke
   first. The smoke must exercise the actual source loader, race construction
   from `hispan`, exact demographic bins, all-ties matching with replacement,
   donor weights, and event-time construction. Save the resulting panel and
   source-key receipt before scaling.
3. **National panel gate.** Run the same code on all states and preserve the
   author provenance classes (original CPS, relabeled CPS, true ACS), matching
   weights, donor roles, `doiy`, `t_es_lw`, and source identifiers. Bridge the
   resulting ACS rows to extract27 housing records on the proven four-part key;
   report matched/unmatched rows by source category, year, sample, state,
   gender, and event time. Require one-to-one source-key joins and
   source-household clusters `(source_origin,YEAR,SAMPLE,SERIAL)`.
4. **Housing gate.** Attach the three raw outcomes, write code-validity and
   missingness tables before fitting, and require positive outcome-valid support
   in every requested pooled gender/event cell. Keep the pooled national fit as
   the primary estimate; state/regional tables are heterogeneity diagnostics,
   not unweighted averages of state estimates.
5. **Fit and receipts.** Estimate rooms and ownership separately with the
   existing adapter, saving raw coefficients, full covariance, matching weights,
   source-household cluster counts, +3/-1 contrasts, and a final manifest before
   rendering figures. Stop the stage if any source, key, sample, weight, or
   support fingerprint differs.

The bounded smoke should use the exact author national loader and matcher on a
manageable state slice, then attach a small extract27 slice containing matched,
unmatched, unknown-code, and repeated-household cases. It must exercise source
loading, `hispan` race construction, matching, four-part key bridge, coding,
support gate, and checkpoint writing. Only after that receipt is reviewed
should the national allocation be submitted. The historical README reports the
author's complete pipeline at roughly 63 hours on a 1.5 TB-RAM server; the
national rematch therefore cannot be budgeted as the earlier 4-CPU/64-GB
Northeast fit. A safe staged estimate is 8 CPUs/128 GB for the smoke and,
conditional on measured peak memory, one 16--32 CPU allocation with 512 GB and
an overnight (12--16 hour) wall limit for national clean/match/bridge/fit.
The driver must write stage checkpoints and peak-memory/time receipts so the
lead can adjust the full allocation before submission. Use a unique output
directory, progress heartbeats, per-stage/per-outcome checkpoints, and
fail-stop behavior; no automatic retry. The separate full-sample PSID run
remains owned by its assigned worker and must not be duplicated here.

The decision gate before launch is narrow: review the bounded national smoke,
its measured memory/time, and the full-key bridge receipt. If they pass, launch
the one national ACS clean/match/bridge/housing allocation under the measured
resource request. Until then, the valid status is “national specification and
source paths mapped; national run not submitted.”

The first implementation artifact is
[`run_national_acs_source_stage.R`](run_national_acs_source_stage.R), with the
45-minute 8-CPU/128-GB launcher
[`run_national_acs_source_stage.sbatch`](run_national_acs_source_stage.sbatch).
It loads ACS and CPS sequentially, writes one raw partition per requested
`STATEFIP`, preserves original column order, records `hispan`/`BIRTHQTR` status,
and can partition a narrow housing packet by the same state. Its receipt is
explicitly `PARTITION_SMOKE_READY` with `no_full_matching=true`; it is not a
claim that the vendor matcher ran. A local synthetic end-to-end
`partition_smoke` passed with ACS/CPS state partitions and housing packet
partition receipts. The next implementation step after lead review is to run
the unchanged vendor cleaner/matcher once per partition, carry the source key
and donor provenance through the existing v6 adapter, and prove protected
column equality on the existing New England/Vermont checkpoint before any
national fit.
