# Ownership-target follow-up: what was retained, and what is actually compared

Read-only follow-up, 9 September 2026. The author has reaffirmed the maintained18-year child-maturation mapping; this note does not propose reopening it. No target, weight, parameter, calibration or source was changed.

## The current calibrated number and its date

Selected external case10 has ownership-gap model value **16.10614646 percentage points**, against retained target **16.766167 percentage points**. Its recorded gap is −0.66002054pp, weight14229.590956, loss contribution0.61987956. The case receipt itself stores the model value, case_id10 and total loss23.791955301663187. The summary SHA `19022673b0f03576d7a211e5d2b79c6e664a998b103a36dab34016a15782a8e0` matches the live tax transition contract.

This ownership statistic is evaluated on the **2023 historical transition cross-section**, not a new steady state: selected summary `target_measurements.remaining_targets` says so explicitly. The scalar empirical target is nevertheless pooled ACS2012–2023. Therefore “2023 is on the modeled path” and “the target is not2023-only” are both true. The CSV's candidate label task_001 is the within-case candidate label, while the enclosing selected receipt is case10. Source: [selected summary](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_simple_fertility_overnight_20260907a/morning_review/selected/summary.json), [selected case receipt](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_simple_fertility_overnight_20260907a/morning_review/selected/case_receipt.json).

Provenance check resolved: the local CSV uses LF line endings whereas the original uses CRLF. Serializing its exact eight columns and all rows with CRLF reproduces the receipt SHA `484cbafeafe254d5de3d98f6a9835174474689e10de3f7e39aae378518bd405d` exactly. All twelve model values equal the receipt values exactly and contributions sum to23.791955301663187. No columns or scientific values changed; this is solely newline normalization.

## Definitions and empirical alternatives already available

The retained empirical16.766167pp is ownership among heads whose **oldest co-resident own child is under four**, minus ownership among heads with **no co-resident own child**. The latter includes empty nesters. The model instead compares **any current dependent-child household** with **never-parent households**, using ages30–55 under its existing age mapping. A model variable called `newparent` does not make it a recent-birth group. Sources: [code/data/mms_center_periphery/audit_ownership_targets.R:67](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/audit_ownership_targets.R:67), [code/model/intergen_eqscale_seq_optimized/solver.py:6310](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_fertility_nest_compute_20260907a/code/model/intergen_eqscale_seq_optimized/solver.py:6310), [code/model/intergen_eqscale_seq_optimized/solver.py:6485](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_fertility_nest_compute_20260907a/code/model/intergen_eqscale_seq_optimized/solver.py:6485).

The same-sample all-parent alternative was already present in the original ACS audit; it was not necessary to invent or estimate a new model to find it. A bounded one-thread reread of the110MB pinned cache reproduced both old rows and extracted the under18 variant below. Cache SHA was reverified; no9.92GB raw file was read. All rows use household heads, HHWT, MMS42-metro geography, positive rooms and owner/renter observations, DUE structures3:10, exact empirical head ages30–55, either sex, with NCHILD0 as the common control. They are descriptive comparisons, **not newly selected targets**.

| Empirical parent definition | Pooled2012–2023 gap (pp) | 2023-only gap (pp) |
|---|---:|---:|
| Oldest co-resident child under4 — retained target | 16.766167 | 21.303955 |
| Any co-resident own child, any age | 15.211698 | 18.739014 |
| At least one co-resident own child under18 | 15.015514 | 18.975823 |

The pooled all-parent value15.21169765pp is independently in [code/data/mms_center_periphery/output_ownership_audit/acs_ownership_window_targets.csv:1](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_ownership_audit/acs_ownership_window_targets.csv:1). These alternatives still do not replicate the model's never-parent control: ACS NCHILD0 is absence of co-resident children, not absence of any previous birth. Nor does a snapshot under18 flag reproduce the model's18-year-mean geometric dependency process exactly.

## Chronology: awareness is documented; explicit group-approximation acceptance was not found

- **May28:** the generated ACS ownership audit already labels the16.766pp target `newparent_minus_nochildren_30_55`; its same-sample table also reports all-parent15.212pp. [code/data/mms_center_periphery/output_ownership_audit/ACS_MMS_OWNERSHIP_TARGET_AUDIT.md:3](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_ownership_audit/ACS_MMS_OWNERSHIP_TARGET_AUDIT.md:3).
- **July2:** audit A10 explicitly flags both broad model parenthood versus recent empirical parenthood and never-parent versus no-child controls. It proposes clarification or remapping; it does not document author adoption. The July3 lead correction withdraws a different parity/top-coding finding, not this ownership finding. [output/model/fable_size_mapping_audit_20260701/CODE_AUDIT_REPORT.md:94](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/fable_size_mapping_audit_20260701/CODE_AUDIT_REPORT.md:94).
- **August5:** the author-directed independent child-count maturation repair is explicit; that establishes authority for the18-year-mean process. It does not itself record acceptance of the ownership target's mismatched groups. See the August5 author-directed section of CALIBRATION_STATUS.md and [code/model/intergen_eqscale_seq_optimized/e5_maturation_repair.py:1](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/tmp/e5f_fertility_nest_compute_20260907a/code/model/intergen_eqscale_seq_optimized/e5_maturation_repair.py:1).
- **September7:** the report label was changed to parent/nonparent, with model parents defined as any child at home and nonparents as no previous birth. The canonical status explicitly left empirical reconciliation outstanding. Neither the active decision ledger nor the targeted earlier decision ledger contains an explicit switch of the retained scalar to15.212pp, or an explicit acceptance of this group approximation. This is a scoped absence of documentation, not proof that the author never discussed it elsewhere.

This ownership target is an **ACS cross-sectional gap**, not the PSID Sun–Abraham first-birth room event study. The latter is the separate first-birth housing-response target0.720246rooms and remains unchanged. Similar “new parent” language does not connect their estimators or identifying samples.

## Measurement routes for an author decision, without choosing one here

1. Keep the recent-parent empirical object and define a matching model measurement of child recency; the current count-only state cannot recover the oldest child's exact age. A birth-event response is a different object and should not silently substitute for this cross-section.
2. Keep the current dependent-child economic mechanism and explicitly align both empirical parent and control groups to a documented current-child comparison. The existing15.212/15.016pp diagnostic rows help assess that option, but simply replacing the scalar leaves the never-parent/control mismatch unresolved.
3. Deliberately retain the approximation and state it as such, after the author considers its identifying role. That is different from having verified a matched target. Any target redefinition requires the complete target/weight fingerprint and identification review; nothing here authorizes it.

The separate choice between pooled2012–2023 and2023-only empirical moments also needs to be explicit. The2023-only diagnostic should not silently replace the retained pooled target just because the modeled endpoint is2023.

Extraction source SHA: `f62a77b0c57351a60997ff3d3b1d2f4ce8676458f9ce9337df76778387439408`. Aggregate CSV SHA: `544b4519ef5c77a041fa8055c32926bdc1af584c9a171d17a66913a636a04c91`. Both pooled previously published rows reproduced within1e-8 and1e-12 respectively. No microdata are included in outputs.
