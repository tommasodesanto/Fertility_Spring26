# ACS date diagnostic — saved ownership aggregates and pinned room cache

Generated 2026-09-09 for lead review. **Diagnostic only: no empirical target, group definition, model, weight, uncertainty or target fingerprint changed.** Ownership rows use saved aggregates only. Under a separately authorized follow-up, room rows use one single-process R extraction from the110MB pinned cache, with one data.table thread. No model, raw-data processing,9.92GB raw-file hashing, bootstrap or package installation was run.

`acs_date_comparison.csv` covers the four active pooled ACS rows. All four original-group pooled-versus2023 comparisons are now available. Ownership and rooms samples remain distinct. No2023 empirical standard errors are supplied.

## Available ownership comparisons

Both periods preserve MMS42-metro household heads of either sex, exact ages30–55, HHWT, owner/renter observations with positive rooms, and the original DUE structures UNITSSTR3:10. The family gap remains oldest co-resident own child under4 (NCHILD>0, observed ELDCH!=99, ELDCH<4) minus NCHILD0, including empty nesters. No all-parent gap is substituted.

The family gap is copied from the previously saved pooled/2023 diagnostic. Overall ownership is recovered from its exhaustive NCHILD>0 and NCHILD0 partition using

    ownership = (parent_HHWT * parent_ownership + nochild_HHWT * nochild_ownership) / (parent_HHWT + nochild_HHWT).

The saved all-housing validation establishes zero unknown NCHILD observations for both prime-age periods; the DUE/head subset inherits that completeness. The pooled reconstructed partition exactly matches the authoritative ownership table's1806067 records and207284892 HHWT total. Its ownership rate reproduces the original pooled rate within1e-14. Source rates are printed decimals, so the reconstructed overall rate has corresponding last-decimal rounding. The2023 partition contains243735 records and29443094 HHWT.

The active rounded scalars remain0.16766167 and0.575472. Their active scales are synthetic5% scales, not empirical standard errors; they have not been recalculated for2023. No2023 uncertainty estimate is available here. The rooms rows preserve their original synthetic weights and separately list the already audited pooled bootstrap SEs; those SEs are not used in the objective.

## Completed room comparison: exact source and formulas

Extraction script: `extract_rooms_dates.R`; unrounded R checks and cache receipt: `rooms_date_extraction_receipt.txt`; group-specific counts, HHWT sums and mean rooms: `rooms_date_extraction.csv`. Runtime7.421seconds, one cache read, one data.table thread. Locale/package-build warnings were nonfatal; R exited0.

The exact authoritative masks come from `code/data/moment_standard_errors/build_active_acs_room_target_receipt.R:131–150`; cached `parent_u18` semantics come from `build_moment_bootstrap_se.R:427`.

- Mean rooms: `age>=18 & age<=85`; sum(HHWT*ROOMS)/sum(HHWT).
- Family gap: `age>=30 & age<=55 & parent_u18==TRUE & NCHILD>=1`; weighted mean rooms among NCHILD>=3 minus weighted mean among NCHILD1–2. Cached parent_u18 is verified equal to `NCHILD>0 & !is.na(YNGCH) & YNGCH!=99 & YNGCH<18`. NCHILD includes co-resident own children of any age.
- Both year2023 comparisons apply `year==2023` before these unchanged masks. No ownership DUE-structure restriction, sex restriction or new child definition was added. Base cache already imposes household head PERNUM1, owner/renter, positive rooms/HHWT and matched MMS geography.

Both pooled points and sample counts/weight totals reproduce the authoritative receipt before any annual output is written. Mean rooms absolute point difference is2.66453525910038e-15, family-gap difference2.22044604925031e-16, each against1e-10 tolerance. The2023 mean-rooms sample has558463 records and61250563 HHWT; the family-gap sample has122545 records and14521384 HHWT, partitioned into29625 large-family records/3740723 HHWT and92920 small-family records/10780661 HHWT. Both periods retain42metros.

The single read input is `code/data/moment_standard_errors/cache/acs_analysis_samples.rds`, object `targets`,4103889records. Its110425398byte SHA was recomputed and equals `0eae9aaed4e1d3b9655235be967378c1d389f29ec9b64f005a811c2f0ced7df0`. The original point/receipt builders were inspected, not executed.

## Verification and source hashes

All small-file SHA256 values below were recomputed during this pass. The preexisting ownership comparison CSV also matched the SHA already reported in the earlier ownership follow-up. Every derived ownership value was recomputed and compared to its pooled reference where available.

- `output/model/e5f_simple_fertility_tax_transition_20260908a/households_children_review/ownership_target_comparison.csv`: `544b4519ef5c77a041fa8055c32926bdc1af584c9a171d17a66913a636a04c91`
- `output/model/e5f_simple_fertility_tax_transition_20260908a/households_children_review/ownership_target_comparison.R`: `f62a77b0c57351a60997ff3d3b1d2f4ce8676458f9ce9337df76778387439408`
- `output/model/e5f_simple_fertility_tax_transition_20260908a/households_children_review/acs_validation.csv`: `b46b15b54ee9eed81adba7f70b17fbde96484f4ace6adbece4e1aaa283aa00fa`
- `output/model/e5f_simple_fertility_tax_transition_20260908a/households_children_review/acs_validation_receipt.json`: `6b6d670b6b36827fcee393d774cc520d26ec072d7885e33eda9d341211ac213f`
- `code/data/mms_center_periphery/output_ownership_audit/acs_ownership_window_targets.csv`: `25fe9b3cea20531bb0cae5e5594279ca6ec0c3ce5f52b176e6a222241b97a507`
- `code/data/mms_center_periphery/audit_ownership_targets.R`: `7d91b53de1f7a25df399ab30e0f65b0d95b6c94cf66d903c82ce581cd57077c5`
- `code/data/moment_standard_errors/output_active_acs_room_target_receipt_20260817/target_receipt.csv`: `c37c100a1ef3ed8c6e6d2449d4c7ab2c37d5e8a5c29da68a5e6f34e596e97fd3`
- `code/data/moment_standard_errors/output_active_acs_room_target_receipt_20260817/provenance.csv`: `4fe50bd66675d33b90990b1b036b3d9d91ec9b2f08e7dec906b762e74450ea56`
- `code/data/moment_standard_errors/build_active_acs_room_target_receipt.R`: `e4307c0ff108b83ff12c41f808d5ad51ce65e555972cdd1315e25f692c5476d2`
- `code/data/moment_standard_errors/build_moment_bootstrap_se.R`: `7731b89c94dcdc2c1dab2408e10bd209e8ee2992d155cbefb782c4fa940c9248`
- `tmp/e5f_matched_pf/code/model/intergen_eqscale_seq_optimized/e5_target_provenance.csv`: `38e6507fdc81ae264d54ab02d9fb6c824114a3da85cfc67c901eca819bf6e78b`
- `output/model/e5f_matched_pf_20260909a/meeting_receipts/acs_date_diagnostic/acs_date_comparison.csv`: `189c0cf939428b9642a87542cccbf439e51172a0d5ba21ed3b2f058300008139`
- `output/model/e5f_matched_pf_20260909a/meeting_receipts/acs_date_diagnostic/extract_rooms_dates.R`: `d478a8796628f748d8fcbfb49036b8c740b48704e075aa41cc467f1638c53055`
- `output/model/e5f_matched_pf_20260909a/meeting_receipts/acs_date_diagnostic/rooms_date_extraction.csv`: `2223a7469b1d8a34664a0129f17609e7189601cbfd472f73e27e519709a6b78f`
- `output/model/e5f_matched_pf_20260909a/meeting_receipts/acs_date_diagnostic/rooms_date_extraction_receipt.txt`: `5df44ef7c441a65e0178bbfe7454a9804f126f127d5b406a3659d8d4a5b84b23`
