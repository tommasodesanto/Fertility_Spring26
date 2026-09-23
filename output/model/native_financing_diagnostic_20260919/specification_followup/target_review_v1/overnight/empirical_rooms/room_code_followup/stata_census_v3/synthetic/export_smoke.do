clear all
set more off
set processors 1
local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/room_code_followup/stata_census_v3/synthetic"
input ID rooms_raw_aligned year rooms_source_year K never_treated first_child_year REL
1 0 1990 1989 -1 0 1991 1
2 0 1990 1989 3 0 1987 2
3 98 1990 1989 . 1 . 1
4 99 1990 1989 . 1 . 1
5 4 1990 1989 1 0 1989 1
end
* No personal IDs are exported. Event time/cohort 999999/0 identify confirmed
* childless controls; REL 1/2 mean reference person/spouse-partner.
preserve
    keep if inlist(rooms_raw_aligned, 0, 98, 99)
    gen int audit_event_time = K
    replace audit_event_time = 999999 if never_treated
    gen int audit_cohort = first_child_year
    replace audit_cohort = 0 if never_treated
    gen byte audit_reporter_role = REL
    egen byte audit_person_tag = tag(rooms_raw_aligned year rooms_source_year ///
        audit_event_time audit_reporter_role audit_cohort ID)
    collapse (count) observations=ID ///
        (sum) people=audit_person_tag, ///
        by(rooms_raw_aligned year rooms_source_year audit_event_time ///
           audit_reporter_role audit_cohort)
    rename rooms_raw_aligned room_code
    rename year aligned_interview_year
    rename audit_event_time event_time
    rename audit_reporter_role reporter_role_code
    rename audit_cohort first_birth_cohort_code
    order room_code aligned_interview_year rooms_source_year event_time ///
        reporter_role_code first_birth_cohort_code observations people
    export delimited using "`outdir'/suspect_code_cells.csv", replace
restore

* Anonymized counts and receipt: prefit and singleton-pruned totals, plus
* overall suspect-code totals to compare with the exact replay diagnostics.
preserve
    keep if inlist(rooms_raw_aligned, 0, 98, 99)
    contract rooms_raw_aligned
    rename rooms_raw_aligned room_code
    rename _freq retained_rows
    export delimited using "`outdir'/suspect_code_totals.csv", replace
restore

di "FULL_CENSUS_EXPORT_SMOKE_PASSED"
