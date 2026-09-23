clear all
set more off
set processors 4
version 17.0

* Read-only data census: clone of the approved rooms-builder prefix through K.
* This stops before eventstudyinteract. It uses only local PSID data and writes
* only to this task's separate stata_census folder.
local project "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local source "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`project'/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/room_code_followup/stata_census_v3"
cap mkdir "`outdir'"
capture confirm file "`outdir'/suspect_code_cells.csv"
if !_rc {
    di as error "Refusing to overwrite prior census output: `outdir'/suspect_code_cells.csv"
    exit 602
}
capture log close _all
log using "`outdir'/room_code_census.log", replace text

* Original first-birth rooms builder: source import, one-observed-row lag,
* person-level biological first birth and history, household FID count,
* year-specific room missing codes, eligibility, and household dedup.
use ID FID HHID year AGEREP EDUYEAR SEX DEATHYEAR REL CURRENT RELCHIREP ///
    RELCHI1TYPE RELCHI1BYEAR RELCHI2TYPE RELCHI2BYEAR ///
    RELCHI3TYPE RELCHI3BYEAR RELCHI4TYPE RELCHI4BYEAR ///
    RELCHI5TYPE RELCHI5BYEAR RELCHI6TYPE RELCHI6BYEAR ///
    RELCHI7TYPE RELCHI7BYEAR RELCHI8TYPE RELCHI8BYEAR ///
    RELCHI9TYPE RELCHI9BYEAR RELCHI10TYPE RELCHI10BYEAR ///
    RELCHI11TYPE RELCHI11BYEAR RELCHI12TYPE RELCHI12BYEAR ///
    RELCHI13TYPE RELCHI13BYEAR RELCHI14TYPE RELCHI14BYEAR ///
    RELCHI15TYPE RELCHI15BYEAR RELCHI16TYPE RELCHI16BYEAR ///
    RELCHI17TYPE RELCHI17BYEAR RELCHI18TYPE RELCHI18BYEAR ///
    RELCHI19TYPE RELCHI19BYEAR RELCHI20TYPE RELCHI20BYEAR ///
    ACTUALROOMS_ IW using "`source'", clear

sort ID year
by ID: gen double rooms = ACTUALROOMS_[_n-1] if _n > 1
by ID: gen double rooms_source_year = year[_n-1] if _n > 1
gen double rooms_raw_aligned = rooms

gen double first_biological_birth_candidate = .
forvalues child = 1/20 {
    replace first_biological_birth_candidate = RELCHI`child'BYEAR ///
        if RELCHI`child'TYPE == 1 & !missing(RELCHI`child'BYEAR) & ///
        (missing(first_biological_birth_candidate) | ///
         RELCHI`child'BYEAR < first_biological_birth_candidate)
}
bysort ID: egen double first_child_year = min(first_biological_birth_candidate)
bysort ID: egen double max_children_reported = max(RELCHIREP)
drop first_biological_birth_candidate

egen byte current_fid_tag = tag(HHID year FID) ///
    if CURRENT == 1 & !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0
bysort HHID year: egen int n_current_fids = total(current_fid_tag)
drop current_fid_tag

replace rooms = . if year <= 1984 & rooms == 9
replace rooms = . if inrange(year, 1985, 1993) & rooms == 99
replace rooms = . if year >= 1994 & inlist(rooms, 98, 99)

drop if year > DEATHYEAR
keep if CURRENT == 1
keep if !missing(AGEREP) & AGEREP >= 18
keep if !missing(IW) & IW > 0
keep if SEX == 2 & inlist(REL, 1, 2)
keep if !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0
keep if n_current_fids == 1

gen byte household_priority = REL != 1
sort HHID year household_priority ID
by HHID year: gen int household_women_before_dedup = _N
by HHID year: keep if _n == 1
isid HHID year

gen double rooms_alignment_gap_years = year - rooms_source_year
replace rooms = . if !inlist(rooms_alignment_gap_years, 1, 2)
drop rooms_alignment_gap_years
drop if missing(rooms) | missing(EDUYEAR)

bysort ID: egen double year_entry = min(year)
drop if !missing(first_child_year) & first_child_year < year_entry

gen byte untimed_known_parent = ///
    missing(first_child_year) & !missing(max_children_reported) & ///
    max_children_reported > 0
gen byte unknown_child_hist = ///
    missing(first_child_year) & missing(max_children_reported)
drop if untimed_known_parent | unknown_child_hist
gen byte never_treated = missing(first_child_year) & max_children_reported == 0
assert never_treated == 1 if missing(first_child_year)
gen double K = year - first_child_year

quietly count
local prefit_rows = r(N)
egen byte prefit_id_tag = tag(ID)
quietly count if prefit_id_tag
local prefit_people = r(N)
quietly count if never_treated
local prefit_control_rows = r(N)
egen byte prefit_control_tag = tag(ID) if never_treated
quietly count if prefit_control_tag
local prefit_controls = r(N)
drop prefit_id_tag prefit_control_tag

* Reconstruct only the documented iterative person/year FE singleton removal.
* No regression or estimator command is called in this file.
local more = 1
local prune_iter = 0
local removed_person_ids = 0
local removed_years = 0
while `more' {
    sort ID year
    by ID: gen long audit_n_id = _N
    bysort year: gen long audit_n_year = _N
    quietly count if audit_n_id == 1
    local one_id_rows = r(N)
    quietly count if audit_n_year == 1
    local one_year_rows = r(N)
    if `one_id_rows' == 0 & `one_year_rows' == 0 {
        local more = 0
        drop audit_n_id audit_n_year
    }
    else {
        egen byte audit_id_tag = tag(ID) if audit_n_id == 1
        quietly count if audit_id_tag
        local one_ids = r(N)
        egen byte audit_year_tag = tag(year) if audit_n_year == 1
        quietly count if audit_year_tag
        local one_years = r(N)
        local removed_person_ids = `removed_person_ids' + `one_ids'
        local removed_years = `removed_years' + `one_years'
        local prune_iter = `prune_iter' + 1
        drop if audit_n_id == 1 | audit_n_year == 1
        drop audit_n_id audit_n_year audit_id_tag audit_year_tag
    }
}

quietly count
local retained_rows = r(N)
egen byte retained_id_tag = tag(ID)
quietly count if retained_id_tag
local retained_people = r(N)
quietly count if never_treated
local retained_control_rows = r(N)
egen byte retained_control_tag = tag(ID) if never_treated
quietly count if retained_control_tag
local retained_controls = r(N)
drop retained_id_tag retained_control_tag

assert `retained_rows' == 49457
assert `retained_people' == 4112
quietly count if rooms_raw_aligned == 0
assert r(N) == 100
quietly count if rooms_raw_aligned == 98
assert r(N) == 1
quietly count if rooms_raw_aligned == 99
assert r(N) == 0

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

clear
set obs 1
gen long prefit_rows = `prefit_rows'
gen long prefit_people = `prefit_people'
gen long prefit_control_rows = `prefit_control_rows'
gen long prefit_confirmed_controls = `prefit_controls'
gen long singleton_pruning_iterations = `prune_iter'
gen long singleton_person_ids_removed = `removed_person_ids'
gen long singleton_years_removed = `removed_years'
gen long singleton_pruned_rows = `retained_rows'
gen long singleton_pruned_people = `retained_people'
gen long singleton_pruned_control_rows = `retained_control_rows'
gen long singleton_pruned_confirmed_controls = `retained_controls'
gen str50 sample_status = "aggregate-equivalent; row-level e(sample) marker unavailable"
export delimited using "`outdir'/sample_receipt.csv", replace
log close _all
