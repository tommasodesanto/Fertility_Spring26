clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/sa_moving_debug_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/sa_moving_sample_audit_v1.log", replace text

cap which eventstudyinteract
cap which reghdfe

tempname fh
file open `fh' using "`outdir'/sa_moving_sample_audit_v1.csv", write replace
file write `fh' "sample,outcome,n_obs,n_ids,n_cohorts,min_cohort,max_cohort,n_eventvars,n_interactions" _n

foreach sample in own_original moving_1984 moving_all_years moving_window {
    use ID year AGEREP EDUYEAR SEX RELCHI1BYEAR MOVEDFREF_ WHYMOVED1_ DEATHYEAR ACTUALROOMS_ HOMEOWN IW using "`dta'", clear

    replace HOMEOWN = 0 if HOMEOWN == 2
    replace HOMEOWN = . if HOMEOWN == 3
    drop if year > DEATHYEAR
    drop if AGEREP < 18

    if "`sample'" == "moving_1984" keep if inrange(year, 1984, 2019)
    if "`sample'" == "moving_window" keep if inrange(year, 1984, 2019)

    rename MOVEDFREF_ movedthisyear
    replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
    replace movedthisyear = 0 if movedthisyear == 5
    rename HOMEOWN own
    rename ACTUALROOMS_ rooms

    gen move_any = movedthisyear
    gen moved_for_size = .
    replace moved_for_size = 1 if movedthisyear == 1 & WHYMOVED1_ == 3 & year > 1974
    replace moved_for_size = 0 if movedthisyear == 0 & year > 1974
    replace moved_for_size = 0 if movedthisyear == 1 & WHYMOVED1_ != 3 & !missing(WHYMOVED1_) & year > 1974

    gen first_child_year = RELCHI1BYEAR
    drop if missing(first_child_year)
    bysort ID: egen year_entry = min(year)
    drop if first_child_year < year_entry
    drop year_entry

    gen K = year - first_child_year
    if "`sample'" == "moving_window" keep if inrange(K, -7, 11)

    quietly summarize first_child_year
    gen lastcohort = first_child_year == r(max)

    cap drop L*event F*event
    forvalues l = 0/10 {
        gen L`l'event = K == `l'
    }
    gen L11event = (K > 10 & K != .)
    forvalues l = 1/5 {
        gen F`l'event = K == -`l'
    }
    gen F7event = (K < -6 & K != .)
    drop F2event

    foreach outcome in own rooms move_any moved_for_size {
        quietly count if !missing(`outcome', first_child_year, year, AGEREP, EDUYEAR, ID)
        local n_obs = r(N)
        capture drop tag_id
        egen tag_id = tag(ID) if !missing(`outcome', first_child_year, year, AGEREP, EDUYEAR, ID)
        quietly count if tag_id == 1
        local n_ids = r(N)
        drop tag_id
        quietly levelsof first_child_year if lastcohort == 0 & !missing(`outcome', first_child_year, year, AGEREP, EDUYEAR, ID), local(cohorts)
        local n_cohorts : word count `cohorts'
        local min_cohort : word 1 of `cohorts'
        local max_cohort : word `n_cohorts' of `cohorts'
        local n_eventvars = 18
        local n_interactions = `n_cohorts' * `n_eventvars'
        file write `fh' "`sample',`outcome',`n_obs',`n_ids',`n_cohorts',`min_cohort',`max_cohort',`n_eventvars',`n_interactions'" _n
    }
}

file close `fh'

di as result "Audit written: `outdir'/sa_moving_sample_audit_v1.csv"

log close _all
