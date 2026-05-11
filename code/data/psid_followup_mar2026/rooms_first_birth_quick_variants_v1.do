clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/rooms_first_birth_quick_variants_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/rooms_first_birth_quick_variants_v1.log", replace text

global weight IW

tempname posth
postfile `posth' ///
    str24 variant str40 description ///
    double mean_rooms_pre mean_rooms_post3 mean_rooms_post5 ///
    double mean_change_post3 mean_change_post5 ///
    double N_post3 N_post5 using "`outdir'/rooms_first_birth_quick_variants_v1.dta", replace

local variants "all drop_twins censor_at_second no_second_by3 no_second_by5"

foreach variant in `variants' {
    di as text "Running quick first-birth rooms variant: `variant'"
    use ID year AGEREP DEATHYEAR RELCHI1BYEAR RELCHI2BYEAR ACTUALROOMS_ ${weight} using "`dta'", clear

    drop if year > DEATHYEAR
    drop if AGEREP < 18

    bysort ID: egen first_child_year = min(RELCHI1BYEAR)
    bysort ID: egen second_child_year = min(RELCHI2BYEAR)
    drop if missing(first_child_year)

    bysort ID: egen year_entry = min(year)
    drop if first_child_year < year_entry
    drop year_entry

    rename ACTUALROOMS_ rooms
    xtset ID year

    local desc "Baseline first-birth quick moment"

    if inlist("`variant'", "drop_twins", "censor_at_second", "no_second_by3", "no_second_by5") {
        drop if !missing(second_child_year) & second_child_year == first_child_year
        local desc "Exclude multiple births at first birth"
    }
    if "`variant'" == "censor_at_second" {
        drop if !missing(second_child_year) & year >= second_child_year
        local desc "Exclude multiples and censor observations at second birth"
    }
    if "`variant'" == "no_second_by3" {
        drop if !missing(second_child_year) & second_child_year <= first_child_year + 3
        local desc "Exclude multiples and require no second birth by +3"
    }
    if "`variant'" == "no_second_by5" {
        drop if !missing(second_child_year) & second_child_year <= first_child_year + 5
        local desc "Exclude multiples and require no second birth by +5"
    }

    gen K = year - first_child_year

    preserve
        keep if K == -2
        keep ID first_child_year rooms ${weight}
        rename rooms rooms_pre
        rename ${weight} iw_pre
        save "`outdir'/tmp_rooms_pre_`variant'.dta", replace
    restore

    preserve
        keep if inrange(K, 0, 3)
        bysort ID first_child_year: egen rooms_post3 = mean(rooms)
        bysort ID first_child_year (year): keep if _n == 1
        keep ID first_child_year rooms_post3
        merge 1:1 ID first_child_year using "`outdir'/tmp_rooms_pre_`variant'.dta", keep(match) nogen
        gen rooms_change_post3 = rooms_post3 - rooms_pre if !missing(rooms_post3, rooms_pre)
        keep if !missing(rooms_change_post3)
        quietly summarize rooms_pre [aw=iw_pre]
        local mean_rooms_pre = r(mean)
        quietly summarize rooms_post3 [aw=iw_pre]
        local mean_rooms_post3 = r(mean)
        quietly summarize rooms_change_post3 [aw=iw_pre]
        local mean_change_post3 = r(mean)
        quietly count
        local N_post3 = r(N)
    restore

    preserve
        keep if inrange(K, 0, 5)
        bysort ID first_child_year: egen rooms_post5 = mean(rooms)
        bysort ID first_child_year (year): keep if _n == 1
        keep ID first_child_year rooms_post5
        merge 1:1 ID first_child_year using "`outdir'/tmp_rooms_pre_`variant'.dta", keep(match) nogen
        gen rooms_change_post5 = rooms_post5 - rooms_pre if !missing(rooms_post5, rooms_pre)
        keep if !missing(rooms_change_post5)
        quietly summarize rooms_post5 [aw=iw_pre]
        local mean_rooms_post5 = r(mean)
        quietly summarize rooms_change_post5 [aw=iw_pre]
        local mean_change_post5 = r(mean)
        quietly count
        local N_post5 = r(N)
    restore

    erase "`outdir'/tmp_rooms_pre_`variant'.dta"

    post `posth' ("`variant'") ("`desc'") ///
        (`mean_rooms_pre') (`mean_rooms_post3') (`mean_rooms_post5') ///
        (`mean_change_post3') (`mean_change_post5') ///
        (`N_post3') (`N_post5')
}

postclose `posth'

use "`outdir'/rooms_first_birth_quick_variants_v1.dta", clear
order variant description mean_rooms_pre mean_rooms_post3 mean_rooms_post5 mean_change_post3 mean_change_post5 N_post3 N_post5
sort variant
export delimited using "`outdir'/rooms_first_birth_quick_variants_v1.csv", replace

di as result "Quick first-birth rooms variants complete."
di as result "Summary: `outdir'/rooms_first_birth_quick_variants_v1.csv"

log close _all
