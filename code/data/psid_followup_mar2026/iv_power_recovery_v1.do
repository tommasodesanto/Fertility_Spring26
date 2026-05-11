clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/iv_power_recovery_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/iv_power_recovery_v1.log", replace text

global weight IW

di as text "Step 1/5: Build parent-level fertility instruments (twins + same-sex)"
use ID year SEX using "`dta'", clear
drop if missing(ID, SEX)
bysort ID (year): keep if _n == 1
rename ID child_id
rename SEX child_sex
save "`outdir'/child_sex_map_tmp.dta", replace

preserve
    use "`outdir'/child_sex_map_tmp.dta", clear
    rename child_id child1_id
    save "`outdir'/child1_sex_map_tmp.dta", replace
restore

preserve
    use "`outdir'/child_sex_map_tmp.dta", clear
    rename child_id child2_id
    save "`outdir'/child2_sex_map_tmp.dta", replace
restore

use ID year RELCHIREP RELCHINUM RELCHI1ID RELCHI2ID RELCHI1BYEAR RELCHI2BYEAR ///
    RELCHI3BYEAR RELCHI4BYEAR RELCHI5BYEAR RELCHI6BYEAR RELCHI7BYEAR RELCHI8BYEAR ///
    RELCHI9BYEAR RELCHI10BYEAR RELCHI11BYEAR RELCHI12BYEAR RELCHI13BYEAR RELCHI14BYEAR ///
    RELCHI15BYEAR RELCHI16BYEAR RELCHI17BYEAR RELCHI18BYEAR RELCHI19BYEAR RELCHI20BYEAR ///
    using "`dta'", clear

bysort ID: egen last_nonmiss_year = max(cond(!missing(RELCHI1BYEAR), year, .))
keep if year == last_nonmiss_year & !missing(last_nonmiss_year)
drop last_nonmiss_year

forvalues k = 1/20 {
    gen sameyr`k' = (RELCHI`k'BYEAR == RELCHI1BYEAR) if !missing(RELCHI`k'BYEAR, RELCHI1BYEAR)
    replace sameyr`k' = 0 if missing(sameyr`k')
}
egen n_sameyear_firstbirth = rowtotal(sameyr1-sameyr20)
gen twin_firstbirth = (n_sameyear_firstbirth >= 2) if !missing(n_sameyear_firstbirth)
drop sameyr1-sameyr20

rename RELCHI1ID child1_id
rename RELCHI2ID child2_id

merge m:1 child1_id using "`outdir'/child1_sex_map_tmp.dta", keep(master match) nogen
rename child_sex child1_sex

merge m:1 child2_id using "`outdir'/child2_sex_map_tmp.dta", keep(master match) nogen
rename child_sex child2_sex

gen same_sex_first2 = (child1_sex == child2_sex) if !missing(child1_sex, child2_sex)
gen two_plus = (RELCHIREP >= 2) if !missing(RELCHIREP)
gen three_plus = (RELCHIREP >= 3) if !missing(RELCHIREP)

keep ID RELCHI1BYEAR RELCHI2BYEAR RELCHIREP RELCHINUM twin_firstbirth same_sex_first2 ///
    child1_sex child2_sex two_plus three_plus
rename RELCHI1BYEAR first_child_year
save "`outdir'/parent_instruments_v1.dta", replace

erase "`outdir'/child_sex_map_tmp.dta"
erase "`outdir'/child1_sex_map_tmp.dta"
erase "`outdir'/child2_sex_map_tmp.dta"

di as text "Step 2/5: Build outcomes with nearest pre-birth baseline and +0..+5 window"
use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR HOMEOWN MOVEDFREF_ WHYMOVED1_ ///
    ACTUALROOMS_ ${weight} using "`dta'", clear

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18
drop if missing(RELCHI1BYEAR)
keep if inrange(year, 1984, 2019)

rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
replace movedthisyear = 0 if movedthisyear == 5

rename HOMEOWN own
rename ACTUALROOMS_ rooms
gen first_child_year = RELCHI1BYEAR

bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

xtset ID year
gen K = year - first_child_year

* Baseline = nearest observed year before first birth (not necessarily t=-2).
bysort ID: egen baseline_year = max(cond(K < 0, year, .))
gen has_baseline = !missing(baseline_year)

preserve
    keep if has_baseline == 1 & year == baseline_year
    keep ID first_child_year baseline_year AGEREP EDUYEAR SEX own rooms ${weight}
    rename AGEREP age_pre
    rename EDUYEAR eduyear_pre
    rename SEX sex_pre
    rename own own_pre
    rename rooms rooms_pre
    rename ${weight} iw_pre
    save "`outdir'/tmp_baseline_nearestpre_v1.dta", replace
restore

* Outcomes in years 0..5 after first birth.
gen own_post5_i = (own == 1) if inrange(K, 0, 5)
gen move_post5_i = (movedthisyear == 1) if inrange(K, 0, 5)
gen moved_to_own_post5_i = (own == 1 & L.own == 0 & movedthisyear == 1) if inrange(K, 0, 5)
gen rooms_post5_i = rooms if inrange(K, 0, 5)

* "Moved for size" = moved and reason code 3.
gen moved_for_size_post5_i = (movedthisyear == 1 & WHYMOVED1_ == 3) if inrange(K, 0, 5)
replace moved_for_size_post5_i = 0 if inrange(K, 0, 5) & movedthisyear == 0
replace moved_for_size_post5_i = 0 if inrange(K, 0, 5) & movedthisyear == 1 & WHYMOVED1_ != 3 & !missing(WHYMOVED1_)

bysort ID: egen own_post5 = max(own_post5_i)
bysort ID: egen move_post5 = max(move_post5_i)
bysort ID: egen moved_to_own_post5 = max(moved_to_own_post5_i)
bysort ID: egen rooms_post5 = mean(rooms_post5_i)
bysort ID: egen moved_for_size_post5 = max(moved_for_size_post5_i)

bysort ID (year): keep if _n == 1
keep ID first_child_year own_post5 move_post5 moved_to_own_post5 rooms_post5 moved_for_size_post5

merge 1:1 ID first_child_year using "`outdir'/tmp_baseline_nearestpre_v1.dta", keep(match) nogen
merge 1:1 ID first_child_year using "`outdir'/parent_instruments_v1.dta", keep(match) nogen

erase "`outdir'/tmp_baseline_nearestpre_v1.dta"

gen female_pre = (sex_pre == 2)
gen rooms_change_post5 = rooms_post5 - rooms_pre if !missing(rooms_post5, rooms_pre)

drop if missing(age_pre, eduyear_pre, female_pre, own_pre, iw_pre)

save "`outdir'/analysis_sample_power_v1.dta", replace

di as text "Step 3/5: Estimate twins and same-sex (all-tenure) for ownership outcomes"
estimates clear

* Twins design: two_plus instrumented by twin_firstbirth.
preserve
    keep if !missing(two_plus, twin_firstbirth, own_post5, moved_to_own_post5, moved_for_size_post5)

    reg two_plus twin_firstbirth c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store fs_twins
    test twin_firstbirth
    scalar F_twins = r(F)

    reg own_post5 two_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store ols_twins_own
    ivregress 2sls own_post5 (two_plus = twin_firstbirth) c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store iv_twins_own

    reg moved_to_own_post5 two_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store ols_twins_moveown
    ivregress 2sls moved_to_own_post5 (two_plus = twin_firstbirth) c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store iv_twins_moveown

    reg moved_for_size_post5 two_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store ols_twins_msize
    ivregress 2sls moved_for_size_post5 (two_plus = twin_firstbirth) c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store iv_twins_msize
restore

* Same-sex design: three_plus instrumented by same_sex_first2.
preserve
    keep if !missing(three_plus, same_sex_first2, own_post5, moved_to_own_post5, moved_for_size_post5)

    reg three_plus same_sex_first2 c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store fs_samesex
    test same_sex_first2
    scalar F_samesex = r(F)

    reg own_post5 three_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store ols_samesex_own
    ivregress 2sls own_post5 (three_plus = same_sex_first2) c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store iv_samesex_own

    reg moved_to_own_post5 three_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store ols_samesex_moveown
    ivregress 2sls moved_to_own_post5 (three_plus = same_sex_first2) c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store iv_samesex_moveown

    reg moved_for_size_post5 three_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store ols_samesex_msize
    ivregress 2sls moved_for_size_post5 (three_plus = same_sex_first2) c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store iv_samesex_msize
restore

capture which esttab
if _rc == 0 {
    esttab fs_twins fs_samesex using "`outdir'/iv_first_stage_power_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("First stage: twins->2+" "First stage: same-sex->3+") ///
        stats(N r2, fmt(0 3) labels("Obs." "R-squared"))

    esttab ols_twins_own iv_twins_own ols_samesex_own iv_samesex_own ///
        using "`outdir'/iv_own_post5_ols_iv_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("OLS twins" "IV twins" "OLS same-sex" "IV same-sex") ///
        stats(N, fmt(0) labels("Obs."))

    esttab ols_twins_moveown iv_twins_moveown ols_samesex_moveown iv_samesex_moveown ///
        using "`outdir'/iv_moved_to_own_post5_ols_iv_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("OLS twins" "IV twins" "OLS same-sex" "IV same-sex") ///
        stats(N, fmt(0) labels("Obs."))

    esttab ols_twins_msize iv_twins_msize ols_samesex_msize iv_samesex_msize ///
        using "`outdir'/iv_moved_for_size_post5_ols_iv_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("OLS twins" "IV twins" "OLS same-sex" "IV same-sex") ///
        stats(N, fmt(0) labels("Obs."))
}

di as text "Step 4/5: Export diagnostics and sample moments"
preserve
    keep own_post5 moved_to_own_post5 moved_for_size_post5 two_plus three_plus twin_firstbirth same_sex_first2 own_pre
    collapse (mean) own_post5 moved_to_own_post5 moved_for_size_post5 two_plus three_plus twin_firstbirth same_sex_first2 own_pre
    export delimited using "`outdir'/iv_power_design_summary_v1.csv", replace
restore

tempname fh
file open `fh' using "`outdir'/iv_power_fstats_v1.csv", write replace
file write `fh' "stat,value" _n
local f_tw : display %10.4f scalar(F_twins)
local f_ss : display %10.4f scalar(F_samesex)
file write `fh' "F_twins_first_stage,`f_tw'" _n
file write `fh' "F_samesex_first_stage,`f_ss'" _n
file close `fh'

di as text "Step 5/5: Save a compact N summary"
preserve
    keep ID own_post5 moved_to_own_post5 moved_for_size_post5 two_plus three_plus twin_firstbirth same_sex_first2
    gen sample_twins = !missing(own_post5, moved_to_own_post5, moved_for_size_post5, two_plus, twin_firstbirth)
    gen sample_samesex = !missing(own_post5, moved_to_own_post5, moved_for_size_post5, three_plus, same_sex_first2)
    collapse (sum) sample_twins sample_samesex
    export delimited using "`outdir'/iv_power_ns_v1.csv", replace
restore

di as result "IV power-recovery script complete."
di as result "Sample file: `outdir'/analysis_sample_power_v1.dta"
di as result "Outcome tables: own/moved-to-own/moved-for-size post5"

log close _all
