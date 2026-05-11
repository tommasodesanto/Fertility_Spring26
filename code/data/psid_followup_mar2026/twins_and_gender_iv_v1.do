clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/twins_gender_iv_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/twins_and_gender_iv_v1.log", replace text

global weight IW

di as text "Step 1/5: Build person-level child sex map"
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

di as text "Step 2/5: Build parent-level fertility instrument dataset"
use ID year RELCHINUM RELCHIREP RELCHI1ID RELCHI2ID RELCHI1BYEAR RELCHI2BYEAR ///
    RELCHI3BYEAR RELCHI4BYEAR RELCHI5BYEAR RELCHI6BYEAR RELCHI7BYEAR RELCHI8BYEAR ///
    RELCHI9BYEAR RELCHI10BYEAR RELCHI11BYEAR RELCHI12BYEAR RELCHI13BYEAR RELCHI14BYEAR ///
    RELCHI15BYEAR RELCHI16BYEAR RELCHI17BYEAR RELCHI18BYEAR RELCHI19BYEAR RELCHI20BYEAR ///
    using "`dta'", clear

* Keep the last observed child roster for each parent.
bysort ID: egen last_nonmiss_year = max(cond(!missing(RELCHI1BYEAR), year, .))
keep if year == last_nonmiss_year & !missing(last_nonmiss_year)
drop last_nonmiss_year

* Twins-at-first-birth proxy: >=2 children with same birth year as first child.
forvalues k = 1/20 {
    gen sameyr`k' = (RELCHI`k'BYEAR == RELCHI1BYEAR) if !missing(RELCHI`k'BYEAR, RELCHI1BYEAR)
    replace sameyr`k' = 0 if missing(sameyr`k')
}
egen n_sameyear_firstbirth = rowtotal(sameyr1-sameyr20)
gen twin_firstbirth = (n_sameyear_firstbirth >= 2) if !missing(n_sameyear_firstbirth)

drop sameyr1-sameyr20

* Same-sex first two children instrument.
rename RELCHI1ID child1_id
rename RELCHI2ID child2_id

merge m:1 child1_id using "`outdir'/child1_sex_map_tmp.dta", keep(master match) nogen
rename child_sex child1_sex

merge m:1 child2_id using "`outdir'/child2_sex_map_tmp.dta", keep(master match) nogen
rename child_sex child2_sex

gen same_sex_first2 = (child1_sex == child2_sex) if !missing(child1_sex, child2_sex)

* Treatment variables for the two IV designs.
gen two_plus = (RELCHIREP >= 2) if !missing(RELCHIREP)
gen three_plus = (RELCHIREP >= 3) if !missing(RELCHIREP)

keep ID RELCHIREP RELCHINUM RELCHI1BYEAR RELCHI2BYEAR twin_firstbirth same_sex_first2 ///
     child1_sex child2_sex two_plus three_plus
save "`outdir'/parent_instruments_v1.dta", replace

di as text "Step 3/5: Build event outcomes + baseline controls (t=-2, +0..+3)"
use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR HOMEOWN MOVEDFREF_ ACTUALROOMS_ ///
    NETWORTH2R ${weight} using "`dta'", clear

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18
drop if missing(RELCHI1BYEAR)
keep if inrange(year, 1984, 2019)

rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
replace movedthisyear = 0 if movedthisyear == 5

rename ACTUALROOMS_ rooms
rename HOMEOWN own
gen first_child_year = RELCHI1BYEAR

bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

xtset ID year
gen K = year - first_child_year

* Baseline (t=-2).
preserve
    keep if K == -2
    keep ID first_child_year AGEREP EDUYEAR SEX own rooms NETWORTH2R ${weight}
    rename AGEREP age_pre
    rename EDUYEAR eduyear_pre
    rename SEX sex_pre
    rename own own_pre
    rename rooms rooms_pre
    rename NETWORTH2R networth2r_pre
    rename ${weight} iw_pre
    save "`outdir'/tmp_baseline_tminus2.dta", replace
restore

* Outcomes in first three post-birth years.
gen own_post3_i = (own == 1) & inrange(K, 0, 3)
gen move_post3_i = (movedthisyear == 1) & inrange(K, 0, 3)
gen moved_to_own_post3_i = (own == 1 & L.own == 0 & movedthisyear == 1) & inrange(K, 0, 3)
gen rooms_post3_i = rooms if inrange(K, 0, 3)

bysort ID: egen own_post3 = max(own_post3_i)
bysort ID: egen move_post3 = max(move_post3_i)
bysort ID: egen moved_to_own_post3 = max(moved_to_own_post3_i)
bysort ID: egen rooms_post3 = mean(rooms_post3_i)

bysort ID (year): keep if _n == 1
keep ID first_child_year own_post3 move_post3 moved_to_own_post3 rooms_post3

merge 1:1 ID first_child_year using "`outdir'/tmp_baseline_tminus2.dta", keep(match) nogen
merge 1:1 ID using "`outdir'/parent_instruments_v1.dta", keep(match) nogen

erase "`outdir'/tmp_baseline_tminus2.dta"
erase "`outdir'/child_sex_map_tmp.dta"
erase "`outdir'/child1_sex_map_tmp.dta"
erase "`outdir'/child2_sex_map_tmp.dta"

gen rooms_change_post3 = rooms_post3 - rooms_pre if !missing(rooms_post3, rooms_pre)
gen female_pre = (sex_pre == 2)
gen asinh_nw2_pre = asinh(networth2r_pre)

* Main sample for ownership-transition outcomes.
keep if own_pre == 0
drop if missing(age_pre, eduyear_pre, female_pre, iw_pre)

save "`outdir'/analysis_sample_iv_v1.dta", replace

di as text "Step 4/5: Estimate twins-IV and same-sex-IV designs"
estimates clear

* A) Twins at first birth -> 2+ children
preserve
    keep if !missing(two_plus, twin_firstbirth, own_post3, moved_to_own_post3)

    reg two_plus twin_firstbirth c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store fs_twins

    reg own_post3 two_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store ols_twins_own

    ivregress 2sls own_post3 (two_plus = twin_firstbirth) c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store iv_twins_own

    reg moved_to_own_post3 two_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store ols_twins_moveown

    ivregress 2sls moved_to_own_post3 (two_plus = twin_firstbirth) c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store iv_twins_moveown
restore

* B) Same-sex first two -> 3+ children
preserve
    keep if !missing(three_plus, same_sex_first2, own_post3, moved_to_own_post3)

    reg three_plus same_sex_first2 c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store fs_samesex

    reg own_post3 three_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store ols_samesex_own

    ivregress 2sls own_post3 (three_plus = same_sex_first2) c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store iv_samesex_own

    reg moved_to_own_post3 three_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store ols_samesex_moveown

    ivregress 2sls moved_to_own_post3 (three_plus = same_sex_first2) c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
    estimates store iv_samesex_moveown
restore

capture which esttab
if _rc == 0 {
    esttab fs_twins fs_samesex using "`outdir'/iv_first_stage_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("First stage: twins->2+" "First stage: same-sex->3+") ///
        stats(N r2, fmt(0 3) labels("Obs." "R-squared"))

    esttab ols_twins_own iv_twins_own ols_samesex_own iv_samesex_own ///
        using "`outdir'/iv_vs_ols_own_post3_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("OLS twins design" "IV twins design" "OLS same-sex design" "IV same-sex design") ///
        stats(N, fmt(0) labels("Obs."))

    esttab ols_twins_moveown iv_twins_moveown ols_samesex_moveown iv_samesex_moveown ///
        using "`outdir'/iv_vs_ols_moved_to_own_post3_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("OLS twins design" "IV twins design" "OLS same-sex design" "IV same-sex design") ///
        stats(N, fmt(0) labels("Obs."))
}

di as text "Step 5/5: Export compact sample diagnostics"
preserve
    keep twin_firstbirth same_sex_first2 two_plus three_plus own_post3 moved_to_own_post3
    collapse (mean) twin_firstbirth same_sex_first2 two_plus three_plus own_post3 moved_to_own_post3
    export delimited using "`outdir'/iv_design_summary_v1.csv", replace
restore

di as result "Twins + same-sex IV script complete."
di as result "Main sample: `outdir'/analysis_sample_iv_v1.dta"
di as result "First stage table: `outdir'/iv_first_stage_v1.tex"
di as result "Own outcome table: `outdir'/iv_vs_ols_own_post3_v1.tex"
di as result "Moved-to-own table: `outdir'/iv_vs_ols_moved_to_own_post3_v1.tex"
di as result "Design summary: `outdir'/iv_design_summary_v1.csv"

log close _all
