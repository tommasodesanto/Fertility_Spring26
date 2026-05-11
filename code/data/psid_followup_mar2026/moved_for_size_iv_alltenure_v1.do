clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/moved_for_size_iv_alltenure_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/moved_for_size_iv_alltenure_v1.log", replace text

global weight IW

di as text "Step 1/4: Build fertility instruments (twins + same-sex)"
use ID year SEX RELCHIREP RELCHI1ID RELCHI2ID RELCHI1BYEAR RELCHI2BYEAR ///
    RELCHI3BYEAR RELCHI4BYEAR RELCHI5BYEAR RELCHI6BYEAR RELCHI7BYEAR RELCHI8BYEAR ///
    RELCHI9BYEAR RELCHI10BYEAR RELCHI11BYEAR RELCHI12BYEAR RELCHI13BYEAR RELCHI14BYEAR ///
    RELCHI15BYEAR RELCHI16BYEAR RELCHI17BYEAR RELCHI18BYEAR RELCHI19BYEAR RELCHI20BYEAR ///
    using "`dta'", clear

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

preserve
    use ID year RELCHIREP RELCHI1ID RELCHI2ID RELCHI1BYEAR RELCHI2BYEAR ///
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

    keep ID RELCHI1BYEAR twin_firstbirth same_sex_first2 child1_sex child2_sex two_plus three_plus
    rename RELCHI1BYEAR first_child_year
    save "`outdir'/parent_instruments_v1.dta", replace
restore

erase "`outdir'/child_sex_map_tmp.dta"
erase "`outdir'/child1_sex_map_tmp.dta"
erase "`outdir'/child2_sex_map_tmp.dta"

di as text "Step 2/4: Build moved-for-size outcome in t=0..3 around first birth"
use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR HOMEOWN WHYMOVED1_ ${weight} ///
    using "`dta'", clear

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18
drop if missing(RELCHI1BYEAR)
keep if inrange(year, 1984, 2019)

gen first_child_year = RELCHI1BYEAR
bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

gen moved_for_size = 1 if WHYMOVED1_ == 3 & WHYMOVED1_ != . & year > 1974
replace moved_for_size = 0 if WHYMOVED1_ != 3 & WHYMOVED1_ != . & year > 1974

xtset ID year
gen K = year - first_child_year

* Baseline controls at t=-2.
preserve
    keep if K == -2
    keep ID first_child_year AGEREP EDUYEAR SEX HOMEOWN ${weight}
    rename AGEREP age_pre
    rename EDUYEAR eduyear_pre
    rename SEX sex_pre
    rename HOMEOWN own_pre
    rename ${weight} iw_pre
    save "`outdir'/tmp_baseline_tminus2.dta", replace
restore

* Main outcome (all tenures): ever moved for size in years 0..3.
gen moved_for_size_post3_i = (moved_for_size == 1) if inrange(K, 0, 3)
bysort ID: egen moved_for_size_post3 = max(moved_for_size_post3_i)

bysort ID (year): keep if _n == 1
keep ID first_child_year moved_for_size_post3

merge 1:1 ID first_child_year using "`outdir'/tmp_baseline_tminus2.dta", keep(match) nogen
merge 1:1 ID first_child_year using "`outdir'/parent_instruments_v1.dta", keep(match) nogen
erase "`outdir'/tmp_baseline_tminus2.dta"

gen female_pre = (sex_pre == 2)
drop if missing(age_pre, eduyear_pre, female_pre, own_pre, iw_pre)

save "`outdir'/analysis_sample_moved_for_size_v1.dta", replace

di as text "Step 3/4: Estimate twins and same-sex designs"
estimates clear

* A) Twins IV: treatment is two_plus.
preserve
    keep if !missing(moved_for_size_post3, two_plus, twin_firstbirth)

    reg two_plus twin_firstbirth c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store fs_twins
    test twin_firstbirth
    scalar F_twins = r(F)

    reg moved_for_size_post3 two_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store ols_twins

    ivregress 2sls moved_for_size_post3 (two_plus = twin_firstbirth) c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store iv_twins
restore

* B) Same-sex IV: treatment is three_plus.
preserve
    keep if !missing(moved_for_size_post3, three_plus, same_sex_first2)

    reg three_plus same_sex_first2 c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store fs_samesex
    test same_sex_first2
    scalar F_samesex = r(F)

    reg moved_for_size_post3 three_plus c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store ols_samesex

    ivregress 2sls moved_for_size_post3 (three_plus = same_sex_first2) c.age_pre i.eduyear_pre i.first_child_year i.female_pre i.own_pre [aw=iw_pre], vce(robust)
    estimates store iv_samesex
restore

capture which esttab
if _rc == 0 {
    esttab fs_twins fs_samesex using "`outdir'/moved_for_size_first_stage_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("First stage: twins->2+" "First stage: same-sex->3+") ///
        stats(N r2, fmt(0 3) labels("Obs." "R-squared"))

    esttab ols_twins iv_twins ols_samesex iv_samesex ///
        using "`outdir'/moved_for_size_second_stage_ols_iv_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("OLS twins" "IV twins" "OLS same-sex" "IV same-sex") ///
        stats(N, fmt(0) labels("Obs."))
}

di as text "Step 4/4: Export compact diagnostics"
preserve
    keep moved_for_size_post3 twin_firstbirth same_sex_first2 two_plus three_plus own_pre
    collapse (mean) moved_for_size_post3 twin_firstbirth same_sex_first2 two_plus three_plus own_pre
    export delimited using "`outdir'/moved_for_size_design_summary_v1.csv", replace
restore

tempname fh
file open `fh' using "`outdir'/moved_for_size_fstats_v1.csv", write replace
file write `fh' "stat,value" _n
local f_tw : display %10.4f scalar(F_twins)
local f_ss : display %10.4f scalar(F_samesex)
file write `fh' "F_twins_first_stage,`f_tw'" _n
file write `fh' "F_samesex_first_stage,`f_ss'" _n
file close `fh'

di as result "Moved-for-size IV analysis complete."
di as result "Output dir: `outdir'"

log close _all
