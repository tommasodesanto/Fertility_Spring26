clear all
set more off
version 17.0

local base "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026"
local indata "`base'/output/twins_gender_iv_v1/analysis_sample_iv_v1.dta"
local outroot "`base'/output"
local outdir "`outroot'/secondborn_gender_v1"

cap mkdir "`outroot'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/secondborn_gender_v1.log", replace text

di as text "Loading IV analysis sample from twins/gender pipeline"
use "`indata'", clear

* Keep complete cases for the second-born-sex design.
keep if !missing(child1_sex, child2_sex, three_plus, own_post3, moved_to_own_post3, ///
    age_pre, eduyear_pre, first_child_year, female_pre, iw_pre)

gen mixed_first2 = (child1_sex != child2_sex)
gen same_first2 = (child1_sex == child2_sex)
gen second_girl = (child2_sex == 2)

estimates clear

di as text "First-stage diagnostics: mixed-sex and second-child-girl"
reg three_plus mixed_first2 i.child1_sex c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
estimates store fs_mixed
test mixed_first2
scalar F_mixed = r(F)

reg three_plus second_girl c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre] if child1_sex == 1, vce(robust)
estimates store fs_secondgirl_firstboy

reg three_plus second_girl c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre] if child1_sex == 2, vce(robust)
estimates store fs_secondgirl_firstgirl

di as text "Ownership outcomes: OLS vs IV (instrument = mixed first-two sex)"
reg own_post3 three_plus i.child1_sex c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
estimates store ols_own

ivregress 2sls own_post3 (three_plus = mixed_first2) i.child1_sex c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
estimates store iv_own

reg moved_to_own_post3 three_plus i.child1_sex c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
estimates store ols_moveown

ivregress 2sls moved_to_own_post3 (three_plus = mixed_first2) i.child1_sex c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
estimates store iv_moveown

capture which esttab
if _rc == 0 {
    esttab fs_mixed fs_secondgirl_firstboy fs_secondgirl_firstgirl ///
        using "`outdir'/secondborn_first_stage_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("three+ on mixed-sex first2" "three+ on second-girl | first boy" "three+ on second-girl | first girl") ///
        stats(N r2, fmt(0 3) labels("Obs." "R-squared"))

    esttab ols_own iv_own using "`outdir'/secondborn_own_post3_ols_iv_v1.tex", replace ///
        b(3) se(3) label mtitles("OLS" "IV (mixed-sex instrument)") ///
        stats(N, fmt(0) labels("Obs."))

    esttab ols_moveown iv_moveown using "`outdir'/secondborn_moved_to_own_post3_ols_iv_v1.tex", replace ///
        b(3) se(3) label mtitles("OLS" "IV (mixed-sex instrument)") ///
        stats(N, fmt(0) labels("Obs."))
}

preserve
    keep mixed_first2 same_first2 second_girl three_plus own_post3 moved_to_own_post3
    collapse (mean) mixed_first2 same_first2 second_girl three_plus own_post3 moved_to_own_post3
    export delimited using "`outdir'/secondborn_design_summary_v1.csv", replace
restore

tempname fh
file open `fh' using "`outdir'/secondborn_design_teststats_v1.csv", write replace
file write `fh' "stat,value" _n
local f_mixed : display %10.4f scalar(F_mixed)
file write `fh' "F_mixed_first_stage,`f_mixed'" _n
file close `fh'

di as result "Second-born gender design complete."
di as result "Output dir: `outdir'"
log close _all
