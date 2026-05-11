clear all
set more off
version 17.0

args variant
if "`variant'" == "" local variant "all"

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/sa_rooms_second_birth_one_variant_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/sa_rooms_second_birth_`variant'_v1.log", replace text

global weight IW
global ref_sa 2
global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

cap which eventstudyinteract
if _rc {
    di as error "eventstudyinteract is not installed. Aborting."
    exit 198
}

cap which svmat2
if _rc {
    di as error "svmat2 is not installed. Aborting."
    exit 198
}

use "`dta'", clear
keep ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR RELCHI2BYEAR RELCHI3BYEAR ACTUALROOMS_ ${weight}

drop if year > DEATHYEAR
drop if AGEREP < 18

bysort ID: egen first_child_year = min(RELCHI1BYEAR)
bysort ID: egen second_child_year = min(RELCHI2BYEAR)
bysort ID: egen third_child_year = min(RELCHI3BYEAR)
drop if missing(second_child_year)

bysort ID: egen year_entry = min(year)
drop if second_child_year < year_entry
drop year_entry

rename ACTUALROOMS_ rooms
xtset ID year

local desc "Baseline second-birth event study"
if inlist("`variant'", "drop_second_twins", "censor_at_third", "no_third_by3", "no_third_by5", "censor_at_third_gap3", "no_third_by3_gap3", "censor_at_third_gap5", "no_third_by3_gap5", "no_third_by5_gap5") {
    drop if !missing(first_child_year) & second_child_year == first_child_year
    drop if !missing(third_child_year) & third_child_year == second_child_year
    local desc "Exclude multiple births at second birth"
}
if "`variant'" == "censor_at_third" {
    drop if !missing(third_child_year) & year >= third_child_year
    local desc "Exclude multiples and censor observations at third birth"
}
if "`variant'" == "no_third_by3" {
    drop if !missing(third_child_year) & third_child_year <= second_child_year + 3
    local desc "Exclude multiples and require no third birth by +3"
}
if "`variant'" == "no_third_by5" {
    drop if !missing(third_child_year) & third_child_year <= second_child_year + 5
    local desc "Exclude multiples and require no third birth by +5"
}
if "`variant'" == "censor_at_third_gap3" {
    drop if missing(first_child_year)
    drop if second_child_year - first_child_year < 3
    drop if !missing(third_child_year) & year >= third_child_year
    local desc "Censor at third birth and require first birth at least 3 years earlier"
}
if "`variant'" == "no_third_by3_gap3" {
    drop if missing(first_child_year)
    drop if second_child_year - first_child_year < 3
    drop if !missing(third_child_year) & third_child_year <= second_child_year + 3
    local desc "Require no third birth by +3 and first birth at least 3 years earlier"
}
if "`variant'" == "censor_at_third_gap5" {
    drop if missing(first_child_year)
    drop if second_child_year - first_child_year < 5
    drop if !missing(third_child_year) & year >= third_child_year
    local desc "Censor at third birth and require first birth at least 5 years earlier"
}
if "`variant'" == "no_third_by3_gap5" {
    drop if missing(first_child_year)
    drop if second_child_year - first_child_year < 5
    drop if !missing(third_child_year) & third_child_year <= second_child_year + 3
    local desc "Require no third birth by +3 and first birth at least 5 years earlier"
}
if "`variant'" == "no_third_by5_gap5" {
    drop if missing(first_child_year)
    drop if second_child_year - first_child_year < 5
    drop if !missing(third_child_year) & third_child_year <= second_child_year + 5
    local desc "Require no third birth by +5 and first birth at least 5 years earlier"
}

gen K = year - second_child_year
sum second_child_year
gen lastcohort = second_child_year == r(max)

cap drop L*event F*event
forvalues l = 0/10 {
    gen L`l'event = K == `l'
}
gen L11event = (K > 10 & K != .)
forvalues l = 1/5 {
    gen F`l'event = K == -`l'
}
gen F7event = (K < -6 & K != .)
drop F${ref_sa}event

eventstudyinteract rooms L*event F*event, ///
    vce(cluster ID) absorb(year) cohort(second_child_year) control_cohort(lastcohort) ///
    covariates(i.AGEREP i.EDUYEAR)

quietly summarize rooms if e(sample) & K == -${ref_sa}
local pre_event_mean = r(mean)
quietly count if e(sample)
local sample_obs = r(N)
egen tag_id = tag(ID) if e(sample)
quietly count if tag_id == 1
local sample_ids = r(N)
drop tag_id

matrix b = e(b_iw)
matrix var = e(V_iw)
matrix vardiag = vecdiag(var)
matrix combine = b \ vardiag
matrix rownames combine = b variance
matrix rooms_second_birth = combine'
matrix drop combine

clear
svmat2 rooms_second_birth, names(col) rnames(coeff)
gen variant = "`variant'"
gen description = "`desc'"
gen pre_event_mean = `pre_event_mean'
gen sample_obs = `sample_obs'
gen sample_ids = `sample_ids'
gen se = sqrt(variance)
drop variance

replace b = . if b == 0 & se == 0
replace se = . if b == . & se == 0
drop if coeff == "_cons"
gen relative_time = subinstr(coeff, "event", "", .)
replace relative_time = subinstr(relative_time, "L", "", .)
replace relative_time = subinstr(relative_time, "F", "-", .)
replace relative_time = subinstr(relative_time, "o.", "", .)
destring relative_time, replace

assert relative_time != -${ref_sa}
assert relative_time != .
set obs `=_N+1'
replace variant = "`variant'" if missing(variant)
replace description = "`desc'" if missing(description)
replace pre_event_mean = `pre_event_mean' if missing(pre_event_mean)
replace sample_obs = `sample_obs' if missing(sample_obs)
replace sample_ids = `sample_ids' if missing(sample_ids)
replace relative_time = -${ref_sa} if missing(relative_time)
replace b = 0 if relative_time == -${ref_sa}
replace se = 0 if relative_time == -${ref_sa}
gen ci_lo = b - 1.96*se
gen ci_hi = b + 1.96*se
sort relative_time
order variant description coeff relative_time b se ci_lo ci_hi pre_event_mean sample_obs sample_ids
save "`outdir'/rooms_s_c_y_`variant'_estimates.dta", replace

preserve
    keep if !missing(b)
    twoway ///
        (rarea ci_lo ci_hi relative_time, fcolor(dkgreen%10) lcolor(dkgreen%10)) ///
        (connected b relative_time, mcolor(dkgreen) msymbol(circle) lcolor(dkgreen) lwidth(medthick) lpattern(dash)), ///
        ytitle("Rooms") xtitle("Years relative to second birth") ///
        xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
        yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
        title("`desc'") ///
        $graphregion_slides legend(off)
    graph export "`outdir'/rooms_s_c_y_`variant'.png", replace
restore

quietly summarize b if relative_time == 3
local coef_p3 = r(mean)
quietly summarize se if relative_time == 3
local se_p3 = r(mean)
quietly summarize b if relative_time == 5
local coef_p5 = r(mean)
quietly summarize se if relative_time == 5
local se_p5 = r(mean)

clear
set obs 1
gen variant = "`variant'"
gen description = "`desc'"
gen pre_event_mean = `pre_event_mean'
gen coef_p3 = `coef_p3'
gen se_p3 = `se_p3'
gen coef_p5 = `coef_p5'
gen se_p5 = `se_p5'
gen sample_obs = `sample_obs'
gen sample_ids = `sample_ids'
save "`outdir'/rooms_s_c_y_`variant'_summary.dta", replace
export delimited using "`outdir'/rooms_s_c_y_`variant'_summary.csv", replace

di as result "Second-birth rooms variant complete: `variant'"
di as result "Summary: `outdir'/rooms_s_c_y_`variant'_summary.csv"
di as result "Estimates: `outdir'/rooms_s_c_y_`variant'_estimates.dta"

log close _all
