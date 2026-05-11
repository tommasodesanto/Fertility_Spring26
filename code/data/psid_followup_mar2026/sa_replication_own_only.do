clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/sa_replication"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/sa_replication_own_only.log", replace text

global weight IW
global ref_sa 2
global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

di as text "Loading PSID panel..."
use "`dta'", clear

* Keep the same core fields used by the original SA script.
keep ID year AGEREP EDUYEAR SEX RELCHI1BYEAR MOVEDFREF_ DEATHYEAR ACTUALROOMS_ HOMEOWN ${weight}

* Match original tenure coding.
replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18

rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
replace movedthisyear = 0 if movedthisyear == 5

gen first_child_year = RELCHI1BYEAR
drop if missing(first_child_year)

* Exclude households with first birth before panel entry.
bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

rename ACTUALROOMS_ rooms
rename HOMEOWN own

xtset ID year

gen K = year - first_child_year
sum first_child_year
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
drop F${ref_sa}event

cap which eventstudyinteract
if _rc {
    di as error "eventstudyinteract is not installed. Aborting replication."
    exit 198
}

cap which svmat2
if _rc {
    di as error "svmat2 is not installed. Aborting replication."
    exit 198
}

di as text "Running Sun-Abraham event-study for outcome: own"
eventstudyinteract own L*event F*event, ///
    vce(cluster ID) absorb(year) cohort(first_child_year) control_cohort(lastcohort) ///
    covariates(i.AGEREP i.EDUYEAR)

summ own if e(sample) & K == -${ref_sa}
local prebirth_mean = r(mean)

matrix b = e(b_iw)
matrix var = e(V_iw)
matrix vardiag = vecdiag(var)
matrix combine = b \ vardiag
matrix rownames combine = b variance
matrix own_f_c_y_all_repl = combine'
matrix drop combine

preserve
    clear
    svmat2 own_f_c_y_all_repl, names(col) rnames(coeff)
    gen prebirth_mean = `prebirth_mean'
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
    replace relative_time = -${ref_sa} if relative_time == .
    replace b = 0 if relative_time == -${ref_sa}
    replace se = 0 if relative_time == -${ref_sa}

    gen ci_lo = b - 1.96*se
    gen ci_hi = b + 1.96*se

    sort relative_time
    order coeff relative_time b se ci_lo ci_hi prebirth_mean
    save "`outdir'/own_f_c_y_all_repl_estimates.dta", replace
restore

use "`outdir'/own_f_c_y_all_repl_estimates.dta", clear
drop if b == .

twoway ///
    (rarea ci_lo ci_hi relative_time, fcolor(dkgreen%10) lcolor(dkgreen%10)) ///
    (connected b relative_time, mcolor(dkgreen) msymbol(circle) lcolor(dkgreen) lwidth(medthick) lpattern(dash)), ///
    ytitle("Homeownership") xtitle("Years relative to childbirth") ///
    xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
    yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
    $graphregion_slides legend(off)

graph export "`outdir'/own_f_c_y_all_repl.png", replace

capture confirm file "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Tables/newtests/noIDFE/own_f_c_y_all_estimates.dta"
if _rc == 0 {
    preserve
        rename b b_repl
        keep relative_time b_repl
        merge 1:1 relative_time using "/Users/tommasodesanto/Desktop/Projects/Fertility/Outputs/Tables/newtests/noIDFE/own_f_c_y_all_estimates.dta", nogen keep(match)
        rename b b_original
        gen b_diff = b_repl - b_original
        save "`outdir'/own_f_c_y_all_repl_vs_original.dta", replace
    restore
}

di as result "Replication complete."
di as result "Graph: `outdir'/own_f_c_y_all_repl.png"
di as result "Estimates: `outdir'/own_f_c_y_all_repl_estimates.dta"

log close _all
