clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/sa_moved_for_size_first_birth_window_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/sa_moved_for_size_first_birth_window_v1.log", replace text

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

di as text "Loading PSID panel for Sun-Abraham moved-for-size event study..."
timer clear
timer on 1

use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR MOVEDFREF_ WHYMOVED1_ ${weight} using "`dta'", clear

drop if year > DEATHYEAR
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
replace movedthisyear = 0 if movedthisyear == 5

gen moved_for_size = .
replace moved_for_size = 1 if movedthisyear == 1 & WHYMOVED1_ == 3 & year > 1974
replace moved_for_size = 0 if movedthisyear == 0 & year > 1974
replace moved_for_size = 0 if movedthisyear == 1 & WHYMOVED1_ != 3 & !missing(WHYMOVED1_) & year > 1974

gen first_child_year = RELCHI1BYEAR
drop if missing(first_child_year)

* Exclude first births observed before panel entry.
bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

gen K = year - first_child_year

* Finite event window. The broader tail-bin script was too slow locally; this is
* the clean window needed for the seminar figure.
keep if inrange(K, -5, 5)
drop if missing(moved_for_size, first_child_year, year, AGEREP, EDUYEAR, ID)

sum first_child_year
gen lastcohort = first_child_year == r(max)

cap drop L*event F*event
forvalues l = 0/5 {
    gen L`l'event = K == `l'
}
forvalues l = 1/5 {
    gen F`l'event = K == -`l'
}
drop F${ref_sa}event

quietly count
local input_obs = r(N)
egen tag_id = tag(ID)
quietly count if tag_id == 1
local input_ids = r(N)
drop tag_id
levelsof first_child_year if lastcohort == 0, local(cohorts)
local n_cohorts : word count `cohorts'
di as text "SA input: obs=`input_obs', ids=`input_ids', non-control cohorts=`n_cohorts'"

timer off 1
timer on 2

di as text "Running Sun-Abraham eventstudyinteract for moved_for_size..."
eventstudyinteract moved_for_size L*event F*event, ///
    vce(cluster ID) absorb(year) cohort(first_child_year) control_cohort(lastcohort) ///
    covariates(i.AGEREP i.EDUYEAR)

timer off 2

quietly summarize moved_for_size if e(sample) & K == -${ref_sa}
local pre_event_mean = r(mean)
quietly summarize moved_for_size if e(sample) & inrange(K, 0, 3)
local post03_raw_mean = r(mean)
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
matrix moved_for_size_sa = combine'
matrix drop combine

clear
svmat2 moved_for_size_sa, names(col) rnames(coeff)
gen outcome = "moved_for_size"
gen estimator = "Sun-Abraham IW"
gen pre_event_mean = `pre_event_mean'
gen post03_raw_mean = `post03_raw_mean'
gen sample_obs = `sample_obs'
gen sample_ids = `sample_ids'
gen input_obs = `input_obs'
gen input_ids = `input_ids'
gen n_cohorts = `n_cohorts'
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
replace outcome = "moved_for_size" if missing(outcome)
replace estimator = "Sun-Abraham IW" if missing(estimator)
replace pre_event_mean = `pre_event_mean' if missing(pre_event_mean)
replace post03_raw_mean = `post03_raw_mean' if missing(post03_raw_mean)
replace sample_obs = `sample_obs' if missing(sample_obs)
replace sample_ids = `sample_ids' if missing(sample_ids)
replace input_obs = `input_obs' if missing(input_obs)
replace input_ids = `input_ids' if missing(input_ids)
replace n_cohorts = `n_cohorts' if missing(n_cohorts)
replace relative_time = -${ref_sa} if missing(relative_time)
replace b = 0 if relative_time == -${ref_sa}
replace se = 0 if relative_time == -${ref_sa}
gen ci_lo = b - 1.96 * se
gen ci_hi = b + 1.96 * se
gen b_pct = 100 * b
gen ci_lo_pct = 100 * ci_lo
gen ci_hi_pct = 100 * ci_hi

sort relative_time
order estimator outcome coeff relative_time b se ci_lo ci_hi pre_event_mean post03_raw_mean sample_obs sample_ids input_obs input_ids n_cohorts
save "`outdir'/moved_for_size_sa_first_birth_window_estimates_v1.dta", replace
export delimited using "`outdir'/moved_for_size_sa_first_birth_window_estimates_v1.csv", replace

twoway ///
    (rarea ci_lo_pct ci_hi_pct relative_time, fcolor(maroon%12) lcolor(maroon%12)) ///
    (connected b_pct relative_time, mcolor(maroon) msymbol(circle) msize(medlarge) ///
        lcolor(maroon) lwidth(medthick)), ///
    ytitle("Effect on move-for-size probability, pp") ///
    xtitle("Years relative to first birth") ///
    xlabel(-5(1)5) ///
    xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
    yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
    ylabel(-4(2)6, angle(horizontal)) ///
    title("Move-for-Size Around First Birth") ///
    note("Sun-Abraham interaction-weighted event study; K=-2 omitted. 95% CI clustered by household.", size(vsmall)) ///
    $graphregion_slides legend(off)

graph export "`outdir'/moved_for_size_sa_first_birth_window_v1.png", replace width(1600)
graph export "`outdir'/moved_for_size_sa_first_birth_window_v1.pdf", replace

timer list

di as result "Sun-Abraham moved-for-size event study complete."
di as result "Output dir: `outdir'"

log close _all
