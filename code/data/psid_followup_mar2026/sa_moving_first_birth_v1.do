clear all
set more off
version 17.0

* Full tail-binned Sun-Abraham moving event-study attempt.
* Status, 2026-04-24: superseded for seminar purposes by
* sa_moving_first_birth_window_v1.do, which runs the finite K=-5..5
* Sun-Abraham object with standard errors for move_any, moved_for_size,
* and moved_to_own. The older full-tail version is left here for reference.

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/sa_moving_first_birth_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/sa_moving_first_birth_v1.log", replace text

global weight IW
global ref_sa 2
global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

cap which eventstudyinteract
if _rc {
    di as error "eventstudyinteract is not installed. Aborting."
    exit 198
}

di as text "Loading PSID panel for first-birth moving event studies..."
use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR HOMEOWN ///
    MOVEDFREF_ WHYMOVED1_ ${weight} using "`dta'", clear

drop if year > DEATHYEAR
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3
rename HOMEOWN own

rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
replace movedthisyear = 0 if movedthisyear == 5

bysort ID: egen first_child_year = min(RELCHI1BYEAR)
drop if missing(first_child_year)

* Exclude first births before panel entry, matching the existing SA scripts.
bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

xtset ID year

gen move_any = movedthisyear

gen moved_for_size = .
replace moved_for_size = 1 if movedthisyear == 1 & WHYMOVED1_ == 3 & year > 1974
replace moved_for_size = 0 if movedthisyear == 0 & year > 1974
replace moved_for_size = 0 if movedthisyear == 1 & WHYMOVED1_ != 3 & !missing(WHYMOVED1_) & year > 1974

gen change_to_own = (own == 1 & L.own == 0) if !missing(own, L.own)
gen moved_to_own = (movedthisyear == 1 & change_to_own == 1) if !missing(movedthisyear, change_to_own)

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

tempfile all_estimates summary status
tempname ph sh sth

postfile `ph' str32 outcome str32 coeff double relative_time b se ci_lo ci_hi ///
    double pre_event_mean post03_raw_mean sample_obs sample_ids using "`all_estimates'", replace

postfile `sh' str32 outcome double pre_event_mean post03_raw_mean ///
    double coef_p0 se_p0 coef_p1 se_p1 coef_p2 se_p2 coef_p3 se_p3 coef_p5 se_p5 ///
    double sample_obs sample_ids using "`summary'", replace

postfile `sth' str32 outcome int rc str120 message using "`status'", replace

local outcomes "move_any moved_for_size moved_to_own"

foreach y of local outcomes {
    local ytitle "`y'"
    if "`y'" == "move_any" local ytitle "Any move"
    if "`y'" == "moved_for_size" local ytitle "Moved for size"
    if "`y'" == "moved_to_own" local ytitle "Moved to ownership"

    di as text "Running Sun-Abraham event study for outcome: `y'"

    capture noisily eventstudyinteract `y' L*event F*event, ///
        vce(cluster ID) absorb(year) cohort(first_child_year) control_cohort(lastcohort) ///
        covariates(i.AGEREP i.EDUYEAR)

    if _rc {
        local rc = _rc
        post `sth' ("`y'") (`rc') ("eventstudyinteract failed")
        di as error "Event study failed for `y' with rc=`rc'. Continuing."
        continue
    }

    quietly summarize `y' if e(sample) & K == -${ref_sa}
    local pre_event_mean = r(mean)
    quietly summarize `y' if e(sample) & inrange(K, 0, 3)
    local post03_raw_mean = r(mean)
    quietly count if e(sample)
    local sample_obs = r(N)
    capture drop tag_id
    egen tag_id = tag(ID) if e(sample)
    quietly count if tag_id == 1
    local sample_ids = r(N)
    drop tag_id

    matrix bmat = e(b_iw)
    matrix vmat = e(V_iw)
    local cn : colnames bmat

    foreach rt in 0 1 2 3 5 {
        local coef_p`rt' = .
        local se_p`rt' = .
    }

    local j = 0
    foreach c of local cn {
        local ++j
        local bb = bmat[1, `j']
        local vv = vmat[`j', `j']
        local ss = sqrt(`vv')

        local cclean "`c'"
        local cclean : subinstr local cclean "event" "", all
        local cclean : subinstr local cclean "L" "", all
        local cclean : subinstr local cclean "F" "-", all
        local cclean : subinstr local cclean "o." "", all
        local cclean : subinstr local cclean "b." "", all
        capture local reltime = real("`cclean'")
        if missing(`reltime') continue

        if `bb' == 0 & `ss' == 0 continue

        local lo = `bb' - 1.96 * `ss'
        local hi = `bb' + 1.96 * `ss'
        post `ph' ("`y'") ("`c'") (`reltime') (`bb') (`ss') (`lo') (`hi') ///
            (`pre_event_mean') (`post03_raw_mean') (`sample_obs') (`sample_ids')

        if inlist(`reltime', 0, 1, 2, 3, 5) {
            local coef_p`reltime' = `bb'
            local se_p`reltime' = `ss'
        }
    }

    post `ph' ("`y'") ("F${ref_sa}event") (-${ref_sa}) (0) (0) (0) (0) ///
        (`pre_event_mean') (`post03_raw_mean') (`sample_obs') (`sample_ids')

    post `sh' ("`y'") (`pre_event_mean') (`post03_raw_mean') ///
        (`coef_p0') (`se_p0') (`coef_p1') (`se_p1') (`coef_p2') (`se_p2') ///
        (`coef_p3') (`se_p3') (`coef_p5') (`se_p5') (`sample_obs') (`sample_ids')

    post `sth' ("`y'") (0) ("completed")
}

postclose `ph'
postclose `sh'
postclose `sth'

use "`all_estimates'", clear
sort outcome relative_time
save "`outdir'/moving_first_birth_eventstudy_all_estimates_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_eventstudy_all_estimates_v1.csv", replace

foreach y of local outcomes {
    local ytitle "`y'"
    if "`y'" == "move_any" local ytitle "Any move"
    if "`y'" == "moved_for_size" local ytitle "Moved for size"
    if "`y'" == "moved_to_own" local ytitle "Moved to ownership"

    preserve
        keep if outcome == "`y'"
        count
        if r(N) > 0 {
            sort relative_time
            export delimited using "`outdir'/`y'_first_birth_eventstudy_estimates_v1.csv", replace
            save "`outdir'/`y'_first_birth_eventstudy_estimates_v1.dta", replace

            keep if !missing(b)
            twoway ///
                (rarea ci_lo ci_hi relative_time, fcolor(navy%10) lcolor(navy%10)) ///
                (connected b relative_time, mcolor(navy) msymbol(circle) lcolor(navy) lwidth(medthick) lpattern(dash)), ///
                ytitle("`ytitle'") xtitle("Years relative to first birth") ///
                xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
                yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
                title("`ytitle' around first birth") ///
                $graphregion_slides legend(off)
            graph export "`outdir'/`y'_first_birth_eventstudy_v1.png", replace
        }
    restore
}

use "`summary'", clear
save "`outdir'/moving_first_birth_eventstudy_summary_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_eventstudy_summary_v1.csv", replace

use "`status'", clear
save "`outdir'/moving_first_birth_eventstudy_status_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_eventstudy_status_v1.csv", replace

di as result "First-birth moving event studies complete."
di as result "Output dir: `outdir'"

log close _all
