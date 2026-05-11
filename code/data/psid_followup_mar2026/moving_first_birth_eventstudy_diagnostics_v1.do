clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/moving_first_birth_eventstudy_diagnostics_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/moving_first_birth_eventstudy_diagnostics_v1.log", replace text

global weight IW
global ref_sa 2
global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

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

local eventvars F7event F5event F4event F3event F1event ///
    L0event L1event L2event L3event L4event L5event ///
    L6event L7event L8event L9event L10event L11event

local outcomes "move_any moved_for_size moved_to_own"

tempfile raw pooled
tempname rh ph

postfile `rh' str32 outcome double relative_time raw_mean raw_se raw_n using "`raw'", replace
postfile `ph' str32 outcome str32 coeff double relative_time b se ci_lo ci_hi pvalue ///
    double pre_event_mean post03_raw_mean sample_obs sample_ids using "`pooled'", replace

foreach y of local outcomes {
    forvalues kk = -5/10 {
        quietly summarize `y' if K == `kk'
        post `rh' ("`y'") (`kk') (r(mean)) (r(sd) / sqrt(r(N))) (r(N))
    }

    regress `y' `eventvars' i.year i.AGEREP i.EDUYEAR, vce(cluster ID)

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

    foreach ev of local eventvars {
        local reltime .
        if substr("`ev'", 1, 1) == "L" {
            local rel = subinstr("`ev'", "L", "", .)
            local rel = subinstr("`rel'", "event", "", .)
            local reltime = real("`rel'")
        }
        if substr("`ev'", 1, 1) == "F" {
            local rel = subinstr("`ev'", "F", "", .)
            local rel = subinstr("`rel'", "event", "", .)
            local reltime = -real("`rel'")
        }

        local bb = _b[`ev']
        local ss = _se[`ev']
        local lo = `bb' - 1.96 * `ss'
        local hi = `bb' + 1.96 * `ss'
        local pval = 2 * ttail(e(df_r), abs(`bb' / `ss'))
        post `ph' ("`y'") ("`ev'") (`reltime') (`bb') (`ss') (`lo') (`hi') (`pval') ///
            (`pre_event_mean') (`post03_raw_mean') (`sample_obs') (`sample_ids')
    }

    post `ph' ("`y'") ("F${ref_sa}event") (-${ref_sa}) (0) (0) (0) (0) (.) ///
        (`pre_event_mean') (`post03_raw_mean') (`sample_obs') (`sample_ids')
}

postclose `rh'
postclose `ph'

use "`raw'", clear
sort outcome relative_time
save "`outdir'/moving_first_birth_raw_means_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_raw_means_v1.csv", replace

foreach y of local outcomes {
    local ytitle "`y'"
    if "`y'" == "move_any" local ytitle "Any move"
    if "`y'" == "moved_for_size" local ytitle "Moved for size"
    if "`y'" == "moved_to_own" local ytitle "Moved to ownership"

    preserve
        keep if outcome == "`y'"
        gen ci_lo = raw_mean - 1.96 * raw_se
        gen ci_hi = raw_mean + 1.96 * raw_se
        twoway ///
            (rarea ci_lo ci_hi relative_time, fcolor(maroon%10) lcolor(maroon%10)) ///
            (connected raw_mean relative_time, mcolor(maroon) msymbol(circle) lcolor(maroon) lwidth(medthick)), ///
            ytitle("`ytitle' raw mean") xtitle("Years relative to first birth") ///
            xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
            $graphregion_slides legend(off) title("Raw mean: `ytitle'")
        graph export "`outdir'/`y'_raw_mean_v1.png", replace
    restore
}

use "`pooled'", clear
sort outcome relative_time
save "`outdir'/moving_first_birth_pooled_controls_estimates_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_pooled_controls_estimates_v1.csv", replace

foreach y of local outcomes {
    local ytitle "`y'"
    if "`y'" == "move_any" local ytitle "Any move"
    if "`y'" == "moved_for_size" local ytitle "Moved for size"
    if "`y'" == "moved_to_own" local ytitle "Moved to ownership"

    preserve
        keep if outcome == "`y'"
        sort relative_time
        twoway ///
            (rarea ci_lo ci_hi relative_time, fcolor(dkgreen%10) lcolor(dkgreen%10)) ///
            (connected b relative_time, mcolor(dkgreen) msymbol(circle) lcolor(dkgreen) lwidth(medthick) lpattern(dash)), ///
            ytitle("`ytitle'") xtitle("Years relative to first birth") ///
            xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
            yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
            $graphregion_slides legend(off) title("Age/year controls: `ytitle'")
        graph export "`outdir'/`y'_pooled_controls_v1.png", replace
    restore
}

di as result "Moving event-study diagnostics complete."
di as result "Output dir: `outdir'"

log close _all
