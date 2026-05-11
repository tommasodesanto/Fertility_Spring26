clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/moving_first_birth_window_summary_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/moving_first_birth_window_summary_v1.log", replace text

global weight IW

use ID year AGEREP DEATHYEAR RELCHI1BYEAR HOMEOWN MOVEDFREF_ WHYMOVED1_ ${weight} using "`dta'", clear

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

tempname fh
file open `fh' using "`outdir'/moving_first_birth_window_summary_v1.csv", write replace
file write `fh' "outcome,window,mean,n_ids,n_obs" _n

foreach y in move_any moved_for_size moved_to_own {
    foreach win in pre5pre3 pre2 pre1to1 post0to3 post0to5 {
        preserve
            if "`win'" == "pre5pre3" keep if inrange(K, -5, -3)
            if "`win'" == "pre2" keep if K == -2
            if "`win'" == "pre1to1" keep if inrange(K, -1, 1)
            if "`win'" == "post0to3" keep if inrange(K, 0, 3)
            if "`win'" == "post0to5" keep if inrange(K, 0, 5)

            keep if !missing(`y')
            bysort ID: egen ever_y = max(`y')
            bysort ID: keep if _n == 1
            quietly summarize ever_y [aw=${weight}]
            local mean = r(mean)
            quietly count
            local n_ids = r(N)
            local n_obs = r(N)
            file write `fh' "`y',`win',`mean',`n_ids',`n_obs'" _n
        restore
    }
}

file close `fh'

import delimited using "`outdir'/moving_first_birth_window_summary_v1.csv", clear

gen window_order = .
replace window_order = 1 if window == "pre5pre3"
replace window_order = 2 if window == "pre2"
replace window_order = 3 if window == "pre1to1"
replace window_order = 4 if window == "post0to3"
replace window_order = 5 if window == "post0to5"

gen window_label = ""
replace window_label = "K=-5..-3" if window == "pre5pre3"
replace window_label = "K=-2" if window == "pre2"
replace window_label = "K=-1..1" if window == "pre1to1"
replace window_label = "K=0..3" if window == "post0to3"
replace window_label = "K=0..5" if window == "post0to5"

foreach y in move_any moved_for_size moved_to_own {
    preserve
        keep if outcome == "`y'"
        sort window_order
        gen rate_pct = 100 * mean

        local ytitle "`y'"
        if "`y'" == "move_any" local ytitle "Any move"
        if "`y'" == "moved_for_size" local ytitle "Moved for size"
        if "`y'" == "moved_to_own" local ytitle "Moved to ownership"

        graph bar rate_pct, over(window_label, sort(window_order) label(angle(35))) ///
            ytitle("Percent of first-birth parents") ///
            title("`ytitle' around first birth") ///
            bar(1, fcolor(navy%70) lcolor(navy)) ///
            blabel(bar, format(%4.1f)) ///
            graphregion(fcolor("250 250 250") lcolor("250 250 250")) ///
            plotregion(fcolor("250 250 250") lcolor("250 250 250"))
        graph export "`outdir'/`y'_window_rates_v1.png", replace
    restore
}

di as result "Window summary complete."
di as result "Output: `outdir'/moving_first_birth_window_summary_v1.csv"

log close _all
