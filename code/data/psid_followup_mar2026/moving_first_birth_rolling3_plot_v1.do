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
log using "`outdir'/moving_first_birth_rolling3_plot_v1.log", replace text

global weight IW
global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

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

tempfile rolling3
tempname ph
postfile `ph' str32 outcome double center start end str12 window_label mean n_ids using "`rolling3'", replace

foreach y in move_any moved_for_size moved_to_own {
    forvalues center = -4/4 {
        local start = `center' - 1
        local end = `center' + 1
        local label "K=`start'..`end'"

        preserve
            keep if inrange(K, `start', `end')
            keep if !missing(`y')
            bysort ID: egen ever_y = max(`y')
            bysort ID: keep if _n == 1
            quietly summarize ever_y [aw=${weight}]
            local mean = r(mean)
            quietly count
            local n_ids = r(N)
            post `ph' ("`y'") (`center') (`start') (`end') ("`label'") (`mean') (`n_ids')
        restore
    }
}

postclose `ph'

use "`rolling3'", clear
gen rate_pct = 100 * mean
sort outcome center
save "`outdir'/moving_first_birth_rolling3_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_rolling3_v1.csv", replace

twoway ///
    (connected rate_pct center if outcome == "move_any", ///
        mcolor(navy) msymbol(circle) msize(medlarge) lcolor(navy) lwidth(medthick)), ///
    ytitle("Percent with any move") ///
    xtitle("") ///
    xlabel(-4 "K=-5..-3" -3 "K=-4..-2" -2 "K=-3..-1" -1 "K=-2..0" ///
        0 "K=-1..1" 1 "K=0..2" 2 "K=1..3" 3 "K=2..4" 4 "K=3..5", angle(35)) ///
    xline(0, lwidth(thin) lpattern(dash) lcolor(gray)) ///
    ylabel(0(20)80, angle(horizontal)) ///
    title("A. Any move") ///
    $graphregion_slides legend(off)
graph save "`outdir'/moving_first_birth_rolling3_any_v1.gph", replace

twoway ///
    (connected rate_pct center if outcome == "moved_for_size", ///
        mcolor(maroon) msymbol(circle) msize(medlarge) lcolor(maroon) lwidth(medthick)) ///
    (connected rate_pct center if outcome == "moved_to_own", ///
        mcolor(dkgreen) msymbol(square) msize(medlarge) lcolor(dkgreen) lwidth(medthick)), ///
    ytitle("Percent with move type") ///
    xtitle("Rolling 3-year window around first birth") ///
    xlabel(-4 "K=-5..-3" -3 "K=-4..-2" -2 "K=-3..-1" -1 "K=-2..0" ///
        0 "K=-1..1" 1 "K=0..2" 2 "K=1..3" 3 "K=2..4" 4 "K=3..5", angle(35)) ///
    xline(0, lwidth(thin) lpattern(dash) lcolor(gray)) ///
    ylabel(0(10)30, angle(horizontal)) ///
    title("B. Housing-relevant moves") ///
    $graphregion_slides ///
    legend(order(1 "Moved for size" 2 "Moved to ownership") cols(2) region(lcolor(none)))
graph save "`outdir'/moving_first_birth_rolling3_types_v1.gph", replace

graph combine "`outdir'/moving_first_birth_rolling3_any_v1.gph" ///
    "`outdir'/moving_first_birth_rolling3_types_v1.gph", ///
    cols(1) imargin(2 2 2 2) ///
    title("Moving Around First Birth", size(medlarge)) ///
    note("Weighted share of eventual first-birth parents with at least one event in each rolling 3-year window.", size(vsmall)) ///
    graphregion(fcolor("$bg") lcolor("$bg"))

graph export "`outdir'/moving_first_birth_rolling3_combined_v1.png", replace width(1600)
graph export "`outdir'/moving_first_birth_rolling3_combined_v1.pdf", replace

di as result "Rolling 3-year plot complete."
di as result "PNG: `outdir'/moving_first_birth_rolling3_combined_v1.png"

log close _all
