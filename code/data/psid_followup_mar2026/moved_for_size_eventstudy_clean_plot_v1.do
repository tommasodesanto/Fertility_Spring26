clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local outdir "`root'/code/data/psid_followup_mar2026/output/moving_first_birth_eventstudy_diagnostics_v1"
local infile "`outdir'/moving_first_birth_pooled_controls_estimates_v1.csv"

capture log close _all
log using "`outdir'/moved_for_size_eventstudy_clean_plot_v1.log", replace text

global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

import delimited using "`infile'", clear
keep if outcome == "moved_for_size"
sort relative_time
gen b_pct = 100 * b
gen ci_lo_pct = 100 * ci_lo
gen ci_hi_pct = 100 * ci_hi

twoway ///
    (rarea ci_lo_pct ci_hi_pct relative_time, fcolor(maroon%12) lcolor(maroon%12)) ///
    (connected b_pct relative_time, mcolor(maroon) msymbol(circle) msize(medlarge) ///
        lcolor(maroon) lwidth(medthick) lpattern(solid)), ///
    ytitle("Effect on move-for-size probability, pp") ///
    xtitle("Years relative to first birth") ///
    xlabel(-7 "-7 or earlier" -5 "-5" -4 "-4" -3 "-3" -2 "-2" -1 "-1" ///
        0 "0" 1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6" 7 "7" ///
        8 "8" 9 "9" 10 "10" 11 "11+", angle(0)) ///
    xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
    yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
    ylabel(-4(2)6, angle(horizontal)) ///
    title("Move-for-Size Around First Birth") ///
    note("Event-time regression with year, age, and education controls; K=-2 omitted. 95% CI clustered by household.", size(vsmall)) ///
    $graphregion_slides legend(off)

graph export "`outdir'/moved_for_size_eventstudy_clean_v1.png", replace width(1600)
graph export "`outdir'/moved_for_size_eventstudy_clean_v1.pdf", replace

di as result "Clean moved-for-size event-study plot written:"
di as result "`outdir'/moved_for_size_eventstudy_clean_v1.pdf"

log close _all
