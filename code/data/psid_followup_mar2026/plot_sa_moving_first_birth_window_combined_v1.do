clear all
set more off
version 17.0

local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_moving_first_birth_window_v1"

global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

import delimited using "`outdir'/moving_first_birth_sa_window_all_estimates_v1.csv", clear

foreach y in move_any moved_for_size moved_to_own {
    preserve
        keep if outcome == "`y'"
        sort relative_time

        local ytitle = label[1]
        local ymin = -6
        local ymax = 3
        if "`y'" == "moved_for_size" {
            local ymin = -4
            local ymax = 5
        }
        if "`y'" == "moved_to_own" {
            local ymin = -2
            local ymax = 2
        }

        twoway ///
            (rarea ci_lo_pct ci_hi_pct relative_time, fcolor(navy%10) lcolor(navy%10)) ///
            (connected b_pct relative_time, mcolor(navy) msymbol(circle) msize(medium) ///
                lcolor(navy) lwidth(medthick)), ///
            ytitle("pp") xtitle("") ///
            xlabel(-5(1)5) ///
            ylabel(`ymin'(2)`ymax', angle(horizontal)) ///
            yscale(range(`ymin' `ymax')) ///
            xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
            yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
            title("`ytitle'", size(medsmall)) ///
            $graphregion_slides legend(off) name(g_`y', replace)
    restore
}

graph combine g_move_any g_moved_for_size g_moved_to_own, ///
    cols(1) imargin(tiny) ///
    title("Moving Around First Birth", size(medium)) ///
    note("Sun-Abraham interaction-weighted event studies; K=-2 omitted. 95% CI clustered by household.", size(vsmall)) ///
    graphregion(fcolor("$bg") lcolor("$bg"))

graph export "`outdir'/moving_sa_first_birth_window_combined_v1.png", replace width(1600)
graph export "`outdir'/moving_sa_first_birth_window_combined_v1.pdf", replace

di as result "Combined Sun-Abraham moving plot complete."
di as result "Output: `outdir'/moving_sa_first_birth_window_combined_v1.png"
