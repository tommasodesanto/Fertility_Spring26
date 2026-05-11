clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/wealth_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/wealth_regressions_v1.log", replace text

global weight IW

di as text "Loading PSID wealth and housing variables..."
use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR HOMEOWN MOVEDFREF_ ACTUALROOMS_ ///
    INCFAMR EARNINDR NETWORTHR NETWORTH2R NETWORTH3R HOMEEQUITYR HOMEMORTOTR ///
    WLTHSAVETOTR WLTHFUNDTOTR WLTHODEBTOTR ${weight} using "`dta'", clear

* Keep the same tenure coding as the baseline PSID scripts.
replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
replace movedthisyear = 0 if movedthisyear == 5

rename ACTUALROOMS_ rooms
rename HOMEOWN own
gen first_child_year = RELCHI1BYEAR
drop if missing(first_child_year)

* Exclude first births observed before panel entry.
bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

xtset ID year
gen K = year - first_child_year

* Baseline characteristics at t = -2 (pre-birth).
preserve
    keep if K == -2
    keep ID first_child_year AGEREP EDUYEAR SEX own rooms INCFAMR EARNINDR ///
         NETWORTHR NETWORTH2R NETWORTH3R HOMEEQUITYR HOMEMORTOTR ///
         WLTHSAVETOTR WLTHFUNDTOTR WLTHODEBTOTR ${weight}

    rename AGEREP age_pre
    rename EDUYEAR eduyear_pre
    rename SEX sex_pre
    rename own own_pre
    rename rooms rooms_pre
    rename INCFAMR incfamr_pre
    rename EARNINDR earnindr_pre
    rename NETWORTHR networthr_pre
    rename NETWORTH2R networth2r_pre
    rename NETWORTH3R networth3r_pre
    rename HOMEEQUITYR homeequityr_pre
    rename HOMEMORTOTR homemortotr_pre
    rename WLTHSAVETOTR wlthsavetotr_pre
    rename WLTHFUNDTOTR wlthfundtotr_pre
    rename WLTHODEBTOTR wlthodebtotr_pre
    rename ${weight} iw_pre

    save "`outdir'/tmp_baseline_tminus2.dta", replace
restore

* Post-birth outcomes in the first four years after childbirth.
gen own_post3_i = (own == 1) & inrange(K, 0, 3)
gen move_post3_i = (movedthisyear == 1) & inrange(K, 0, 3)
gen moved_to_own_post3_i = (own == 1 & L.own == 0 & movedthisyear == 1) & inrange(K, 0, 3)
gen rooms_post3_i = rooms if inrange(K, 0, 3)

bysort ID: egen own_post3 = max(own_post3_i)
bysort ID: egen move_post3 = max(move_post3_i)
bysort ID: egen moved_to_own_post3 = max(moved_to_own_post3_i)
bysort ID: egen rooms_post3 = mean(rooms_post3_i)

bysort ID (year): keep if _n == 1
keep ID first_child_year own_post3 move_post3 moved_to_own_post3 rooms_post3

merge 1:1 ID first_child_year using "`outdir'/tmp_baseline_tminus2.dta", keep(match) nogen
erase "`outdir'/tmp_baseline_tminus2.dta"

gen rooms_change_post3 = rooms_post3 - rooms_pre if !missing(rooms_post3, rooms_pre)
gen female_pre = (sex_pre == 2)

* Focus on households that were renters before childbirth.
keep if own_pre == 0
drop if missing(networth2r_pre, age_pre, eduyear_pre, female_pre, own_post3, iw_pre)

* Wealth transformations robust to zeros/negatives.
gen asinh_nw_pre = asinh(networthr_pre)
gen asinh_nw2_pre = asinh(networth2r_pre)
gen asinh_sav_pre = asinh(wlthsavetotr_pre)
gen asinh_fund_pre = asinh(wlthfundtotr_pre)
gen asinh_odebt_pre = asinh(wlthodebtotr_pre)

xtile nw2_tercile = networth2r_pre, n(3)

save "`outdir'/wealth_event_sample_pre_renters_v1.dta", replace

di as text "Running wealth regressions..."
estimates clear

reg own_post3 c.asinh_nw2_pre c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
estimates store r1

reg own_post3 c.asinh_sav_pre c.asinh_fund_pre c.asinh_odebt_pre c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
estimates store r2

reg moved_to_own_post3 c.asinh_nw2_pre c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre], vce(robust)
estimates store r3

reg rooms_change_post3 c.asinh_nw2_pre c.age_pre i.eduyear_pre i.first_child_year i.female_pre [aw=iw_pre] if !missing(rooms_change_post3), vce(robust)
estimates store r4

capture which esttab
if _rc == 0 {
    esttab r1 r2 r3 r4 using "`outdir'/wealth_transition_regressions_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("Own by +3y" "Own by +3y (components)" "Moved-to-own by +3y" "Rooms change by +3y") ///
        stats(N r2, fmt(0 3) labels("Obs." "R-squared"))
}

preserve
    collapse (mean) own_post3 moved_to_own_post3 move_post3 rooms_change_post3 [aw=iw_pre], by(nw2_tercile)
    export delimited using "`outdir'/wealth_tercile_outcomes_v1.csv", replace
restore

di as result "Wealth regressions complete."
di as result "Sample: `outdir'/wealth_event_sample_pre_renters_v1.dta"
di as result "Regression table: `outdir'/wealth_transition_regressions_v1.tex (if esttab installed)"
di as result "Tercile summary: `outdir'/wealth_tercile_outcomes_v1.csv"

log close _all
