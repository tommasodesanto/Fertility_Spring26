clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/rooms_second_birth_quick_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/rooms_second_birth_quick_moment_v1.log", replace text

global weight IW

di as text "Loading PSID panel for quick second-birth rooms moment..."
use ID year AGEREP DEATHYEAR RELCHI1BYEAR RELCHI2BYEAR ACTUALROOMS_ ${weight} using "`dta'", clear

drop if year > DEATHYEAR
drop if AGEREP < 18

gen first_child_year  = RELCHI1BYEAR
gen second_child_year = RELCHI2BYEAR

drop if missing(second_child_year)
drop if !missing(first_child_year) & second_child_year == first_child_year

bysort ID: egen year_entry = min(year)
drop if second_child_year < year_entry
drop year_entry

rename ACTUALROOMS_ rooms

xtset ID year
gen K = year - second_child_year

preserve
    keep if K == -2
    keep ID second_child_year rooms ${weight}
    rename rooms rooms_pre
    rename ${weight} iw_pre
    save "`outdir'/tmp_rooms_pre_v1.dta", replace
restore

keep if inrange(K, 0, 3)
bysort ID second_child_year: egen rooms_post = mean(rooms)
bysort ID second_child_year: egen iw_post = mean(${weight})
bysort ID second_child_year (year): keep if _n == 1
keep ID second_child_year rooms_post iw_post

merge 1:1 ID second_child_year using "`outdir'/tmp_rooms_pre_v1.dta", keep(match) nogen
erase "`outdir'/tmp_rooms_pre_v1.dta"

gen rooms_change_post3 = rooms_post - rooms_pre if !missing(rooms_post, rooms_pre)
keep if !missing(rooms_change_post3)

quietly summarize rooms_pre [aw=iw_pre]
local mean_pre = r(mean)

quietly summarize rooms_post [aw=iw_pre]
local mean_post = r(mean)

quietly summarize rooms_change_post3 [aw=iw_pre]
local mean_change = r(mean)

quietly centile rooms_change_post3, centile(50)
local p50_change = r(c_1)

quietly count
local N = r(N)

clear
set obs 1
gen mean_rooms_pre = `mean_pre'
gen mean_rooms_post = `mean_post'
gen mean_rooms_change_post3 = `mean_change'
gen median_rooms_change_post3 = `p50_change'
gen N = `N'

save "`outdir'/rooms_second_birth_quick_moment_v1.dta", replace
export delimited using "`outdir'/rooms_second_birth_quick_moment_v1.csv", replace

di as result "Quick second-birth rooms moment complete."
di as result "Mean rooms, t=-2: `mean_pre'"
di as result "Mean rooms, t=0..3: `mean_post'"
di as result "Mean change: `mean_change'"
di as result "Median change: `p50_change'"
di as result "Observations: `N'"

log close _all
