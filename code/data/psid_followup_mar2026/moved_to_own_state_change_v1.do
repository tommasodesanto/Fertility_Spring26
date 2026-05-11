clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/moved_to_own_state_change_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/moved_to_own_state_change_v1.log", replace text

global weight IW

use ID year DEATHYEAR AGEREP RELCHI1BYEAR HOMEOWN MOVEDFREF_ GEOSTATE STATEFIPS_ ${weight} using "`dta'", clear

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
replace movedthisyear = 0 if movedthisyear == 5

rename HOMEOWN own

xtset ID year
gen change_own = (own == 1 & L.own == 0) if !missing(own, L.own)
gen moved_to_own = (movedthisyear == 1 & change_own == 1) if !missing(movedthisyear, change_own)

gen interstate_geo = (GEOSTATE != L.GEOSTATE) if !missing(GEOSTATE, L.GEOSTATE)
gen interstate_fips = (STATEFIPS_ != L.STATEFIPS_) if !missing(STATEFIPS_, L.STATEFIPS_)

gen first_child_year = RELCHI1BYEAR
gen K = year - first_child_year
gen moved_to_own_post3 = (moved_to_own == 1 & inrange(K, 0, 3)) if !missing(moved_to_own, K)
gen moved_to_own_post5 = (moved_to_own == 1 & inrange(K, 0, 5)) if !missing(moved_to_own, K)

tempname fh
file open `fh' using "`outdir'/moved_to_own_state_change_summary_v1.csv", write replace
file write `fh' "sample,measure,value" _n

quietly count if moved_to_own == 1
local n_mto = r(N)
file write `fh' "all,moved_to_own_events,`n_mto'" _n

quietly count if moved_to_own == 1 & interstate_geo == 1
local n_mto_geo_inter = r(N)
file write `fh' "all,moved_to_own_interstate_geo_events,`n_mto_geo_inter'" _n

quietly summarize interstate_geo if moved_to_own == 1
local share_geo = r(mean)
file write `fh' "all,share_interstate_geo_unweighted,`share_geo'" _n

quietly summarize interstate_geo [aw=${weight}] if moved_to_own == 1
local share_geo_w = r(mean)
file write `fh' "all,share_interstate_geo_weighted,`share_geo_w'" _n

quietly summarize interstate_fips if moved_to_own == 1
local share_fips = r(mean)
file write `fh' "all,share_interstate_fips_unweighted,`share_fips'" _n

quietly summarize interstate_fips [aw=${weight}] if moved_to_own == 1
local share_fips_w = r(mean)
file write `fh' "all,share_interstate_fips_weighted,`share_fips_w'" _n

quietly count if moved_to_own_post3 == 1
local n_post3 = r(N)
file write `fh' "post3,moved_to_own_events,`n_post3'" _n

quietly summarize interstate_geo if moved_to_own_post3 == 1
local share_post3_geo = r(mean)
file write `fh' "post3,share_interstate_geo_unweighted,`share_post3_geo'" _n

quietly summarize interstate_geo [aw=${weight}] if moved_to_own_post3 == 1
local share_post3_geo_w = r(mean)
file write `fh' "post3,share_interstate_geo_weighted,`share_post3_geo_w'" _n

quietly count if moved_to_own_post5 == 1
local n_post5 = r(N)
file write `fh' "post5,moved_to_own_events,`n_post5'" _n

quietly summarize interstate_geo if moved_to_own_post5 == 1
local share_post5_geo = r(mean)
file write `fh' "post5,share_interstate_geo_unweighted,`share_post5_geo'" _n

quietly summarize interstate_geo [aw=${weight}] if moved_to_own_post5 == 1
local share_post5_geo_w = r(mean)
file write `fh' "post5,share_interstate_geo_weighted,`share_post5_geo_w'" _n

file close `fh'

preserve
    keep if moved_to_own == 1
    keep ID year first_child_year K interstate_geo interstate_fips ${weight}
    save "`outdir'/moved_to_own_event_rows_v1.dta", replace
restore

di as result "Completed moved-to-own interstate check."
di as result "Summary CSV: `outdir'/moved_to_own_state_change_summary_v1.csv"
di as result "Event-level rows: `outdir'/moved_to_own_event_rows_v1.dta"

log close _all
