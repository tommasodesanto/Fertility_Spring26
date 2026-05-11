clear all
set more off
version 17.0

local in "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/moved_to_own_state_change_v1/moved_to_own_event_rows_v1.dta"
local out "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/moved_to_own_state_change_v1"

capture log close _all
log using "`out'/moved_to_own_state_change_compare_v1.log", replace text

use "`in'", clear

di as text "Counts and means"
count if interstate_geo < .
di as text "N nonmissing interstate_geo = " r(N)
sum interstate_geo

count if interstate_fips < .
di as text "N nonmissing interstate_fips = " r(N)
sum interstate_fips

di as text "Agreement where both present"
count if interstate_geo < . & interstate_fips < .
di as text "N both = " r(N)
gen agree = (interstate_geo == interstate_fips) if interstate_geo < . & interstate_fips < .
sum agree

di as text "Mismatch table (both present)"
tab interstate_geo interstate_fips if interstate_geo < . & interstate_fips < .

di as text "Post-birth windows"
sum interstate_geo if inrange(K, 0, 3)
sum interstate_geo if inrange(K, 0, 5)

log close _all
