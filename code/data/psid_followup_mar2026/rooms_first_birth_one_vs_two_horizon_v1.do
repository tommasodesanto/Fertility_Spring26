clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/rooms_first_birth_one_vs_two_horizon_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/rooms_first_birth_one_vs_two_horizon_v1.log", replace text

global weight IW

use "`dta'", clear
keep ID year AGEREP DEATHYEAR RELCHI1BYEAR RELCHI2BYEAR ACTUALROOMS_ ${weight}

drop if year > DEATHYEAR
drop if AGEREP < 18

bysort ID: egen first_child_year = min(RELCHI1BYEAR)
bysort ID: egen second_child_year = min(RELCHI2BYEAR)
drop if missing(first_child_year)

bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

* Drop multiple births at first birth.
drop if !missing(second_child_year) & second_child_year == first_child_year

rename ACTUALROOMS_ rooms
gen K = year - first_child_year

tempfile base
preserve
    keep if K == -2
    keep ID rooms ${weight}
    rename rooms rooms_pre
    rename ${weight} wt_pre
    save `base'
restore

tempname posth
postfile `posth' str24 horizon_group double horizon years_post mean_pre mean_post mean_change N_ids using "`outdir'/rooms_first_birth_one_vs_two_horizon_v1.dta", replace

foreach h in 2 3 {
    foreach grp in one_kid two_plus {
        preserve
            if "`grp'" == "one_kid" {
                keep if missing(second_child_year) | second_child_year > first_child_year + `h'
            }
            else {
                keep if !missing(second_child_year) & second_child_year <= first_child_year + `h'
            }

            keep if inlist(K, -2, `h')
            merge m:1 ID using `base', keep(match) nogen
            keep if K == `h'

            quietly summarize rooms_pre [aw=wt_pre]
            local m_pre = r(mean)
            quietly summarize rooms [aw=${weight}]
            local m_post = r(mean)
            quietly count if !missing(ID)
            local n_ids = r(N)
            local d = `m_post' - `m_pre'
            post `posth' ("`grp'_by`h'") (`h') (`h') (`m_pre') (`m_post') (`d') (`n_ids')
        restore
    }
}

postclose `posth'

use "`outdir'/rooms_first_birth_one_vs_two_horizon_v1.dta", clear
export delimited using "`outdir'/rooms_first_birth_one_vs_two_horizon_v1.csv", replace
list, clean noobs

di as result "Summary: `outdir'/rooms_first_birth_one_vs_two_horizon_v1.csv"

log close _all
