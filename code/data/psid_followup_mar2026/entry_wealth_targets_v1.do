clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/entry_wealth_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/entry_wealth_targets_v1.log", replace text

global weight IW

di as text "Loading PSID wealth data for candidate entry-wealth targets..."
use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHIREP RELCHI1BYEAR HOMEOWN ///
    INCFAMR EARNINDR NETWORTHR NETWORTH2R NETWORTH3R HOMEEQUITYR HOMEMORTOTR ///
    WLTHSAVETOTR WLTHFUNDTOTR WLTHODEBTOTR ${weight} using "`dta'", clear

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

gen own = HOMEOWN
gen n_children = RELCHIREP
gen childless = (n_children == 0) if !missing(n_children)
gen first_child_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_child_year) | year < first_child_year)

gen young_25_35 = inrange(AGEREP, 25, 35)
gen young_25_30 = inrange(AGEREP, 25, 30)

* Candidate liquid-wealth objects.
gen liquid_nw_ex_housing = NETWORTH2R
gen liquid_savings       = WLTHSAVETOTR

gen liq_nw_to_earn = liquid_nw_ex_housing / EARNINDR if EARNINDR > 1000 & !missing(liquid_nw_ex_housing, EARNINDR)
gen liq_nw_to_inc  = liquid_nw_ex_housing / INCFAMR if INCFAMR > 1000 & !missing(liquid_nw_ex_housing, INCFAMR)
gen sav_to_earn    = liquid_savings / EARNINDR if EARNINDR > 1000 & !missing(liquid_savings, EARNINDR)

tempname handle
postfile `handle' str32 sample str24 moment double mean p50 N using ///
    "`outdir'/entry_wealth_candidate_targets_v1.dta", replace

* Stata's centile does not accept analytic weights, so we report
* weighted means and unweighted medians for candidate target construction.

local moments liquid_nw_ex_housing liquid_savings liq_nw_to_earn liq_nw_to_inc sav_to_earn

local sample1_name young_all_25_35
local sample1_if   young_25_35 == 1

local sample2_name young_childless_25_35
local sample2_if   young_25_35 == 1 & childless == 1

local sample3_name young_childless_rent_25_35
local sample3_if   young_25_35 == 1 & childless == 1 & own == 0

local sample4_name young_prebirth_25_35
local sample4_if   young_25_35 == 1 & pre_first_birth == 1

local sample5_name young_prebirth_rent_25_35
local sample5_if   young_25_35 == 1 & pre_first_birth == 1 & own == 0

local sample6_name young_all_25_30
local sample6_if   young_25_30 == 1

local sample7_name young_childless_25_30
local sample7_if   young_25_30 == 1 & childless == 1

local sample8_name young_childless_rent_25_30
local sample8_if   young_25_30 == 1 & childless == 1 & own == 0

forvalues s = 1/8 {
    local sname = "`sample`s'_name'"
    local sif   = "`sample`s'_if'"

    foreach m of local moments {
        quietly count if `sif' & !missing(`m')
        local N = r(N)

        if `N' > 0 {
            quietly summarize `m' [aw=${weight}] if `sif' & !missing(`m')
            local mean = r(mean)

            quietly centile `m' if `sif' & !missing(`m'), centile(50)
            local p50 = r(c_1)

            post `handle' ("`sname'") ("`m'") (`mean') (`p50') (`N')
        }
    }
}

postclose `handle'

use "`outdir'/entry_wealth_candidate_targets_v1.dta", clear
order sample moment mean p50 N
sort sample moment
export delimited using "`outdir'/entry_wealth_candidate_targets_v1.csv", replace

preserve
    keep if inlist(sample, "young_childless_25_35", "young_childless_rent_25_35", ///
        "young_prebirth_25_35", "young_prebirth_rent_25_35")
    export delimited using "`outdir'/entry_wealth_candidate_targets_focus_v1.csv", replace
restore

di as result "Entry-wealth candidate target extraction complete."
di as result "Main output: `outdir'/entry_wealth_candidate_targets_v1.csv"
di as result "Focused output: `outdir'/entry_wealth_candidate_targets_focus_v1.csv"

log close _all
