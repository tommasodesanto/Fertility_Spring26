clear all
set more off

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
use ID year AGEREP RELCHIREP RELCHI1BYEAR INCFAMR NETWORTH2R IW using "`dta'", clear

drop if missing(AGEREP) | AGEREP < 18
gen childless = (RELCHIREP == 0) if !missing(RELCHIREP)
gen ever_parent = !missing(RELCHI1BYEAR)
gen liq_nw_to_inc = NETWORTH2R / INCFAMR if INCFAMR > 1000 & !missing(NETWORTH2R, INCFAMR)

preserve
    keep if inrange(AGEREP, 25, 30)
    keep if childless == 1
    collapse (mean) liq_nw_to_inc_young = liq_nw_to_inc [pweight=IW], by(ID)
    keep if !missing(liq_nw_to_inc_young)
    di "wealth_young sample:"
    count
    summ liq_nw_to_inc_young, detail
    tempfile wealth_young
    save `wealth_young', replace
restore

preserve
    keep if inrange(AGEREP, 35, 40)
    di "n_children distribution at 35-40 (raw):"
    tab RELCHIREP, missing
    collapse (mean) n_children_35 = RELCHIREP (max) ever_parent_35 = ever_parent ///
        (mean) weight_35 = IW, by(ID)
    di "fertility_35 sample:"
    count
    summ n_children_35 ever_parent_35, detail
    tempfile fertility_35
    save `fertility_35', replace
restore

use `wealth_young', clear
merge 1:1 ID using `fertility_35', keep(match) nogen
di "Merged sample:"
count
summ n_children_35 ever_parent_35 liq_nw_to_inc_young, detail

di "Tabulation of n_children_35:"
tab n_children_35, missing
