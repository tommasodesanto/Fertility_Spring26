clear all
set more off

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
use ID year AGEREP RELCHIREP RELCHI1BYEAR IW using "`dta'", clear

keep if inrange(AGEREP, 20, 50)

* Distribution of RELCHIREP
di as text "RELCHIREP distribution (all):"
tabulate RELCHIREP, missing

di as text "RELCHIREP at ages 35-40:"
preserve
    keep if inrange(AGEREP, 35, 40)
    tabulate RELCHIREP, missing
restore

di as text "RELCHI1BYEAR distribution (head of sample):"
summ RELCHI1BYEAR
tabulate RELCHI1BYEAR if inrange(AGEREP,35,40) & !missing(RELCHI1BYEAR), missing

di as text "How many with first-birth year observed at 35-40 vs missing?"
gen has_first_birth = !missing(RELCHI1BYEAR)
tab has_first_birth if inrange(AGEREP,35,40)

di as text "RELCHIREP stats by whether RELCHI1BYEAR present"
summ RELCHIREP if inrange(AGEREP,35,40) & has_first_birth == 1
summ RELCHIREP if inrange(AGEREP,35,40) & has_first_birth == 0
