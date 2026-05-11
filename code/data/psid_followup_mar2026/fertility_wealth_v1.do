clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"

cap mkdir "`outdir'"
capture log close _all
log using "`outdir'/fertility_wealth_v1.log", replace text

global weight IW

di as text "Loading PSID for fertility-wealth descriptives..."
use ID year AGEREP SEX EDUYEAR DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR EARNINDR NETWORTH2R HOMEOWN ${weight} using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

gen n_children = RELCHIREP
gen childless = (n_children == 0) if !missing(n_children)
gen first_birth_year = RELCHI1BYEAR
gen ever_parent = !missing(first_birth_year)

gen liq_nw_to_inc = NETWORTH2R / INCFAMR if INCFAMR > 1000 & !missing(NETWORTH2R, INCFAMR)

* Restrict attention to people we actually observe at young and at 35+ ages,
* so wealth-at-young and fertility-at-35 can both be measured for them.

* ------------------------------------------------------------------
* (A) Completed fertility at age 35-40 by wealth tercile at 25-30
* ------------------------------------------------------------------

preserve
    keep if inrange(AGEREP, 25, 30)
    keep if childless == 1
    collapse (mean) liq_nw_to_inc_young = liq_nw_to_inc [pweight=${weight}], by(ID)
    tempfile wealth_young
    save `wealth_young', replace
restore

preserve
    keep if inrange(AGEREP, 35, 40)
    collapse (mean) n_children_35 = n_children (max) ever_parent_35 = ever_parent ///
        (mean) weight_35 = ${weight}, by(ID)
    tempfile fertility_35
    save `fertility_35', replace
restore

use `wealth_young', clear
merge 1:1 ID using `fertility_35', keep(match) nogen
keep if !missing(liq_nw_to_inc_young, n_children_35)

* Drop extreme tails
_pctile liq_nw_to_inc_young [pw=weight_35], p(1 99)
drop if liq_nw_to_inc_young < r(r1) | liq_nw_to_inc_young > r(r2)

xtile wealth_tercile = liq_nw_to_inc_young [pw=weight_35], n(3)
xtile wealth_decile  = liq_nw_to_inc_young [pw=weight_35], n(10)

* Tercile summary
preserve
    collapse (mean) children_at_35 = n_children_35 ever_parent_at_35 = ever_parent_35 ///
        (p50) median_wealth = liq_nw_to_inc_young (count) n_obs = ID [pweight=weight_35], ///
        by(wealth_tercile)
    list
    export delimited using "`outdir'/fertility_by_wealth_tercile.csv", replace
restore

* Decile summary for plotting
preserve
    collapse (mean) children_at_35 = n_children_35 ever_parent_at_35 = ever_parent_35 ///
        (p50) median_wealth = liq_nw_to_inc_young [pweight=weight_35], by(wealth_decile)
    list
    export delimited using "`outdir'/fertility_by_wealth_decile.csv", replace

    * Plot: children at 35 by wealth decile
    twoway (connected children_at_35 wealth_decile, mcolor(navy) lcolor(navy) lwidth(medthick)) ///
        , xtitle("Wealth decile at ages 25-30 (liquid NW / family income)") ///
          ytitle("Children ever born by age 35-40") ///
          xlabel(1(1)10) ///
          title("Completed fertility by entry-wealth") ///
          graphregion(color(white)) ///
          scheme(s1color)
    graph export "`outdir'/fertility_by_wealth_decile.png", replace width(1600)

    * Ever-parent version
    twoway (connected ever_parent_at_35 wealth_decile, mcolor(maroon) lcolor(maroon) lwidth(medthick)) ///
        , xtitle("Wealth decile at ages 25-30 (liquid NW / family income)") ///
          ytitle("Share who are parents by age 35-40") ///
          xlabel(1(1)10) ///
          title("Probability of parenthood by entry-wealth") ///
          graphregion(color(white)) ///
          scheme(s1color)
    graph export "`outdir'/parenthood_by_wealth_decile.png", replace width(1600)
restore

* ------------------------------------------------------------------
* (B) Annual first-birth hazard by wealth tercile (among childless)
* ------------------------------------------------------------------

use ID year AGEREP SEX RELCHIREP RELCHI1BYEAR INCFAMR NETWORTH2R DEATHYEAR ${weight} using "`dta'", clear
drop if year > DEATHYEAR
drop if missing(AGEREP) | !inrange(AGEREP, 22, 40)
keep if inrange(year, 1984, 2019)
keep if SEX == 2  /* women — cleaner first-birth timing */

gen n_children = RELCHIREP
gen first_birth_year = RELCHI1BYEAR
gen had_first = !missing(first_birth_year) & first_birth_year <= year
gen childless = (n_children == 0) if !missing(n_children)

* First-birth event dummy: this year is the first-birth year
gen first_birth_event = !missing(first_birth_year) & first_birth_year == year

* Risk set: childless at start of period (haven't had first birth yet)
gen pre_first_birth = missing(first_birth_year) | year < first_birth_year

gen liq_nw_to_inc = NETWORTH2R / INCFAMR if INCFAMR > 1000 & !missing(NETWORTH2R, INCFAMR)

* Panel-consistent lagged wealth using prior observation
sort ID year
by ID: gen liq_nw_lag = liq_nw_to_inc[_n-1] if ID == ID[_n-1]

* Keep risk set only: childless at start of year
keep if pre_first_birth == 1
keep if !missing(liq_nw_lag)

* Trim extreme tails (persistent outliers in ratios)
_pctile liq_nw_lag [pw=${weight}], p(1 99)
drop if liq_nw_lag < r(r1) | liq_nw_lag > r(r2)

xtile wealth_tercile = liq_nw_lag [pw=${weight}], n(3)

* Also compute a quintile version
xtile wealth_quintile = liq_nw_lag [pw=${weight}], n(5)

* Hazard by age * wealth tercile
preserve
    collapse (mean) first_birth_rate = first_birth_event [pweight=${weight}], by(AGEREP wealth_tercile)
    export delimited using "`outdir'/first_birth_hazard_by_age_tercile.csv", replace

    twoway (connected first_birth_rate AGEREP if wealth_tercile == 1, mcolor(cranberry) lcolor(cranberry) lwidth(medthick)) ///
           (connected first_birth_rate AGEREP if wealth_tercile == 2, mcolor(dkgreen) lcolor(dkgreen) lwidth(medthick)) ///
           (connected first_birth_rate AGEREP if wealth_tercile == 3, mcolor(navy) lcolor(navy) lwidth(medthick)) ///
        , legend(order(1 "Bottom wealth tercile" 2 "Middle wealth tercile" 3 "Top wealth tercile") position(6) rows(1)) ///
          xtitle("Age") ytitle("Annual P(first birth | childless)") ///
          title("First-birth hazard by lagged wealth") ///
          graphregion(color(white)) ///
          scheme(s1color)
    graph export "`outdir'/first_birth_hazard_by_age_tercile.png", replace width(1600)
restore

* Summary across ages, by tercile
preserve
    collapse (mean) first_birth_rate = first_birth_event ///
        (p50) median_lag_wealth = liq_nw_lag ///
        (count) n_obs = ID [pweight=${weight}], by(wealth_tercile)
    list
    export delimited using "`outdir'/first_birth_hazard_by_tercile.csv", replace
restore

* Also by quintile (for decile-ish picture)
preserve
    collapse (mean) first_birth_rate = first_birth_event ///
        (p50) median_lag_wealth = liq_nw_lag ///
        (count) n_obs = ID [pweight=${weight}], by(wealth_quintile)
    list
    export delimited using "`outdir'/first_birth_hazard_by_quintile.csv", replace
restore

log close
