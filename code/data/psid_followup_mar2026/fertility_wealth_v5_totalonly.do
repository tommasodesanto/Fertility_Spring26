clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v5.log", replace text

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen total_nw_to_inc = NETWORTHR / INCFAMR if INCFAMR > 1000 & !missing(NETWORTHR, INCFAMR)

* Young-entry total wealth, pre-birth
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if pre_first_birth == 1
    keep if !missing(total_nw_to_inc)
    collapse (mean) total_nw_young = total_nw_to_inc ///
        (mean) weight_young = IW, by(ID)
    tempfile wy
    save `wy', replace
restore

* First-birth outcome
preserve
    collapse (min) min_year = year (min) min_age = AGEREP ///
        (first) first_birth_year_b = first_birth_year, by(ID)
    gen birth_year = min_year - min_age
    gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
    keep ID birth_year age_at_first_birth
    tempfile fert
    save `fert', replace
restore

use `wy', clear
merge 1:1 ID using `fert', keep(match) nogen

_pctile total_nw_young [pw=weight_young], p(1 99)
drop if total_nw_young < r(r1) | total_nw_young > r(r2)

xtile wealth_tercile = total_nw_young [pw=weight_young], n(3)
xtile wealth_decile  = total_nw_young [pw=weight_young], n(10)

gen observed_by_35 = (2019 - birth_year) >= 35
keep if observed_by_35 == 1
gen parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)

preserve
    collapse (mean) parent_by_35 (p50) median_wealth = total_nw_young ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_tercile)
    list
    export delimited using "`outdir'/fertility_by_totalwealth_tercile_v5.csv", replace
restore

preserve
    collapse (mean) parent_by_35 (p50) median_wealth = total_nw_young ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_decile)
    list
    export delimited using "`outdir'/fertility_by_totalwealth_decile_v5.csv", replace

    twoway (connected parent_by_35 wealth_decile, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(medlarge) lwidth(medthick)) ///
        , xtitle("Total net worth / income decile at ages 25-30") ///
          ytitle("Share parents by age 35") ///
          xlabel(1(1)10) ylabel(,angle(h)) ///
          title("Parenthood by 35 by total entry wealth (PSID 1984-2019)") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/fertility_by_totalwealth_decile_v5.png", replace width(1600)
restore

log close
