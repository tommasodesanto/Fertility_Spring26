clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
cap mkdir "`outdir'"
capture log close _all
log using "`outdir'/fertility_wealth_v2.log", replace text

use ID year AGEREP SEX EDUYEAR DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR EARNINDR NETWORTH2R HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen liq_nw_to_inc = NETWORTH2R / INCFAMR if INCFAMR > 1000 & !missing(NETWORTH2R, INCFAMR)

* ------------------------------------------------------------------
* Young-entry wealth: average liquid NW/income at ages 25-30, BEFORE first birth
* ------------------------------------------------------------------
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if pre_first_birth == 1
    keep if !missing(liq_nw_to_inc)
    collapse (mean) liq_nw_to_inc_young = liq_nw_to_inc ///
        (mean) weight_young = IW ///
        (first) first_birth_year_p = first_birth_year, by(ID)
    di "Young-entry wealth sample:"
    count
    summ liq_nw_to_inc_young, detail
    tempfile wy
    save `wy', replace
restore

* ------------------------------------------------------------------
* Fertility outcome: did they have a first birth by age 35? 40?
* Computed from RELCHI1BYEAR — no need to observe them at 35+
* ------------------------------------------------------------------
preserve
    * Birth cohort from earliest observation
    collapse (min) min_year = year (min) min_age = AGEREP ///
        (first) first_birth_year_b = first_birth_year, by(ID)
    gen birth_year = min_year - min_age
    gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
    gen ever_parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)
    gen ever_parent_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
    gen ever_parent_any   = (!missing(first_birth_year_b))
    keep ID birth_year age_at_first_birth ever_parent_by_35 ever_parent_by_40 ever_parent_any
    tempfile fert
    save `fert', replace
restore

use `wy', clear
merge 1:1 ID using `fert', keep(match) nogen

di "Merged sample:"
count
summ liq_nw_to_inc_young age_at_first_birth ever_parent_by_35 ever_parent_by_40 ever_parent_any

* Keep those where we can actually observe outcome — i.e., person's birth cohort is old enough
* that by 2019 they reached 35 or 40 years of age.
gen observed_by_35 = (2019 - birth_year) >= 35
gen observed_by_40 = (2019 - birth_year) >= 40

* ------------------------------------------------------------------
* Trim wealth outliers
* ------------------------------------------------------------------
_pctile liq_nw_to_inc_young [pw=weight_young], p(1 99)
drop if liq_nw_to_inc_young < r(r1) | liq_nw_to_inc_young > r(r2)

xtile wealth_tercile = liq_nw_to_inc_young [pw=weight_young], n(3)
xtile wealth_decile  = liq_nw_to_inc_young [pw=weight_young], n(10)

* ------------------------------------------------------------------
* Parenthood by 35 (and 40) by wealth tercile
* ------------------------------------------------------------------
preserve
    keep if observed_by_35 == 1
    collapse (mean) parent_by_35 = ever_parent_by_35 ///
        (p50) median_age_fb = age_at_first_birth ///
        (p50) median_wealth = liq_nw_to_inc_young ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_tercile)
    list
    export delimited using "`outdir'/fertility_by_wealth_tercile_v2.csv", replace
restore

preserve
    keep if observed_by_40 == 1
    collapse (mean) parent_by_40 = ever_parent_by_40 ///
        (p50) median_age_fb = age_at_first_birth ///
        (p50) median_wealth = liq_nw_to_inc_young ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_tercile)
    list
    export delimited using "`outdir'/fertility_by_wealth_tercile40_v2.csv", replace
restore

* Decile version
preserve
    keep if observed_by_35 == 1
    collapse (mean) parent_by_35 = ever_parent_by_35 ///
        (p50) median_wealth = liq_nw_to_inc_young [pweight=weight_young], by(wealth_decile)
    list
    export delimited using "`outdir'/fertility_by_wealth_decile_v2.csv", replace

    twoway (connected parent_by_35 wealth_decile, mcolor(navy) lcolor(navy) ///
            msymbol(O) msize(medium) lwidth(medthick)) ///
        , xtitle("Liquid NW/income decile at ages 25-30 (pre-birth)") ///
          ytitle("Share who are parents by age 35") ///
          xlabel(1(1)10) ///
          title("Parenthood by 35 by entry wealth (PSID 1984-2019)") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/parent_by_35_wealth_decile_v2.png", replace width(1600)
restore

* ------------------------------------------------------------------
* Age at first birth by wealth tercile (among parents)
* ------------------------------------------------------------------
preserve
    keep if ever_parent_any == 1 & !missing(age_at_first_birth)
    collapse (mean) mean_age_fb = age_at_first_birth ///
        (p50) median_age_fb = age_at_first_birth ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_tercile)
    list
    export delimited using "`outdir'/age_fb_by_wealth_tercile_v2.csv", replace
restore

log close
