clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
cap mkdir "`outdir'"
capture log close _all
log using "`outdir'/fertility_wealth_v3.log", replace text

use ID year AGEREP SEX EDUYEAR DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR EARNINDR NETWORTH2R HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen liq_nw_to_inc = NETWORTH2R / INCFAMR if INCFAMR > 1000 & !missing(NETWORTH2R, INCFAMR)
gen own = HOMEOWN

* ------------------------------------------------------------------
* Per-person wealth AND modal tenure at ages 25-30, pre-birth
* ------------------------------------------------------------------
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if pre_first_birth == 1
    keep if !missing(liq_nw_to_inc)
    collapse (mean) liq_nw_to_inc_young = liq_nw_to_inc ///
        (mean) own_rate_young = own ///
        (mean) weight_young = IW, by(ID)
    gen owner_young = (own_rate_young >= 0.5) if !missing(own_rate_young)
    tempfile wy
    save `wy', replace
restore

* Fertility outcome file: age at first birth
preserve
    collapse (min) min_year = year (min) min_age = AGEREP ///
        (first) first_birth_year_b = first_birth_year, by(ID)
    gen birth_year = min_year - min_age
    gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
    gen ever_parent = !missing(age_at_first_birth)
    keep ID birth_year age_at_first_birth ever_parent
    tempfile fert
    save `fert', replace
restore

use `wy', clear
merge 1:1 ID using `fert', keep(match) nogen

* Trim wealth outliers
_pctile liq_nw_to_inc_young [pw=weight_young], p(1 99)
drop if liq_nw_to_inc_young < r(r1) | liq_nw_to_inc_young > r(r2)

xtile wealth_tercile = liq_nw_to_inc_young [pw=weight_young], n(3)

* Observed-by-age flags (birth cohort old enough by 2019)
gen observed_by_40 = (2019 - birth_year) >= 40

* ------------------------------------------------------------------
* (1) By-tenure: parenthood by 35 and by 40, wealth tercile
* ------------------------------------------------------------------
preserve
    keep if !missing(owner_young)
    gen parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)
    gen parent_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
    gen tenure_label = cond(owner_young == 1, "Owner", "Renter")

    collapse (mean) parent_by_35 parent_by_40 ///
        (count) n_obs = ID [pweight=weight_young], by(tenure_label wealth_tercile)
    list
    export delimited using "`outdir'/fertility_by_tenure_wealth_v3.csv", replace

    twoway (connected parent_by_35 wealth_tercile if tenure_label=="Owner", ///
            mcolor(navy) lcolor(navy) msymbol(O) lwidth(medthick)) ///
           (connected parent_by_35 wealth_tercile if tenure_label=="Renter", ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) lwidth(medthick)) ///
        , legend(order(1 "Owner at 25-30" 2 "Renter at 25-30") position(6) rows(1)) ///
          xtitle("Liquid NW/income tercile at ages 25-30") ///
          ytitle("Share parents by age 35") ///
          xlabel(1 "Low" 2 "Mid" 3 "High") ///
          title("Parenthood by 35: by tenure and entry wealth") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/fertility_by_tenure_wealth_v3.png", replace width(1600)
restore

* ------------------------------------------------------------------
* (2) Cumulative parenthood by age, by wealth tercile
* ------------------------------------------------------------------
preserve
    * Build long panel: one row per (ID, age_grid) where age_grid=22..40
    tempfile core
    keep ID wealth_tercile liq_nw_to_inc_young weight_young birth_year age_at_first_birth ever_parent
    save `core', replace

    clear
    set obs 19
    gen age_grid = 21 + _n  /* 22..40 */
    tempfile ages
    save `ages', replace

    use `core', clear
    cross using `ages'
    gen parent_at_age = (!missing(age_at_first_birth) & age_at_first_birth <= age_grid)
    * Restrict to ages where the cohort is old enough to be observed
    gen observable = (2019 - birth_year) >= age_grid
    keep if observable == 1

    collapse (mean) parent_at_age = parent_at_age [pweight=weight_young], ///
        by(age_grid wealth_tercile)
    list, noobs
    export delimited using "`outdir'/cum_parenthood_by_age_tercile_v3.csv", replace

    twoway (connected parent_at_age age_grid if wealth_tercile==1, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(O) lwidth(medthick)) ///
           (connected parent_at_age age_grid if wealth_tercile==2, ///
            mcolor(dkgreen) lcolor(dkgreen) msymbol(D) lwidth(medthick)) ///
           (connected parent_at_age age_grid if wealth_tercile==3, ///
            mcolor(navy) lcolor(navy) msymbol(S) lwidth(medthick)) ///
        , legend(order(1 "Low wealth" 2 "Mid wealth" 3 "High wealth") position(6) rows(1)) ///
          xtitle("Age") ytitle("Share who are parents by age") ///
          xlabel(22(2)40) ylabel(0(0.1)0.7, angle(h)) ///
          title("Age profile of parenthood by entry wealth") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/cum_parenthood_by_age_tercile_v3.png", replace width(1600)
restore

log close
