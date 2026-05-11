clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v8.log", replace text

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen total_nw_to_inc = NETWORTHR / INCFAMR if INCFAMR > 1000 & !missing(NETWORTHR, INCFAMR)
gen own = HOMEOWN

* ------------------------------------------------------------------
* Sample: renters at 25-30, pre-birth, with wealth observed
* ------------------------------------------------------------------
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if pre_first_birth == 1
    keep if !missing(own)
    collapse (mean) own_rate_2530 = own ///
        (mean) total_nw_2530 = total_nw_to_inc ///
        (mean) weight_young = IW, by(ID)
    gen renter_2530 = (own_rate_2530 < 0.5)
    keep if renter_2530 == 1
    keep if !missing(total_nw_2530)
    tempfile r2530
    save `r2530', replace
restore

* ------------------------------------------------------------------
* Tenure by age 31-35: did they transition to ownership?
* ------------------------------------------------------------------
preserve
    keep if inrange(AGEREP, 31, 35)
    keep if !missing(own)
    collapse (mean) own_rate_3135 = own (mean) weight_3135 = IW, by(ID)
    gen owner_by_35 = (own_rate_3135 >= 0.5) if !missing(own_rate_3135)
    keep ID owner_by_35 own_rate_3135
    tempfile t3135
    save `t3135', replace
restore

* ------------------------------------------------------------------
* Fertility outcome
* ------------------------------------------------------------------
preserve
    collapse (min) min_year = year (min) min_age = AGEREP ///
        (first) first_birth_year_b = first_birth_year, by(ID)
    gen birth_year = min_year - min_age
    gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
    keep ID birth_year age_at_first_birth
    tempfile fert
    save `fert', replace
restore

use `r2530', clear
merge 1:1 ID using `t3135', keep(match master) nogen
merge 1:1 ID using `fert', keep(match) nogen

gen observed_by_45 = (2019 - birth_year) >= 45
keep if observed_by_45 == 1
keep if !missing(owner_by_35)

gen parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)

* Trim wealth tails
_pctile total_nw_2530 [pw=weight_young], p(1 99)
drop if total_nw_2530 < r(r1) | total_nw_2530 > r(r2)

* ------------------------------------------------------------------
* Key analysis: wealth ventiles among initial renters
* Plot both ownership transition AND parenthood by wealth bin
* ------------------------------------------------------------------
xtile wealth_ventile = total_nw_2530 [pw=weight_young], n(20)

preserve
    collapse (mean) became_owner = owner_by_35 ///
        (mean) parent_by_45 ///
        (p50) median_wealth = total_nw_2530 ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_ventile)
    list
    export delimited using "`outdir'/margin_by_ventile_v8.csv", replace

    * Dual-axis plot: ownership transition + parenthood
    twoway (connected became_owner wealth_ventile, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(medlarge) lwidth(medthick) ///
            yaxis(1)) ///
           (connected parent_by_45 wealth_ventile, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(medlarge) lwidth(medthick) ///
            yaxis(1)) ///
        , legend(order(1 "Became owner by 35" 2 "Parent by 45") ///
                 position(6) rows(1)) ///
          xtitle("Net worth / income ventile at 25-30 (among renters)") ///
          ytitle("Share", axis(1)) ///
          xlabel(1(1)20) ylabel(0(0.1)1, angle(h)) ///
          title("The wealth-ownership-fertility margin") ///
          subtitle("Among renters at 25-30, by entry wealth") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/margin_ventile_v8.png", replace width(1600)
restore

* ------------------------------------------------------------------
* Decile version (more power per bin)
* ------------------------------------------------------------------
xtile wealth_decile = total_nw_2530 [pw=weight_young], n(10)

preserve
    collapse (mean) became_owner = owner_by_35 ///
        (mean) parent_by_45 ///
        (p50) median_wealth = total_nw_2530 ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_decile)
    list
    export delimited using "`outdir'/margin_by_decile_v8.csv", replace

    twoway (connected became_owner wealth_decile, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(large) lwidth(medthick)) ///
           (connected parent_by_45 wealth_decile, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(large) lwidth(medthick)) ///
        , legend(order(1 "Became owner by 35" 2 "Parent by 45") ///
                 position(6) rows(1)) ///
          xtitle("Net worth / income decile at 25-30 (among renters)") ///
          ytitle("Share") ///
          xlabel(1(1)10) ylabel(0(0.1)1, angle(h)) ///
          title("Wealth drives ownership AND fertility") ///
          subtitle("Among renters at 25-30, by entry wealth") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/margin_decile_v8.png", replace width(1600)
restore

* ------------------------------------------------------------------
* Quintile version (maximum power)
* ------------------------------------------------------------------
xtile wealth_quintile = total_nw_2530 [pw=weight_young], n(5)

preserve
    collapse (mean) became_owner = owner_by_35 ///
        (mean) parent_by_45 ///
        (p50) median_wealth = total_nw_2530 ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_quintile)
    list
    export delimited using "`outdir'/margin_by_quintile_v8.csv", replace

    twoway (connected became_owner wealth_quintile, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(vlarge) lwidth(thick)) ///
           (connected parent_by_45 wealth_quintile, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 "Became owner by 35" 2 "Parent by 45") ///
                 position(6) rows(1)) ///
          xtitle("Net worth / income quintile at 25-30 (among renters)") ///
          ytitle("Share") ///
          xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5") ///
          ylabel(0(0.1)1, angle(h)) ///
          title("Wealth drives ownership AND fertility") ///
          subtitle("Among renters at 25-30, by entry wealth") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/margin_quintile_v8.png", replace width(1600)
restore

* ------------------------------------------------------------------
* Correlation between the two outcomes at individual level
* ------------------------------------------------------------------
di _n "=== Individual-level correlation ==="
corr became_owner parent_by_45 [aw=weight_young] if !missing(owner_by_35)

di _n "=== Probit: parent_by_45 on wealth + became_owner ==="
probit parent_by_45 total_nw_2530 owner_by_35 [pw=weight_young], robust
margins, dydx(total_nw_2530 owner_by_35) post

di _n "=== How much of the wealth-fertility gradient works through ownership? ==="
di "Running mediation-style: (1) fertility on wealth, (2) fertility on wealth + ownership"

probit parent_by_45 total_nw_2530 [pw=weight_young], robust
est store no_own
margins, dydx(total_nw_2530) post
di "Wealth effect WITHOUT controlling for ownership:"
mat list e(b)

probit parent_by_45 total_nw_2530 owner_by_35 [pw=weight_young], robust
est store with_own
margins, dydx(total_nw_2530) post
di "Wealth effect WITH controlling for ownership:"
mat list e(b)

log close
