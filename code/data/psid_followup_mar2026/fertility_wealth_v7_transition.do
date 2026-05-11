clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v7.log", replace text

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
* Sample: renters at 25-30, pre-birth, with wealth observed.
* ------------------------------------------------------------------
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if pre_first_birth == 1
    keep if !missing(own)
    collapse (mean) own_rate_2530 = own ///
        (mean) total_nw_2530 = total_nw_to_inc ///
        (mean) weight_young = IW, by(ID)
    * Classify as renter at 25-30 if own_rate < 0.5 (mostly renting)
    gen renter_2530 = (own_rate_2530 < 0.5)
    keep if renter_2530 == 1
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
* Fertility outcome: age at first birth, childless by 45
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

gen childless_by_45 = (missing(age_at_first_birth) | age_at_first_birth > 45)

* Trim wealth tails
_pctile total_nw_2530 [pw=weight_young], p(1 99)
drop if total_nw_2530 < r(r1) | total_nw_2530 > r(r2)

* ------------------------------------------------------------------
* Headline: childless by 45, by transition status
* ------------------------------------------------------------------
preserve
    gen transition_label = cond(owner_by_35 == 1, ///
        "Renter 25-30 → Owner by 35", "Renter 25-30 → Still renter by 35")
    collapse (mean) childless_by_45 (p50) median_wealth = total_nw_2530 ///
        (count) n_obs = ID [pweight=weight_young], by(transition_label)
    list
    export delimited using "`outdir'/childless_by_tenure_transition_v7.csv", replace
restore

* Same split by wealth tercile at 25-30
xtile wealth_tercile_renter = total_nw_2530 [pw=weight_young], n(3)

preserve
    gen transition_label = cond(owner_by_35 == 1, "Became Owner", "Stayed Renter")
    collapse (mean) childless_by_45 ///
        (count) n_obs = ID [pweight=weight_young], ///
        by(transition_label wealth_tercile_renter)
    list
    export delimited using "`outdir'/childless_by_transition_wealth_v7.csv", replace

    twoway (connected childless_by_45 wealth_tercile_renter if transition_label=="Became Owner", ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(large) lwidth(medthick)) ///
           (connected childless_by_45 wealth_tercile_renter if transition_label=="Stayed Renter", ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(large) lwidth(medthick)) ///
        , legend(order(1 "Renter @ 25-30 → Owner by 35" 2 "Renter @ 25-30 → Still renter by 35") ///
                 position(6) rows(2)) ///
          xtitle("Wealth tercile at 25-30 (within renters)") ///
          ytitle("Share childless by age 45") ///
          xlabel(1 "Low" 2 "Mid" 3 "High") ///
          ylabel(, angle(h)) ///
          title("Clearing the ownership hurdle predicts having kids") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/childless_by_transition_wealth_v7.png", replace width(1600)
restore

* Plot 2: simple bar chart — just two bars by transition status
preserve
    collapse (mean) childless_by_45 ///
        (count) n_obs = ID [pweight=weight_young], by(owner_by_35)
    gen transition = "Became owner by 35" if owner_by_35 == 1
    replace transition = "Stayed renter by 35" if owner_by_35 == 0
    list

    graph bar (asis) childless_by_45, over(transition) ///
        blabel(bar, format(%4.1f) position(outside)) ///
        ytitle("Share childless by 45 (among renters at 25-30)") ///
        ylabel(, angle(h)) ///
        bar(1, fcolor(navy) lcolor(navy)) ///
        title("Who remains childless among initial renters?") ///
        graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/childless_by_transition_v7_bar.png", replace width(1600)
restore

log close
