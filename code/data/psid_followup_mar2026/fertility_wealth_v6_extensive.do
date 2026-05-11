clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v6.log", replace text

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

* Young-entry pre-birth wealth and tenure
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if pre_first_birth == 1
    keep if !missing(total_nw_to_inc)
    collapse (mean) total_nw_young = total_nw_to_inc ///
        (mean) own_rate_young = own ///
        (mean) weight_young = IW, by(ID)
    gen owner_young = (own_rate_young >= 0.5) if !missing(own_rate_young)
    tempfile wy
    save `wy', replace
restore

* Outcome: age at first birth (from RELCHI1BYEAR)
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

* Only keep those observable up to age 45 (cohort old enough by 2019)
gen observed_by_45 = (2019 - birth_year) >= 45
keep if observed_by_45 == 1

* Outcome 1: childless by 45 (i.e., age_at_first_birth > 45 OR missing)
gen childless_by_45 = (missing(age_at_first_birth) | age_at_first_birth > 45)

* Outcome 2: ever-parent (by 45) — complement
gen parent_by_45 = 1 - childless_by_45

* --------------------------
* By wealth tercile
* --------------------------
preserve
    collapse (mean) childless_by_45 parent_by_45 ///
        (mean) own_rate_young ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_tercile)
    list
    export delimited using "`outdir'/childless_by_wealth_tercile_v6.csv", replace
restore

* --------------------------
* By tenure × wealth tercile
* --------------------------
preserve
    keep if !missing(owner_young)
    gen tenure_label = cond(owner_young == 1, "Owner", "Renter")
    collapse (mean) childless_by_45 ///
        (count) n_obs = ID [pweight=weight_young], by(tenure_label wealth_tercile)
    list
    export delimited using "`outdir'/childless_by_tenure_wealth_v6.csv", replace

    twoway (connected childless_by_45 wealth_tercile if tenure_label=="Owner", ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(medlarge) lwidth(medthick)) ///
           (connected childless_by_45 wealth_tercile if tenure_label=="Renter", ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(medlarge) lwidth(medthick)) ///
        , legend(order(1 "Owner at 25-30" 2 "Renter at 25-30") position(6) rows(1)) ///
          xtitle("Total net worth / income tercile at 25-30") ///
          ytitle("Share childless by age 45") ///
          xlabel(1 "Low" 2 "Mid" 3 "High") ///
          ylabel(, angle(h)) ///
          title("Childlessness by 45: by tenure and entry wealth") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/childless_by_tenure_wealth_v6.png", replace width(1600)
restore

* --------------------------
* Also: decile version, total wealth only
* --------------------------
xtile wealth_decile = total_nw_young [pw=weight_young], n(10)
preserve
    collapse (mean) childless_by_45 (p50) median_wealth = total_nw_young ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_decile)
    list
    export delimited using "`outdir'/childless_by_wealth_decile_v6.csv", replace

    twoway (connected childless_by_45 wealth_decile, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(medlarge) lwidth(medthick)) ///
        , xtitle("Total net worth / income decile at ages 25-30") ///
          ytitle("Share childless by age 45") ///
          xlabel(1(1)10) ylabel(, angle(h)) ///
          title("Childlessness by 45 by total entry wealth") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/childless_by_wealth_decile_v6.png", replace width(1600)
restore

log close
