clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v4.log", replace text

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR EARNINDR NETWORTHR NETWORTH2R NETWORTH3R HOMEEQUITYR HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)

* Two wealth measures
gen total_nw_to_inc = NETWORTHR  / INCFAMR if INCFAMR > 1000 & !missing(NETWORTHR, INCFAMR)
gen liq_nw_to_inc   = NETWORTH2R / INCFAMR if INCFAMR > 1000 & !missing(NETWORTH2R, INCFAMR)
gen own = HOMEOWN

* Young-entry wealth (pre-birth, 25-30)
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if pre_first_birth == 1
    keep if !missing(total_nw_to_inc, liq_nw_to_inc)
    collapse (mean) total_nw_young = total_nw_to_inc ///
        (mean) liq_nw_young = liq_nw_to_inc ///
        (mean) own_rate_young = own ///
        (mean) weight_young = IW, by(ID)
    gen owner_young = (own_rate_young >= 0.5) if !missing(own_rate_young)
    tempfile wy
    save `wy', replace
restore

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

* Trim
_pctile total_nw_young [pw=weight_young], p(1 99)
drop if total_nw_young < r(r1) | total_nw_young > r(r2)

xtile wealth_tercile = total_nw_young [pw=weight_young], n(3)

* Observable flag
gen observed_by_35 = (2019 - birth_year) >= 35
keep if observed_by_35 == 1

gen parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)

* (A) Overall: parenthood by 35 by TOTAL NW tercile
preserve
    collapse (mean) parent_by_35 (mean) own_rate_young ///
        (p50) median_total_nw = total_nw_young (p50) median_liq_nw = liq_nw_young ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_tercile)
    list
    export delimited using "`outdir'/fertility_by_total_wealth_tercile_v4.csv", replace
restore

* (B) By tenure × total-NW tercile
preserve
    keep if !missing(owner_young)
    gen tenure_label = cond(owner_young == 1, "Owner", "Renter")
    collapse (mean) parent_by_35 (count) n_obs = ID [pweight=weight_young], ///
        by(tenure_label wealth_tercile)
    list
    export delimited using "`outdir'/fertility_by_tenure_totalwealth_v4.csv", replace

    twoway (connected parent_by_35 wealth_tercile if tenure_label=="Owner", ///
            mcolor(navy) lcolor(navy) msymbol(O) lwidth(medthick)) ///
           (connected parent_by_35 wealth_tercile if tenure_label=="Renter", ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) lwidth(medthick)) ///
        , legend(order(1 "Owner at 25-30" 2 "Renter at 25-30") position(6) rows(1)) ///
          xtitle("Total net worth / income tercile at 25-30") ///
          ytitle("Share parents by age 35") ///
          xlabel(1 "Low" 2 "Mid" 3 "High") ///
          title("Parenthood by 35: by tenure and TOTAL wealth (incl. housing)") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/fertility_by_tenure_totalwealth_v4.png", replace width(1600)
restore

log close
