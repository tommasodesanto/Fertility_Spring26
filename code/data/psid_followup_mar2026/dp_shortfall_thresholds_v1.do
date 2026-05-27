clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/dp_shortfall_v1"
capture log close _all
log using "`outdir'/dp_shortfall_thresholds_v1.log", replace text

* ==================================================================
* DP-SHORTFALL ROBUSTNESS ACROSS THRESHOLDS
*
* The original Object B uses 20% as the threshold. But low-income
* buyers in reality can use FHA (3.5% min) or low-DP conventional
* (5%). So the relevant question is: does the shortfall-childless
* link survive at lower, more-realistic thresholds for low-income
* buyers?
*
* Also: split by income tercile and check whether low-income renters
* are even FHA-short (wealth < 3.5% of typical first-buyer home).
* ==================================================================

* Same target price build as dp_shortfall_v1.do
use ID year AGEREP HOMEOWN HOMEVALUER HOMEMORTOTR YEARMOVED_ IW ///
    using "`dta'", clear
gen byte own = .
replace own = 1 if HOMEOWN == 1
replace own = 0 if HOMEOWN == 2
keep if inrange(year, 1984, 2019)
keep if !missing(AGEREP) & inrange(AGEREP, 18, 65)
sort ID year
by ID: gen any_renter_before = sum(own == 0)
by ID: gen byte first_purchase = own == 1 & any_renter_before[_n-1] >= 1 & ///
    sum(own == 1) == 1
keep if first_purchase == 1
keep if inrange(YEARMOVED_, year - 2, year)
keep if HOMEVALUER > 0 & !missing(HOMEVALUER)
gen year_bin = .
replace year_bin = 1 if inrange(year, 1984, 1994)
replace year_bin = 2 if inrange(year, 1995, 2004)
replace year_bin = 3 if inrange(year, 2005, 2010)
replace year_bin = 4 if inrange(year, 2011, 2019)
collapse (p50) p_target = HOMEVALUER [aweight=IW], by(year_bin)
sum p_target if !missing(p_target)
local fill = r(p50)
replace p_target = `fill' if missing(p_target) | p_target == 0
tempfile target
save `target', replace

* Renter panel at 25-30
use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN IW ///
    using "`dta'", clear
drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)
gen byte renter = HOMEOWN == 2
keep if inrange(AGEREP, 25, 30)
keep if renter == 1
keep if !missing(NETWORTHR) & !missing(INCFAMR)
gen liquid_nw = cond(missing(NETWORTH2R), NETWORTHR, NETWORTH2R)
gen year_bin = .
replace year_bin = 1 if inrange(year, 1984, 1994)
replace year_bin = 2 if inrange(year, 1995, 2004)
replace year_bin = 3 if inrange(year, 2005, 2010)
replace year_bin = 4 if inrange(year, 2011, 2019)
merge m:1 year_bin using `target', keep(master match) nogen

* Three thresholds: 20% (GSE), 10% (conv. PMI), 3.5% (FHA)
gen byte sh_200 = liquid_nw < 0.20  * p_target if !missing(liquid_nw) & !missing(p_target)
gen byte sh_100 = liquid_nw < 0.10  * p_target if !missing(liquid_nw) & !missing(p_target)
gen byte sh_050 = liquid_nw < 0.05  * p_target if !missing(liquid_nw) & !missing(p_target)
gen byte sh_035 = liquid_nw < 0.035 * p_target if !missing(liquid_nw) & !missing(p_target)

* Collapse to one obs per ID
preserve
    collapse (mean) liquid_nw_mean = liquid_nw ///
             (mean) inc_mean = INCFAMR ///
             (mean) p_target_mean = p_target ///
             (mean) sh_200_m = sh_200 (mean) sh_100_m = sh_100 ///
             (mean) sh_050_m = sh_050 (mean) sh_035_m = sh_035 ///
             (mean) weight_young = IW, by(ID)
    foreach t in 200 100 050 035 {
        gen byte sh_`t' = sh_`t'_m > 0.5
    }
    tempfile entry
    save `entry', replace
restore

* Fertility outcomes
use ID year AGEREP RELCHI1BYEAR using "`dta'", clear
drop if missing(AGEREP)
collapse (min) min_year = year (min) min_age = AGEREP ///
    (first) first_birth_year_b = RELCHI1BYEAR, by(ID)
gen birth_year = min_year - min_age
gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
keep ID birth_year age_at_first_birth
tempfile fert
save `fert', replace

use `entry', clear
merge 1:1 ID using `fert', keep(match) nogen
gen byte observed_by_45 = (2019 - birth_year) >= 45
gen byte parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
gen byte childless_by_45 = 1 - parent_by_45

* ------------------------------------------------------------------
* HEADLINE: shortfall rates by income tercile, for each threshold
* ------------------------------------------------------------------
keep if observed_by_45 == 1
keep if !missing(inc_mean)
xtile inc_tercile = inc_mean [pweight=weight_young], n(3)
label define inct 1 "Bottom inc tercile" 2 "Middle" 3 "Top"
label values inc_tercile inct

di _n "=== Share of renters who are short at each threshold, by income tercile ==="
table inc_tercile [aweight=weight_young], stat(mean sh_200 sh_100 sh_050 sh_035) stat(frequency)

di _n "=== Childlessness gap at each threshold (within income tercile) ==="
foreach t in 200 100 050 035 {
    di _n "--- Threshold: `t' (basis points of price; 200=20%, 035=3.5%) ---"
    table inc_tercile sh_`t' [aweight=weight_young], stat(mean childless_by_45) stat(frequency)
}

* ------------------------------------------------------------------
* Regression: childless ~ shortfall at each threshold + log income
* ------------------------------------------------------------------
gen log_inc = log(inc_mean) if inc_mean > 0
gen log_nw  = log(max(liquid_nw_mean, 1))

foreach t in 200 100 050 035 {
    di _n "--- Childless by 45 ~ sh_`t' + log income + log NW ---"
    reg childless_by_45 sh_`t' log_inc log_nw [pweight=weight_young] if !missing(log_inc), robust
}

* ------------------------------------------------------------------
* Export shortfall rates by tercile-threshold for plotting
* ------------------------------------------------------------------
preserve
    collapse (mean) sh_200 sh_100 sh_050 sh_035 ///
             (count) n_obs = ID [aweight=weight_young], by(inc_tercile)
    list
    export delimited using "`outdir'/shortfall_rates_by_threshold_inc_tercile_v1.csv", replace
restore

* Childless gap by threshold-tercile
preserve
    keep if !missing(sh_200)
    foreach t in 200 100 050 035 {
        preserve
            collapse (mean) childless_by_45 (count) n_obs = ID ///
                [aweight=weight_young], by(inc_tercile sh_`t')
            gen threshold = `t'
            tempfile gap_`t'
            save `gap_`t'', replace
        restore
    }
    use `gap_200', clear
    foreach t in 100 050 035 {
        gen byte sh = sh_200 if "`t'" == "200"
        drop sh
        append using `gap_`t''
    }
restore

* Cleaner export: long-format table of childlessness by tercile x threshold x shortfall
clear
gen str10 threshold = ""
gen byte tercile = .
gen byte short = .
gen double childless = .
gen long n = .
local row = 0
foreach t in 200 100 050 035 {
    use `entry', clear
    merge 1:1 ID using `fert', keep(match) nogen
    keep if (2019 - birth_year) >= 45
    keep if !missing(inc_mean)
    gen byte parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
    gen byte childless_by_45 = 1 - parent_by_45
    xtile inc_tercile = inc_mean [pweight=IW], n(3)
    collapse (mean) childless = childless_by_45 (count) n = ID ///
        [aweight=IW], by(inc_tercile sh_`t')
    gen str10 threshold = "`t'"
    rename sh_`t' short
    rename inc_tercile tercile
    if `row' == 0 {
        save "`outdir'/childless_gap_by_threshold_tercile_v1.dta", replace
    }
    else {
        append using "`outdir'/childless_gap_by_threshold_tercile_v1.dta"
        save "`outdir'/childless_gap_by_threshold_tercile_v1.dta", replace
    }
    local ++row
}
use "`outdir'/childless_gap_by_threshold_tercile_v1.dta", clear
sort threshold tercile short
list
export delimited using "`outdir'/childless_gap_by_threshold_tercile_v1.csv", replace

log close
