clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/dp_shortfall_with_controls_v1"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/dp_shortfall_with_controls_v1.log", replace text

* ==================================================================
* Wealth-gate regression with demographic controls
* Sample, shortfall definition, and target-price construction match
* dp_shortfall_v1.do; this version layers in education, race,
* marital status, employment, occupation broad category, year FE.
* ==================================================================

* Target price by year-bin (replicate from dp_shortfall_v1)
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

* Renter entry sample with demographic controls
use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN IW ///
    EDUYEAR RACE FAMMARRIED EMPWORK OCC2010C1M ///
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

gen byte sh_035 = liquid_nw < 0.035 * p_target if !missing(liquid_nw) & !missing(p_target)
gen byte sh_200 = liquid_nw < 0.20  * p_target if !missing(liquid_nw) & !missing(p_target)

* Broad occupation: bucket OCC2010 codes (4-digit 2010 SOC); simple major-group bin
gen occ_major = .
replace occ_major = 1 if inrange(OCC2010C1M, 10, 950)      // Mgmt/Prof
replace occ_major = 2 if inrange(OCC2010C1M, 1000, 3550)   // STEM/Healthcare
replace occ_major = 3 if inrange(OCC2010C1M, 3600, 4650)   // Service
replace occ_major = 4 if inrange(OCC2010C1M, 4700, 5940)   // Sales/Office
replace occ_major = 5 if inrange(OCC2010C1M, 6000, 7630)   // Production/Construction
replace occ_major = 6 if inrange(OCC2010C1M, 7700, 9750)   // Transport/Material moving
replace occ_major = 0 if missing(OCC2010C1M) | OCC2010C1M == 0 | OCC2010C1M >= 9800  // unemp/N/A
label define occ_major 0 "None/NA" 1 "Mgmt" 2 "STEM" 3 "Service" 4 "Sales/Office" 5 "Production" 6 "Transport"
label values occ_major occ_major

* Collapse to one obs per ID with means over the 25-30 window
collapse (mean) liquid_nw_mean = liquid_nw ///
         (mean) inc_mean = INCFAMR ///
         (mean) sh_035_m = sh_035 ///
         (mean) sh_200_m = sh_200 ///
         (mean) eduyear_mean = EDUYEAR ///
         (mean) fammarried_mean = FAMMARRIED ///
         (mean) empwork_mean = EMPWORK ///
         (mean) race_mode = RACE ///
         (mean) occ_mode = occ_major ///
         (mean) sex_mode = SEX ///
         (mean) year_mean = year ///
         (mean) weight_young = IW, by(ID)
gen byte sh_035 = sh_035_m > 0.5 if !missing(sh_035_m)
gen byte sh_200 = sh_200_m > 0.5 if !missing(sh_200_m)
tempfile entry
save `entry', replace

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

keep if observed_by_45 == 1
keep if !missing(sh_035) & !missing(inc_mean) & inc_mean > 0

gen log_inc = log(inc_mean)
gen log_nw  = log(max(liquid_nw_mean, 1))

* Bucket year_mean for year FE
gen year_fe = round(year_mean)
* Race buckets (modal value); convert to discrete bins for indicator FE
gen race_bin = round(race_mode)
* Occ buckets
gen occ_bin = round(occ_mode)

di _n "=== Sample size ==="
count

* ------------------------------------------------------------------
* Progressive specification: just shortfall; add controls one set
* at a time so we see how the coefficient moves.
* ------------------------------------------------------------------
di _n "----- (1) Bare shortfall -----"
reg childless_by_45 sh_035 [pweight=weight_young], robust

di _n "----- (2) + log income + log NW -----"
reg childless_by_45 sh_035 log_inc log_nw [pweight=weight_young], robust

di _n "----- (3) + education + sex -----"
reg childless_by_45 sh_035 log_inc log_nw eduyear_mean sex_mode [pweight=weight_young], robust

di _n "----- (4) + marital status + employment -----"
reg childless_by_45 sh_035 log_inc log_nw eduyear_mean sex_mode ///
    fammarried_mean empwork_mean [pweight=weight_young], robust

di _n "----- (5) + race FE + occupation FE -----"
reg childless_by_45 sh_035 log_inc log_nw eduyear_mean sex_mode ///
    fammarried_mean empwork_mean i.race_bin i.occ_bin [pweight=weight_young], robust

di _n "----- (6) full: + year FE -----"
reg childless_by_45 sh_035 log_inc log_nw eduyear_mean sex_mode ///
    fammarried_mean empwork_mean i.race_bin i.occ_bin i.year_fe [pweight=weight_young], robust

* Also: at 20% threshold
di _n "----- (7) 20% shortfall + full controls -----"
reg childless_by_45 sh_200 log_inc log_nw eduyear_mean sex_mode ///
    fammarried_mean empwork_mean i.race_bin i.occ_bin i.year_fe [pweight=weight_young], robust

log close
