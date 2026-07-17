clear all
set more off
version 17.0

* ==================================================================
* Future-parent down-payment moments among young childless renters
*
* Object:
*   For PSID respondents observed as renters and still childless at ages
*   25-30, measure liquid wealth and down-payment distance before a first
*   birth, then split by eventual parenthood horizons.
*
* Main outputs:
*   output/future_parent_dp_moments_v1/future_parent_dp_entry_sample_v1.csv
*   output/future_parent_dp_moments_v1/future_parent_dp_moments_v1.csv
*   output/future_parent_dp_moments_v1/future_parent_dp_shortfall_bins_v1.csv
* ==================================================================

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local repo "`root'/Fertility_Spring26"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`repo'/code/data/psid_followup_mar2026/output/future_parent_dp_moments_v1"

cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/future_parent_dp_moments_v1.log", replace text

global weight IW

* ------------------------------------------------------------------
* STEP 1: target first-time-buyer home price by broad year bin
* ------------------------------------------------------------------
use ID year AGEREP HOMEOWN HOMEVALUER YEARMOVED_ ${weight} using "`dta'", clear

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

collapse (p50) p_target = HOMEVALUER ///
         (count) n_buyers = ID [aweight=${weight}], by(year_bin)

sum p_target if !missing(p_target), detail
local fill_price = r(p50)

tempfile target
save `target', replace

* ------------------------------------------------------------------
* STEP 2: childless-renter entry observations at ages 25-30
* ------------------------------------------------------------------
use ID year AGEREP SEX DEATHYEAR RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN ${weight} using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

gen byte renter = HOMEOWN == 2
gen byte childless_year = missing(RELCHI1BYEAR) | RELCHI1BYEAR > year

keep if inrange(AGEREP, 25, 30)
keep if renter == 1
keep if childless_year == 1

gen liquid_nw = cond(!missing(NETWORTH2R), NETWORTH2R, NETWORTHR)
keep if !missing(liquid_nw) & !missing(INCFAMR)
keep if !missing(${weight}) & ${weight} > 0

gen year_bin = .
replace year_bin = 1 if inrange(year, 1984, 1994)
replace year_bin = 2 if inrange(year, 1995, 2004)
replace year_bin = 3 if inrange(year, 2005, 2010)
replace year_bin = 4 if inrange(year, 2011, 2019)

merge m:1 year_bin using `target', keep(master match) nogen
replace p_target = `fill_price' if missing(p_target) | p_target <= 0

gen req_dp20 = 0.20 * p_target
gen req_dp035 = 0.035 * p_target

collapse (count) n_entry_years = year ///
         (mean) age_entry = AGEREP ///
         (mean) sex_entry = SEX ///
         (mean) liquid_nw_mean = liquid_nw ///
         (mean) income_mean = INCFAMR ///
         (mean) p_target_mean = p_target ///
         (mean) req_dp20_mean = req_dp20 ///
         (mean) req_dp035_mean = req_dp035 ///
         (mean) weight_entry = ${weight}, by(ID)

tempfile entry
save `entry', replace

* ------------------------------------------------------------------
* STEP 3: eventual fertility timing and observation horizons
* ------------------------------------------------------------------
use ID year AGEREP RELCHI1BYEAR using "`dta'", clear
drop if missing(AGEREP)
keep if inrange(year, 1984, 2019)

collapse (min) min_year = year (min) min_age = AGEREP ///
         (max) first_birth_year = RELCHI1BYEAR, by(ID)

gen birth_year = min_year - min_age
gen age_at_first_birth = first_birth_year - birth_year if !missing(first_birth_year)

foreach a in 35 40 45 {
    gen byte observed_by_`a' = (2019 - birth_year) >= `a'
    gen byte parent_by_`a' = (!missing(age_at_first_birth) & age_at_first_birth <= `a')
    gen byte childless_by_`a' = 1 - parent_by_`a'
}

keep ID birth_year age_at_first_birth observed_by_* parent_by_* childless_by_*

tempfile fert
save `fert', replace

use `entry', clear
merge 1:1 ID using `fert', keep(match) nogen

order ID n_entry_years weight_entry age_entry sex_entry liquid_nw_mean income_mean ///
    p_target_mean req_dp20_mean req_dp035_mean birth_year age_at_first_birth ///
    observed_by_35 parent_by_35 childless_by_35 ///
    observed_by_40 parent_by_40 childless_by_40 ///
    observed_by_45 parent_by_45 childless_by_45

export delimited using "`outdir'/future_parent_dp_entry_sample_v1.csv", replace

local py "/Users/tommasodesanto/miniconda3/bin/python"
local builder "`repo'/code/data/psid_followup_mar2026/build_future_parent_dp_moments_v1.py"

di as text "Building future-parent down-payment moment CSVs..."
shell "`py'" "`builder'"

capture confirm file "`outdir'/future_parent_dp_moments_v1.csv"
if _rc {
    di as error "Python builder did not produce future_parent_dp_moments_v1.csv."
    exit 601
}

di as result "Future-parent down-payment diagnostic complete."
di as result "Entry sample: `outdir'/future_parent_dp_entry_sample_v1.csv"
di as result "Moments:      `outdir'/future_parent_dp_moments_v1.csv"
di as result "Bins:         `outdir'/future_parent_dp_shortfall_bins_v1.csv"

log close _all
