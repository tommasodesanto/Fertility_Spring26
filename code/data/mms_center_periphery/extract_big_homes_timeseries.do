clear all
set more off
version 17.0

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/extract27.dta"
local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_big_homes_timeseries"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/extract_big_homes_timeseries.log", replace text

* ==================================================================
* Time series of intergenerational allocation of large homes.
* Sample: heads of owner-occupied units with 3+ bedrooms, ACS 2012-2023.
*
* Output: cells by (year, age_bin, hh_type) with weighted shares,
* so we can plot:
*   (a) hh_type composition over time (1-2 adults vs with-kids vs 3+ no kids)
*   (b) within 1-2 adults, age composition over time
*   (c) the Redfin-style fixed-cohort (generation) composition over time
* ==================================================================

use year age perwt gq pernum ownershp bedrooms serial nchild yngch ///
    using "`dta'", clear

di _n "=== Raw rows ==="
count

keep if year >= 2012 & year <= 2023
keep if inlist(gq, 1, 2)
bysort year serial: gen numprec = _N
keep if pernum == 1
keep if ownershp == 1

* Bedrooms count from IPUMS coding
gen bedrooms_count = .
replace bedrooms_count = 0 if bedrooms == 1
replace bedrooms_count = bedrooms - 2 if bedrooms >= 2 & bedrooms <= 13
keep if !missing(bedrooms_count)
keep if bedrooms_count >= 3
keep if age >= 22 & age <= 95

di _n "=== After filter: heads, owner-occupied, 3+ bedroom, 2012-2023 ==="
count

* Generation (birth-year fixed)
gen birth_year = year - age
gen str20 generation = ""
replace generation = "Silent"      if inrange(birth_year, 1928, 1945)
replace generation = "Boomer"      if inrange(birth_year, 1946, 1964)
replace generation = "Gen X"       if inrange(birth_year, 1965, 1980)
replace generation = "Millennial"  if inrange(birth_year, 1981, 1996)
replace generation = "Gen Z"       if inrange(birth_year, 1997, 2012)
drop if generation == ""

* Household composition (minor children: yngch < 18)
gen byte has_minor_children = (nchild > 0 & yngch < 18 & yngch != 99) if !missing(nchild)
gen str30 hh_type = ""
replace hh_type = "1-2 adults"             if has_minor_children == 0 & numprec <= 2
replace hh_type = "3+ people, no children" if has_minor_children == 0 & numprec >= 3
replace hh_type = "with own children"      if has_minor_children == 1
drop if hh_type == ""

* Age bin
gen str10 age_bin = ""
replace age_bin = "22-39" if inrange(age, 22, 39)
replace age_bin = "40-59" if inrange(age, 40, 59)
replace age_bin = "60+"   if age >= 60
drop if age_bin == ""

* --------------------------------------------------------------
* (1) hh_type x year cells -> share within year
* --------------------------------------------------------------
preserve
    collapse (sum) cell_weight = perwt, by(year hh_type)
    bysort year: egen yr_total = total(cell_weight)
    gen share = cell_weight / yr_total
    sort year hh_type
    list, sepby(year) noobs abbrev(25)
    export delimited using "`outdir'/big_homes_hhtype_by_year.csv", replace
restore

* --------------------------------------------------------------
* (2) age_bin x year cells WITHIN 1-2 adults -> who are the empty nesters?
* --------------------------------------------------------------
preserve
    keep if hh_type == "1-2 adults"
    collapse (sum) cell_weight = perwt, by(year age_bin)
    bysort year: egen yr_total = total(cell_weight)
    gen share_within_1_2adults = cell_weight / yr_total
    sort year age_bin
    list, sepby(year) noobs abbrev(20)
    export delimited using "`outdir'/big_homes_emptynester_by_age_year.csv", replace
restore

* --------------------------------------------------------------
* (3) age_bin x hh_type x year for the full breakdown
* --------------------------------------------------------------
preserve
    collapse (sum) cell_weight = perwt, by(year age_bin hh_type)
    bysort year: egen yr_total = total(cell_weight)
    gen share = cell_weight / yr_total
    sort year hh_type age_bin
    export delimited using "`outdir'/big_homes_full_breakdown_by_year.csv", replace
restore

* --------------------------------------------------------------
* (4) Generation (fixed cohort) x hh_type x year
*     Tracks "the Redfin cells" across years
* --------------------------------------------------------------
preserve
    collapse (sum) cell_weight = perwt, by(year generation hh_type)
    bysort year: egen yr_total = total(cell_weight)
    gen share = cell_weight / yr_total
    sort year generation hh_type
    export delimited using "`outdir'/big_homes_generation_by_year.csv", replace
restore

di _n "=== Done ==="

log close
