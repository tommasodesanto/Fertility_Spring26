clear all
set more off
version 17.0

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/extract27.dta"
local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_redfin_big_homes_replica"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/extract_redfin_big_homes_replica.log", replace text

* ==================================================================
* Replicate Redfin chart: share of 3+ BEDROOM OWNER-OCCUPIED homes
* held by (generation x household composition) cells.
* Redfin's headline: "Boomer 1-2 adults own 28%; Millennial w/ kids own 16%".
*
* Generations (used in 2023; we use ACS 2022-2023 to align):
*   Silent:     birth 1928-1945
*   Boomer:     birth 1946-1964
*   Gen X:      birth 1965-1980
*   Millennial: birth 1981-1996
*   Gen Z:      birth 1997-2012
*
* Household composition:
*   "1-2 adults" = no own children present, household size <= 2
*   "Households w/ own children" = nchild > 0
*   "3+ people, no children" = no own children, household size >= 3
* ==================================================================

use year age perwt gq pernum ownershp rooms bedrooms serial nchild yngch ///
    using "`dta'", clear

di _n "=== Raw rows ==="
count

* Use ACS 2022-2023 for generational labels to align with Redfin
keep if year >= 2022 & year <= 2023
keep if inlist(gq, 1, 2)

* Compute household size from serial + pernum (number of records per HH)
bysort year serial: gen numprec = _N

keep if pernum == 1  // heads only
keep if ownershp == 1  // owner-occupied only

* Bedrooms in IPUMS coded: 0 = N/A; 1 = "no bedrooms" (efficiency);
* 2 = 1 bedroom; 3 = 2; 4 = 3; ... 13 = 12+. Convert to count of bedrooms.
gen bedrooms_count = .
replace bedrooms_count = 0 if bedrooms == 1
replace bedrooms_count = bedrooms - 2 if bedrooms >= 2 & bedrooms <= 13
keep if !missing(bedrooms_count)
keep if bedrooms_count >= 3  // 3+ bedrooms only

di _n "=== After heads + owner + 3+ bedroom filter ==="
count

* Generation labels (birth year = year - age)
gen birth_year = year - age
gen str20 generation = ""
replace generation = "Silent"      if inrange(birth_year, 1928, 1945)
replace generation = "Boomer"      if inrange(birth_year, 1946, 1964)
replace generation = "Gen X"       if inrange(birth_year, 1965, 1980)
replace generation = "Millennial"  if inrange(birth_year, 1981, 1996)
replace generation = "Gen Z"       if inrange(birth_year, 1997, 2012)
drop if generation == ""

* Household composition. Redfin's "with own children" means MINOR (under-18) children.
* IPUMS nchild counts all own children including adults; restrict via yngch < 18.
* yngch = age of youngest own child; 99 = N/A.
gen byte has_minor_children = (nchild > 0 & yngch < 18 & yngch != 99) if !missing(nchild)
gen str30 hh_type = ""
replace hh_type = "1-2 adults"             if has_minor_children == 0 & numprec <= 2
replace hh_type = "3+ people, no children" if has_minor_children == 0 & numprec >= 3
replace hh_type = "with own children"      if has_minor_children == 1
drop if hh_type == ""

di _n "=== Final sample ==="
count

* ------------------------------------------------------------------
* Cell counts: generation x hh_type, owner-occupied 3+ bedrooms
* ------------------------------------------------------------------
preserve
    collapse (sum) cell_weight = perwt, by(generation hh_type)
    egen total_weight = total(cell_weight)
    gen share = cell_weight / total_weight
    gsort -share
    list, sepby(generation) noobs abbrev(25)
    export delimited using "`outdir'/redfin_big_homes_replica.csv", replace
restore

* Also: by individual generation (for header bar charts)
preserve
    collapse (sum) gen_weight = perwt, by(generation)
    egen total = total(gen_weight)
    gen gen_share = gen_weight / total
    list, noobs abbrev(20)
    export delimited using "`outdir'/redfin_big_homes_by_generation.csv", replace
restore

* And by hh_type only
preserve
    collapse (sum) typ_weight = perwt, by(hh_type)
    egen total = total(typ_weight)
    gen typ_share = typ_weight / total
    list, noobs abbrev(30)
    export delimited using "`outdir'/redfin_big_homes_by_hhtype.csv", replace
restore

di _n "=== Done. Output dir: `outdir' ==="

log close
