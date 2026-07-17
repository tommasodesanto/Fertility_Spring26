clear all
set more off
version 17.0

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/extract27.dta"
local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_big_homes_timeseries"
capture mkdir "`outdir'"

capture log close _all
log using "`outdir'/extract_excess_bedrooms_by_age.log", replace text

* ==================================================================
* Age-conditional space-held measure, to net out mechanical cohort
* aging from the "intergenerational allocation" time series.
*
* For owner-occupied heads, compute:
*   excess_bedrooms     = bedrooms - household size  (spare rooms)
*   bedrooms_per_person = bedrooms / household size
* by (year, 5-year age band). If a given age band holds more spare
* space over time, that is real intensification net of cohort flow.
*
* Also tab the FULL year range in extract27 up front, to answer
* whether the file supports a window before 2012.
* ==================================================================

use year age perwt gq pernum ownershp bedrooms serial nchild yngch ///
    using "`dta'", clear

di _n "=== Full year range available in extract27 (before any filter) ==="
tab year

keep if year >= 2000 & year <= 2023
keep if inlist(gq, 1, 2)

* Household size = person records within (year, serial)
bysort year serial: gen hhsize = _N

keep if pernum == 1
keep if ownershp == 1

* Bedrooms count from IPUMS coding (1 = no bedrooms; >=2 -> code-2)
gen bedrooms_count = .
replace bedrooms_count = 0 if bedrooms == 1
replace bedrooms_count = bedrooms - 2 if bedrooms >= 2 & bedrooms <= 13
keep if !missing(bedrooms_count)

keep if age >= 25 & age <= 89

gen excess_bedrooms     = bedrooms_count - hhsize
gen bedrooms_per_person = bedrooms_count / hhsize

* 5-year age bands keyed by lower edge (25,30,...,85)
gen age5 = 5 * floor(age / 5)

gen one = 1

di _n "=== Owner-occupied heads per year after filters ==="
tab year

* (1) Age x year cells: weighted-mean space measures + cell sizes
preserve
    collapse (mean) excess_bedrooms bedrooms_per_person bedrooms_count hhsize ///
             (rawsum) pop = perwt (count) ncell = one [pw = perwt], ///
             by(year age5)
    sort year age5
    export delimited using "`outdir'/excess_bedrooms_by_age_year.csv", replace
restore

* (2) Headline: overall and 60+ owner heads, weighted mean excess bedrooms by year
preserve
    gen byte old = age >= 60
    collapse (mean) excess_bedrooms bedrooms_per_person bedrooms_count hhsize ///
             (rawsum) pop = perwt (count) ncell = one [pw = perwt], ///
             by(year old)
    sort year old
    list, sepby(year) noobs abbrev(20)
    export delimited using "`outdir'/excess_bedrooms_overall_by_year.csv", replace
restore

di _n "=== Done ==="
log close
