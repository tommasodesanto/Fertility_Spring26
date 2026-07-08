clear all
set more off
version 17.0

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/extract27.dta"
local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_misallocation_table"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/extract_misallocation_table.log", replace text

* ==================================================================
* ACS household-head family-housing misallocation table
*
* Template matched:
*   code/data/mms_center_periphery/extract_heads_only_age_rooms.do
*   Source extract: extract27.dta
*   Filters: ACS years >= 2012, group quarters in {1,2}, household
*            heads only (pernum == 1), owners/renters only, valid rooms.
*
* Age support:
*   The requested table needs bands {18-29, 30-45, 46-64, 65+}.
*   Relative to the template's age-profile script, the lower age cutoff is
*   extended from 22 to 18 so the first requested band is literal. The
*   template's upper cap age <= 85 is retained for consistency with that
*   heads-only profile.
*
* Weight convention:
*   The script checks whether hhwt is available in extract27.dta. If hhwt
*   exists, weighted counts are household-weighted housing-unit counts.
*   If hhwt is unavailable, weighted counts fall back to perwt after
*   restricting to household heads, matching the template convention.
*
* MMS note:
*   The template .do/.log apply no explicit MMS lookup, metro, CBSA, PUMA,
*   or matched-metro filter. This script deliberately mirrors that template
*   rather than adding a new geography merge.
* ==================================================================

quietly describe using "`dta'", simple
local vars "`r(varlist)'"
local has_hhwt = strpos(" `vars' ", " hhwt ") > 0

local weight_var "perwt"
local weight_note "perwt restricted to household heads"
local load_weight_vars "perwt"
if `has_hhwt' {
    local weight_var "hhwt"
    local weight_note "hhwt"
    local load_weight_vars "perwt hhwt"
}

di _n "=== Source extract ==="
di "`dta'"
di "Weight variable selected: `weight_var' (`weight_note')"

use year statefip age `load_weight_vars' gq pernum ownershp rooms nchild ///
    using "`dta'", clear

di _n "=== Raw rows ==="
count

* Sample filter, mirroring extract_heads_only_age_rooms.do except for the
* requested 18-29 lower-band support.
keep if year >= 2012
keep if age >= 18 & age <= 85
keep if inlist(gq, 1, 2)
keep if pernum == 1  // HEAD ONLY

* Recode tenure and keep valid housing observations.
gen byte owner = ownershp == 1 if !missing(ownershp)
gen byte renter = ownershp == 2 if !missing(ownershp)
keep if owner == 1 | renter == 1
keep if !missing(rooms) & rooms > 0 & rooms < 99
keep if !missing(nchild)
keep if !missing(`weight_var') & `weight_var' > 0

di _n "=== After heads-only + requested table filters ==="
count
local final_n = r(N)

quietly summarize year, meanonly
local year_min = r(min)
local year_max = r(max)

di "ACS years in final sample: `year_min'-`year_max'"

* Requested table dimensions.
gen byte age_band_id = .
replace age_band_id = 1 if inrange(age, 18, 29)
replace age_band_id = 2 if inrange(age, 30, 45)
replace age_band_id = 3 if inrange(age, 46, 64)
replace age_band_id = 4 if age >= 65
label define ageband 1 "18-29" 2 "30-45" 3 "46-64" 4 "65+"
label values age_band_id ageband
drop if missing(age_band_id)

gen str6 tenure = cond(owner == 1, "owner", "renter")
gen byte rooms_ge_6 = rooms >= 6 if !missing(rooms)
gen byte children_present = nchild > 0 if !missing(nchild)

tempfile cells skeleton
preserve
    collapse (count) n_households = age ///
             (sum) weighted_housing_units = `weight_var', ///
             by(age_band_id tenure rooms_ge_6 children_present)
    decode age_band_id, gen(age_band)
    save `cells', replace
restore

* Build an explicit 4 x 2 x 2 x 2 skeleton so zero cells would be emitted.
preserve
    clear
    set obs 32
    gen byte age_band_id = .
    gen str6 tenure = ""
    gen byte rooms_ge_6 = .
    gen byte children_present = .
    label define ageband 1 "18-29" 2 "30-45" 3 "46-64" 4 "65+"

    local row = 0
    forvalues a = 1/4 {
        foreach t in owner renter {
            forvalues r = 0/1 {
                forvalues c = 0/1 {
                    local ++row
                    replace age_band_id = `a' in `row'
                    replace tenure = "`t'" in `row'
                    replace rooms_ge_6 = `r' in `row'
                    replace children_present = `c' in `row'
                }
            }
        }
    }
    label values age_band_id ageband
    decode age_band_id, gen(age_band)
    save `skeleton', replace
restore

use `skeleton', clear
merge 1:1 age_band_id tenure rooms_ge_6 children_present using `cells', nogen
replace n_households = 0 if missing(n_households)
replace weighted_housing_units = 0 if missing(weighted_housing_units)

gen str20 weight_var = "`weight_var'"
gen long source_observations = `final_n'
gen int acs_year_min = `year_min'
gen int acs_year_max = `year_max'
sort age_band_id tenure rooms_ge_6 children_present

order age_band_id age_band tenure rooms_ge_6 children_present ///
      weighted_housing_units n_households weight_var source_observations ///
      acs_year_min acs_year_max

di _n "=== Cell-count diagnostics ==="
summarize n_households, meanonly
di "Total source household-head observations in cross-tab: `final_n'"
di "Minimum unweighted cell count: " r(min)
di "Maximum unweighted cell count: " r(max)
di "No cells suppressed by this script."

export delimited using "`outdir'/acs_misallocation_cells.csv", replace nolabel
di _n "Saved `outdir'/acs_misallocation_cells.csv"

log close
