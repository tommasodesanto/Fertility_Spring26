clear all
set more off
version 17.0

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/extract27.dta"
local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_heads_only_age_rooms"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/extract_heads_only_age_rooms.log", replace text

* ==================================================================
* HEADS-ONLY ACS AGE PROFILES + ROOMS x AGE CELLS
*
* Restrict to household heads (pernum == 1) so the model-vs-data
* comparison is apples-to-apples: model agents are household heads.
* ==================================================================

* Load only needed columns. Use column-select to avoid loading 9.2GB.
use year statefip age perwt gq pernum ownershp rooms ///
    using "`dta'", clear

di _n "=== Raw rows ==="
count

* Sample filter
keep if year >= 2012
keep if age >= 22 & age <= 85
keep if inlist(gq, 1, 2)
keep if pernum == 1  // HEAD ONLY

* Recode tenure
gen byte owner = ownershp == 1 if !missing(ownershp)
gen byte renter = ownershp == 2 if !missing(ownershp)
keep if owner == 1 | renter == 1
keep if !missing(rooms) & rooms > 0 & rooms < 99

di _n "=== After heads-only + age filter ==="
count

* ------------------------------------------------------------------
* (1) Age profile: owner_rate and mean rooms by integer age
* ------------------------------------------------------------------
preserve
    collapse (mean) owner_rate = owner mean_rooms = rooms ///
             (sum) pop_weight = perwt [aweight=perwt], by(age)
    export delimited using "`outdir'/acs_heads_age_profile.csv", replace
    di _n "Saved acs_heads_age_profile.csv"
restore

* ------------------------------------------------------------------
* (2) Rooms x age cells: for the inverted view (rooms on x, age stack)
* Bin rooms into {<=4, 5, 6, 7-8, 9-10, 11+} matching housing_size_fit
* Bin ages into 10-year groups
* ------------------------------------------------------------------
gen rooms_bin = .
replace rooms_bin = 1 if rooms <= 4
replace rooms_bin = 2 if rooms == 5
replace rooms_bin = 3 if rooms == 6
replace rooms_bin = 4 if rooms >= 7 & rooms <= 8
replace rooms_bin = 5 if rooms >= 9 & rooms <= 10
replace rooms_bin = 6 if rooms >= 11
label define rb 1 "<=4" 2 "5" 3 "6" 4 "7-8" 5 "9-10" 6 "11+"
label values rooms_bin rb

gen age_bin = .
replace age_bin = 1 if inrange(age, 22, 29)
replace age_bin = 2 if inrange(age, 30, 39)
replace age_bin = 3 if inrange(age, 40, 49)
replace age_bin = 4 if inrange(age, 50, 59)
replace age_bin = 5 if inrange(age, 60, 69)
replace age_bin = 6 if inrange(age, 70, 85)
label define ab 1 "22-29" 2 "30-39" 3 "40-49" 4 "50-59" 5 "60-69" 6 "70+"
label values age_bin ab

* All-tenure cells
preserve
    collapse (sum) cell_weight = perwt, by(rooms_bin age_bin)
    decode rooms_bin, gen(rooms_label)
    decode age_bin, gen(age_label)
    export delimited using "`outdir'/acs_heads_rooms_x_age_all.csv", replace
restore

* By tenure
preserve
    gen tenure_label = cond(owner == 1, "owner", "renter")
    collapse (sum) cell_weight = perwt, by(tenure_label rooms_bin age_bin)
    decode rooms_bin, gen(rooms_label)
    decode age_bin, gen(age_label)
    export delimited using "`outdir'/acs_heads_rooms_x_age_by_tenure.csv", replace
restore

di _n "=== Done. Output dir: `outdir' ==="

log close
