clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local extract "`root'/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
local lookup2010 "`root'/code/data/mms_center_periphery/data/puma_mms_lookup_2010.csv"
local lookup2020 "`root'/code/data/mms_center_periphery/data/puma_mms_lookup_2020.csv"
local outdir "`root'/april_26_discrete_time/output/acs_renter_price_elasticity_v1"

cap mkdir "`outdir'"
capture log close _all
log using "`outdir'/acs_renter_price_elasticity_v1.log", replace text

tempfile lookup base cellrent

import delimited using "`lookup2010'", clear varnames(1) stringcols(_all)
gen lookup_period = "pre2022"
destring statefip puma cbsacode, replace force
keep lookup_period statefip puma cbsacode mms_location
tempfile pre
save `pre', replace

import delimited using "`lookup2020'", clear varnames(1) stringcols(_all)
gen lookup_period = "post2021"
destring statefip puma cbsacode, replace force
keep lookup_period statefip puma cbsacode mms_location
append using `pre'
rename cbsacode met2013
replace mms_location = "center" if mms_location == "middle"
keep if inlist(mms_location, "center", "periphery")
save `lookup', replace

use year statefip puma met2013 age sex race perwt gq ownershp rooms rent hhincome nchild ///
    using "`extract'", clear

keep if year >= 2012 & year <= 2023
keep if age >= 22 & age <= 45
keep if inlist(gq, 1, 2)

gen lookup_period = cond(year <= 2021, "pre2022", "post2021")
gen renter = ownershp == 2 if !missing(ownershp)
gen any_parent = nchild > 0 if !missing(nchild)
gen childless = nchild == 0 if !missing(nchild)
replace any_parent = 0 if childless == 1 & missing(any_parent)
gen ln_income = ln(hhincome) if hhincome > 1000 & !missing(hhincome)
gen ln_rooms = ln(rooms) if rooms > 0 & !missing(rooms)
gen rent_per_room = rent / rooms if renter == 1 & rent > 0 & rooms > 0

merge m:1 lookup_period statefip puma met2013 using `lookup', keep(match) nogen
drop if missing(mms_location)

quietly _pctile hhincome [aw=perwt] if hhincome > 1000 & !missing(hhincome), p(1 99)
local q1 = r(r1)
local q99 = r(r2)
keep if missing(hhincome) | (hhincome >= `q1' & hhincome <= `q99')

gen loc_center = mms_location == "center"
save `base', replace

use `base', clear
keep if renter == 1 & rent_per_room > 0 & !missing(met2013, year, mms_location)
collapse (mean) local_rpr=rent_per_room (count) cell_n=rent_per_room [aw=perwt], by(year met2013 mms_location)
keep if cell_n >= 50
save `cellrent', replace

use `base', clear
keep if renter == 1
merge m:1 year met2013 mms_location using `cellrent', keep(match) nogen
drop if missing(ln_rooms, ln_income, any_parent, local_rpr)
gen ln_local_rpr = ln(local_rpr) if local_rpr > 0

reg ln_rooms c.ln_local_rpr##i.any_parent c.ln_income i.sex i.race i.age i.year i.met2013 i.loc_center ///
    [aw=perwt], vce(robust)

lincom ln_local_rpr
local child = r(estimate)
local child_se = r(se)
local child_p = r(p)
lincom 1.any_parent#c.ln_local_rpr
local diff = r(estimate)
local diff_se = r(se)
local diff_p = r(p)
lincom ln_local_rpr + 1.any_parent#c.ln_local_rpr
local parent = r(estimate)
local parent_se = r(se)
local parent_p = r(p)
local n = e(N)

tempname fh
file open `fh' using "`outdir'/acs_renter_price_elasticity_summary_v1.csv", write replace
file write `fh' "sample,outcome,childless_price_slope,childless_price_slope_se,childless_price_slope_p,parent_diff,parent_diff_se,parent_diff_p,parent_price_slope,parent_price_slope_se,parent_price_slope_p,n_obs" _n
file write `fh' "acs_renters,ln_rooms_on_ln_local_rpr,`child',`child_se',`child_p',`diff',`diff_se',`diff_p',`parent',`parent_se',`parent_p',`n'" _n
file close `fh'

di as result "ACS renter price elasticity complete."
di as result "Summary CSV: `outdir'/acs_renter_price_elasticity_summary_v1.csv"

log close _all
