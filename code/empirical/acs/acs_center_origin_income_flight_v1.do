clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local extract27 "`root'/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
local extract28 "`root'/code/data/Spatial_aggregate_withmicrodata/raw_data/extract28_origin_geo_2012_2023_age22_45_households.csv"
local lookup2010 "`root'/code/data/mms_center_periphery/data/puma_mms_lookup_2010.csv"
local lookup2020 "`root'/code/data/mms_center_periphery/data/puma_mms_lookup_2020.csv"
local originbridge "`root'/code/data/mms_center_periphery/data/migpuma_mms_origin_bridge.csv"
local outdir "`root'/april_26_discrete_time/output/acs_center_origin_income_flight_v1"

cap mkdir "`outdir'"
capture log close _all
log using "`outdir'/acs_center_origin_income_flight_v1.log", replace text

tempfile lookup origin_geo origin_bridge base within_origin tertiles

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

import delimited using "`extract28'", clear varnames(1) stringcols(_all)
destring year sample serial pernum migplac1 migmetro1, replace force
keep year sample serial pernum migplac1 migmetro1
save `origin_geo', replace

import delimited using "`originbridge'", clear varnames(1) stringcols(_all)
destring migplac1 migpuma1 cbsacode center_share periphery_share middle_share, replace force
rename cbsacode migmet131
rename period lookup_period
keep lookup_period migplac1 migpuma1 migmet131 center_share periphery_share middle_share
save `origin_bridge', replace

use year sample serial pernum statefip puma migpuma1 metro met2013 migmet131 ///
    age sex race perwt gq hhincome nchild eldch migrate1 using "`extract27'", clear

keep if year >= 2012 & year <= 2023
keep if age >= 22 & age <= 45
keep if inlist(gq, 1, 2)

gen lookup_period = cond(year <= 2021, "pre2022", "post2021")
gen has_children = nchild > 0 if !missing(nchild)
gen childless = nchild == 0 if !missing(nchild)
gen newparent = eldch < 4 & eldch != 99 & nchild > 0 if !missing(eldch, nchild)
replace newparent = 0 if childless == 1 & missing(newparent)
gen moved1y = migrate1 >= 2 if !missing(migrate1)
gen same_msa = moved1y == 1 & met2013 > 0 & migmet131 > 0 & met2013 == migmet131 if !missing(moved1y, met2013, migmet131)
gen ln_income = ln(hhincome) if hhincome > 1000 & !missing(hhincome)
gen parent_compare = ""
replace parent_compare = "Non-Parents" if childless == 1
replace parent_compare = "New Parents" if newparent == 1

merge 1:1 year sample serial pernum using `origin_geo', keep(match master) nogen
merge m:1 lookup_period statefip puma met2013 using `lookup', keep(match) nogen
merge m:1 lookup_period migplac1 migpuma1 migmet131 using `origin_bridge', ///
    keep(master match) nogen

replace center_share = 0 if missing(center_share)
replace periphery_share = 0 if missing(periphery_share)
replace middle_share = 0 if missing(middle_share)
replace center_share = center_share + middle_share
gen origin_cp_share = center_share + periphery_share
gen origin_center_weight = perwt * center_share

drop if missing(ln_income, parent_compare, same_msa, mms_location)

quietly _pctile hhincome [aw=perwt] if hhincome > 1000 & !missing(hhincome), p(1 99)
local q1 = r(r1)
local q99 = r(r2)
keep if missing(hhincome) | (hhincome >= `q1' & hhincome <= `q99')

gen periphery_dest = mms_location == "periphery"

keep if same_msa == 1
keep if origin_cp_share > 0
keep if origin_center_weight > 0

save `within_origin', replace

quietly _pctile ln_income [aw=origin_center_weight], p(25 75)
local p25 = r(r1)
local p75 = r(r2)

reg periphery_dest c.ln_income##i.newparent i.sex i.race i.age i.year [aw=origin_center_weight], vce(robust)

lincom 1.newparent + `p25' * 1.newparent#c.ln_income
local low = r(estimate)
local low_se = r(se)
local low_p = r(p)
lincom 1.newparent + `p75' * 1.newparent#c.ln_income
local high = r(estimate)
local high_se = r(se)
local high_p = r(p)
local inter = _b[1.newparent#c.ln_income]
local inter_se = _se[1.newparent#c.ln_income]
local n = e(N)

tempname fh
file open `fh' using "`outdir'/acs_center_origin_income_flight_regression_v1.csv", write replace
file write `fh' "outcome,income_p25,income_p75,newparent_income_interaction,interaction_se,newparent_effect_p25,newparent_effect_p25_se,newparent_effect_p25_p,newparent_effect_p75,newparent_effect_p75_se,newparent_effect_p75_p,n_obs" _n
file write `fh' "periphery_dest_center_origin_mass,`p25',`p75',`inter',`inter_se',`low',`low_se',`low_p',`high',`high_se',`high_p',`n'" _n
file close `fh'

use `within_origin', clear
xtile income_tercile = ln_income [aw=origin_center_weight], n(3)
collapse (mean) periphery_share=periphery_dest (count) n_obs=periphery_dest [aw=origin_center_weight], by(parent_compare income_tercile)
export delimited using "`outdir'/acs_center_origin_income_flight_terciles_v1.csv", replace

di as result "ACS center-origin income flight test complete."
di as result "Regression CSV: `outdir'/acs_center_origin_income_flight_regression_v1.csv"

log close _all
