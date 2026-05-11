clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local extract "`root'/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
local outdir "`root'/april_26_discrete_time/output/acs_renter_elasticity_v1"

cap mkdir "`outdir'"
capture log close _all
log using "`outdir'/acs_renter_elasticity_v1.log", replace text

use year met2013 age sex race perwt gq ownershp rooms rent hhincome nchild ///
    using "`extract'", clear

keep if year >= 2012 & year <= 2023
keep if age >= 22 & age <= 45
keep if inlist(gq, 1, 2)

gen renter = ownershp == 2 if !missing(ownershp)
gen any_parent = nchild > 0 if !missing(nchild)
gen childless = nchild == 0 if !missing(nchild)
gen ln_income = ln(hhincome) if hhincome > 1000 & !missing(hhincome)
gen ln_rooms = ln(rooms) if rooms > 0 & !missing(rooms)
gen ln_rent = ln(rent) if rent > 0 & !missing(rent)

keep if renter == 1
keep if met2013 > 0

quietly _pctile hhincome [aw=perwt] if hhincome > 1000 & !missing(hhincome), p(1 99)
local q1 = r(r1)
local q99 = r(r2)
keep if missing(hhincome) | (hhincome >= `q1' & hhincome <= `q99')

tempname fh
file open `fh' using "`outdir'/acs_renter_elasticity_summary_v1.csv", write replace
file write `fh' "sample,outcome,childless_slope,childless_slope_se,childless_slope_p,parent_diff,parent_diff_se,parent_diff_p,parent_slope,parent_slope_se,parent_slope_p,n_obs" _n

reg ln_rooms c.ln_income##i.any_parent i.sex i.race i.age i.year i.met2013 [aw=perwt] ///
    if !missing(ln_rooms, ln_income, any_parent), vce(robust)
lincom ln_income
local r_child = r(estimate)
local r_child_se = r(se)
local r_child_p = r(p)
lincom 1.any_parent#c.ln_income
local r_diff = r(estimate)
local r_diff_se = r(se)
local r_diff_p = r(p)
lincom ln_income + 1.any_parent#c.ln_income
local r_parent = r(estimate)
local r_parent_se = r(se)
local r_parent_p = r(p)
local rn = e(N)
file write `fh' "acs_renters,ln_rooms,`r_child',`r_child_se',`r_child_p',`r_diff',`r_diff_se',`r_diff_p',`r_parent',`r_parent_se',`r_parent_p',`rn'" _n

reg ln_rent c.ln_income##i.any_parent i.sex i.race i.age i.year i.met2013 [aw=perwt] ///
    if !missing(ln_rent, ln_income, any_parent), vce(robust)
lincom ln_income
local t_child = r(estimate)
local t_child_se = r(se)
local t_child_p = r(p)
lincom 1.any_parent#c.ln_income
local t_diff = r(estimate)
local t_diff_se = r(se)
local t_diff_p = r(p)
lincom ln_income + 1.any_parent#c.ln_income
local t_parent = r(estimate)
local t_parent_se = r(se)
local t_parent_p = r(p)
local tn = e(N)
file write `fh' "acs_renters,ln_rent,`t_child',`t_child_se',`t_child_p',`t_diff',`t_diff_se',`t_diff_p',`t_parent',`t_parent_se',`t_parent_p',`tn'" _n

file close `fh'

di as result "ACS renter-only elasticity complete."
di as result "Summary CSV: `outdir'/acs_renter_elasticity_summary_v1.csv"

log close _all
