clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local extract "`root'/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta"
local lookup2010 "`root'/code/data/mms_center_periphery/data/puma_mms_lookup_2010.csv"
local lookup2020 "`root'/code/data/mms_center_periphery/data/puma_mms_lookup_2020.csv"
local outdir "`root'/april_26_discrete_time/output/acs_income_crossmetro_v1"

cap mkdir "`outdir'"
capture log close _all
log using "`outdir'/acs_income_crossmetro_v1.log", replace text

tempfile lookup extract_base restricted_np within_cbsa rentgap centergap destgap

import delimited using "`lookup2010'", clear varnames(1) stringcols(_all)
gen lookup_period = "pre2022"
destring statefip puma cbsacode, replace force
keep lookup_period statefip puma cbsacode cbsatitle mms_location
tempfile lookup_pre
save `lookup_pre', replace

import delimited using "`lookup2020'", clear varnames(1) stringcols(_all)
gen lookup_period = "post2021"
destring statefip puma cbsacode, replace force
keep lookup_period statefip puma cbsacode cbsatitle mms_location
append using `lookup_pre'
rename cbsacode met2013
replace mms_location = "center" if mms_location == "middle"
keep if inlist(mms_location, "center", "periphery")
save `lookup', replace

use year statefip puma metro met2013 migmet131 age sex race perwt gq ownershp ///
    rooms rent hhincome nchild nchlt5 eldch migrate1 using "`extract'", clear

keep if year >= 2012 & year <= 2023
keep if age >= 22 & age <= 45
keep if inlist(gq, 1, 2)

gen lookup_period = cond(year <= 2021, "pre2022", "post2021")

gen owner = ownershp == 1 if !missing(ownershp)
gen renter = ownershp == 2 if !missing(ownershp)
gen has_children = nchild > 0 if !missing(nchild)
gen childless = nchild == 0 if !missing(nchild)
gen newparent = eldch < 4 & eldch != 99 & nchild > 0 if !missing(eldch, nchild)
gen any_parent = has_children
replace newparent = 0 if childless == 1 & missing(newparent)
replace any_parent = 0 if childless == 1 & missing(any_parent)
gen moved1y = migrate1 >= 2 if !missing(migrate1)
gen same_msa = moved1y == 1 & met2013 > 0 & migmet131 > 0 & met2013 == migmet131 if !missing(moved1y, met2013, migmet131)

gen ln_income = ln(hhincome) if hhincome > 1000 & !missing(hhincome)
gen ln_rooms = ln(rooms) if rooms > 0 & !missing(rooms)
gen rent_per_room = rent / rooms if renter == 1 & rent > 0 & rooms > 0

merge m:1 lookup_period statefip puma met2013 using `lookup', keep(match) nogen

drop if missing(mms_location, ln_income)
gen dest_label = cond(mms_location == "center", "Center", "Periphery")
gen center_live = mms_location == "center"
gen center_dest = center_live

quietly _pctile hhincome [aw=perwt] if hhincome > 1000 & !missing(hhincome), p(1 99)
local q1 = r(r1)
local q99 = r(r2)
keep if missing(hhincome) | (hhincome >= `q1' & hhincome <= `q99')

gen parent_compare = ""
replace parent_compare = "Non-Parents" if childless == 1
replace parent_compare = "New Parents" if newparent == 1

save `extract_base', replace

use `extract_base', clear
keep if inlist(parent_compare, "Non-Parents", "New Parents")
save `restricted_np', replace

quietly _pctile ln_income [aw=perwt], p(25 75)
local p25 = r(r1)
local p75 = r(r2)

reg center_live c.ln_income##i.newparent i.sex i.race i.age i.year [aw=perwt], vce(robust)
lincom 1.newparent + `p25' * 1.newparent#c.ln_income
local clow = r(estimate)
local clow_se = r(se)
local clow_p = r(p)
lincom 1.newparent + `p75' * 1.newparent#c.ln_income
local chigh = r(estimate)
local chigh_se = r(se)
local chigh_p = r(p)
local cint = _b[1.newparent#c.ln_income]
local cint_se = _se[1.newparent#c.ln_income]
local cn = e(N)

use `restricted_np', clear
keep if same_msa == 1
save `within_cbsa', replace

reg center_dest c.ln_income##i.newparent i.sex i.race i.age i.year [aw=perwt], vce(robust)
lincom 1.newparent + `p25' * 1.newparent#c.ln_income
local dlow = r(estimate)
local dlow_se = r(se)
local dlow_p = r(p)
lincom 1.newparent + `p75' * 1.newparent#c.ln_income
local dhigh = r(estimate)
local dhigh_se = r(se)
local dhigh_p = r(p)
local dint = _b[1.newparent#c.ln_income]
local dint_se = _se[1.newparent#c.ln_income]
local dn = e(N)

use `extract_base', clear
reg ln_rooms c.ln_income##i.any_parent i.sex i.race i.age i.year [aw=perwt] if !missing(ln_rooms), vce(robust)
lincom ln_income
local r_child = r(estimate)
local r_child_se = r(se)
local r_child_p = r(p)
lincom ln_income + 1.any_parent#c.ln_income
local r_parent = r(estimate)
local r_parent_se = r(se)
local r_parent_p = r(p)
local r_diff = _b[1.any_parent#c.ln_income]
local r_diff_se = _se[1.any_parent#c.ln_income]
local rn = e(N)

reg owner c.ln_income##i.any_parent i.sex i.race i.age i.year [aw=perwt], vce(robust)
lincom ln_income
local o_child = r(estimate)
local o_child_se = r(se)
local o_child_p = r(p)
lincom ln_income + 1.any_parent#c.ln_income
local o_parent = r(estimate)
local o_parent_se = r(se)
local o_parent_p = r(p)
local o_diff = _b[1.any_parent#c.ln_income]
local o_diff_se = _se[1.any_parent#c.ln_income]
local on = e(N)

tempname fh1
file open `fh1' using "`outdir'/acs_income_parent_penalty_v1.csv", write replace
file write `fh1' "outcome,income_p25,income_p75,parent_income_interaction,interaction_se,parent_effect_p25,parent_effect_p25_se,parent_effect_p25_p,parent_effect_p75,parent_effect_p75_se,parent_effect_p75_p,n" _n
file write `fh1' "center_live_newparent_income,`p25',`p75',`cint',`cint_se',`clow',`clow_se',`clow_p',`chigh',`chigh_se',`chigh_p',`cn'" _n
file write `fh1' "center_dest_newparent_income,`p25',`p75',`dint',`dint_se',`dlow',`dlow_se',`dlow_p',`dhigh',`dhigh_se',`dhigh_p',`dn'" _n
file close `fh1'

tempname fh2
file open `fh2' using "`outdir'/acs_income_housing_slopes_v1.csv", write replace
file write `fh2' "outcome,childless_slope,childless_slope_se,childless_slope_p,parent_interaction,parent_interaction_se,parent_slope,parent_slope_se,parent_slope_p,n" _n
file write `fh2' "ln_rooms_anyparent_income,`r_child',`r_child_se',`r_child_p',`r_diff',`r_diff_se',`r_parent',`r_parent_se',`r_parent_p',`rn'" _n
file write `fh2' "own_anyparent_income,`o_child',`o_child_se',`o_child_p',`o_diff',`o_diff_se',`o_parent',`o_parent_se',`o_parent_p',`on'" _n
file close `fh2'

use `restricted_np', clear
xtile income_tercile = ln_income [aw=perwt], n(3)
collapse (mean) share=center_live (count) n_obs=center_live [aw=perwt], by(parent_compare income_tercile)
gen sample = "all_mms"
tempfile tertile_all
save `tertile_all', replace

use `within_cbsa', clear
xtile income_tercile = ln_income [aw=perwt], n(3)
collapse (mean) share=center_dest (count) n_obs=center_dest [aw=perwt], by(parent_compare income_tercile)
gen sample = "within_cbsa_movers"
append using `tertile_all'
order sample parent_compare income_tercile share n_obs
export delimited using "`outdir'/acs_income_proxy_terciles_v1.csv", replace

use `extract_base', clear
keep if renter == 1 & rent_per_room > 0 & !missing(cbsatitle, met2013, dest_label)
collapse (mean) mean_rpr=rent_per_room (count) rent_obs=rent_per_room [aw=perwt], by(cbsatitle met2013 dest_label)
reshape wide mean_rpr rent_obs, i(cbsatitle met2013) j(dest_label) string
gen rent_gap_ratio = mean_rprCenter / mean_rprPeriphery
gen log_rent_gap = ln(rent_gap_ratio)
save `rentgap', replace

use `restricted_np', clear
collapse (mean) center_share=center_live (count) obs=center_live [aw=perwt], by(cbsatitle met2013 parent_compare)
gen parent_key = cond(parent_compare == "Non-Parents", "np", "newp")
drop parent_compare
reshape wide center_share obs, i(cbsatitle met2013) j(parent_key) string
gen gap_pp = 100 * (center_sharenp - center_sharenewp)
gen metro_obs = obsnp + obsnewp
gen np_obs = obsnp
gen newparent_obs = obsnewp
merge 1:1 cbsatitle met2013 using `rentgap', keep(match) nogen
keep if np_obs >= 200 & newparent_obs >= 50 & rent_obsCenter >= 100 & rent_obsPeriphery >= 100 & !missing(log_rent_gap)
save `centergap', replace
export delimited using "`outdir'/acs_cross_metro_center_gap_v1.csv", replace

reg gap_pp log_rent_gap [aw=metro_obs], vce(robust)
local m1_b = _b[log_rent_gap]
local m1_se = _se[log_rent_gap]
test log_rent_gap = 0
local m1_p = r(p)
local m1_n = e(N)

twoway ///
    (scatter gap_pp log_rent_gap [aw=metro_obs], mcolor(navy%70) msymbol(circle_hollow)) ///
    (lfit gap_pp log_rent_gap [aw=metro_obs], lcolor(maroon) lwidth(medthick)), ///
    ytitle("Non-parent minus new-parent center share (pp)") ///
    xtitle("Log center/periphery rent-per-room gap") ///
    title("Cross-metro sorting gap vs rent gap")
graph export "`outdir'/acs_cross_metro_center_gap_v1.png", replace

use `within_cbsa', clear
collapse (mean) center_dest_share=center_dest (count) obs=center_dest [aw=perwt], by(cbsatitle met2013 parent_compare)
gen parent_key = cond(parent_compare == "Non-Parents", "np", "newp")
drop parent_compare
reshape wide center_dest_share obs, i(cbsatitle met2013) j(parent_key) string
gen gap_pp = 100 * (center_dest_sharenp - center_dest_sharenewp)
gen mover_obs = obsnp + obsnewp
gen np_obs = obsnp
gen newparent_obs = obsnewp
merge 1:1 cbsatitle met2013 using `rentgap', keep(match) nogen
keep if np_obs >= 50 & newparent_obs >= 25 & rent_obsCenter >= 100 & rent_obsPeriphery >= 100 & !missing(log_rent_gap)
save `destgap', replace
export delimited using "`outdir'/acs_cross_metro_dest_gap_v1.csv", replace

reg gap_pp log_rent_gap [aw=mover_obs], vce(robust)
local m2_b = _b[log_rent_gap]
local m2_se = _se[log_rent_gap]
test log_rent_gap = 0
local m2_p = r(p)
local m2_n = e(N)

twoway ///
    (scatter gap_pp log_rent_gap [aw=mover_obs], mcolor(forest_green%70) msymbol(circle_hollow)) ///
    (lfit gap_pp log_rent_gap [aw=mover_obs], lcolor(maroon) lwidth(medthick)), ///
    ytitle("Non-parent minus new-parent center destination gap (pp)") ///
    xtitle("Log center/periphery rent-per-room gap") ///
    title("Within-CBSA mover sorting gap vs rent gap")
graph export "`outdir'/acs_cross_metro_dest_gap_v1.png", replace

tempname fh3
file open `fh3' using "`outdir'/acs_cross_metro_regressions_v1.csv", write replace
file write `fh3' "outcome,slope,slope_se,p_value,n_metros" _n
file write `fh3' "center_residence_gap_pp,`m1_b',`m1_se',`m1_p',`m1_n'" _n
file write `fh3' "within_cbsa_center_dest_gap_pp,`m2_b',`m2_se',`m2_p',`m2_n'" _n
file close `fh3'

di as result "ACS income and cross-metro tests complete."
di as result "`outdir'/acs_income_parent_penalty_v1.csv"
di as result "`outdir'/acs_income_housing_slopes_v1.csv"
di as result "`outdir'/acs_cross_metro_regressions_v1.csv"

log close _all
