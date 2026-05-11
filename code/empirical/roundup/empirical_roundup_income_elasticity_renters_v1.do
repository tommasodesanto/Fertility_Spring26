clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "`root'/Fertility_Spring26/april_26_discrete_time/output"
local outdir "`out_root'/empirical_roundup_income_elasticity_renters_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/empirical_roundup_income_elasticity_renters_v1.log", replace text

global weight IW

use ID year AGEREP DEATHYEAR RELCHIREP HOMEOWN ACTUALROOMS_ INCFAMR ${weight} ///
    using "`dta'", clear

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

drop if year > DEATHYEAR
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

rename HOMEOWN own
rename ACTUALROOMS_ rooms

gen any_kid = (RELCHIREP > 0) if !missing(RELCHIREP)
gen ln_income = log(INCFAMR) if INCFAMR > 1000 & !missing(INCFAMR)
gen ln_rooms = log(rooms) if rooms > 0 & !missing(rooms)

keep if own == 0
drop if missing(ID, year, AGEREP, any_kid, ln_income, ln_rooms)

xtset ID year

xtreg ln_rooms c.ln_income##i.any_kid i.year i.AGEREP, fe vce(cluster ID)
estimates store fe_renter_rooms

lincom ln_income
scalar b_childless = r(estimate)
scalar se_childless = r(se)
scalar p_childless = r(p)

lincom 1.any_kid#c.ln_income
scalar b_diff = r(estimate)
scalar se_diff = r(se)
scalar p_diff = r(p)

lincom ln_income + 1.any_kid#c.ln_income
scalar b_parent = r(estimate)
scalar se_parent = r(se)
scalar p_parent = r(p)

tempvar tag_samp
gen `tag_samp' = e(sample)
egen tag_id = tag(ID) if `tag_samp' == 1
quietly count if tag_id == 1
scalar n_ids = r(N)
drop tag_id `tag_samp'

quietly count if e(sample)
scalar n_obs = r(N)

capture which esttab
if _rc == 0 {
    esttab fe_renter_rooms using "`outdir'/income_elasticity_renters_model_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("PSID renters: log rooms FE") ///
        stats(N, fmt(0) labels("Obs."))
}

tempname fh
file open `fh' using "`outdir'/income_elasticity_renters_summary_v1.csv", write replace
file write `fh' "sample,outcome,childless_effect,childless_se,childless_p,parent_diff,parent_diff_se,parent_diff_p,parent_effect,parent_effect_se,parent_effect_p,n_obs,n_ids" _n

local b1 : display %12.6f scalar(b_childless)
local s1 : display %12.6f scalar(se_childless)
local p1 : display %12.6f scalar(p_childless)
local b2 : display %12.6f scalar(b_diff)
local s2 : display %12.6f scalar(se_diff)
local p2 : display %12.6f scalar(p_diff)
local b3 : display %12.6f scalar(b_parent)
local s3 : display %12.6f scalar(se_parent)
local p3 : display %12.6f scalar(p_parent)
local n1 : display %12.0f scalar(n_obs)
local n2 : display %12.0f scalar(n_ids)
file write `fh' "psid_renters,ln_rooms,`b1',`s1',`p1',`b2',`s2',`p2',`b3',`s3',`p3',`n1',`n2'" _n
file close `fh'

di as result "Renter-only PSID elasticity complete."
di as result "Summary CSV: `outdir'/income_elasticity_renters_summary_v1.csv"

log close _all
