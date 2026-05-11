clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "`root'/Fertility_Spring26/april_26_discrete_time/output"
local outdir "`out_root'/empirical_roundup_income_elasticity_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/empirical_roundup_income_elasticity_v1.log", replace text

global weight IW

use ID year AGEREP DEATHYEAR RELCHIREP RELCHI1BYEAR HOMEOWN ACTUALROOMS_ INCFAMR ${weight} ///
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

drop if missing(ID, year, AGEREP, any_kid, ln_income, own)

xtset ID year

estimates clear

di as text "Estimating FE rooms elasticity by parental status"
xtreg ln_rooms c.ln_income##i.any_kid i.year i.AGEREP if !missing(ln_rooms), fe vce(cluster ID)
estimates store fe_rooms

lincom ln_income
scalar b_rooms_childless = r(estimate)
scalar se_rooms_childless = r(se)
scalar p_rooms_childless = r(p)

lincom 1.any_kid#c.ln_income
scalar b_rooms_diff = r(estimate)
scalar se_rooms_diff = r(se)
scalar p_rooms_diff = r(p)

lincom ln_income + 1.any_kid#c.ln_income
scalar b_rooms_parent = r(estimate)
scalar se_rooms_parent = r(se)
scalar p_rooms_parent = r(p)

tempvar tag_rooms
gen `tag_rooms' = e(sample)
egen tag_id_rooms = tag(ID) if `tag_rooms' == 1
quietly count if tag_id_rooms == 1
scalar n_ids_rooms = r(N)
drop tag_id_rooms `tag_rooms'

di as text "Estimating FE ownership-income semi-elasticity by parental status"
xtreg own c.ln_income##i.any_kid i.year i.AGEREP, fe vce(cluster ID)
estimates store fe_own

lincom ln_income
scalar b_own_childless = r(estimate)
scalar se_own_childless = r(se)
scalar p_own_childless = r(p)

lincom 1.any_kid#c.ln_income
scalar b_own_diff = r(estimate)
scalar se_own_diff = r(se)
scalar p_own_diff = r(p)

lincom ln_income + 1.any_kid#c.ln_income
scalar b_own_parent = r(estimate)
scalar se_own_parent = r(se)
scalar p_own_parent = r(p)

tempvar tag_own
gen `tag_own' = e(sample)
egen tag_id_own = tag(ID) if `tag_own' == 1
quietly count if tag_id_own == 1
scalar n_ids_own = r(N)
drop tag_id_own `tag_own'

capture which esttab
if _rc == 0 {
    esttab fe_rooms fe_own using "`outdir'/income_elasticity_models_v1.tex", replace ///
        b(3) se(3) label ///
        mtitles("Log rooms FE" "Ownership FE") ///
        stats(N, fmt(0) labels("Obs."))
}

tempname fh
file open `fh' using "`outdir'/income_elasticity_summary_v1.csv", write replace
file write `fh' "outcome,childless_effect,childless_se,childless_p,parent_diff,parent_diff_se,parent_diff_p,parent_effect,parent_effect_se,parent_effect_p,n_ids" _n

local b1 : display %10.6f scalar(b_rooms_childless)
local s1 : display %10.6f scalar(se_rooms_childless)
local p1 : display %10.6f scalar(p_rooms_childless)
local b2 : display %10.6f scalar(b_rooms_diff)
local s2 : display %10.6f scalar(se_rooms_diff)
local p2 : display %10.6f scalar(p_rooms_diff)
local b3 : display %10.6f scalar(b_rooms_parent)
local s3 : display %10.6f scalar(se_rooms_parent)
local p3 : display %10.6f scalar(p_rooms_parent)
local n1 : display %10.0f scalar(n_ids_rooms)
file write `fh' "ln_rooms,`b1',`s1',`p1',`b2',`s2',`p2',`b3',`s3',`p3',`n1'" _n

local c1 : display %10.6f scalar(b_own_childless)
local c2 : display %10.6f scalar(se_own_childless)
local c3 : display %10.6f scalar(p_own_childless)
local c4 : display %10.6f scalar(b_own_diff)
local c5 : display %10.6f scalar(se_own_diff)
local c6 : display %10.6f scalar(p_own_diff)
local c7 : display %10.6f scalar(b_own_parent)
local c8 : display %10.6f scalar(se_own_parent)
local c9 : display %10.6f scalar(p_own_parent)
local n2 : display %10.0f scalar(n_ids_own)
file write `fh' "own,`c1',`c2',`c3',`c4',`c5',`c6',`c7',`c8',`c9',`n2'" _n
file close `fh'

di as result "Income-elasticity regressions complete."
di as result "Summary CSV: `outdir'/income_elasticity_summary_v1.csv"

log close _all
