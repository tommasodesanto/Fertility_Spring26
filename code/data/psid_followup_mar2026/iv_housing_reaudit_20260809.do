clear all
set more off
version 17.0

local project "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`project'/code/data/psid_followup_mar2026/output/iv_housing_reaudit_20260809"

capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/iv_housing_reaudit_20260809.log", replace text

tempfile sex_map sex1_map sex2_map parents panel first_baseline first_outcomes first_sample second_baseline second_outcomes second_sample

di as text "Step 1/7: Build child-sex map"
use ID year SEX using "`dta'", clear
drop if missing(ID, SEX)
bysort ID (year): keep if _n == 1
rename ID child_id
rename SEX child_sex
save `sex_map', replace

use `sex_map', clear
rename child_id child1_id
rename child_sex child1_sex
save `sex1_map', replace

use `sex_map', clear
rename child_id child2_id
rename child_sex child2_sex
save `sex2_map', replace

di as text "Step 2/7: Build one final fertility roster per parent"
use ID year RELCHIREP RELCHINUM RELCHI1ID RELCHI2ID RELCHI1BYEAR RELCHI2BYEAR ///
    RELCHI3BYEAR RELCHI4BYEAR RELCHI5BYEAR RELCHI6BYEAR RELCHI7BYEAR RELCHI8BYEAR ///
    RELCHI9BYEAR RELCHI10BYEAR RELCHI11BYEAR RELCHI12BYEAR RELCHI13BYEAR RELCHI14BYEAR ///
    RELCHI15BYEAR RELCHI16BYEAR RELCHI17BYEAR RELCHI18BYEAR RELCHI19BYEAR RELCHI20BYEAR ///
    using "`dta'", clear

bysort ID: egen roster_year = max(cond(!missing(RELCHI1BYEAR), year, .))
keep if year == roster_year & !missing(roster_year)
bysort ID (year): keep if _n == 1

rename RELCHI1ID child1_id
rename RELCHI2ID child2_id
merge m:1 child1_id using `sex1_map', keep(master match) nogen
merge m:1 child2_id using `sex2_map', keep(master match) nogen

rename RELCHI1BYEAR birth1
rename RELCHI2BYEAR birth2
rename RELCHI3BYEAR birth3

forvalues k = 1/20 {
    local src = cond(`k' == 1, "birth1", cond(`k' == 2, "birth2", cond(`k' == 3, "birth3", "RELCHI`k'BYEAR")))
    gen byte same_first_`k' = (`src' == birth1) if !missing(`src', birth1)
    replace same_first_`k' = 0 if missing(same_first_`k')
}
egen n_sameyear_first = rowtotal(same_first_1-same_first_20)
gen byte twin_first_proxy = n_sameyear_first >= 2 if !missing(birth1)
drop same_first_1-same_first_20

gen byte same_sex_first2 = child1_sex == child2_sex if !missing(child1_sex, child2_sex)
gen byte first_child_girl = child1_sex == 2 if !missing(child1_sex)

* Horizon-aligned treatments. These avoid using births that occur after the housing window.
gen byte two_plus_by5 = !missing(birth2) & birth2 <= birth1 + 5 if !missing(birth1)
gen byte three_plus_by5_after_second = !missing(birth3) & birth3 <= birth2 + 5 if !missing(birth2)

* The sex-composition design is restricted to singleton first and second births.
gen byte singleton_first2 = birth2 > birth1 if !missing(birth1, birth2)
gen byte singleton_second = birth3 != birth2 if !missing(birth2)
replace singleton_second = 1 if missing(birth3) & !missing(birth2)

keep ID child1_id child2_id child1_sex child2_sex first_child_girl birth1 birth2 birth3 ///
    roster_year RELCHIREP RELCHINUM twin_first_proxy same_sex_first2 ///
    two_plus_by5 three_plus_by5_after_second singleton_first2 singleton_second
isid ID
save `parents', replace

di as text "Step 3/7: Build person-year housing panel with explicit missingness"
use ID year AGEREP EDUYEAR SEX DEATHYEAR HOMEOWN MOVEDFREF_ WHYMOVED1_ ACTUALROOMS_ IW ///
    using "`dta'", clear

drop if year > DEATHYEAR & !missing(DEATHYEAR)
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

rename HOMEOWN own
replace own = 0 if own == 2
replace own = . if own == 3

rename MOVEDFREF_ moved_since_interview
replace moved_since_interview = . if inlist(moved_since_interview, 8, 9)
replace moved_since_interview = 0 if moved_since_interview == 5
rename ACTUALROOMS_ rooms
gen byte invalid_rooms_code = !missing(rooms) & !inrange(rooms, 1, 20)
replace rooms = . if !inrange(rooms, 1, 20)

sort ID year
by ID: gen prior_observed_own = own[_n-1]
by ID: gen prior_observed_year = year[_n-1]

gen byte moved_to_own = .
replace moved_to_own = (own == 1 & prior_observed_own == 0 & moved_since_interview == 1) ///
    if !missing(own, prior_observed_own, moved_since_interview)

gen byte moved_for_size = .
replace moved_for_size = 0 if moved_since_interview == 0
replace moved_for_size = (WHYMOVED1_ == 3) if moved_since_interview == 1 & !missing(WHYMOVED1_)

merge m:1 ID using `parents', keep(match) nogen
save `panel', replace

di as text "Step 4/7: Construct first-birth outcomes for the twins design"
use `panel', clear
gen k1 = year - birth1

bysort ID: egen first_baseline_year = max(cond(k1 < 0 & !missing(AGEREP, EDUYEAR, SEX, own, IW), year, .))
preserve
    keep if year == first_baseline_year & !missing(first_baseline_year)
    keep ID first_baseline_year AGEREP EDUYEAR SEX own rooms IW
    rename AGEREP age_pre
    rename EDUYEAR eduyear_pre
    rename SEX sex_pre
    rename own own_pre
    rename rooms rooms_pre
    rename IW iw_pre
    save `first_baseline', replace
restore

gen own_first5_i = own if inrange(k1, 0, 5)
gen moveown_first5_i = moved_to_own if inrange(k1, 0, 5)
gen msize_first5_i = moved_for_size if inrange(k1, 0, 5)
gen rooms_first5_i = rooms if inrange(k1, 0, 5)
gen observed_first5_i = k1 if inrange(k1, 0, 5)

bysort ID: egen own_post5 = max(own_first5_i)
bysort ID: egen moved_to_own_post5 = max(moveown_first5_i)
bysort ID: egen moved_for_size_post5 = max(msize_first5_i)
bysort ID: egen rooms_post5 = mean(rooms_first5_i)
bysort ID: egen latest_post_k = max(observed_first5_i)

bysort ID (year): keep if _n == 1
keep ID child1_id birth1 twin_first_proxy two_plus_by5 own_post5 moved_to_own_post5 ///
    moved_for_size_post5 rooms_post5 latest_post_k first_child_girl
save `first_outcomes', replace

use `first_outcomes', clear
merge 1:1 ID using `first_baseline', keep(match) nogen
gen female_pre = sex_pre == 2 if !missing(sex_pre)
gen baseline_gap = birth1 - first_baseline_year
gen rooms_change_post5 = rooms_post5 - rooms_pre if !missing(rooms_post5, rooms_pre)

* Mother-level main sample, with a nearby pre-birth baseline and late-window observation.
keep if female_pre == 1
keep if birth1 <= 2014
keep if inrange(baseline_gap, 1, 4)
keep if latest_post_k >= 4
drop if missing(age_pre, eduyear_pre, own_pre, iw_pre, twin_first_proxy, two_plus_by5)
isid ID
save `first_sample', replace
save "`outdir'/twins_firstbirth_clean_sample.dta", replace

di as text "Step 5/7: Construct second-birth outcomes for the same-sex design"
use `panel', clear
drop if missing(birth2)
gen k2 = year - birth2

bysort ID: egen second_baseline_year = max(cond(k2 < 0 & !missing(AGEREP, EDUYEAR, SEX, own, IW), year, .))
preserve
    keep if year == second_baseline_year & !missing(second_baseline_year)
    keep ID second_baseline_year AGEREP EDUYEAR SEX own rooms IW
    rename AGEREP age_pre
    rename EDUYEAR eduyear_pre
    rename SEX sex_pre
    rename own own_pre
    rename rooms rooms_pre
    rename IW iw_pre
    save `second_baseline', replace
restore

gen own_second5_i = own if inrange(k2, 0, 5)
gen moveown_second5_i = moved_to_own if inrange(k2, 0, 5)
gen msize_second5_i = moved_for_size if inrange(k2, 0, 5)
gen rooms_second5_i = rooms if inrange(k2, 0, 5)
gen observed_second5_i = k2 if inrange(k2, 0, 5)

bysort ID: egen own_post5 = max(own_second5_i)
bysort ID: egen moved_to_own_post5 = max(moveown_second5_i)
bysort ID: egen moved_for_size_post5 = max(msize_second5_i)
bysort ID: egen rooms_post5 = mean(rooms_second5_i)
bysort ID: egen latest_post_k = max(observed_second5_i)

bysort ID (year): keep if _n == 1
keep ID child1_id birth1 birth2 birth3 same_sex_first2 first_child_girl ///
    three_plus_by5_after_second singleton_first2 singleton_second own_post5 ///
    moved_to_own_post5 moved_for_size_post5 rooms_post5 latest_post_k
save `second_outcomes', replace

use `second_outcomes', clear
merge 1:1 ID using `second_baseline', keep(match) nogen
gen female_pre = sex_pre == 2 if !missing(sex_pre)
gen baseline_gap = birth2 - second_baseline_year
gen rooms_change_post5 = rooms_post5 - rooms_pre if !missing(rooms_post5, rooms_pre)

* Align the design at the second birth and exclude same-year multiple-birth proxies.
keep if female_pre == 1
keep if birth2 <= 2014
keep if singleton_first2 == 1 & singleton_second == 1
keep if inrange(baseline_gap, 1, 4)
keep if latest_post_k >= 4
drop if missing(age_pre, eduyear_pre, own_pre, iw_pre, same_sex_first2, ///
    first_child_girl, three_plus_by5_after_second)
isid ID
save `second_sample', replace
save "`outdir'/samesex_secondbirth_clean_sample.dta", replace

di as text "Step 6/7: Estimate first stages, reduced forms, OLS, and 2SLS"
tempname estimates samplepost balancepost
postfile `estimates' str12 design str32 specification str28 outcome str24 estimand ///
    double b se p ci_low ci_high N first_stage_F using "`outdir'/iv_housing_reaudit_estimates.dta", replace
postfile `samplepost' str12 design str32 specification str28 outcome ///
    double N instrument_ones instrument_share treatment_share outcome_mean using "`outdir'/iv_housing_reaudit_samples.dta", replace
postfile `balancepost' str12 design str24 covariate double b se p N using "`outdir'/iv_housing_reaudit_balance.dta", replace

* Twins at first birth. Housing-space outcomes use all baseline tenures;
* ownership-transition outcomes use baseline renters.
use `first_sample', clear
local controls "c.age_pre##c.age_pre i.eduyear_pre i.birth1 i.own_pre c.rooms_pre"
foreach y in moved_for_size_post5 rooms_change_post5 own_post5 moved_to_own_post5 {
    preserve
        if inlist("`y'", "own_post5", "moved_to_own_post5") {
            keep if own_pre == 0
        }
        keep if !missing(`y', rooms_pre)

        quietly count
        local n_sample = r(N)
        quietly count if twin_first_proxy == 1
        local n_z1 = r(N)
        quietly summarize twin_first_proxy [aw=iw_pre], meanonly
        local zbar = r(mean)
        quietly summarize two_plus_by5 [aw=iw_pre], meanonly
        local dbar = r(mean)
        quietly summarize `y' [aw=iw_pre], meanonly
        local ybar = r(mean)
        post `samplepost' ("twins") ("mother_firstbirth_post5") ("`y'") ///
            (`n_sample') (`n_z1') (`zbar') (`dbar') (`ybar')

        quietly reg two_plus_by5 twin_first_proxy `controls' [aw=iw_pre], vce(robust)
        quietly test twin_first_proxy
        local firstF = r(F)
        local bb = _b[twin_first_proxy]
        local ss = _se[twin_first_proxy]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("twins") ("mother_firstbirth_post5") ("`y'") ("first_stage") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF')

        quietly reg `y' twin_first_proxy `controls' [aw=iw_pre], vce(robust)
        local bb = _b[twin_first_proxy]
        local ss = _se[twin_first_proxy]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("twins") ("mother_firstbirth_post5") ("`y'") ("reduced_form") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF')

        quietly reg `y' two_plus_by5 `controls' [aw=iw_pre], vce(robust)
        local bb = _b[two_plus_by5]
        local ss = _se[two_plus_by5]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("twins") ("mother_firstbirth_post5") ("`y'") ("ols") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF')

        quietly ivregress 2sls `y' (two_plus_by5 = twin_first_proxy) `controls' [aw=iw_pre], vce(robust)
        local bb = _b[two_plus_by5]
        local ss = _se[two_plus_by5]
        local pp = 2 * normal(-abs(`bb'/`ss'))
        local cc = invnormal(.975)
        post `estimates' ("twins") ("mother_firstbirth_post5") ("`y'") ("iv_2sls") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF')

        * Unweighted sensitivity: useful because there are few same-birth-year events.
        quietly reg two_plus_by5 twin_first_proxy `controls', vce(robust)
        quietly test twin_first_proxy
        local firstF_u = r(F)
        local bb = _b[twin_first_proxy]
        local ss = _se[twin_first_proxy]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("twins") ("mother_firstbirth_post5") ("`y'") ("first_stage_unweighted") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF_u')

        quietly reg `y' twin_first_proxy `controls', vce(robust)
        local bb = _b[twin_first_proxy]
        local ss = _se[twin_first_proxy]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("twins") ("mother_firstbirth_post5") ("`y'") ("reduced_form_unweighted") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF_u')
    restore
}

* Covariate balance for the twins instrument, conditional on birth cohort.
foreach x in age_pre eduyear_pre own_pre rooms_pre first_child_girl {
    quietly reg `x' twin_first_proxy i.birth1 [aw=iw_pre], vce(robust)
    local bb = _b[twin_first_proxy]
    local ss = _se[twin_first_proxy]
    local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
    post `balancepost' ("twins") ("`x'") (`bb') (`ss') (`pp') (e(N))
}

* Same sex of the first two children, indexed at the second birth.
use `second_sample', clear
local controls "i.first_child_girl c.age_pre##c.age_pre i.eduyear_pre i.birth2 i.own_pre c.rooms_pre"
foreach y in moved_for_size_post5 rooms_change_post5 own_post5 moved_to_own_post5 {
    preserve
        if inlist("`y'", "own_post5", "moved_to_own_post5") {
            keep if own_pre == 0
        }
        keep if !missing(`y', rooms_pre)

        quietly count
        local n_sample = r(N)
        quietly count if same_sex_first2 == 1
        local n_z1 = r(N)
        quietly summarize same_sex_first2 [aw=iw_pre], meanonly
        local zbar = r(mean)
        quietly summarize three_plus_by5_after_second [aw=iw_pre], meanonly
        local dbar = r(mean)
        quietly summarize `y' [aw=iw_pre], meanonly
        local ybar = r(mean)
        post `samplepost' ("same_sex") ("mother_secondbirth_post5") ("`y'") ///
            (`n_sample') (`n_z1') (`zbar') (`dbar') (`ybar')

        quietly reg three_plus_by5_after_second same_sex_first2 `controls' [aw=iw_pre], vce(robust)
        quietly test same_sex_first2
        local firstF = r(F)
        local bb = _b[same_sex_first2]
        local ss = _se[same_sex_first2]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("same_sex") ("mother_secondbirth_post5") ("`y'") ("first_stage") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF')

        quietly reg `y' same_sex_first2 `controls' [aw=iw_pre], vce(robust)
        local bb = _b[same_sex_first2]
        local ss = _se[same_sex_first2]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("same_sex") ("mother_secondbirth_post5") ("`y'") ("reduced_form") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF')

        quietly reg `y' three_plus_by5_after_second `controls' [aw=iw_pre], vce(robust)
        local bb = _b[three_plus_by5_after_second]
        local ss = _se[three_plus_by5_after_second]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("same_sex") ("mother_secondbirth_post5") ("`y'") ("ols") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF')

        quietly ivregress 2sls `y' (three_plus_by5_after_second = same_sex_first2) ///
            `controls' [aw=iw_pre], vce(robust)
        local bb = _b[three_plus_by5_after_second]
        local ss = _se[three_plus_by5_after_second]
        local pp = 2 * normal(-abs(`bb'/`ss'))
        local cc = invnormal(.975)
        post `estimates' ("same_sex") ("mother_secondbirth_post5") ("`y'") ("iv_2sls") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF')

        quietly reg three_plus_by5_after_second same_sex_first2 `controls', vce(robust)
        quietly test same_sex_first2
        local firstF_u = r(F)
        local bb = _b[same_sex_first2]
        local ss = _se[same_sex_first2]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("same_sex") ("mother_secondbirth_post5") ("`y'") ("first_stage_unweighted") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF_u')

        quietly reg `y' same_sex_first2 `controls', vce(robust)
        local bb = _b[same_sex_first2]
        local ss = _se[same_sex_first2]
        local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
        local cc = invttail(e(df_r), .025)
        post `estimates' ("same_sex") ("mother_secondbirth_post5") ("`y'") ("reduced_form_unweighted") ///
            (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`firstF_u')
    restore
}

* Covariate balance for sex composition, conditional on first-child sex and cohort.
foreach x in age_pre eduyear_pre own_pre rooms_pre {
    quietly reg `x' same_sex_first2 i.first_child_girl i.birth2 [aw=iw_pre], vce(robust)
    local bb = _b[same_sex_first2]
    local ss = _se[same_sex_first2]
    local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
    post `balancepost' ("same_sex") ("`x'") (`bb') (`ss') (`pp') (e(N))
}

postclose `estimates'
postclose `samplepost'
postclose `balancepost'

use "`outdir'/iv_housing_reaudit_estimates.dta", clear
sort design outcome estimand
export delimited using "`outdir'/iv_housing_reaudit_estimates.csv", replace

use "`outdir'/iv_housing_reaudit_samples.dta", clear
sort design outcome
export delimited using "`outdir'/iv_housing_reaudit_samples.csv", replace

use "`outdir'/iv_housing_reaudit_balance.dta", clear
sort design covariate
export delimited using "`outdir'/iv_housing_reaudit_balance.csv", replace

di as text "Step 7/7: Save construction audit"
use `panel', clear
gen k1 = year - birth1
xtset ID year

tempname constructionpost
postfile `constructionpost' str48 metric double value denominator using ///
    "`outdir'/iv_housing_reaudit_construction.dta", replace

quietly count if invalid_rooms_code == 1
local invalid_rooms_n = r(N)
quietly count
post `constructionpost' ("invalid_room_codes_person_years") (`invalid_rooms_n') (r(N))

* Original calendar-lag transition rule versus the corrected previous-observation rule.
gen byte legacy_moveown_i = (own == 1 & L.own == 0 & moved_since_interview == 1) if inrange(k1, 0, 5)
gen byte clean_moveown_i = moved_to_own if inrange(k1, 0, 5)
bysort ID: egen legacy_moveown = max(legacy_moveown_i)
bysort ID: egen clean_moveown = max(clean_moveown_i)

* The original +3 ownership code evaluated rows outside the window as zero.
gen byte legacy_own3_i = (own == 1) & inrange(k1, 0, 3)
gen byte clean_own3_i = own if inrange(k1, 0, 3)
bysort ID: egen legacy_own3 = max(legacy_own3_i)
bysort ID: egen clean_own3 = max(clean_own3_i)

* The original all-tenure +3 script used the reason code without confirming a move.
gen byte reason_only_msize3_i = (WHYMOVED1_ == 3) if inrange(k1, 0, 3) & !missing(WHYMOVED1_)
gen byte clean_msize3_i = moved_for_size if inrange(k1, 0, 3)
bysort ID: egen reason_only_msize3 = max(reason_only_msize3_i)
bysort ID: egen clean_msize3 = max(clean_msize3_i)

bysort ID (year): keep if _n == 1

quietly count if !missing(legacy_moveown, clean_moveown)
local transition_n = r(N)
quietly count if legacy_moveown == 1 & !missing(clean_moveown)
post `constructionpost' ("calendar_lag_move_to_own_positive") (r(N)) (`transition_n')
quietly count if clean_moveown == 1 & !missing(legacy_moveown)
post `constructionpost' ("observed_lag_move_to_own_positive") (r(N)) (`transition_n')
quietly count if legacy_moveown != clean_moveown & !missing(legacy_moveown, clean_moveown)
post `constructionpost' ("move_to_own_classification_changed") (r(N)) (`transition_n')

quietly count if missing(clean_own3) & legacy_own3 == 0
local false_zero_n = r(N)
quietly count if missing(clean_own3)
post `constructionpost' ("legacy_own_post3_false_zeros") (`false_zero_n') (r(N))

quietly count if !missing(reason_only_msize3, clean_msize3)
local msize_n = r(N)
quietly count if reason_only_msize3 == 1 & !missing(clean_msize3)
post `constructionpost' ("reason_only_moved_for_size_positive") (r(N)) (`msize_n')
quietly count if clean_msize3 == 1 & !missing(reason_only_msize3)
post `constructionpost' ("move_confirmed_moved_for_size_positive") (r(N)) (`msize_n')
quietly count if reason_only_msize3 != clean_msize3 & !missing(reason_only_msize3, clean_msize3)
post `constructionpost' ("moved_for_size_classification_changed") (r(N)) (`msize_n')

postclose `constructionpost'
use "`outdir'/iv_housing_reaudit_construction.dta", clear
export delimited using "`outdir'/iv_housing_reaudit_construction.csv", replace

di as result "IV housing re-audit complete: `outdir'"
log close _all
