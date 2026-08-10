clear all
set more off
version 17.0

local project "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`project'/code/data/psid_followup_mar2026/output/iv_housing_postonly_reaudit_20260810"

capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/iv_housing_postonly_reaudit_20260810.log", replace text

tempfile sex_map sex1_map sex2_map parents panel twins_sample samesex_sample

di as text "Step 1/5: Build child-sex map and final fertility roster"
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

use ID year RELCHI1ID RELCHI2ID RELCHI1BYEAR RELCHI2BYEAR ///
    RELCHI3BYEAR RELCHI4BYEAR RELCHI5BYEAR RELCHI6BYEAR RELCHI7BYEAR ///
    RELCHI8BYEAR RELCHI9BYEAR RELCHI10BYEAR RELCHI11BYEAR RELCHI12BYEAR ///
    RELCHI13BYEAR RELCHI14BYEAR RELCHI15BYEAR RELCHI16BYEAR RELCHI17BYEAR ///
    RELCHI18BYEAR RELCHI19BYEAR RELCHI20BYEAR using "`dta'", clear

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

forvalues j = 1/20 {
    local source = cond(`j' == 1, "birth1", cond(`j' == 2, "birth2", cond(`j' == 3, "birth3", "RELCHI`j'BYEAR")))
    gen byte same_first_`j' = (`source' == birth1) if !missing(`source', birth1)
    replace same_first_`j' = 0 if missing(same_first_`j')
}
egen n_sameyear_first = rowtotal(same_first_1-same_first_20)
gen byte twin_first_proxy = n_sameyear_first >= 2 if !missing(birth1)
drop same_first_1-same_first_20

gen byte same_sex_first2 = child1_sex == child2_sex if !missing(child1_sex, child2_sex)
gen byte first_child_girl = child1_sex == 2 if !missing(child1_sex)
gen byte singleton_first2 = birth2 > birth1 if !missing(birth1, birth2)
gen byte singleton_second = birth3 != birth2 if !missing(birth2)
replace singleton_second = 1 if missing(birth3) & !missing(birth2)
gen byte two_plus_by5 = !missing(birth2) & birth2 <= birth1 + 5 if !missing(birth1)
gen byte three_plus_by5 = !missing(birth3) & birth3 <= birth2 + 5 if !missing(birth2)

keep ID birth1 birth2 birth3 twin_first_proxy same_sex_first2 first_child_girl ///
    singleton_first2 singleton_second two_plus_by5 three_plus_by5
isid ID
save `parents', replace

di as text "Step 2/5: Build clean person-year housing panel"
use ID year AGEREP EDUYEAR SEX RACE DEATHYEAR HOMEOWN MOVEDFREF_ WHYMOVED1_ ///
    ACTUALROOMS_ IW using "`dta'", clear

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
replace rooms = . if !inrange(rooms, 1, 20)

sort ID year
by ID: gen prior_observed_own = own[_n-1]
gen byte moved_to_own = .
replace moved_to_own = (own == 1 & prior_observed_own == 0 & moved_since_interview == 1) ///
    if !missing(own, prior_observed_own, moved_since_interview)

gen byte moved_for_size = .
replace moved_for_size = 0 if moved_since_interview == 0
replace moved_for_size = (WHYMOVED1_ == 3) if moved_since_interview == 1 & !missing(WHYMOVED1_)

merge m:1 ID using `parents', keep(match) nogen
bysort ID: egen race_mode = mode(RACE), minmode
replace race_mode = 0 if missing(race_mode) | inlist(race_mode, 8, 9)
save `panel', replace

di as text "Step 3/5: Construct post-only first-birth twins sample"
use `panel', clear
keep if SEX == 2 & !missing(birth1, twin_first_proxy, two_plus_by5)
gen event_time = year - birth1
gen age_at_event_i = AGEREP - event_time if !missing(AGEREP, event_time)
bysort ID: egen age_at_event = median(age_at_event_i)
gen age_event_round = round(age_at_event)

gen rooms_i = rooms if inrange(event_time, 0, 5)
gen own_i = own if inrange(event_time, 0, 5)
gen moved_size_i = moved_for_size if inrange(event_time, 0, 5)
gen moved_own_i = moved_to_own if inrange(event_time, 0, 5)
gen observed_post_i = event_time if inrange(event_time, 0, 5)
gen end_weight_year_i = year if inrange(event_time, 4, 5) & !missing(IW)

bysort ID: egen rooms_mean_post5 = mean(rooms_i)
bysort ID: egen own_mean_post5 = mean(own_i)
bysort ID: egen moved_for_size_post5 = max(moved_size_i)
bysort ID: egen moved_to_own_post5 = max(moved_own_i)
bysort ID: egen latest_post_k = max(observed_post_i)
bysort ID: egen end_weight_year = max(end_weight_year_i)
bysort ID: egen iw_end = max(cond(year == end_weight_year, IW, .))

bysort ID (year): keep if _n == 1
keep if birth1 <= 2014
keep if latest_post_k >= 4
keep if inrange(age_event_round, 15, 50)
keep if !missing(iw_end) & iw_end > 0
gen str12 design = "twins"
gen byte z = twin_first_proxy
gen byte treatment = two_plus_by5
gen event_year = birth1
keep ID design z treatment event_year age_at_event age_event_round race_mode first_child_girl iw_end ///
    rooms_mean_post5 own_mean_post5 moved_for_size_post5 moved_to_own_post5
isid ID
save `twins_sample', replace
save "`outdir'/twins_postonly_sample.dta", replace

di as text "Step 4/5: Construct post-only second-birth same-sex sample"
use `panel', clear
keep if SEX == 2 & !missing(birth2, same_sex_first2, first_child_girl, three_plus_by5)
keep if singleton_first2 == 1 & singleton_second == 1
gen event_time = year - birth2
gen age_at_event_i = AGEREP - event_time if !missing(AGEREP, event_time)
bysort ID: egen age_at_event = median(age_at_event_i)
gen age_event_round = round(age_at_event)

gen rooms_i = rooms if inrange(event_time, 0, 5)
gen own_i = own if inrange(event_time, 0, 5)
gen moved_size_i = moved_for_size if inrange(event_time, 0, 5)
gen moved_own_i = moved_to_own if inrange(event_time, 0, 5)
gen observed_post_i = event_time if inrange(event_time, 0, 5)
gen end_weight_year_i = year if inrange(event_time, 4, 5) & !missing(IW)

bysort ID: egen rooms_mean_post5 = mean(rooms_i)
bysort ID: egen own_mean_post5 = mean(own_i)
bysort ID: egen moved_for_size_post5 = max(moved_size_i)
bysort ID: egen moved_to_own_post5 = max(moved_own_i)
bysort ID: egen latest_post_k = max(observed_post_i)
bysort ID: egen end_weight_year = max(end_weight_year_i)
bysort ID: egen iw_end = max(cond(year == end_weight_year, IW, .))

bysort ID (year): keep if _n == 1
keep if birth2 <= 2014
keep if latest_post_k >= 4
keep if inrange(age_event_round, 15, 50)
keep if !missing(iw_end) & iw_end > 0
gen str12 design = "same_sex"
gen byte z = same_sex_first2
gen byte treatment = three_plus_by5
gen event_year = birth2
keep ID design z treatment event_year age_at_event age_event_round race_mode first_child_girl iw_end ///
    rooms_mean_post5 own_mean_post5 moved_for_size_post5 moved_to_own_post5
isid ID
save `samesex_sample', replace
save "`outdir'/samesex_postonly_sample.dta", replace

di as text "Step 5/5: Estimate first stages, reduced forms, OLS, and 2SLS"
tempname estimates samplepost balancepost
postfile `estimates' str12 design str28 specification str28 outcome str10 weighting ///
    str18 estimand double b se p ci_low ci_high N first_stage_F ///
    using "`outdir'/iv_housing_postonly_estimates.dta", replace
postfile `samplepost' str12 design str28 outcome double N z1 z_share treatment_share outcome_mean ///
    using "`outdir'/iv_housing_postonly_samples.dta", replace
postfile `balancepost' str12 design str20 covariate str10 weighting double b se p N ///
    using "`outdir'/iv_housing_postonly_balance.dta", replace

foreach d in twins same_sex {
    if "`d'" == "twins" use `twins_sample', clear
    if "`d'" == "same_sex" use `samesex_sample', clear

    foreach y in rooms_mean_post5 own_mean_post5 moved_for_size_post5 moved_to_own_post5 {
        preserve
            keep if !missing(`y')

            quietly count
            local n_sample = r(N)
            quietly count if z == 1
            local n_z1 = r(N)
            quietly summarize z [aw=iw_end], meanonly
            local zbar = r(mean)
            quietly summarize treatment [aw=iw_end], meanonly
            local dbar = r(mean)
            quietly summarize `y' [aw=iw_end], meanonly
            local ybar = r(mean)
            post `samplepost' ("`d'") ("`y'") (`n_sample') (`n_z1') (`zbar') (`dbar') (`ybar')

            foreach controlspec in age_fe_race age_quadratic_race {
                local controls "i.age_event_round i.event_year i.race_mode"
                if "`controlspec'" == "age_quadratic_race" ///
                    local controls "c.age_at_event##c.age_at_event i.event_year i.race_mode"
                if "`d'" == "same_sex" local controls "`controls' i.first_child_girl"

                foreach w in weighted unweighted {
                    local wt ""
                    if "`w'" == "weighted" local wt "[aw=iw_end]"

                    quietly reg treatment z `controls' `wt', vce(robust)
                    quietly test z
                    local firstF = r(F)
                    local bb = _b[z]
                    local fs_b = `bb'
                    local ss = _se[z]
                    local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
                    local cc = invttail(e(df_r), .025)
                    post `estimates' ("`d'") ("`controlspec'") ("`y'") ("`w'") ///
                        ("first_stage") (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') ///
                        (e(N)) (`firstF')

                    quietly reg `y' z `controls' `wt', vce(robust)
                    local bb = _b[z]
                    local rf_b = `bb'
                    local ss = _se[z]
                    local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
                    local cc = invttail(e(df_r), .025)
                    post `estimates' ("`d'") ("`controlspec'") ("`y'") ("`w'") ///
                        ("reduced_form") (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') ///
                        (e(N)) (`firstF')

                    quietly reg `y' treatment `controls' `wt', vce(robust)
                    local bb = _b[treatment]
                    local ss = _se[treatment]
                    local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
                    local cc = invttail(e(df_r), .025)
                    post `estimates' ("`d'") ("`controlspec'") ("`y'") ("`w'") ///
                        ("ols") (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') ///
                        (e(N)) (`firstF')

                    quietly ivregress 2sls `y' (treatment = z) `controls' `wt', vce(robust)
                    local bb = _b[treatment]
                    local ss = _se[treatment]
                    assert abs(`bb' - (`rf_b'/`fs_b')) < 1e-6
                    local pp = 2 * normal(-abs(`bb'/`ss'))
                    local cc = invnormal(.975)
                    post `estimates' ("`d'") ("`controlspec'") ("`y'") ("`w'") ///
                        ("iv_2sls") (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') ///
                        (e(N)) (`firstF')
                }
            }
        restore
    }

    foreach x in age_at_event first_child_girl {
        if "`d'" == "same_sex" & "`x'" == "first_child_girl" continue
        foreach w in weighted unweighted {
            local wt ""
            if "`w'" == "weighted" local wt "[aw=iw_end]"
            quietly reg `x' z i.event_year `wt', vce(robust)
            local bb = _b[z]
            local ss = _se[z]
            local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
            post `balancepost' ("`d'") ("`x'") ("`w'") (`bb') (`ss') (`pp') (e(N))
        }
    }
}

postclose `estimates'
postclose `samplepost'
postclose `balancepost'

use "`outdir'/iv_housing_postonly_estimates.dta", clear
assert _N == 128
assert !missing(b, se, p, N, first_stage_F)
sort design outcome weighting estimand
export delimited using "`outdir'/iv_housing_postonly_estimates.csv", replace

use "`outdir'/iv_housing_postonly_samples.dta", clear
assert _N == 8
assert !missing(N, z1, z_share, treatment_share, outcome_mean)
sort design outcome
export delimited using "`outdir'/iv_housing_postonly_samples.csv", replace

use "`outdir'/iv_housing_postonly_balance.dta", clear
assert _N == 6
assert !missing(b, se, p, N)
sort design covariate weighting
export delimited using "`outdir'/iv_housing_postonly_balance.csv", replace

di as result "Post-only fertility-IV diagnostic complete."
di as result "Output directory: `outdir'"
log close _all
