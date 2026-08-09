clear all
set more off
version 17.0

local project "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`project'/code/data/psid_followup_mar2026/output/iv_housing_panel_reaudit_20260809"

capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/iv_housing_panel_reaudit_20260809.log", replace text

capture which reghdfe
if _rc {
    di as error "reghdfe is required."
    exit 198
}
capture which xtivreg2
if _rc {
    di as error "xtivreg2 is required."
    exit 198
}

tempfile sex_map sex1_map sex2_map parents panel twins_panel samesex_panel stacked

di as text "Step 1/6: Build child-sex map and final fertility roster"
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

keep ID birth1 birth2 birth3 twin_first_proxy same_sex_first2 first_child_girl ///
    singleton_first2 singleton_second
isid ID
save `parents', replace

di as text "Step 2/6: Build clean person-year housing panel"
use ID year AGEREP EDUYEAR SEX DEATHYEAR HOMEOWN MOVEDFREF_ WHYMOVED1_ ///
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
by ID: gen prior_observed_year = year[_n-1]
gen interview_gap = year - prior_observed_year if !missing(prior_observed_year)

gen byte moved_to_own = .
replace moved_to_own = (own == 1 & prior_observed_own == 0 & moved_since_interview == 1) ///
    if !missing(own, prior_observed_own, moved_since_interview)

gen byte moved_for_size = .
replace moved_for_size = 0 if moved_since_interview == 0
replace moved_for_size = (WHYMOVED1_ == 3) if moved_since_interview == 1 & !missing(WHYMOVED1_)

merge m:1 ID using `parents', keep(match) nogen
save `panel', replace

di as text "Step 3/6: Build first-birth twins panel"
use `panel', clear
keep if SEX == 2 & !missing(birth1, twin_first_proxy)
gen event_time = year - birth1
keep if inrange(event_time, -4, 5)
gen str12 design = "twins"
gen byte z = twin_first_proxy
gen byte additional_child = !missing(birth2) & birth2 <= year
gen event_year = birth1

bysort ID: egen baseline_year = max(cond(event_time < 0 & ///
    !missing(AGEREP, EDUYEAR, own, IW), year, .))
bysort ID: egen age_pre = max(cond(year == baseline_year, AGEREP, .))
bysort ID: egen eduyear_pre = max(cond(year == baseline_year, EDUYEAR, .))
bysort ID: egen own_pre = max(cond(year == baseline_year, own, .))
bysort ID: egen iw_pre = max(cond(year == baseline_year, IW, .))
bysort ID: egen rooms_baseline_year = max(cond(event_time < 0 & !missing(rooms), year, .))
bysort ID: egen rooms_pre = max(cond(year == rooms_baseline_year, rooms, .))
save `twins_panel', replace

di as text "Step 4/6: Build second-birth same-sex panel"
use `panel', clear
keep if SEX == 2 & !missing(birth2, same_sex_first2, first_child_girl)
keep if singleton_first2 == 1 & singleton_second == 1
gen event_time = year - birth2
keep if inrange(event_time, -4, 5)
gen str12 design = "same_sex"
gen byte z = same_sex_first2
gen byte additional_child = !missing(birth3) & birth3 <= year
gen event_year = birth2

bysort ID: egen baseline_year = max(cond(event_time < 0 & ///
    !missing(AGEREP, EDUYEAR, own, IW), year, .))
bysort ID: egen age_pre = max(cond(year == baseline_year, AGEREP, .))
bysort ID: egen eduyear_pre = max(cond(year == baseline_year, EDUYEAR, .))
bysort ID: egen own_pre = max(cond(year == baseline_year, own, .))
bysort ID: egen iw_pre = max(cond(year == baseline_year, IW, .))
bysort ID: egen rooms_baseline_year = max(cond(event_time < 0 & !missing(rooms), year, .))
bysort ID: egen rooms_pre = max(cond(year == rooms_baseline_year, rooms, .))
save `samesex_panel', replace

use `twins_panel', clear
append using `samesex_panel'

gen byte post = event_time >= 0
gen byte z_post = z * post
gen age_pre_post = age_pre * post
gen age2_pre_post = age_pre^2 * post
gen eduyear_pre_post = eduyear_pre * post
gen own_pre_post = own_pre * post
gen rooms_pre_post = rooms_pre * post
gen first_child_girl_post = first_child_girl * post
gen byte gap_one_year = interview_gap == 1 & !missing(interview_gap)
gen byte gap_three_plus = interview_gap >= 3 & !missing(interview_gap)
gen byte gap_missing = missing(interview_gap)

foreach suffix in m4 m3 m1 p0 p1 p2 p3 p4 p5 {
    local kval = .
    if "`suffix'" == "m4" local kval = -4
    if "`suffix'" == "m3" local kval = -3
    if "`suffix'" == "m1" local kval = -1
    if "`suffix'" == "p0" local kval = 0
    if "`suffix'" == "p1" local kval = 1
    if "`suffix'" == "p2" local kval = 2
    if "`suffix'" == "p3" local kval = 3
    if "`suffix'" == "p4" local kval = 4
    if "`suffix'" == "p5" local kval = 5
    gen byte ev_`suffix' = event_time == `kval'
    gen byte zev_`suffix' = z * ev_`suffix'
}
quietly tabulate year, generate(year_fe_)
drop year_fe_1
unab year_controls : year_fe_*

label var z_post "Instrument x post-birth"
label var additional_child "Additional-child status at interview"
save `stacked', replace
save "`outdir'/iv_housing_panel_analysis_sample.dta", replace

di as text "Step 5/6: Estimate average panel reduced forms and IV"
tempname estimates samplepost eventpost testpost
postfile `estimates' str12 design str28 specification str24 outcome str10 weighting ///
    str22 estimand double b se p ci_low ci_high N_obs N_ids z1_ids first_stage_F ///
    using "`outdir'/iv_housing_panel_estimates.dta", replace
postfile `samplepost' str12 design str24 outcome double N_obs N_ids z1_ids ///
    pre_obs post_obs using "`outdir'/iv_housing_panel_samples.dta", replace
postfile `eventpost' str12 design str24 outcome double event_time b se p ci_low ci_high ///
    N_obs N_ids z1_ids using "`outdir'/iv_housing_panel_eventstudy.dta", replace
postfile `testpost' str12 design str24 outcome str16 test double statistic df_num df_den p ///
    using "`outdir'/iv_housing_panel_tests.dta", replace

local common_events "ev_m4 ev_m3 ev_m1 ev_p0 ev_p1 ev_p2 ev_p3 ev_p4 ev_p5"
local z_events "zev_m4 zev_m3 zev_m1 zev_p0 zev_p1 zev_p2 zev_p3 zev_p4 zev_p5"

foreach d in twins same_sex {
    foreach y in rooms moved_for_size own moved_to_own {
        use `stacked', clear
        keep if design == "`d'"
        if inlist("`y'", "own", "moved_to_own") keep if own_pre == 0
        keep if !missing(`y', additional_child, z, age_pre, eduyear_pre, own_pre, iw_pre)
        keep if iw_pre > 0
        if "`y'" == "rooms" keep if !missing(rooms_pre)

        bysort ID: egen has_pre = max(event_time < 0)
        bysort ID: egen has_post = max(event_time >= 0)
        keep if has_pre == 1 & has_post == 1
        drop has_pre has_post

        quietly count
        local n_obs = r(N)
        quietly count if event_time < 0
        local pre_obs = r(N)
        quietly count if event_time >= 0
        local post_obs = r(N)
        egen byte id_tag = tag(ID)
        quietly count if id_tag == 1
        local n_ids = r(N)
        quietly count if id_tag == 1 & z == 1
        local n_z1 = r(N)
        drop id_tag
        post `samplepost' ("`d'") ("`y'") (`n_obs') (`n_ids') (`n_z1') (`pre_obs') (`post_obs')

        local controls "`common_events' age_pre_post age2_pre_post eduyear_pre_post own_pre_post"
        if "`d'" == "same_sex" local controls "`controls' first_child_girl_post"
        if "`y'" == "rooms" local controls "`controls' rooms_pre_post"
        if inlist("`y'", "moved_for_size", "moved_to_own") ///
            local controls "`controls' gap_one_year gap_three_plus gap_missing"

        xtset ID year
        foreach w in weighted unweighted {
            local wt ""
            if "`w'" == "weighted" local wt "[aw=iw_pre]"

            quietly reghdfe additional_child z_post `controls' `wt', absorb(ID year) vce(cluster ID)
            quietly test z_post
            local firstF = r(F)
            local bb = _b[z_post]
            local fs_b = `bb'
            local ss = _se[z_post]
            local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
            local cc = invttail(e(df_r), .025)
            post `estimates' ("`d'") ("panel_average_idfe") ("`y'") ("`w'") ("first_stage") ///
                (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') ///
                (e(N)) (`n_ids') (`n_z1') (`firstF')

            quietly reghdfe `y' z_post `controls' `wt', absorb(ID year) vce(cluster ID)
            local bb = _b[z_post]
            local rf_b = `bb'
            local ss = _se[z_post]
            local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
            local cc = invttail(e(df_r), .025)
            post `estimates' ("`d'") ("panel_average_idfe") ("`y'") ("`w'") ("reduced_form") ///
                (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') ///
                (e(N)) (`n_ids') (`n_z1') (`firstF')

            quietly reghdfe `y' additional_child `controls' `wt', absorb(ID year) vce(cluster ID)
            local bb = _b[additional_child]
            local ss = _se[additional_child]
            local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
            local cc = invttail(e(df_r), .025)
            post `estimates' ("`d'") ("panel_average_idfe") ("`y'") ("`w'") ("ols_fe") ///
                (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') ///
                (e(N)) (`n_ids') (`n_z1') (`firstF')

            quietly xtivreg2 `y' `controls' `year_controls' (additional_child = z_post) `wt', ///
                fe cluster(ID)
            local bb = _b[additional_child]
            local ss = _se[additional_child]
            assert abs(`bb' - (`rf_b'/`fs_b')) < 1e-6
            local pp = 2 * normal(-abs(`bb'/`ss'))
            local cc = invnormal(.975)
            post `estimates' ("`d'") ("panel_average_idfe") ("`y'") ("`w'") ("iv_2sls") ///
                (`bb') (`ss') (`pp') (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') ///
                (e(N)) (`n_ids') (`n_z1') (`firstF')
        }

        * Weighted dynamic reduced form, relative to event time -2.
        local dynamic_controls "age_pre_post age2_pre_post eduyear_pre_post own_pre_post"
        if "`d'" == "same_sex" local dynamic_controls "`dynamic_controls' first_child_girl_post"
        if "`y'" == "rooms" local dynamic_controls "`dynamic_controls' rooms_pre_post"
        if inlist("`y'", "moved_for_size", "moved_to_own") ///
            local dynamic_controls "`dynamic_controls' gap_one_year gap_three_plus gap_missing"
        quietly reghdfe `y' `z_events' `common_events' `dynamic_controls' [aw=iw_pre], ///
            absorb(ID year) vce(cluster ID)

        foreach suffix in m4 m3 m1 p0 p1 p2 p3 p4 p5 {
            local kval = .
            if "`suffix'" == "m4" local kval = -4
            if "`suffix'" == "m3" local kval = -3
            if "`suffix'" == "m1" local kval = -1
            if "`suffix'" == "p0" local kval = 0
            if "`suffix'" == "p1" local kval = 1
            if "`suffix'" == "p2" local kval = 2
            if "`suffix'" == "p3" local kval = 3
            if "`suffix'" == "p4" local kval = 4
            if "`suffix'" == "p5" local kval = 5
            local bb = _b[zev_`suffix']
            local ss = _se[zev_`suffix']
            local pp = 2 * ttail(e(df_r), abs(`bb'/`ss'))
            local cc = invttail(e(df_r), .025)
            post `eventpost' ("`d'") ("`y'") (`kval') (`bb') (`ss') (`pp') ///
                (`bb'-`cc'*`ss') (`bb'+`cc'*`ss') (e(N)) (`n_ids') (`n_z1')
        }
        post `eventpost' ("`d'") ("`y'") (-2) (0) (0) (.) (0) (0) ///
            (e(N)) (`n_ids') (`n_z1')

        quietly test zev_m4 zev_m3 zev_m1
        post `testpost' ("`d'") ("`y'") ("pre_joint") (r(F)) (r(df)) (r(df_r)) (r(p))
        quietly test zev_p0 zev_p1 zev_p2 zev_p3 zev_p4 zev_p5
        post `testpost' ("`d'") ("`y'") ("post_joint") (r(F)) (r(df)) (r(df_r)) (r(p))
    }
}

postclose `estimates'
postclose `samplepost'
postclose `eventpost'
postclose `testpost'

use "`outdir'/iv_housing_panel_estimates.dta", clear
assert _N == 64
assert !missing(b, se, p, N_obs, N_ids, z1_ids, first_stage_F)
sort design outcome weighting estimand
export delimited using "`outdir'/iv_housing_panel_estimates.csv", replace

use "`outdir'/iv_housing_panel_samples.dta", clear
assert _N == 8
assert !missing(N_obs, N_ids, z1_ids, pre_obs, post_obs)
sort design outcome
export delimited using "`outdir'/iv_housing_panel_samples.csv", replace

use "`outdir'/iv_housing_panel_tests.dta", clear
assert _N == 16
assert !missing(statistic, df_num, df_den, p)
sort design outcome test
export delimited using "`outdir'/iv_housing_panel_tests.csv", replace

di as text "Step 6/6: Export weighted reduced-form event-study paths"
use "`outdir'/iv_housing_panel_eventstudy.dta", clear
assert _N == 80
assert !missing(b, se, N_obs, N_ids, z1_ids)
sort design outcome event_time
export delimited using "`outdir'/iv_housing_panel_eventstudy.csv", replace

foreach d in twins same_sex {
    foreach y in rooms moved_for_size own moved_to_own {
        preserve
            keep if design == "`d'" & outcome == "`y'"
            sort event_time
            twoway ///
                (rcap ci_low ci_high event_time if se > 0, lcolor(navy%35)) ///
                (connected b event_time, mcolor(navy) lcolor(navy) msymbol(circle)), ///
                yline(0, lcolor(gs10)) xline(-1.5, lpattern(dash) lcolor(gs10)) ///
                xtitle("Years relative to event birth") ytitle("Instrument differential") ///
                title("`d': `y' reduced form") legend(off) ///
                graphregion(fcolor(white)) plotregion(fcolor(white))
            graph export "`outdir'/eventstudy_`d'_`y'.png", replace width(1600)
        restore
    }
}

di as result "Panel fertility-IV diagnostic complete."
di as result "Output directory: `outdir'"
log close _all
