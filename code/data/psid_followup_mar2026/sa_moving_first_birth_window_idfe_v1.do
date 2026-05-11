clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/sa_moving_first_birth_window_idfe_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/sa_moving_first_birth_window_idfe_v1.log", replace text

global weight IW
global ref_sa 2
global bg "250 250 250"
global graphregion_slides graphregion(fcolor("$bg") lcolor("$bg")) plotregion(fcolor("$bg") lcolor("$bg"))

cap which eventstudyinteract
if _rc {
    di as error "eventstudyinteract is not installed. Aborting."
    exit 198
}

use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR HOMEOWN ///
    MOVEDFREF_ WHYMOVED1_ ${weight} using "`dta'", clear

drop if year > DEATHYEAR
drop if AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3
rename HOMEOWN own

rename MOVEDFREF_ movedthisyear
replace movedthisyear = . if movedthisyear == 8 | movedthisyear == 9
replace movedthisyear = 0 if movedthisyear == 5

gen first_child_year = RELCHI1BYEAR
drop if missing(first_child_year)

bysort ID: egen year_entry = min(year)
drop if first_child_year < year_entry
drop year_entry

xtset ID year

gen move_any = movedthisyear
gen moved_for_size = .
replace moved_for_size = 1 if movedthisyear == 1 & WHYMOVED1_ == 3 & year > 1974
replace moved_for_size = 0 if movedthisyear == 0 & year > 1974
replace moved_for_size = 0 if movedthisyear == 1 & WHYMOVED1_ != 3 & !missing(WHYMOVED1_) & year > 1974

gen change_to_own = (own == 1 & L.own == 0) if !missing(own, L.own)
gen moved_to_own = (movedthisyear == 1 & change_to_own == 1) if !missing(movedthisyear, change_to_own)

gen K = year - first_child_year
keep if inrange(K, -5, 5)

sum first_child_year
gen lastcohort = first_child_year == r(max)

cap drop L*event F*event
forvalues l = 0/5 {
    gen L`l'event = K == `l'
}
forvalues l = 1/5 {
    gen F`l'event = K == -`l'
}
drop F${ref_sa}event

quietly count
local base_obs = r(N)
egen tag_id = tag(ID)
quietly count if tag_id == 1
local base_ids = r(N)
drop tag_id
levelsof first_child_year if lastcohort == 0, local(cohorts)
local n_cohorts : word count `cohorts'
di as text "Base finite-window sample: obs=`base_obs', ids=`base_ids', non-control cohorts=`n_cohorts'"

tempfile all_estimates summary status
tempname ph sh sth

postfile `ph' str32 outcome str32 label str32 coeff double relative_time b se ci_lo ci_hi ///
    double pre_event_mean post03_raw_mean sample_obs sample_ids base_obs base_ids n_cohorts ///
    using "`all_estimates'", replace

postfile `sh' str32 outcome str32 label double pre_event_mean post03_raw_mean ///
    double coef_m1 se_m1 coef_0 se_0 coef_1 se_1 coef_2 se_2 coef_3 se_3 ///
    double sample_obs sample_ids using "`summary'", replace

postfile `sth' str32 outcome str32 label int rc double runtime_seconds str120 message ///
    using "`status'", replace

local outcomes "move_any moved_for_size moved_to_own"

foreach y of local outcomes {
    local ytitle "`y'"
    if "`y'" == "move_any" local ytitle "Any move"
    if "`y'" == "moved_for_size" local ytitle "Moved for size"
    if "`y'" == "moved_to_own" local ytitle "Moved to ownership"

    quietly count if !missing(`y', first_child_year, year, AGEREP, EDUYEAR, ID)
    local input_obs = r(N)
    capture drop tag_id
    egen tag_id = tag(ID) if !missing(`y', first_child_year, year, AGEREP, EDUYEAR, ID)
    quietly count if tag_id == 1
    local input_ids = r(N)
    drop tag_id

    di as text "Running IDFE Sun-Abraham for `y' (`ytitle'): obs=`input_obs', ids=`input_ids'"

    timer clear 10
    timer on 10
    capture noisily eventstudyinteract `y' L*event F*event, ///
        vce(cluster ID) absorb(ID year) cohort(first_child_year) control_cohort(lastcohort) ///
        covariates(i.AGEREP i.EDUYEAR)
    local rc = _rc
    timer off 10
    quietly timer list 10
    local runtime_seconds = r(t10)

    if `rc' {
        post `sth' ("`y'") ("`ytitle'") (`rc') (`runtime_seconds') ("eventstudyinteract failed")
        di as error "Event study failed for `y' with rc=`rc'. Continuing."
        continue
    }

    quietly summarize `y' if e(sample) & K == -${ref_sa}
    local pre_event_mean = r(mean)
    quietly summarize `y' if e(sample) & inrange(K, 0, 3)
    local post03_raw_mean = r(mean)
    quietly count if e(sample)
    local sample_obs = r(N)
    capture drop tag_id
    egen tag_id = tag(ID) if e(sample)
    quietly count if tag_id == 1
    local sample_ids = r(N)
    drop tag_id

    matrix bmat = e(b_iw)
    matrix vmat = e(V_iw)
    local cn : colnames bmat

    local coef_m1 = .
    local se_m1 = .
    forvalues rt = 0/3 {
        local coef_`rt' = .
        local se_`rt' = .
    }

    local j = 0
    foreach c of local cn {
        local ++j
        local bb = bmat[1, `j']
        local vv = vmat[`j', `j']
        local ss = sqrt(`vv')

        local cclean "`c'"
        local cclean : subinstr local cclean "event" "", all
        local cclean : subinstr local cclean "L" "", all
        local cclean : subinstr local cclean "F" "-", all
        local cclean : subinstr local cclean "o." "", all
        local cclean : subinstr local cclean "b." "", all
        local reltime = real("`cclean'")
        if missing(`reltime') continue
        if `bb' == 0 & `ss' == 0 continue

        local lo = `bb' - 1.96 * `ss'
        local hi = `bb' + 1.96 * `ss'
        post `ph' ("`y'") ("`ytitle'") ("`c'") (`reltime') (`bb') (`ss') (`lo') (`hi') ///
            (`pre_event_mean') (`post03_raw_mean') (`sample_obs') (`sample_ids') ///
            (`base_obs') (`base_ids') (`n_cohorts')

        if `reltime' == -1 {
            local coef_m1 = `bb'
            local se_m1 = `ss'
        }
        forvalues rt = 0/3 {
            if `reltime' == `rt' {
                local coef_`rt' = `bb'
                local se_`rt' = `ss'
            }
        }
    }

    post `ph' ("`y'") ("`ytitle'") ("F${ref_sa}event") (-${ref_sa}) (0) (0) (0) (0) ///
        (`pre_event_mean') (`post03_raw_mean') (`sample_obs') (`sample_ids') ///
        (`base_obs') (`base_ids') (`n_cohorts')

    post `sh' ("`y'") ("`ytitle'") (`pre_event_mean') (`post03_raw_mean') ///
        (`coef_m1') (`se_m1') (`coef_0') (`se_0') (`coef_1') (`se_1') ///
        (`coef_2') (`se_2') (`coef_3') (`se_3') (`sample_obs') (`sample_ids')

    post `sth' ("`y'") ("`ytitle'") (0) (`runtime_seconds') ("completed")
}

postclose `ph'
postclose `sh'
postclose `sth'

use "`all_estimates'", clear
gen b_pct = 100 * b
gen ci_lo_pct = 100 * ci_lo
gen ci_hi_pct = 100 * ci_hi
sort outcome relative_time
save "`outdir'/moving_first_birth_sa_window_idfe_all_estimates_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_sa_window_idfe_all_estimates_v1.csv", replace

foreach y of local outcomes {
    preserve
        keep if outcome == "`y'"
        count
        if r(N) > 0 {
            local ytitle = label[1]
            sort relative_time
            export delimited using "`outdir'/`y'_sa_first_birth_window_idfe_estimates_v1.csv", replace
            save "`outdir'/`y'_sa_first_birth_window_idfe_estimates_v1.dta", replace

            twoway ///
                (rarea ci_lo_pct ci_hi_pct relative_time, fcolor(dkgreen%10) lcolor(dkgreen%10)) ///
                (connected b_pct relative_time, mcolor(dkgreen) msymbol(circle) msize(medlarge) ///
                    lcolor(dkgreen) lwidth(medthick) lpattern(dash)), ///
                ytitle("Effect on `ytitle' probability, pp") ///
                xtitle("Years relative to first birth") ///
                xlabel(-5(1)5) ///
                xline(-1.5, lwidth(thin) lpattern(dash) lcolor(gray)) ///
                yline(0, lwidth(vthin) lpattern(solid) lcolor(gray)) ///
                title("`ytitle' Around First Birth") ///
                note("Sun-Abraham IW; household FE + year FE; K=-2 omitted. 95% CI clustered by household.", size(vsmall)) ///
                $graphregion_slides legend(off)
            graph export "`outdir'/`y'_sa_first_birth_window_idfe_v1.png", replace width(1600)
            graph export "`outdir'/`y'_sa_first_birth_window_idfe_v1.pdf", replace
        }
    restore
}

use "`summary'", clear
save "`outdir'/moving_first_birth_sa_window_idfe_summary_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_sa_window_idfe_summary_v1.csv", replace

use "`status'", clear
save "`outdir'/moving_first_birth_sa_window_idfe_status_v1.dta", replace
export delimited using "`outdir'/moving_first_birth_sa_window_idfe_status_v1.csv", replace

di as result "Finite-window IDFE Sun-Abraham moving event studies complete."
di as result "Output dir: `outdir'"

log close _all
