clear all
set more off
set processors 4
version 17.0

* Corrected first-birth housing-space target.
* Unit: one weighted woman (reference person or spouse/partner) per
* single-family-unit household-year.
* Outcome: ACTUALROOMS_ shifted forward by one observed interview within person,
* because the PSIDSHELF extraction attaches the next-wave rooms item to the
* preceding row. Event-study baseline is K=-2; controls are women whose full
* relationship history confirms no children.

local project "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local source  "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local outroot "`project'/output/model/native_financing_diagnostic_20260919/specification_followup/target_review_v1/overnight/empirical_rooms/covariance_replay"
local outdir  "`outroot'/output"

cap mkdir "`outroot'"
cap mkdir "`outdir'"
capture confirm file "`outdir'/target_receipt.csv"
if !_rc {
    di as error "Refusing to overwrite completed target receipt: `outdir'/target_receipt.csv"
    exit 602
}

capture log close _all
log using "`outdir'/sa_rooms_first_birth_household_aligned_v1.log", replace text
timer clear 1
timer on 1

cap which eventstudyinteract
if _rc {
    di as error "eventstudyinteract is not installed. Aborting."
    exit 198
}
cap which svmat2
if _rc {
    di as error "svmat2 is not installed. Aborting."
    exit 198
}

di as text "Loading the PSID shelf and aligning ACTUALROOMS_ to its interview..."
use ID FID HHID year AGEREP EDUYEAR SEX DEATHYEAR REL CURRENT RELCHIREP ///
    RELCHI1TYPE RELCHI1BYEAR RELCHI2TYPE RELCHI2BYEAR ///
    RELCHI3TYPE RELCHI3BYEAR RELCHI4TYPE RELCHI4BYEAR ///
    RELCHI5TYPE RELCHI5BYEAR RELCHI6TYPE RELCHI6BYEAR ///
    RELCHI7TYPE RELCHI7BYEAR RELCHI8TYPE RELCHI8BYEAR ///
    RELCHI9TYPE RELCHI9BYEAR RELCHI10TYPE RELCHI10BYEAR ///
    RELCHI11TYPE RELCHI11BYEAR RELCHI12TYPE RELCHI12BYEAR ///
    RELCHI13TYPE RELCHI13BYEAR RELCHI14TYPE RELCHI14BYEAR ///
    RELCHI15TYPE RELCHI15BYEAR RELCHI16TYPE RELCHI16BYEAR ///
    RELCHI17TYPE RELCHI17BYEAR RELCHI18TYPE RELCHI18BYEAR ///
    RELCHI19TYPE RELCHI19BYEAR RELCHI20TYPE RELCHI20BYEAR ///
    ACTUALROOMS_ IW using "`source'", clear

sort ID year
by ID: gen double rooms = ACTUALROOMS_[_n-1] if _n > 1
by ID: gen double rooms_source_year = year[_n-1] if _n > 1

* Diagnostic-only flags record raw aligned room codes; they never filter rows.
gen byte room_code_0 = rooms == 0 if !missing(rooms)
gen byte room_code_9 = rooms == 9 if !missing(rooms)
gen byte room_code_98 = rooms == 98 if !missing(rooms)
gen byte room_code_99 = rooms == 99 if !missing(rooms)

* Fertility history is a person-level object and is constructed before any
* sample restriction or household reporter selection. TYPE==1 is biological.
gen double first_biological_birth_candidate = .
forvalues child = 1/20 {
    replace first_biological_birth_candidate = RELCHI`child'BYEAR ///
        if RELCHI`child'TYPE == 1 & !missing(RELCHI`child'BYEAR) & ///
        (missing(first_biological_birth_candidate) | ///
         RELCHI`child'BYEAR < first_biological_birth_candidate)
}
bysort ID: egen double first_child_year = min(first_biological_birth_candidate)
bysort ID: egen double max_children_reported = max(RELCHIREP)
drop first_biological_birth_candidate

* HHID is a physical dwelling and can contain several PSID family units.
* Count current FIDs using every current member before selecting women.
egen byte current_fid_tag = tag(HHID year FID) ///
    if CURRENT == 1 & !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0
bysort HHID year: egen int n_current_fids = total(current_fid_tag)
drop current_fid_tag

* Missing/non-room codes follow the year to which the lagged value is aligned.
replace rooms = . if year <= 1984 & rooms == 9
replace rooms = . if inrange(year, 1985, 1993) & rooms == 99
replace rooms = . if year >= 1994 & inlist(rooms, 98, 99)

drop if year > DEATHYEAR
keep if CURRENT == 1
keep if !missing(AGEREP) & AGEREP >= 18
keep if !missing(IW) & IW > 0
keep if SEX == 2 & inlist(REL, 1, 2)
keep if !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0

egen byte multi_fu_hhyear_tag = tag(HHID year) if n_current_fids > 1
quietly count if multi_fu_hhyear_tag
local excl_multi_hh = r(N)
quietly count if n_current_fids > 1
local excl_multi_women = r(N)
drop multi_fu_hhyear_tag
keep if n_current_fids == 1

* Household outcomes enter once per household-year. Reference women receive
* deterministic priority in the remaining within-family duplicates. Reporter
* selection does not condition on outcome or education availability.
gen byte household_priority = REL != 1
sort HHID year household_priority ID
by HHID year: gen int household_women_before_dedup = _N
by HHID year: keep if _n == 1
isid HHID year

gen double rooms_alignment_gap_years = year - rooms_source_year
quietly count if !missing(rooms) & !inlist(rooms_alignment_gap_years, 1, 2)
local excl_bad_room_gap = r(N)
replace rooms = . if !inlist(rooms_alignment_gap_years, 1, 2)
drop rooms_alignment_gap_years
* Pre-fit code census before the original room/education missingness filter.
local prefit_candidate_rows = _N
quietly count if room_code_0 == 1
local prefit_room0 = r(N)
quietly count if room_code_9 == 1
local prefit_room9 = r(N)
quietly count if room_code_98 == 1
local prefit_room98 = r(N)
quietly count if room_code_99 == 1
local prefit_room99 = r(N)
quietly summarize rooms
local prefit_min_rooms = r(min)

drop if missing(rooms) | missing(EDUYEAR)
local estimation_input_rows = _N
quietly count if room_code_0 == 1
local input_room0 = r(N)
quietly count if room_code_9 == 1
local input_room9 = r(N)
quietly count if room_code_98 == 1
local input_room98 = r(N)
quietly count if room_code_99 == 1
local input_room99 = r(N)
quietly summarize rooms
local input_min_rooms = r(min)

bysort ID: egen double year_entry = min(year)
drop if !missing(first_child_year) & first_child_year < year_entry

gen byte untimed_known_parent = ///
    missing(first_child_year) & !missing(max_children_reported) & ///
    max_children_reported > 0
gen byte unknown_child_hist = ///
    missing(first_child_year) & missing(max_children_reported)
egen byte excl_known_tag = tag(ID) if untimed_known_parent
quietly count if excl_known_tag
local excl_known_ids = r(N)
egen byte excl_unknown_tag = tag(ID) if unknown_child_hist
quietly count if excl_unknown_tag
local excl_unknown_ids = r(N)
drop excl_known_tag excl_unknown_tag
drop if untimed_known_parent | unknown_child_hist
gen byte never_treated = missing(first_child_year) & max_children_reported == 0
assert never_treated == 1 if missing(first_child_year)
gen double K = year - first_child_year

cap drop L*event F*event
forvalues k = 0/10 {
    gen byte L`k'event = K == `k'
}
gen byte L11event = K > 10 & !missing(K)
gen byte F1event = K == -1
gen byte F3event = K == -3
gen byte F4event = K == -4
gen byte F5event = K == -5
gen byte F6event = K == -6
gen byte F7event = K <= -7 & !missing(K)

quietly count
local input_obs = r(N)
egen byte input_id_tag = tag(ID)
quietly count if input_id_tag
local input_ids = r(N)
drop input_id_tag
quietly count if never_treated
local never_obs = r(N)
egen byte treated_id_tag = tag(ID) if !never_treated
quietly count if treated_id_tag
local treated_ids = r(N)
drop treated_id_tag
quietly summarize household_women_before_dedup
local max_women_before_dedup = r(max)

di as text "Running weighted Sun-Abraham event study: obs=`input_obs', IDs=`input_ids', treated IDs=`treated_ids'"
eventstudyinteract rooms L*event F*event [pw=IW], ///
    vce(cluster ID) absorb(ID year) cohort(first_child_year) ///
    control_cohort(never_treated) covariates(i.AGEREP i.EDUYEAR)
quietly count if e(sample) & room_code_0 == 1
local es_room0 = r(N)
quietly count if e(sample) & room_code_9 == 1
local es_room9 = r(N)
quietly count if e(sample) & room_code_98 == 1
local es_room98 = r(N)
quietly count if e(sample) & room_code_99 == 1
local es_room99 = r(N)
quietly summarize rooms if e(sample)
local es_min_rooms = r(min)

quietly count if e(sample)
local sample_obs = r(N)
egen byte sample_id_tag = tag(ID) if e(sample)
quietly count if sample_id_tag
local sample_ids = r(N)
drop sample_id_tag
egen byte sample_hh_tag = tag(HHID year) if e(sample)
quietly count if sample_hh_tag
local sample_household_years = r(N)
drop sample_hh_tag
quietly count if e(sample) & never_treated
local sample_never_obs = r(N)
egen byte sample_never_id_tag = tag(ID) if e(sample) & never_treated
quietly count if sample_never_id_tag
local sample_never_ids = r(N)
drop sample_never_id_tag
quietly summarize rooms [aw=IW] if e(sample) & K == -2
local pre_event_mean = r(mean)

matrix b = e(b_iw)
matrix V = e(V_iw)
matrix b_interact = e(b_interact)
matrix V_interact = e(V_interact)
matrix b_full = e(b)
matrix V_full = e(V)

* e(b_interact) is cohort x event time. e(V_interact) is only its marginal
* variance array: the ado applies diagonal(e(V)) before reshaping it. The full
* interaction covariance is the leading q x q block of e(V), ordered event-time
* outermost and cohort innermost by the ado's bcohort_rel_varlist construction.
local cohortnames : rownames b_interact
local eventnames : colnames b_interact
local rawnames : colnames b_full
local ncohort = rowsof(b_interact)
local nevent = colsof(b_interact)
local n_interactions = `ncohort' * `nevent'
file open bfh using "`outdir'/cohort_interaction_coefficients.csv", write replace
file write bfh "matrix_index,cohort_index,cohort,event_index,event,stata_interaction_name,estimate,marginal_variance" _n
forvalues i = 1/`n_interactions' {
    local event_index = ceil(`i' / `ncohort')
    local cohort_index = mod(`i' - 1, `ncohort') + 1
    local cohort : word `cohort_index' of `cohortnames'
    local event : word `event_index' of `eventnames'
    local rawname : word `i' of `rawnames'
    local beta = b_interact[`cohort_index', `event_index']
    local beta_full = b_full[1, `i']
    local vari = V_interact[`cohort_index', `event_index']
    local vari_full = V_full[`i', `i']
    assert abs(`beta' - `beta_full') <= 1e-14
    assert abs(`vari' - `vari_full') <= 1e-14
    local beta_export : display %24.17g `beta'
    local vari_export : display %24.17g `vari'
    file write bfh `"`i',`cohort_index',`cohort',`event_index',`event',`rawname',`beta_export',`vari_export'"' _n
}
file close bfh
file open vfh using "`outdir'/cohort_interaction_covariance.csv", write replace
file write vfh "row_index,column_index,covariance" _n
forvalues i = 1/`n_interactions' {
    forvalues j = 1/`n_interactions' {
        local cov = V_full[`i', `j']
        local cov_export : display %24.17g `cov'
        file write vfh `"`i',`j',`cov_export'"' _n
    }
}
file close vfh

* Also export the aggregation used by the frozen target, with its full covariance.
local coefnames : colnames b
local ncoefs = colsof(b)
file open iwfh using "`outdir'/iw_coefficients.csv", write replace
file write iwfh "matrix_index,event_coefficient,estimate,variance,standard_error" _n
forvalues i = 1/`ncoefs' {
    local coef : word `i' of `coefnames'
    local beta = b[1, `i']
    local vari = V[`i', `i']
    local serr = sqrt(`vari')
    local beta_export : display %24.17g `beta'
    local vari_export : display %24.17g `vari'
    local serr_export : display %24.17g `serr'
    file write iwfh `"`i',`coef',`beta_export',`vari_export',`serr_export'"' _n
}
file close iwfh
file open iwvfh using "`outdir'/iw_covariance.csv", write replace
file write iwvfh "row_index,row_name,column_index,column_name,covariance" _n
forvalues i = 1/`ncoefs' {
    local rowname : word `i' of `coefnames'
    forvalues j = 1/`ncoefs' {
        local colname : word `j' of `coefnames'
        local cov = V[`i', `j']
        local cov_export : display %24.17g `cov'
        file write iwvfh `"`i',`rowname',`j',`colname',`cov_export'"' _n
    }
}
file close iwvfh
estimates save "`outdir'/eventstudyinteract_replay.ster", replace

gen byte __replay_esample = e(sample)
preserve
    keep if __replay_esample
    gen str24 cohort_group = cond(never_treated, "never_treated", string(first_child_year, "%9.0g"))
    collapse (count) estimation_observations=ID (sum) longitudinal_weight=IW, by(cohort_group K)
    rename K event_time
    sort cohort_group event_time
    export delimited using "`outdir'/estimation_cohort_event_support.csv", replace
restore
drop __replay_esample
matrix variance = vecdiag(V)
matrix combined = b \ variance
matrix rownames combined = b variance
matrix estimates = combined'

* The calibration moment is the four-year post-birth contrast from the closest
* pre-birth observation (k=-1) to k=+3. This removes the already-realized
* anticipatory housing adjustment between k=-2 and k=-1. Its uncertainty uses
* the full Sun--Abraham covariance matrix, not the sum of marginal variances.
local index_p3 = colnumb(b, "L3event")
local index_m1 = colnumb(b, "F1event")
assert `index_p3' > 0 & `index_m1' > 0
local component_p3 = b[1, `index_p3']
local component_m1 = b[1, `index_m1']
local covariance_p3_m1 = V[`index_p3', `index_m1']
local contrast_variance = V[`index_p3', `index_p3'] + ///
    V[`index_m1', `index_m1'] - 2 * `covariance_p3_m1'
assert `contrast_variance' > 0
local target_estimate = `component_p3' - `component_m1'
local target_se = sqrt(`contrast_variance')
* Nonfatal reproduction gates. Numeric tolerances match the original log's
* nine-decimal result display; sample counts must match exactly.
local gate_numeric_tolerance = 1e-9
local gate_target_pass = abs(`target_estimate' - .7202462623815278) <= `gate_numeric_tolerance'
local gate_se_pass = abs(`target_se' - .0852600513385958) <= `gate_numeric_tolerance'
local gate_obs_pass = `sample_obs' == 49457
local gate_ids_pass = `sample_ids' == 4112

preserve
    clear
    svmat2 estimates, names(col) rnames(coefficient_name)
    gen double se = sqrt(variance)
    drop variance
    replace b = . if b == 0 & se == 0
    replace se = . if missing(b) & se == 0
    drop if coefficient_name == "_cons"
    gen str12 event_string = subinstr(coefficient_name, "event", "", .)
    replace event_string = subinstr(event_string, "L", "", .)
    replace event_string = subinstr(event_string, "F", "-", .)
    replace event_string = subinstr(event_string, "o.", "", .)
    destring event_string, gen(relative_time)
    drop event_string
    assert !missing(relative_time)
    set obs `=_N+1'
    replace coefficient_name = "F2event" if missing(coefficient_name)
    replace relative_time = -2 if missing(relative_time)
    replace b = 0 if relative_time == -2
    replace se = 0 if relative_time == -2
    gen double ci_lo = b - 1.96 * se
    gen double ci_hi = b + 1.96 * se
    format b se ci_lo ci_hi %24.17g
    sort relative_time
    order coefficient_name relative_time b se ci_lo ci_hi
    export delimited using "`outdir'/event_study_estimates.csv", replace
    save "`outdir'/event_study_estimates.dta", replace

restore

timer off 1
quietly timer list 1
local runtime_seconds = r(t1)

clear
set obs 1
gen str50 moment = "housing_increment_0to1"
gen double estimate = `target_estimate'
gen double standard_error = `target_se'
gen int contrast_start_time = -1
gen int contrast_end_time = 3
gen int regression_omitted_time = -2
gen double component_l3 = `component_p3'
gen double component_f1 = `component_m1'
gen double covariance_l3_f1 = `covariance_p3_m1'
gen long input_observations = `input_obs'
gen long input_individuals = `input_ids'
gen long treated_individuals = `treated_ids'
gen long never_treated_observations = `never_obs'
gen long estimation_observations = `sample_obs'
gen long estimation_individuals = `sample_ids'
gen long estimation_household_years = `sample_household_years'
gen long est_never_obs = `sample_never_obs'
gen long est_confirmed_never_ids = `sample_never_ids'
gen long excl_multi_fu_hhyears = `excl_multi_hh'
gen long excl_multi_fu_women = `excl_multi_women'
gen long excl_untimed_parent_ids = `excl_known_ids'
gen long excl_unknown_history_ids = `excl_unknown_ids'
gen long excl_bad_room_gap_rows = `excl_bad_room_gap'
gen double pre_event_mean_rooms = `pre_event_mean'
gen int max_women_before_dedup = `max_women_before_dedup'
gen double runtime_seconds = `runtime_seconds'
gen str24 estimator = "eventstudyinteract"
gen str120 sample = "current women age 18+, ref/spouse, positive IW, single-FID dwelling, one woman per HH-year"
gen str100 fixed_effects = "person and survey-year fixed effects; age and education controls"
gen str40 clustering = "individual ID"
gen str40 weighting = "PSID longitudinal pweight IW"
gen str100 control_group = "confirmed zero-child women from full relationship history"
gen str80 rooms_alignment = "ACTUALROOMS_ shifted forward one observed interview within individual"
gen str80 source_file = "PSID/PSIDSHELF_MOBILITY.dta"
gen str80 fertility_timing = "first biological child across RELCHI1-20 TYPE/BYEAR records"
gen byte single_fu_only = 1
gen str32 status = "corrected_primary_target"
export delimited using "`outdir'/target_receipt.csv", replace
save "`outdir'/target_receipt.dta", replace
preserve
    clear
    set obs 1
    gen double target_expected = .7202462623815278
    gen double target_observed = `target_estimate'
    gen double target_abs_tolerance = `gate_numeric_tolerance'
    gen byte target_pass = `gate_target_pass'
    gen double se_expected = .0852600513385958
    gen double se_observed = `target_se'
    gen double se_abs_tolerance = `gate_numeric_tolerance'
    gen byte se_pass = `gate_se_pass'
    gen long observations_expected = 49457
    gen long observations_observed = `sample_obs'
    gen byte observations_pass = `gate_obs_pass'
    gen long people_expected = 4112
    gen long people_observed = `sample_ids'
    gen byte people_pass = `gate_ids_pass'
    export delimited using "`outdir'/reproduction_gates.csv", replace
restore
preserve
    clear
    set obs 3
    gen str40 stage = ""
    replace stage = "before_room_edu_drop" in 1
    replace stage = "post_room_edu_pre_cohort_filters" in 2
    replace stage = "e_sample" in 3
    gen long observations = .
    replace observations = `prefit_candidate_rows' in 1
    replace observations = `estimation_input_rows' in 2
    replace observations = `sample_obs' in 3
    gen long room_code_0_count = .
    replace room_code_0_count = `prefit_room0' in 1
    replace room_code_0_count = `input_room0' in 2
    replace room_code_0_count = `es_room0' in 3
    gen long room_code_9_count = .
    replace room_code_9_count = `prefit_room9' in 1
    replace room_code_9_count = `input_room9' in 2
    replace room_code_9_count = `es_room9' in 3
    gen long room_code_98_count = .
    replace room_code_98_count = `prefit_room98' in 1
    replace room_code_98_count = `input_room98' in 2
    replace room_code_98_count = `es_room98' in 3
    gen long room_code_99_count = .
    replace room_code_99_count = `prefit_room99' in 1
    replace room_code_99_count = `input_room99' in 2
    replace room_code_99_count = `es_room99' in 3
    gen double minimum_nonmissing_rooms = .
    replace minimum_nonmissing_rooms = `prefit_min_rooms' in 1
    replace minimum_nonmissing_rooms = `input_min_rooms' in 2
    replace minimum_nonmissing_rooms = `es_min_rooms' in 3
    export delimited using "`outdir'/room_code_diagnostics.csv", replace
restore

di as result "CORRECTED_FIRST_BIRTH_ROOMS_TARGET estimate=" %12.9f `target_estimate' " se=" %12.9f `target_se'
di as result "Receipt: `outdir'/target_receipt.csv"
di as result "Full estimates: `outdir'/event_study_estimates.csv"
di as result "Runtime seconds: `runtime_seconds'"
log close _all
