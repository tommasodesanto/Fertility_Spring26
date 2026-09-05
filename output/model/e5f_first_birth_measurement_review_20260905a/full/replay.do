clear all
set more off
set processors 8
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
local outroot "`project'/code/data/psid_followup_mar2026/output"
local outdir  "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_first_birth_measurement_review_20260905a/full"

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
drop if missing(rooms) | missing(EDUYEAR)

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


* Supplemental input-support receipts; no estimation changes.
assert `input_obs' == 49872
assert `input_ids' == 4527
assert `treated_ids' == 2486
preserve
    gen byte post1997 = year >= 1997
    gen long audit_one = 1
    collapse (sum) observations=audit_one weight_sum=IW, ///
        by(first_child_year K post1997 never_treated)
    export delimited using "`outdir'/input_cohort_event_support.csv", replace
restore

di as text "Running weighted Sun-Abraham event study: obs=`input_obs', IDs=`input_ids', treated IDs=`treated_ids'"
eventstudyinteract rooms L*event F*event [pw=IW], ///
    vce(cluster ID) absorb(ID year) cohort(first_child_year) ///
    control_cohort(never_treated) covariates(i.AGEREP i.EDUYEAR)

* Capture original estimation objects before any preserve/clear/export.
gen byte audit_estimation_sample = e(sample)
bysort ID: egen byte audit_observed_kminus2 = ///
    max(audit_estimation_sample & K == -2)
foreach matrix_name in b_iw V_iw ff_w b_interact V_interact {
    matrix audit_`matrix_name' = e(`matrix_name')
}

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


* Export captured objects after all primary estimator calculations.
foreach matrix_name in b_iw V_iw ff_w b_interact V_interact {
    matrix audit_matrix = audit_`matrix_name'
    preserve
        clear
        svmat2 audit_matrix, names(col) rnames(row_name)
        export delimited using "`outdir'/`matrix_name'.csv", replace
    restore
}
preserve
    keep if audit_estimation_sample
    gen byte post1997 = year >= 1997
    gen long audit_one = 1
    gen double audit_weighted_rooms = IW * rooms
    gen byte audit_nonpositive_rooms = rooms <= 0
    collapse (sum) observations=audit_one weight_sum=IW ///
        weighted_rooms_sum=audit_weighted_rooms ///
        nonpositive_room_rows=audit_nonpositive_rooms, ///
        by(first_child_year K post1997 never_treated audit_observed_kminus2)
    export delimited using "`outdir'/estimation_cohort_event_support.csv", replace
restore
di as result "FIRST_BIRTH_SUPPORT_OBJECTS_SAVED"

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

di as result "CORRECTED_FIRST_BIRTH_ROOMS_TARGET estimate=" %12.9f `target_estimate' " se=" %12.9f `target_se'
di as result "Receipt: `outdir'/target_receipt.csv"
di as result "Full estimates: `outdir'/event_study_estimates.csv"
di as result "Runtime seconds: `runtime_seconds'"
log close _all
