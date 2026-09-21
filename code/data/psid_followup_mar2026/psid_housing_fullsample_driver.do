clear all
set more off
version 17.0

* Custom full-sample extensions. Exact author arms are staged and executed by
* launch_psid_housing_fullsample_torch.sh in separate Stata processes because
* the author files intentionally begin with clear all.
args source outroot arm variant
if "`arm'" == "" local arm "aligned_first_ownership"
if "`variant'" == "" local variant "all"
capture confirm file "`source'"
if _rc {
    di as error "Missing PSID shelf: `source'"
    exit 601
}
cap mkdir "`outroot'"

cap which eventstudyinteract
if _rc {
    di as error "eventstudyinteract is not installed"
    exit 198
}
cap which svmat2
if _rc {
    di as error "svmat2 is not installed"
    exit 198
}

if "`arm'" == "aligned_first_ownership" {
    local outdir "`outroot'/first_birth_aligned_ownership"
    cap mkdir "`outdir'"
    log using "`outdir'/aligned_first_ownership.log", replace text
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
        HOMEOWN IW using "`source'", clear
    sort ID year
    gen double birth_candidate = .
    forvalues child = 1/20 {
        replace birth_candidate = RELCHI`child'BYEAR if ///
            RELCHI`child'TYPE == 1 & !missing(RELCHI`child'BYEAR) & ///
            (missing(birth_candidate) | RELCHI`child'BYEAR < birth_candidate)
    }
    bysort ID: egen double first_child_year = min(birth_candidate)
    bysort ID: egen double max_children_reported = max(RELCHIREP)
    drop birth_candidate
    egen byte current_fid_tag = tag(HHID year FID) if CURRENT == 1 & ///
        !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0
    bysort HHID year: egen int n_current_fids = total(current_fid_tag)
    drop current_fid_tag
    drop if year > DEATHYEAR
    keep if CURRENT == 1 & !missing(AGEREP) & AGEREP >= 18
    keep if !missing(IW) & IW > 0 & SEX == 2 & inlist(REL,1,2)
    keep if !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0
    keep if n_current_fids == 1
    gen byte household_priority = REL != 1
    sort HHID year household_priority ID
    by HHID year: keep if _n == 1
    isid HHID year
    * Outcome-specific missingness: this arm does not require rooms observed.
    replace HOMEOWN = 0 if HOMEOWN == 2
    replace HOMEOWN = . if HOMEOWN == 3
    gen double own = HOMEOWN
    drop if missing(own)
    bysort ID: egen double year_entry = min(year)
    drop if !missing(first_child_year) & first_child_year < year_entry
    drop if missing(first_child_year) & (max_children_reported > 0 | missing(max_children_reported))
    gen byte never_treated = missing(first_child_year) & max_children_reported == 0
    gen double K = year - first_child_year
    forvalues k = 0/10 {
        gen byte L`k'event = K == `k'
    }
    gen byte L11event = K > 10 & !missing(K)
    forvalues k = 1/6 {
        gen byte F`k'event = K == -`k'
    }
    gen byte F7event = K <= -7 & !missing(K)
    drop F2event
    eventstudyinteract own L*event F*event [pw=IW], ///
        vce(cluster ID) absorb(ID year) cohort(first_child_year) ///
        control_cohort(never_treated) covariates(i.AGEREP i.EDUYEAR)
    matrix b = e(b_iw)
    matrix V = e(V_iw)
    local nobs = e(N)
    local ip = colnumb(b,"L3event")
    local im = colnumb(b,"F1event")
    assert !missing(`ip') & !missing(`im') & `ip' > 0 & `im' > 0
    local vv = V[`ip',`ip'] + V[`im',`im'] - 2*V[`ip',`im']
    assert `vv' >= 0
    local est = b[1,`ip'] - b[1,`im']
    preserve
        clear
        svmat2 b, names(col) rnames(coefficient_name)
        export delimited using "`outdir'/event_study_estimates.csv", replace
    restore
    preserve
        clear
        svmat2 V, names(col) rnames(row_name)
        export delimited using "`outdir'/event_study_covariance.csv", replace
    restore
    clear
    set obs 1
    gen str40 arm = "first_birth_aligned_ownership"
    gen double contrast_l3_minus_f1 = `est'
    gen double contrast_se = sqrt(`vv')
    gen double contrast_ci_lo = contrast_l3_minus_f1 - 1.96*contrast_se
    gen double contrast_ci_hi = contrast_l3_minus_f1 + 1.96*contrast_se
    gen long estimation_observations = `nobs'
    gen str16 weighting = "IW pweight"
    gen str80 sample_note = "corrected HH-year selection; ownership missingness only"
    export delimited using "`outdir'/contrast.csv", replace
    save "`outdir'/contrast.dta", replace
    log close _all
    exit
}

if "`arm'" == "second_ownership" {
    local outdir "`outroot'/second_birth_ownership/`variant'"
    cap mkdir "`outroot'/second_birth_ownership"
    cap mkdir "`outdir'"
    log using "`outdir'/second_ownership.log", replace text
    use ID year AGEREP EDUYEAR SEX DEATHYEAR RELCHI1BYEAR RELCHI2BYEAR ///
        RELCHI3BYEAR HOMEOWN IW using "`source'", clear
    drop if year > DEATHYEAR | AGEREP < 18 | missing(AGEREP)
    bysort ID: egen double first_child_year = min(RELCHI1BYEAR)
    bysort ID: egen double second_child_year = min(RELCHI2BYEAR)
    bysort ID: egen double third_child_year = min(RELCHI3BYEAR)
    drop if missing(first_child_year)
    bysort ID: egen double year_entry = min(year)
    drop if first_child_year < year_entry
    drop if year < first_child_year
    drop if !missing(second_child_year) & second_child_year == first_child_year
    drop if !missing(third_child_year) & third_child_year == second_child_year
    if "`variant'" == "no_third_by3" {
        drop if !missing(second_child_year) & !missing(third_child_year) & third_child_year <= second_child_year + 3
    }
    if "`variant'" == "no_third_by3_gap5" {
        drop if !missing(second_child_year) & second_child_year - first_child_year < 5
        drop if !missing(second_child_year) & !missing(third_child_year) & third_child_year <= second_child_year + 3
    }
    replace HOMEOWN = 0 if HOMEOWN == 2
    replace HOMEOWN = . if HOMEOWN == 3
    gen double own = HOMEOWN
    drop if missing(own)
    gen byte stay_one_control = missing(second_child_year)
    gen double cohort = second_child_year
    gen double K = year - second_child_year
    replace K = . if stay_one_control
    forvalues k = 0/10 {
        gen byte L`k'event = cond(stay_one_control,0,K==`k')
    }
    gen byte L11event = cond(stay_one_control,0,K>10 & !missing(K))
    forvalues k = 1/6 {
        gen byte F`k'event = cond(stay_one_control,0,K==-`k')
    }
    gen byte F7event = cond(stay_one_control,0,K < -6 & !missing(K))
    drop F2event
    * Exact legacy second-birth scripts are unweighted; this extension matches that command.
    eventstudyinteract own L*event F*event, vce(cluster ID) absorb(year) ///
        cohort(cohort) control_cohort(stay_one_control) covariates(i.AGEREP i.EDUYEAR)
    matrix b = e(b_iw)
    matrix V = e(V_iw)
    local nobs = e(N)
    local ip = colnumb(b,"L3event")
    local im = colnumb(b,"F1event")
    assert !missing(`ip') & !missing(`im') & `ip' > 0 & `im' > 0
    local vv = V[`ip',`ip'] + V[`im',`im'] - 2*V[`ip',`im']
    assert `vv' >= 0
    preserve
        clear
        svmat2 b, names(col) rnames(coefficient_name)
        export delimited using "`outdir'/event_study_estimates.csv", replace
    restore
    preserve
        clear
        svmat2 V, names(col) rnames(row_name)
        export delimited using "`outdir'/event_study_covariance.csv", replace
    restore
    clear
    set obs 1
    gen str40 arm = "second_birth_ownership"
    gen str24 variant = "`variant'"
    gen double contrast_l3_minus_f1 = b[1,`ip'] - b[1,`im']
    gen double contrast_se = sqrt(`vv')
    gen double contrast_ci_lo = contrast_l3_minus_f1 - 1.96*contrast_se
    gen double contrast_ci_hi = contrast_l3_minus_f1 + 1.96*contrast_se
    gen long estimation_observations = `nobs'
    gen str32 weighting = "unweighted author-style extension"
    export delimited using "`outdir'/contrast.csv", replace
    save "`outdir'/contrast.dta", replace
    log close _all
    exit
}

di as error "Unknown custom arm: `arm'"
exit 198
