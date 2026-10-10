* Later-birth housing-space event study (October 10, 2026).
* The A2h first-birth design applied to the q-th biological birth, q = 2 or 3.
* Treated: current adults whose q-th biological birth is a separate event from
* the (q-1)-th (no twins across the two orders), observed as an adult, and who
* were reference person or spouse in the -3/-2 baseline window.
* Controls: confirmed exactly q-1 children (q-1 biological records, reported
* number of children q-1 throughout).
* Rows: person-years at or after the (q-1)-th birth year, so the comparison is
* between households already at q-1 children.
* support = full   : the above only.
* support = model  : also requires a valid DOB proxy, observation age cell 0..16,
*                    (q-1)-th birth cell >= 0, q-th birth cell <= 6 and strictly
*                    after the (q-1)-th birth cell (the baseline period is then a
*                    model period with q-1 children); controls need a fertility
*                    report censor age and (q-1)-th birth cell 0..6.
* extra = prevtime : sensitivity adding categorical years since the (q-1)-th
*                    birth (0..14, 15+) as a covariate.
* Estimator, weights, fixed effects, covariates, clustering, windows, cohort
* support and headline are exactly those of sa_rooms_first_birth_v2.do (A2h).
version 17.0
args q support outdir extra
assert inlist(`q',2,3)
assert inlist("`support'","full","model")
capture mkdir "`outdir'"
timer clear 1
timer on 1
local p = `q' - 1
local evvar = cond(`q'==2, "bio_second_year", "bio_third_year")
local prevvar = cond(`q'==2, "bio_first_year", "bio_second_year")
local evj = cond(`q'==2, "b2_j", "b3_j")
local prevj = cond(`q'==2, "b1_j", "b2_j")

keep if current & !missing(rooms, AGEREP, EDUYEAR) & AGEREP >= 18
keep if !missing(iw) & iw > 0
* Rows in the (q-1)-child-or-more state.
keep if !missing(`prevvar') & year >= `prevvar'
gen double f_c_y = `evvar'
gen byte twin_cross = !missing(f_c_y) & f_c_y == `prevvar'
quietly count if twin_cross
local twin_rows = r(N)
drop if twin_cross
drop twin_cross
gen byte control = missing(f_c_y) & relchirep_max == `p' & nbio_max == `p'
drop if missing(f_c_y) & !control
drop if !missing(f_c_y) & f_c_y < year_entry_adult
if "`support'" == "model" {
    gen byte drop_model = missing(model_dob_proxy) | !inrange(model_age_index,0,16) | missing(`prevj') | `prevj' < 0
    replace drop_model = 1 if !control & (missing(`evj') | `evj' > 6 | `evj' <= `prevj')
    replace drop_model = 1 if control & (missing(censor_j) | `prevj' > 6)
    quietly count if drop_model
    local model_dropped_rows = r(N)
    drop if drop_model
    drop drop_model
}
else local model_dropped_rows = 0
gen byte hs_row = inlist(rel,1,2) & inrange(year - f_c_y, -3, -2)
bysort ID: egen byte hs_base = max(hs_row)
keep if control | hs_base == 1
drop hs_row hs_base
gen double K = year - f_c_y
local covs i.AGEREP i.EDUYEAR
if "`extra'" == "prevtime" {
    gen int since_prev = min(year - `prevvar', 15)
    local covs `covs' i.since_prev
}

* Cohort support: every treated cohort must be observed in the baseline window.
gen byte in_reference = inrange(K,-3,-2)
bysort f_c_y: egen long reference_rows = total(in_reference)
gen byte treated = !missing(f_c_y) & control==0
gen byte unsupported = treated & reference_rows==0
quietly count if unsupported
local dropped_rows = r(N)
quietly levelsof f_c_y if unsupported, local(dropped_cohorts) separate(" ")
drop if unsupported
drop in_reference reference_rows unsupported treated

gen byte Wleft = K<=-8 & !missing(K)
gen byte Wm2 = inrange(K,-7,-6)
gen byte Wm1 = inrange(K,-5,-4)
gen byte Wp1 = inrange(K,-1,0)
gen byte Wp2 = inrange(K,1,2)
gen byte Wp3 = inrange(K,3,4)
gen byte Wp4 = inrange(K,5,6)
gen byte Wp5 = inrange(K,7,8)
gen byte Wp6 = inrange(K,9,10)
gen byte Wright = K>=11 & !missing(K)
local dummies Wleft Wm2 Wm1 Wp1 Wp2 Wp3 Wp4 Wp5 Wp6 Wright
gen byte checksum = Wleft+Wm2+Wm1+Wp1+Wp2+Wp3+Wp4+Wp5+Wp6+Wright
assert checksum+inrange(K,-3,-2)==1 if !missing(K)
assert checksum==0 if missing(K)
drop checksum
local headline Wp3

* Cell counts by window for the treated, before estimation.
preserve
    gen long rows = 1
    gen str8 window = "base"
    foreach d of local dummies {
        replace window = "`d'" if `d'
    }
    replace window = "control" if control
    egen byte idtag = tag(ID window)
    collapse (sum) rows persons=idtag, by(window)
    export delimited using "`outdir'/window_cells.csv", replace
restore
quietly count
local input_obs = r(N)
eventstudyinteract rooms `dummies' [pw=iw], ///
    vce(cluster ID) absorb(ID year) cohort(f_c_y) control_cohort(control) ///
    covariates(`covs')

matrix b = e(b_iw)
matrix V = e(V_iw)
local n = e(N)
local clusters = e(N_clust)
local h = colnumb("b","`headline'")
local effect = b[1,`h']
local effect_se = sqrt(V[`h',`h'])
assert !missing(`effect',`effect_se') & `effect_se'>0
quietly summarize rooms [aw=iw] if e(sample) & inrange(K,-3,-2)
local premean = r(mean)
egen byte ctag = tag(ID) if e(sample) & control
quietly count if ctag
local control_ids = r(N)
egen byte ttag = tag(ID) if e(sample) & !control & !missing(f_c_y)
quietly count if ttag
local treated_ids = r(N)
quietly count if e(sample) & !control & inrange(K,3,4)
local headline_rows = r(N)
egen byte htag = tag(ID) if e(sample) & !control & inrange(K,3,4)
quietly count if htag
local headline_ids = r(N)
quietly summarize f_c_y if e(sample) & !control
local cmin = r(min)
local cmax = r(max)
quietly summarize `prevvar' if ttag
local gap_note = ""
gen double spacing = f_c_y - `prevvar' if ttag
quietly summarize spacing, detail
local sp_p50 = r(p50)
local sp_mean = r(mean)
drop ctag ttag htag spacing
preserve
    keep if e(sample)
    sort ID year
    export delimited ID year using "/tmp/later_birth_rooms_20261010/private_keys_q`q'_`support'`extra'.csv", replace
    gen long rows = 1
    collapse (sum) rows, by(f_c_y K control)
    export delimited using "`outdir'/fitted_support.csv", replace
restore
local names : colnames b
tempname points covariance
postfile `points' str20 coefficient double estimate standard_error using "`outdir'/coefficients.dta", replace
postfile `covariance' str20 coefficient_i str20 coefficient_j double covariance using "`outdir'/covariance.dta", replace
forvalues i = 1/`=colsof(b)' {
    local name : word `i' of `names'
    post `points' ("`name'") (b[1,`i']) (sqrt(V[`i',`i']))
    forvalues j = 1/`=colsof(b)' {
        local other : word `j' of `names'
        post `covariance' ("`name'") ("`other'") (V[`i',`j'])
    }
}
postclose `points'
postclose `covariance'
preserve
    use "`outdir'/coefficients.dta", clear
    format estimate standard_error %24.17g
    export delimited using "`outdir'/coefficients.csv", replace
    use "`outdir'/covariance.dta", clear
    format covariance %24.17g
    export delimited using "`outdir'/covariance.csv", replace
restore
erase "`outdir'/coefficients.dta"
erase "`outdir'/covariance.dta"
timer off 1
quietly timer list 1
local seconds = r(t1)
clear
set obs 1
gen int birth_order = `q'
gen str8 support = "`support'"
gen str12 extra = "`extra'"
gen double input_observations = `input_obs'
gen double observations = `n'
gen double clusters = `clusters'
gen double treated_individuals = `treated_ids'
gen double control_individuals = `control_ids'
gen double headline_window_rows = `headline_rows'
gen double headline_window_treated = `headline_ids'
gen str12 headline = "`headline'"
gen double headline_effect = `effect'
gen double headline_se = `effect_se'
gen double reference_mean_rooms = `premean'
gen double spacing_median = `sp_p50'
gen double spacing_mean = `sp_mean'
gen long twin_cross_rows_dropped = `twin_rows'
gen long model_support_rows_dropped = `model_dropped_rows'
gen long dropped_unsupported_rows = `dropped_rows'
gen str200 dropped_cohorts = "`dropped_cohorts'"
gen int cohort_min = `cmin'
gen int cohort_max = `cmax'
gen double runtime_seconds = `seconds'
format headline_effect headline_se reference_mean_rooms %24.17g
export delimited using "`outdir'/fit_receipt.csv", replace
di "LATER_BIRTH_PASS q=`q' support=`support' effect=`effect' se=`effect_se' N=`n' treated=`treated_ids' controls=`control_ids'"
