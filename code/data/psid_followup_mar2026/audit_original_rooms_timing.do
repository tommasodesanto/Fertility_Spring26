* Timing-only comparisons against the author's recognized regression.
* These diagnostics intentionally retain original room codes and controls.
* They cannot promote a replacement empirical target.
version 17.0
args mode outdir
timer clear 1
timer on 1
if "`mode'" == "toy" {
    clear
    set seed 12092026
    set obs 2400
    gen long ID = ceil(_n/20)
    bysort ID: gen int year = 1989+_n
    gen int f_c_y = 1995 + 2*mod(ID,5)
    gen double K = year-f_c_y
    gen byte lastcohort = f_c_y==2003
    gen byte AGEREP = 25+mod(ID,10)
    gen byte EDUYEAR = 12+mod(ID,4)
    gen double rooms = 4+ID/100+0.2*(K>=0)+rnormal()
    sort ID year
    by ID: gen double rooms_aligned = rooms[_n-1]
    gen byte common_rooms = !missing(rooms,rooms_aligned)
    gen byte cohort_has_reference = 1
    gen byte cohort_has_m2 = 1
    tempfile toydata
    save `toydata'
    foreach arm in original_native original_common aligned_common {
        use `toydata', clear
        do "audit_original_rooms_timing.do" `arm' "`outdir'/`arm'"
    }
    di "ORIGINAL_TIMING_TOY_PASS"
    exit
}
capture mkdir "`outdir'"
* Rows unusable in either fit are excluded by the original marksample/markout.
* Remove them before allocating interactions, preserving the full cohort list.
quietly levelsof f_c_y if lastcohort==0, local(cohorts_before)
drop if missing(rooms,AGEREP,EDUYEAR) & missing(rooms_aligned,AGEREP,EDUYEAR)
quietly levelsof f_c_y if lastcohort==0, local(cohorts_after)
assert "`cohorts_before'" == "`cohorts_after'"
if inlist("`mode'","original_common","aligned_common") keep if common_rooms
if "`mode'" == "aligned_common" replace rooms = rooms_aligned
assert inlist("`mode'","original_native","original_common","aligned_common")

* Verbatim event construction and command from the recognized Stata file.
capture drop L*event F*event
forvalues l = 0/10 {
    gen L`l'event = K==`l'
}
gen L11event = (K>10 & K!=.)
forvalues l = 1/5 {
    gen F`l'event = K==-`l'
}
gen F7event = (K<-6 & K!=.)
drop F2event
eventstudyinteract rooms L*event F*event, ///
    vce(cluster ID) absorb(ID year) cohort(f_c_y) control_cohort(lastcohort) ///
    covariates(i.AGEREP i.EDUYEAR)

matrix b = e(b_iw)
matrix V = e(V_iw)
local n = e(N)
local clusters = e(N_clust)
local p3 = colnumb(b,"L3event")
local m1 = colnumb(b,"F1event")
assert !missing(`p3',`m1') & `p3'>0 & `m1'>0
local difference = b[1,`p3']-b[1,`m1']
local difference_se = sqrt(V[`p3',`p3']+V[`m1',`m1']-2*V[`p3',`m1'])
assert !missing(`difference',`difference_se') & `difference_se'>0
quietly summarize rooms if e(sample) & K==-2
local premean = r(mean)
preserve
    keep if e(sample)
    sort ID year
    export delimited ID year using "`outdir'/private_sample_keys.csv", replace
    gen long rows = 1
    collapse (sum) rows, by(f_c_y K lastcohort)
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
timer off 1
quietly timer list 1
local seconds = r(t1)
clear
set obs 1
gen str24 arm = "`mode'"
gen double observations = `n'
gen double clusters = `clusters'
gen double contrast_p3_m1 = `difference'
gen double contrast_se = `difference_se'
gen double premean = `premean'
gen double runtime_seconds = `seconds'
format contrast_p3_m1 contrast_se premean %24.17g
export delimited using "`outdir'/fit_receipt.csv", replace
di "ORIGINAL_TIMING_ARM_PASS `mode'"
