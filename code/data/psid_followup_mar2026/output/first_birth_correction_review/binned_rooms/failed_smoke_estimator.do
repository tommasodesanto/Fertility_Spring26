* Diagnostic only: common observations, original specification, two-year bins.
* Baseline -3/-2; supported cohorts; last cohort used only before treatment.
version 17.0
args mode outdir
timer clear 1
timer on 1
if "`mode'" == "toy" {
    clear
    set seed 12092026
    set obs 4800
    gen long ID = ceil(_n/40)
    bysort ID: gen int year = 1979+_n
    keep if year<=1997 | mod(year,2)==1
    gen int f_c_y = 1989+3*mod(ID,5)
    gen double K = year-f_c_y
    gen byte lastcohort = f_c_y==2001
    gen byte AGEREP = 25+mod(ID,10)
    gen byte EDUYEAR = 12+mod(ID,4)
    gen double rooms = 4+ID/100+0.2*(K>=0)+rnormal()
    sort ID year
    by ID: gen double rooms_aligned = rooms[_n-1]
    gen byte common_rooms = !missing(rooms,rooms_aligned)
    tempfile toydata
    save `toydata'
    foreach arm in original_binned aligned_binned {
        use `toydata', clear
        do "audit_binned_rooms.do" `arm' "`outdir'/`arm'"
    }
    di "BINNED_ROOMS_TOY_PASS"
    exit
}
assert inlist("`mode'","original_binned","aligned_binned")
capture mkdir "`outdir'"
quietly summarize f_c_y if lastcohort==1, meanonly
local control_year = r(min)
assert r(min)==r(max) & !missing(r(min))
keep if year < `control_year'
keep if common_rooms & !missing(rooms,rooms_aligned,AGEREP,EDUYEAR,f_c_y)
* Each displayed two-year bin, including the baseline, must have support.
gen byte bin = floor((K+3)/2)
gen byte supported = 1
forvalues j = -2/3 {
    gen byte here = bin==`j'
    bysort f_c_y: egen long nbin = total(here)
    replace supported = 0 if nbin==0
    drop here nbin
}
preserve
    gen long rows = 1
    collapse (sum) rows (min) supported, by(f_c_y bin lastcohort)
    export delimited using "`outdir'/input_support.csv", replace
restore
keep if supported | lastcohort
if "`mode'" == "aligned_binned" replace rooms = rooms_aligned
gen byte Dm7 = inrange(K,-7,-6)
gen byte Dm5 = inrange(K,-5,-4)
gen byte Dm1 = inrange(K,-1,0)
gen byte Dp1 = inrange(K,1,2)
gen byte Dp3 = inrange(K,3,4)
gen byte Dleft = K<=-8
gen byte Dright = K>=5
assert Dm7+Dm5+Dm1+Dp1+Dp3+Dleft+Dright+inrange(K,-3,-2)==1
eventstudyinteract rooms Dm7 Dm5 Dm1 Dp1 Dp3 Dleft Dright, ///
    vce(cluster ID) absorb(ID year) cohort(f_c_y) control_cohort(lastcohort) ///
    covariates(i.AGEREP i.EDUYEAR)
matrix b = e(b_iw)
matrix V = e(V_iw)
local n = e(N)
local clusters = e(N_clust)
local p3 = colnumb(b,"Dp3")
assert !missing(`p3') & `p3'>0
local effect = b[1,`p3']
local effect_se = sqrt(V[`p3',`p3'])
assert !missing(`effect',`effect_se') & `effect_se'>0
preserve
    keep if e(sample)
    sort ID year
    export delimited ID year using "`outdir'/private_sample_keys.csv", replace
    forvalues j = -2/3 {
        gen byte here = bin==`j'
        bysort f_c_y: egen long nbin = total(here)
        assert nbin>0 if !lastcohort
        drop here nbin
    }
    gen long rows = 1
    collapse (sum) rows, by(f_c_y bin lastcohort)
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
use "`outdir'/coefficients.dta", clear
format estimate standard_error %24.17g
export delimited using "`outdir'/coefficients.csv", replace
use "`outdir'/covariance.dta", clear
format covariance %24.17g
export delimited using "`outdir'/covariance.csv", replace
timer off 1
quietly timer list 1
local seconds = r(t1)
clear
set obs 1
gen str24 arm = "`mode'"
gen double observations = `n'
gen double clusters = `clusters'
gen double post_3_4_vs_pre_3_2 = `effect'
gen double standard_error = `effect_se'
gen double control_birth_year = `control_year'
gen double runtime_seconds = `seconds'
format post_3_4_vs_pre_3_2 standard_error %24.17g
export delimited using "`outdir'/fit_receipt.csv", replace
di "BINNED_ROOMS_ARM_PASS `mode'"
