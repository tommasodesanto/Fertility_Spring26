* Data-only designs on the author's original Sun-Abraham specification.
* Two designs, each with the original and the verified (own-interview) room
* dates: (i) two-year interview windows with a -3/-2 baseline window, and
* (ii) the original annual event time restricted to birth cohorts observed at
* event year -2, read at the even horizons 0, +2, +4, +6, +8.
* Retained from the original specification: all adults, no weights, ID/year
* fixed effects, age and education covariates, ID clustering, last first-birth
* cohort as control_cohort, people without a recorded first birth retained,
* original room codes. Only the room date assignment differs within a design.
* The -2-only design adds the K=-6 indicator so that -2 is the sole omitted
* event year; the window design has no omitted year other than the baseline.
version 17.0
args mode outdir
timer clear 1
timer on 1
if "`mode'" == "toy" {
    capture mkdir "`outdir'"
    clear
    set seed 23092026
    set obs 6000
    gen long ID = ceil(_n/40)
    bysort ID: gen int year = 1979+_n
    keep if year<=1997 | mod(year,2)==1
    gen int f_c_y = 1985+3*mod(ID,7)
    replace f_c_y = . if mod(ID,11)==0
    gen double K = year-f_c_y
    quietly summarize f_c_y
    gen byte lastcohort = f_c_y==r(max)
    gen byte AGEREP = 25+mod(ID,10)
    gen byte EDUYEAR = 12+mod(ID,4)
    gen double rooms = 4+ID/100+0.2*(K>=0 & !missing(K))+rnormal()
    sort ID year
    by ID: gen double rooms_aligned = rooms[_n-1]
    gen byte common_rooms = !missing(rooms,rooms_aligned)
    tempfile toydata
    save `toydata'
    foreach arm in window_original window_aligned m2only_original m2only_aligned {
        use `toydata', clear
        do "audit_rooms_window_designs.do" `arm' "`outdir'/`arm'"
    }
    di "ROOMS_WINDOW_DESIGNS_TOY_PASS"
    exit
}
assert inlist("`mode'","window_original","window_aligned","m2only_original","m2only_aligned")
capture mkdir "`outdir'"
local design = cond(strpos("`mode'","window")>0,"window","m2only")
local assignment = cond(strpos("`mode'","aligned")>0,"aligned","original")

* Observations complete under both date assignments, as in the timing-only
* comparison. Rows without a recorded first birth and the last cohort stay.
keep if common_rooms & !missing(AGEREP,EDUYEAR)
if "`assignment'" == "aligned" replace rooms = rooms_aligned

* Cohort support: every treated cohort must be observed in its reference.
if "`design'" == "window" {
    gen byte in_reference = inrange(K,-3,-2)
}
else {
    gen byte in_reference = K==-2
}
bysort f_c_y: egen long reference_rows = total(in_reference)
gen byte treated = !missing(f_c_y) & lastcohort==0
gen byte unsupported = treated & reference_rows==0
preserve
    gen long rows = 1
    collapse (sum) rows (max) unsupported, by(f_c_y K lastcohort)
    export delimited using "`outdir'/input_support.csv", replace
restore
quietly count if unsupported
local dropped_rows = r(N)
quietly levelsof f_c_y if unsupported, local(dropped_cohorts) separate(" ")
drop if unsupported
drop in_reference reference_rows unsupported treated

if "`design'" == "window" {
    gen byte Dleft = K<=-8 & !missing(K)
    gen byte Dm7 = inrange(K,-7,-6)
    gen byte Dm5 = inrange(K,-5,-4)
    gen byte Dm1 = inrange(K,-1,0)
    gen byte Dp1 = inrange(K,1,2)
    gen byte Dp3 = inrange(K,3,4)
    gen byte Dp5 = inrange(K,5,6)
    gen byte Dp7 = inrange(K,7,8)
    gen byte Dp9 = inrange(K,9,10)
    gen byte Dright = K>=11 & !missing(K)
    local dummies Dleft Dm7 Dm5 Dm1 Dp1 Dp3 Dp5 Dp7 Dp9 Dright
    gen byte checksum = Dleft+Dm7+Dm5+Dm1+Dp1+Dp3+Dp5+Dp7+Dp9+Dright
    assert checksum+inrange(K,-3,-2)==1 if !missing(K)
    assert checksum==0 if missing(K)
    drop checksum
    local headline Dp3
    local reference_condition inrange(K,-3,-2)
}
else {
    forvalues l = 0/10 {
        gen byte L`l'event = K==`l'
    }
    gen byte L11event = (K>10 & K!=.)
    foreach l in 1 3 4 5 6 {
        gen byte F`l'event = K==-`l'
    }
    gen byte F7event = (K<-6 & K!=.)
    local dummies L0event L1event L2event L3event L4event L5event L6event ///
        L7event L8event L9event L10event L11event F1event F3event F4event ///
        F5event F6event F7event
    gen byte checksum = 0
    foreach d of local dummies {
        replace checksum = checksum+`d'
    }
    assert checksum+(K==-2)==1 if !missing(K)
    assert checksum==0 if missing(K)
    drop checksum
    local headline L4event
    local reference_condition K==-2
}

eventstudyinteract rooms `dummies', ///
    vce(cluster ID) absorb(ID year) cohort(f_c_y) control_cohort(lastcohort) ///
    covariates(i.AGEREP i.EDUYEAR)

matrix b = e(b_iw)
matrix V = e(V_iw)
local n = e(N)
local clusters = e(N_clust)
local h = colnumb("b","`headline'")
assert !missing(`h') & `h'>0
local effect = b[1,`h']
local effect_se = sqrt(V[`h',`h'])
assert !missing(`effect',`effect_se') & `effect_se'>0
quietly summarize rooms if e(sample) & `reference_condition'
local premean = r(mean)
preserve
    keep if e(sample)
    sort ID year
    export delimited ID year using "`outdir'/private_sample_keys.csv", replace
    gen byte in_reference = `reference_condition'
    bysort f_c_y: egen long reference_rows = total(in_reference)
    quietly count if !missing(f_c_y) & lastcohort==0 & reference_rows==0
    local fitted_unsupported = r(N)
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
gen str8 design = "`design'"
gen str8 assignment = "`assignment'"
gen double observations = `n'
gen double clusters = `clusters'
gen str12 headline = "`headline'"
gen double headline_effect = `effect'
gen double headline_se = `effect_se'
gen double reference_mean_rooms = `premean'
gen long dropped_unsupported_rows = `dropped_rows'
gen str80 dropped_cohorts = "`dropped_cohorts'"
gen long fitted_unsupported_rows = `fitted_unsupported'
gen double runtime_seconds = `seconds'
format headline_effect headline_se reference_mean_rooms %24.17g
export delimited using "`outdir'/fit_receipt.csv", replace
di "ROOMS_WINDOW_DESIGNS_ARM_PASS `mode'"
