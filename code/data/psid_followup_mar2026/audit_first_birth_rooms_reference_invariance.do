* Diagnostic only: same interacted least-squares regression, two normalizations.
* Called by audit_first_birth_rooms_measurement.py after the primary sample
* builder, or with the small synthetic sample used to exercise this exact loop.
args mode outdir
if "`mode'" == "toy" {
    clear
    set seed 51986
    set obs 2520
    gen long ID = ceil(_n/21)
    bysort ID: gen int year = 1978 + _n - 1
    gen int first_child_year = 1982 + 4*mod(ID,3) if ID <= 90
    gen byte never_treated = ID > 90
    gen int K = year-first_child_year
    drop if first_child_year == 1986 & K == -2
    gen double IW = 1+mod(ID,4)+.01*(year-1978)*mod(ID,3)
    gen int AGEREP = 20+year-1978+mod(ID,5)
    gen int EDUYEAR = 12+mod(ID,4)
    gen double rooms = 4+.01*ID+.03*(year-1978)+.6*(K>=0 & K<.)+rnormal(0,.1)
    forvalues k=0/10 {
        gen byte L`k'event = K==`k'
    }
    gen byte L11event = K>10 & !missing(K)
    foreach k in 1 3 4 5 6 {
        gen byte F`k'event = K==-`k'
    }
    gen byte F7event = K<=-7 & !missing(K)
}
assert !missing(ID, year, AGEREP, EDUYEAR, IW, rooms)
bysort ID: assert first_child_year==first_child_year[1]
egen byte indicator_total = rowtotal(L*event F*event)
assert indicator_total==(K!=-2) if !never_treated
assert indicator_total==0 if never_treated
drop indicator_total
quietly count if first_child_year==1986 & K==-2
assert r(N)==0
quietly count if first_child_year==1986 & K==-1
assert r(N)>0
quietly count if first_child_year==1986 & K==3
assert r(N)>0
preserve
    gen long one=1
    collapse (sum) observations=one weight_sum=IW, by(first_child_year K never_treated)
    export delimited using "`outdir'/input_support.csv", replace
restore

quietly levelsof first_child_year if !never_treated, local(cohorts)
local interactions ""
foreach l of varlist L*event F*event {
    foreach g of local cohorts {
        gen byte sa_`l'_`g' = (first_child_year==`g')*`l'
        local interactions "`interactions' sa_`l'_`g'"
    }
}
local alternative ""
foreach v of local interactions {
    if "`v'"!="sa_F1event_1986" local alternative "`alternative' `v'"
}
tempname summaries
foreach specification in original reference {
    local regressors "`interactions'"
    if "`specification'"=="reference" local regressors "`alternative'"
    di as result "REFERENCE_CHECK_START `specification'"
    quietly reghdfe rooms `regressors' i.AGEREP i.EDUYEAR [pw=IW], ///
        absorb(ID year) vce(cluster ID) residuals(residual_`specification')
    local obs = e(N)
    local clusters = e(N_clust)
    local rss = e(rss)
    gen byte sample_`specification'=e(sample)
    if "`specification'"=="original" & "`mode'"=="full" {
        assert `obs'==49457
        assert `clusters'==4112
    }
    local maxdiff=0
    if "`specification'"=="reference" {
        assert sample_original==sample_reference
        gen double fitted_difference=abs(residual_reference-residual_original)
        quietly summarize fitted_difference if sample_original, meanonly
        local maxdiff=r(max)
        assert `maxdiff'<2e-6
    }
    postfile `summaries' str12 specification double observations clusters rss ///
        max_fitted_difference using "`outdir'/fit_receipt_`specification'.dta", replace
    post `summaries' ("`specification'") (`obs') (`clusters') (`rss') (`maxdiff')
    postclose `summaries'
    tempname coefficients
    postfile `coefficients' int cohort str12 event double coefficient variance ///
        byte explicit_reference using "`outdir'/coefficients_`specification'.dta", replace
    foreach l of varlist L*event F*event {
        foreach g of local cohorts {
            if "`specification'"=="reference" & "`l'"=="F1event" & `g'==1986 {
                post `coefficients' (`g') ("`l'") (0) (0) (1)
            }
            else {
                post `coefficients' (`g') ("`l'") (_b[sa_`l'_`g']) ///
                    (_se[sa_`l'_`g']^2) (0)
            }
        }
    }
    postclose `coefficients'
    preserve
        use "`outdir'/fit_receipt_`specification'.dta", clear
        export delimited using "`outdir'/fit_receipt_`specification'.csv", replace
    restore
    preserve
        use "`outdir'/coefficients_`specification'.dta", clear
        export delimited using "`outdir'/coefficients_`specification'.csv", replace
    restore
    preserve
        keep if sample_`specification'
        gen long one=1
        collapse (sum) observations=one weight_sum=IW, by(first_child_year K never_treated)
        export delimited using "`outdir'/estimation_support_`specification'.csv", replace
    restore
    di as result "REFERENCE_CHECK_DONE `specification' obs=`obs' max_fitted_difference=`maxdiff'"
}
di as result "FIRST_BIRTH_REFERENCE_INVARIANCE_CHECK_COMPLETED"
log close _all
