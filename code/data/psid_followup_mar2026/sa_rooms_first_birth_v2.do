* First-birth housing-space event study, version 2 (September 24, 2026).
* Sun-Abraham interacted-cohort event study in two-year interview windows.
* Arms:
*   H    household design, baseline window -3/-2 (headline)
*   Hb4, Hb6  household design with earlier baseline windows -5/-4 and -7/-6
*   A2   all current adults 18+ with every other household-design choice:
*        IW weights, biological first birth, confirmed-childless controls,
*        entry at first adult observation, codes cleaned (baseline -3/-2)
*   A2h  A2, treated adults who were reference person or spouse in the
*        baseline window (pre-treatment status); all controls retained
*   A2n  A2, treated adults who were NOT reference person or spouse in the
*        baseline window; all controls retained
*   H21  household design, baseline window -2/-1 (sensitivity)
*   A    all-adult design, baseline window -3/-2 (robustness)
* Reconciliation arms (one change from H or A): Hentry = August entry rule
* (birth after first observation as head/spouse); Hctrl = last first-birth
* cohort as control cohort with childless women retained as untreated;
* Hshift, Ashift = shelf rooms shifted one interview (codes cleaned) instead
* of the rebuilt official variables.
* Household design: one woman per single-family-unit household-year (current
* reference person or spouse), PSID longitudinal weight IW, first birth from
* the full biological-child history, confirmed-childless women as the control
* cohort, first birth no earlier than the woman's first adult observation.
* All-adult design: every current adult 18+, unweighted, last first-birth
* cohort as control cohort, first child record as in the original code.
* Rooms: official year-specific PSID variables merged to the interview year;
* codes 9 (through 1984), 99 (1985-1993), 98/99 (1994 on) set to missing.
version 17.0
args mode outdir
timer clear 1
timer on 1
if "`mode'" == "toy" {
    capture mkdir "`outdir'"
    clear
    set seed 24092026
    set obs 8000
    gen long ID = ceil(_n/40)
    bysort ID: gen int year = 1979+_n
    keep if year<=1997 | mod(year,2)==1
    gen double bio_first_year = 1985+3*mod(ID,7)
    replace bio_first_year = . if mod(ID,5)==0
    gen double relchi1_year = bio_first_year
    replace relchi1_year = bio_first_year-1 if mod(ID,13)==0
    gen double K = year-bio_first_year
    gen double rooms = 4+ID/100+0.2*(K>=0 & !missing(K))+rnormal()
    gen byte AGEREP = 25+mod(ID,10)
    gen byte EDUYEAR = 12+mod(ID,4)
    gen double relchirep_max = cond(missing(bio_first_year), cond(mod(ID,10)==0,0,1), 1)
    replace relchirep_max = . if mod(ID,17)==0 & missing(bio_first_year)
    gen byte woman = mod(ID,2)==0
    gen byte current = 1
    gen double rel = 1+mod(ID,3)
    gen double hhid = ceil(ID/2)
    gen double fid = hhid
    gen int n_current_fids = 1+(mod(hhid,9)==0)
    gen double iw = 1+mod(ID,4)
    replace iw = 0 if mod(ID,19)==0
    bysort ID: egen double year_entry_adult = min(year)
    drop K
    gen double rooms_shift = rooms
    tempfile toydata
    save `toydata'
    foreach arm in H H21 A Hentry Hctrl Hshift Ashift Hb4 Hb6 A2 A2h A2n {
        use `toydata', clear
        do "sa_rooms_first_birth_v2.do" `arm' "`outdir'/`arm'"
    }
    di "ROOMS_V2_TOY_PASS"
    exit
}
assert inlist("`mode'","H","H21","A","Hentry","Hctrl","Hshift","Ashift","Hb4","Hb6") | inlist("`mode'","A2","A2h","A2n")
local household = substr("`mode'",1,1) == "H"
local alladult2 = substr("`mode'",1,2) == "A2"
if inlist("`mode'","Hshift","Ashift") {
    replace rooms = rooms_shift
    replace rooms = . if year <= 1984 & rooms == 9
    replace rooms = . if inrange(year,1985,1993) & rooms == 99
    replace rooms = . if year >= 1994 & inlist(rooms,98,99)
}
capture mkdir "`outdir'"
local baseline_lo = cond("`mode'"=="H21", -2, cond("`mode'"=="Hb4", -5, cond("`mode'"=="Hb6", -7, -3)))
local baseline_hi = cond("`mode'"=="H21", -1, cond("`mode'"=="Hb4", -4, cond("`mode'"=="Hb6", -6, -2)))

keep if current & !missing(rooms, AGEREP, EDUYEAR) & AGEREP >= 18
if `alladult2' {
    * All current adults with the household design's other choices.
    keep if !missing(iw) & iw > 0
    gen double f_c_y = bio_first_year
    drop if missing(f_c_y) & !(relchirep_max == 0)
    drop if !missing(f_c_y) & f_c_y < year_entry_adult
    gen byte control = missing(f_c_y) & relchirep_max == 0
    local weightspec "[pw=iw]"
    if "`mode'" != "A2" {
        gen byte hs_row = inlist(rel,1,2) & inrange(year - f_c_y, `baseline_lo', `baseline_hi')
        bysort ID: egen byte hs_base = max(hs_row)
        gen byte base_row = inrange(year - f_c_y, `baseline_lo', `baseline_hi')
        bysort ID: egen byte has_base = max(base_row)
        keep if control | (has_base & hs_base == cond("`mode'" == "A2h", 1, 0))
        drop hs_row hs_base base_row has_base
    }
}
else if !`household' {
    * All adults; first child record; untreated = no recorded first child.
    gen double f_c_y = relchi1_year
    drop if !missing(f_c_y) & f_c_y < year_entry_adult
    quietly summarize f_c_y
    gen byte control = f_c_y == r(max)
    local weightspec ""
}
else {
    keep if woman & inlist(rel,1,2)
    keep if !missing(iw) & iw > 0
    keep if !missing(hhid) & hhid > 0 & !missing(fid) & fid > 0 & n_current_fids == 1
    gen byte priority = rel != 1
    sort hhid year priority ID
    by hhid year: keep if _n == 1
    drop priority
    isid hhid year
    gen double f_c_y = bio_first_year
    drop if missing(f_c_y) & !(relchirep_max == 0)
    if "`mode'" == "Hentry" {
        bysort ID: egen double year_entry_hh = min(year)
        drop if !missing(f_c_y) & f_c_y < year_entry_hh
        drop year_entry_hh
    }
    else {
        drop if !missing(f_c_y) & f_c_y < year_entry_adult
    }
    gen byte control = missing(f_c_y) & relchirep_max == 0
    if "`mode'" == "Hctrl" {
        quietly summarize f_c_y
        replace control = f_c_y == r(max)
    }
    local weightspec "[pw=iw]"
}
gen double K = year - f_c_y

* Cohort support: every treated cohort must be observed in the baseline window.
gen byte in_reference = inrange(K,`baseline_lo',`baseline_hi')
bysort f_c_y: egen long reference_rows = total(in_reference)
gen byte treated = !missing(f_c_y) & control==0
gen byte unsupported = treated & reference_rows==0
preserve
    gen long rows = 1
    collapse (sum) rows (max) unsupported, by(f_c_y K control)
    export delimited using "`outdir'/input_support.csv", replace
restore
quietly count if unsupported
local dropped_rows = r(N)
quietly levelsof f_c_y if unsupported, local(dropped_cohorts) separate(" ")
drop if unsupported
drop in_reference reference_rows unsupported treated

local b = `baseline_hi'
gen byte Wleft = K<=`b'-6 & !missing(K)
gen byte Wm2 = inrange(K,`b'-5,`b'-4)
gen byte Wm1 = inrange(K,`b'-3,`b'-2)
gen byte Wp1 = inrange(K,`b'+1,`b'+2)
gen byte Wp2 = inrange(K,`b'+3,`b'+4)
gen byte Wp3 = inrange(K,`b'+5,`b'+6)
gen byte Wp4 = inrange(K,`b'+7,`b'+8)
gen byte Wp5 = inrange(K,`b'+9,`b'+10)
gen byte Wp6 = inrange(K,`b'+11,`b'+12)
gen byte Wright = K>=`b'+13 & !missing(K)
local dummies Wleft Wm2 Wm1 Wp1 Wp2 Wp3 Wp4 Wp5 Wp6 Wright
gen byte checksum = Wleft+Wm2+Wm1+Wp1+Wp2+Wp3+Wp4+Wp5+Wp6+Wright
assert checksum+inrange(K,`baseline_lo',`baseline_hi')==1 if !missing(K)
assert checksum==0 if missing(K)
drop checksum
* Headline = the window containing +3/+4: Wp((4-b)/2) for even baselines; +2/+3 for the -2/-1 baseline.
local headline = cond("`mode'"=="H21", "Wp2", "Wp" + string((4-`baseline_hi')/2))

quietly count
local input_obs = r(N)
eventstudyinteract rooms `dummies' `weightspec', ///
    vce(cluster ID) absorb(ID year) cohort(f_c_y) control_cohort(control) ///
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
if "`weightspec'" == "" {
    quietly summarize rooms if e(sample) & inrange(K,`baseline_lo',`baseline_hi')
}
else {
    quietly summarize rooms [aw=iw] if e(sample) & inrange(K,`baseline_lo',`baseline_hi')
}
local premean = r(mean)
quietly count if e(sample) & control
local control_rows = r(N)
egen byte ctag = tag(ID) if e(sample) & control
quietly count if ctag
local control_ids = r(N)
egen byte ttag = tag(ID) if e(sample) & !control & !missing(f_c_y)
quietly count if ttag
local treated_ids = r(N)
drop ctag ttag
preserve
    keep if e(sample)
    sort ID year
    export delimited ID year using "`outdir'/private_sample_keys.csv", replace
    gen byte in_reference = inrange(K,`baseline_lo',`baseline_hi')
    bysort f_c_y: egen long reference_rows = total(in_reference)
    quietly count if !missing(f_c_y) & control==0 & reference_rows==0
    local fitted_unsupported = r(N)
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
timer off 1
quietly timer list 1
local seconds = r(t1)
clear
set obs 1
gen str8 arm = "`mode'"
gen int baseline_lo = `baseline_lo'
gen int baseline_hi = `baseline_hi'
gen double input_observations = `input_obs'
gen double observations = `n'
gen double clusters = `clusters'
gen double treated_individuals = `treated_ids'
gen double control_individuals = `control_ids'
gen double control_rows = `control_rows'
gen str12 headline = "`headline'"
gen double headline_effect = `effect'
gen double headline_se = `effect_se'
gen double reference_mean_rooms = `premean'
gen long dropped_unsupported_rows = `dropped_rows'
gen str120 dropped_cohorts = "`dropped_cohorts'"
gen long fitted_unsupported_rows = `fitted_unsupported'
gen str12 weights = cond("`weightspec'"=="","none","pw=IW")
gen double runtime_seconds = `seconds'
format headline_effect headline_se reference_mean_rooms %24.17g
export delimited using "`outdir'/fit_receipt.csv", replace
di "ROOMS_V2_ARM_PASS `mode'"
