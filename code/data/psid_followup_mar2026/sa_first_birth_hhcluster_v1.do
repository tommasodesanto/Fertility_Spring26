* First-birth A2h event study with standard errors clustered by household (Oct 19, 2026).
* Point estimates are those of the person-clustered runs; only vce changes.
*   full support: sa_rooms_first_birth_v2.do A2h logic (rooms) and the v3 arms
*     A2h_own, A2h_moved (A2h filter on each outcome's non-missing rows)
*   common support (rooms_cs): sa_first_birth_model_bridge.do, original clock
* Household cluster: PSID HHID is wave-specific, so a cluster is a (year, HHID)
* pair. Treated adults take the household of their earliest reference-window
* (-3/-2) interview at which they are reference person or spouse; never-parent
* controls take the household of their first interview in the estimation sample.
* Each person keeps one cluster for all years, so persons nest in clusters.
* Panel: sa_first_birth_hhcluster_prepare_v1.do plus first_birth_model_bridge.py
* (scratch, deleted after the run).
version 17.0
args panel outdir
capture mkdir "`outdir'"
tempname C R
postfile `C' str12 run str8 window double estimate se_household using "`outdir'/coef.dta", replace
postfile `R' str12 run double observations persons clusters treated_ids control_ids missing_hh_persons persons_sharing seconds using "`outdir'/receipt.dta", replace
foreach run in rooms_full own moved rooms_cs {
    timer clear 2
    timer on 2
    use "`panel'", clear
    local y = cond("`run'"=="rooms_cs" | "`run'"=="rooms_full", "rooms", "`run'")
    keep if current & !missing(`y', AGEREP, EDUYEAR) & AGEREP >= 18 & !missing(iw) & iw > 0
    gen double f_c_y = bio_first_year
    drop if missing(f_c_y) & !(relchirep_max == 0)
    drop if !missing(f_c_y) & f_c_y < year_entry_adult
    gen byte control = missing(f_c_y) & relchirep_max == 0
    gen byte hs_row = inlist(rel, 1, 2) & inrange(year - f_c_y, -3, -2)
    bysort ID: egen byte hs_base = max(hs_row)
    gen byte base_row = inrange(year - f_c_y, -3, -2)
    bysort ID: egen byte has_base = max(base_row)
    keep if control | (has_base & hs_base == 1)
    drop hs_base base_row has_base
    gen byte treated = !control & !missing(f_c_y)
    if "`run'" == "rooms_cs" {
        gen byte bad_f = treated & (missing(f_j) | f_j < 1 | f_j > fertility_last_index)
        gen byte miss_model = missing(model_boundary_censored, model_boundary_lifelong, model_age_censored)
        drop if bad_f | miss_model
        gen byte original_ref = treated & inrange(year - f_c_y, -3, -2)
        bysort ID: egen byte id_original_ref = max(original_ref)
        drop if treated & !id_original_ref
    }
    else {
        * v2 cohort support: every treated cohort observed in the reference window.
        gen byte in_ref = inrange(year - f_c_y, -3, -2)
        bysort f_c_y: egen long ref_rows = total(in_ref)
        drop if treated & ref_rows == 0
    }
    gen double K = year - f_c_y
    gen byte Wleft = K <= -8 & !missing(K)
    gen byte Wm2 = inrange(K, -7, -6)
    gen byte Wm1 = inrange(K, -5, -4)
    gen byte Wp1 = inrange(K, -1, 0)
    gen byte Wp2 = inrange(K, 1, 2)
    gen byte Wp3 = inrange(K, 3, 4)
    gen byte Wp4 = inrange(K, 5, 6)
    gen byte Wp5 = inrange(K, 7, 8)
    gen byte Wp6 = inrange(K, 9, 10)
    gen byte Wright = K >= 11 & !missing(K)
    local dummies Wleft Wm2 Wm1 Wp1 Wp2 Wp3 Wp4 Wp5 Wp6 Wright
    * Household cluster.
    gen double hhkey = year * 100000 + hhid if !missing(hhid) & hhid > 0
    gen double tref_row = year if treated & hs_row == 1 & !missing(hhkey)
    bysort ID: egen double tref = min(tref_row)
    gen double kref_row = hhkey if treated & year == tref
    bysort ID: egen double hh_t = min(kref_row)
    gen double tfirst_row = year if control & !missing(hhkey)
    bysort ID: egen double tfirst = min(tfirst_row)
    gen double kfirst_row = hhkey if control & year == tfirst
    bysort ID: egen double hh_c = min(kfirst_row)
    gen double hhc = cond(control, hh_c, hh_t)
    egen byte ptag = tag(ID)
    quietly count if ptag & missing(hhc)
    local miss = r(N)
    * Persons without a valid household id keep a cluster of their own.
    replace hhc = -ID if missing(hhc)
    bysort hhc: egen long npers = total(ptag)
    quietly count if ptag & npers > 1
    local sharing = r(N)
    eventstudyinteract `y' `dummies' [pw=iw], vce(cluster hhc) absorb(ID year) ///
        cohort(f_c_y) control_cohort(control) covariates(i.AGEREP i.EDUYEAR)
    matrix b = e(b_iw)
    matrix V = e(V_iw)
    local n = e(N)
    local nc = e(N_clust)
    local names : colnames b
    forvalues i = 1/`=colsof(b)' {
        local w : word `i' of `names'
        post `C' ("`run'") ("`w'") (b[1,`i']) (sqrt(V[`i',`i']))
    }
    egen byte ttag = tag(ID) if e(sample) & treated
    quietly count if ttag
    local tids = r(N)
    egen byte ctag = tag(ID) if e(sample) & control
    quietly count if ctag
    local cids = r(N)
    egen byte etag = tag(ID) if e(sample)
    quietly count if etag
    local pids = r(N)
    timer off 2
    quietly timer list 2
    post `R' ("`run'") (`n') (`pids') (`nc') (`tids') (`cids') (`miss') (`sharing') (r(t2))
    di "HHC_PASS `run' N=`n' clusters=`nc' Wp3=" b[1, colnumb("b","Wp3")]
}
postclose `C'
postclose `R'
use "`outdir'/coef.dta", clear
format estimate se_household %24.17g
export delimited using "`outdir'/coefficients_household.csv", replace
use "`outdir'/receipt.dta", clear
export delimited using "`outdir'/receipt.csv", replace
erase "`outdir'/coef.dta"
erase "`outdir'/receipt.dta"
di "HHCLUSTER_DONE"
