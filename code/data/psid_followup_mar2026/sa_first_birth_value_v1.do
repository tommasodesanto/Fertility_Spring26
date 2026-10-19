* First-birth house value, mortgage and interstate-move event study, v1 (Oct 18, 2026).
* Question: do families get space at first birth by moving to cheaper places,
* so that house value and mortgage barely rise?
* Design copied from the A2h arm of sa_rooms_first_birth_v2.do (rooms target):
*   all current adults 18+, PSID longitudinal weight IW, first biological birth,
*   confirmed-childless never-parents as the control cohort, first birth no earlier
*   than first adult observation, treated adults reference person or spouse in the
*   -3/-2 baseline window; Sun-Abraham IW estimator (eventstudyinteract), person
*   and year fixed effects, age and education indicators, two-year windows,
*   reference window -3/-2, SE clustered by person.
* Sample preparation copied from prepare_first_birth_outcomes_v3.do (official
* year-specific rooms and tenure items merged to the interview year).
* The A2h filter is evaluated on the rooms sample, as in the rooms target; each
* outcome then keeps its own non-missing rows. Arm own_rep re-runs the v3 A2h_own
* arm with the original logic (filter on the outcome sample) as a replication check.
* Shelf HOMEVALUER and HOMEMORTOTR are real 2022 USD at the interview year
* (HOMEVALUER verified equal to official A20 times a year constant, 1999 and 2001).
* HOMEMORTOTR exists only in 1984, 1989, 1994 and 1999-2019.
clear all
set more off
set processors 4
version 17.0
local out "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_first_birth_value_v1"
local cms "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta"
local shelf "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local arms `"`1'"'
if `"`arms'"' == "" local arms "own_rep lnval lnrooms_own lnvpr lnval_preown lnrooms_preown lnvpr_preown anymort lnmort ltv interstate"
timer clear 1
timer on 1
* 1. Official year-specific family items, merged to the interview year (copied from v3).
use ID V102 V592 V1263 V1966 V2565 V3107 V3521 V3937 V4448 V5362 V5862 V6477 V7080 V7671 V8360 V8969 V10432 V11614 V13019 V14122 V15138 V16639 V18070 V19370 V20670 V22425 ER2029 ER5028 ER7028 ER10032 ER13037 ER17040 ER21039 ER25027 ER36027 ER42028 ER47328 ER53028 ER60029 ER66029 ER72029 using "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta", clear
isid ID
rename V102 rooms1968
rename V592 rooms1969
rename V1263 rooms1970
rename V1966 rooms1971
rename V2565 rooms1972
rename V3107 rooms1973
rename V3521 rooms1974
rename V3937 rooms1975
rename V4448 rooms1976
rename V5362 rooms1977
rename V5862 rooms1978
rename V6477 rooms1979
rename V7080 rooms1980
rename V7671 rooms1981
rename V8360 rooms1982
rename V8969 rooms1983
rename V10432 rooms1984
rename V11614 rooms1985
rename V13019 rooms1986
rename V14122 rooms1987
rename V15138 rooms1988
rename V16639 rooms1989
rename V18070 rooms1990
rename V19370 rooms1991
rename V20670 rooms1992
rename V22425 rooms1993
rename ER2029 rooms1994
rename ER5028 rooms1995
rename ER7028 rooms1996
rename ER10032 rooms1997
rename ER13037 rooms1999
rename ER17040 rooms2001
rename ER21039 rooms2003
rename ER25027 rooms2005
rename ER36027 rooms2007
rename ER42028 rooms2009
rename ER47328 rooms2011
rename ER53028 rooms2013
rename ER60029 rooms2015
rename ER66029 rooms2017
rename ER72029 rooms2019
reshape long rooms, i(ID) j(year)
drop if missing(rooms)
tempfile src_rooms
save `src_rooms'
use ID V103 V593 V1264 V1967 V2566 V3108 V3522 V3939 V4450 V5364 V5864 V6479 V7084 V7675 V8364 V8974 V10437 V11618 V13023 V14126 V15140 V16641 V18072 V19372 V20672 V22427 ER2032 ER5031 ER7031 ER10035 ER13040 ER17043 ER21042 ER25028 ER36028 ER42029 ER47329 ER53029 ER60030 ER66030 ER72030 using "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta", clear
isid ID
rename V103 own_raw1968
rename V593 own_raw1969
rename V1264 own_raw1970
rename V1967 own_raw1971
rename V2566 own_raw1972
rename V3108 own_raw1973
rename V3522 own_raw1974
rename V3939 own_raw1975
rename V4450 own_raw1976
rename V5364 own_raw1977
rename V5864 own_raw1978
rename V6479 own_raw1979
rename V7084 own_raw1980
rename V7675 own_raw1981
rename V8364 own_raw1982
rename V8974 own_raw1983
rename V10437 own_raw1984
rename V11618 own_raw1985
rename V13023 own_raw1986
rename V14126 own_raw1987
rename V15140 own_raw1988
rename V16641 own_raw1989
rename V18072 own_raw1990
rename V19372 own_raw1991
rename V20672 own_raw1992
rename V22427 own_raw1993
rename ER2032 own_raw1994
rename ER5031 own_raw1995
rename ER7031 own_raw1996
rename ER10035 own_raw1997
rename ER13040 own_raw1999
rename ER17043 own_raw2001
rename ER21042 own_raw2003
rename ER25028 own_raw2005
rename ER36028 own_raw2007
rename ER42029 own_raw2009
rename ER47329 own_raw2011
rename ER53029 own_raw2013
rename ER60030 own_raw2015
rename ER66030 own_raw2017
rename ER72030 own_raw2019
reshape long own_raw, i(ID) j(year)
drop if missing(own_raw)
tempfile src_own_raw
save `src_own_raw'

* 2. Shelf columns.
use ID year AGEREP EDUYEAR SEX REL CURRENT HHID FID IW DEATHYEAR RELCHIREP RELCHI1BYEAR ACTUALROOMS_ HOMEOWN HOMEVALUER HOMEMORTOTR GEOSTATE RELCHI1TYPE RELCHI2TYPE RELCHI2BYEAR RELCHI3TYPE RELCHI3BYEAR RELCHI4TYPE RELCHI4BYEAR RELCHI5TYPE RELCHI5BYEAR RELCHI6TYPE RELCHI6BYEAR RELCHI7TYPE RELCHI7BYEAR RELCHI8TYPE RELCHI8BYEAR RELCHI9TYPE RELCHI9BYEAR RELCHI10TYPE RELCHI10BYEAR RELCHI11TYPE RELCHI11BYEAR RELCHI12TYPE RELCHI12BYEAR RELCHI13TYPE RELCHI13BYEAR RELCHI14TYPE RELCHI14BYEAR RELCHI15TYPE RELCHI15BYEAR RELCHI16TYPE RELCHI16BYEAR RELCHI17TYPE RELCHI17BYEAR RELCHI18TYPE RELCHI18BYEAR RELCHI19TYPE RELCHI19BYEAR RELCHI20TYPE RELCHI20BYEAR using "`shelf'", clear
sort ID year
by ID: gen double rooms_shift = ACTUALROOMS_[_n-1] if _n > 1 & inlist(year-year[_n-1],1,2)
foreach s in rooms own_raw {
    merge 1:1 ID year using `src_`s'', keep(master match) nogenerate
}
* Verification as in v3: rebuilt rooms vs shifted shelf column; own vs shelf HOMEOWN.
quietly count if CURRENT == 1 & !missing(rooms, rooms_shift) & rooms != rooms_shift
assert r(N) == 0
gen byte own = cond(own_raw == 1, 1, cond(own_raw == 5, 0, .)) if !missing(own_raw)
gen byte own_shelf = cond(HOMEOWN == 1, 1, cond(HOMEOWN == 2, 0, .)) if !missing(HOMEOWN)
quietly count if CURRENT == 1 & !missing(own, own_shelf) & own != own_shelf
assert r(N) == 0
replace rooms = . if year <= 1984 & rooms == 9
replace rooms = . if inrange(year,1985,1993) & rooms == 99
replace rooms = . if year >= 1994 & inlist(rooms,98,99)

* 3. Code and timing checks, 2001 owners: value and mortgage DK/NA codes; state timing.
preserve
    keep if year == 2001 & CURRENT == 1
    keep ID own HOMEVALUER HOMEMORTOTR GEOSTATE
    tempfile chk2001
    save `chk2001'
    use ID ER17004 ER17044 ER17052 using "`cms'", clear
    merge 1:1 ID using `chk2001', keep(using match) nogenerate
    di "CHECK owners 2001"
    count if own == 1
    count if own == 1 & ER17044 >= 9999997 & !missing(ER17044)
    count if own == 1 & ER17044 >= 9999997 & !missing(ER17044) & !missing(HOMEVALUER)
    count if own == 1 & missing(HOMEVALUER)
    count if own == 1 & HOMEVALUER == 0
    count if own == 1 & HOMEMORTOTR == 0
    count if own == 1 & missing(HOMEMORTOTR)
    count if own == 1 & ER17052 >= 9999997 & !missing(ER17052)
    count if own == 1 & ER17052 >= 9999997 & !missing(ER17052) & !missing(HOMEMORTOTR)
    summarize HOMEVALUER HOMEMORTOTR if own == 1, detail
    di "CHECK state timing 2001: GEOSTATE vs official ER17004 (mismatch with the modal map)"
    egen double modal = mode(GEOSTATE), by(ER17004)
    count if !missing(GEOSTATE, ER17004)
    count if !missing(GEOSTATE, ER17004) & GEOSTATE != modal
restore

* 4. Person-level constructions before any restriction (copied from v3).
gen double bio_candidate = .
forvalues c = 1/20 {
    replace bio_candidate = RELCHI`c'BYEAR if RELCHI`c'TYPE == 1 & !missing(RELCHI`c'BYEAR) & (missing(bio_candidate) | RELCHI`c'BYEAR < bio_candidate)
}
bysort ID: egen double bio_first_year = min(bio_candidate)
drop bio_candidate
bysort ID: egen double relchirep_max = max(RELCHIREP)
gen double relchi1_year = RELCHI1BYEAR
egen byte fid_tag = tag(HHID year FID) if CURRENT == 1 & !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0
bysort HHID year: egen int n_current_fids = total(fid_tag)
drop fid_tag
gen byte woman = SEX == 2
gen byte current = CURRENT == 1
gen double rel = REL
gen double hhid = HHID
gen double fid = FID
gen double iw = IW
gen double adult_year = year if CURRENT == 1 & AGEREP >= 18 & !missing(AGEREP) & year <= DEATHYEAR
bysort ID: egen double year_entry_adult = min(adult_year)
drop adult_year
keep if CURRENT == 1
drop if year > DEATHYEAR
* Interstate move: state differs from the previous observed interview, one or two years earlier.
sort ID year
by ID: gen double prev_state = GEOSTATE[_n-1] if _n > 1 & inlist(year - year[_n-1], 1, 2)
gen byte interstate = GEOSTATE != prev_state if !missing(GEOSTATE, prev_state)
gen double value = HOMEVALUER
gen double mort = HOMEMORTOTR
keep ID year AGEREP EDUYEAR rooms own value mort interstate bio_first_year relchirep_max current rel iw year_entry_adult
isid ID year
compress
tempfile base
save `base'

* 5. Event studies, A2h design.
tempname P S
tempfile coeffile sumfile
postfile `P' str16 arm str8 window int k_lo int k_hi double estimate double standard_error using `coeffile', replace
postfile `S' str16 arm str70 definition int status double observations clusters treated_individuals control_individuals ref_mean_outcome ref_mean_level control_mean_outcome headline_effect headline_se dropped_unsupported_rows topcoded_rows runtime_seconds using `sumfile', replace
local ownsamp "own == 1 & value > 0 & !missing(value) & rooms > 0 & !missing(rooms)"
foreach arm of local arms {
    timer clear 2
    timer on 2
    use `base', clear
    local basis = cond("`arm'" == "own_rep", "own", "rooms")
    keep if current & !missing(`basis', AGEREP, EDUYEAR) & AGEREP >= 18
    keep if !missing(iw) & iw > 0
    gen double f_c_y = bio_first_year
    drop if missing(f_c_y) & !(relchirep_max == 0)
    drop if !missing(f_c_y) & f_c_y < year_entry_adult
    gen byte control = missing(f_c_y) & relchirep_max == 0
    gen byte hs_row = inlist(rel,1,2) & inrange(year - f_c_y, -3, -2)
    bysort ID: egen byte hs_base = max(hs_row)
    gen byte base_row = inrange(year - f_c_y, -3, -2)
    bysort ID: egen byte has_base = max(base_row)
    keep if control | (has_base & hs_base == 1)
    drop hs_row hs_base base_row has_base
    gen double K = year - f_c_y
    if strpos("`arm'", "_preown") {
        * Treated must own at every observed baseline (-3/-2) interview.
        gen byte b_own = inrange(K,-3,-2) & own == 1
        gen byte b_not = inrange(K,-3,-2) & own != 1
        bysort ID: egen byte any_b_own = max(b_own)
        bysort ID: egen byte any_b_not = max(b_not)
        keep if control | (any_b_own & !any_b_not)
        drop b_own b_not any_b_own any_b_not
    }
    gen double y = .
    gen double lvl = .
    local tc = 0
    if "`arm'" == "own_rep" {
        replace y = own
        replace lvl = own
        local def "owner; replication of v3 A2h_own"
    }
    if inlist("`arm'", "lnval", "lnval_preown") {
        replace y = ln(value) if `ownsamp'
        replace lvl = value if !missing(y)
        local def "log real house value (2022 USD); owners"
    }
    if inlist("`arm'", "lnrooms_own", "lnrooms_preown") {
        replace y = ln(rooms) if `ownsamp'
        replace lvl = rooms if !missing(y)
        local def "log rooms; owners (same rows as lnval)"
    }
    if inlist("`arm'", "lnvpr", "lnvpr_preown") {
        replace y = ln(value / rooms) if `ownsamp'
        replace lvl = value / rooms if !missing(y)
        local def "log real house value per room; owners (same rows as lnval)"
    }
    if "`arm'" == "anymort" {
        replace y = mort > 0 if `ownsamp' & !missing(mort)
        replace lvl = y
        local def "any mortgage; owners; mortgage years only"
    }
    if "`arm'" == "lnmort" {
        replace y = ln(mort) if `ownsamp' & mort > 0 & !missing(mort)
        replace lvl = mort if !missing(y)
        local def "log real mortgage principal; owners with a mortgage"
    }
    if "`arm'" == "ltv" {
        quietly count if `ownsamp' & !missing(mort) & mort / value > 2
        local tc = r(N)
        replace y = min(mort / value, 2) if `ownsamp' & !missing(mort)
        replace lvl = mort / value if !missing(y)
        local def "mortgage / value, top-coded at 2; owners incl. zero mortgage"
    }
    if "`arm'" == "interstate" {
        replace y = interstate
        replace lvl = interstate
        local def "state differs from previous interview; all tenures"
    }
    if strpos("`arm'", "_preown") local def "`def'; treated owned at -3/-2"
    keep if !missing(y)
    * Cohort support: every treated cohort must be observed in the baseline window.
    gen byte in_reference = inrange(K,-3,-2)
    bysort f_c_y: egen long reference_rows = total(in_reference)
    gen byte unsupported = !missing(f_c_y) & control == 0 & reference_rows == 0
    quietly count if unsupported
    local dropped = r(N)
    drop if unsupported
    drop in_reference reference_rows unsupported
    gen byte Wleft = K <= -8 & !missing(K)
    gen byte Wm2 = inrange(K,-7,-6)
    gen byte Wm1 = inrange(K,-5,-4)
    gen byte Wp1 = inrange(K,-1,0)
    gen byte Wp2 = inrange(K,1,2)
    gen byte Wp3 = inrange(K,3,4)
    gen byte Wp4 = inrange(K,5,6)
    gen byte Wp5 = inrange(K,7,8)
    gen byte Wp6 = inrange(K,9,10)
    gen byte Wright = K >= 11 & !missing(K)
    local dummies Wleft Wm2 Wm1 Wp1 Wp2 Wp3 Wp4 Wp5 Wp6 Wright
    di "ARM `arm': `def'"
    capture noisily eventstudyinteract y `dummies' [pw=iw], ///
        vce(cluster ID) absorb(ID year) cohort(f_c_y) control_cohort(control) ///
        covariates(i.AGEREP i.EDUYEAR)
    local rc = _rc
    timer off 2
    quietly timer list 2
    local secs = r(t2)
    if `rc' {
        di "ARM_FAIL `arm' rc=`rc'"
        post `S' ("`arm'") ("`def'") (`rc') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (`dropped') (`tc') (`secs')
        continue
    }
    matrix b = e(b_iw)
    matrix V = e(V_iw)
    local n = e(N)
    local clusters = e(N_clust)
    local names : colnames b
    forvalues i = 1/`=colsof(b)' {
        local w : word `i' of `names'
        local lo = cond("`w'"=="Wleft", ., cond("`w'"=="Wm2", -7, cond("`w'"=="Wm1", -5, cond("`w'"=="Wp1", -1, cond("`w'"=="Wp2", 1, cond("`w'"=="Wp3", 3, cond("`w'"=="Wp4", 5, cond("`w'"=="Wp5", 7, cond("`w'"=="Wp6", 9, 11)))))))))
        local hi = cond("`w'"=="Wleft", -8, cond("`w'"=="Wright", ., `lo' + 1))
        post `P' ("`arm'") ("`w'") (`lo') (`hi') (b[1,`i']) (sqrt(V[`i',`i']))
    }
    local h = colnumb("b", "Wp3")
    quietly summarize y [aw=iw] if e(sample) & inrange(K,-3,-2)
    local rm = r(mean)
    quietly summarize lvl [aw=iw] if e(sample) & inrange(K,-3,-2)
    local rl = r(mean)
    quietly summarize y [aw=iw] if e(sample) & control
    local cm = r(mean)
    egen byte ctag = tag(ID) if e(sample) & control
    quietly count if ctag
    local cids = r(N)
    egen byte ttag = tag(ID) if e(sample) & !control & !missing(f_c_y)
    quietly count if ttag
    local tids = r(N)
    post `S' ("`arm'") ("`def'") (0) (`n') (`clusters') (`tids') (`cids') (`rm') (`rl') (`cm') (b[1,`h']) (sqrt(V[`h',`h'])) (`dropped') (`tc') (`secs')
    di "ARM_PASS `arm' N=`n' Wp3=" b[1,`h'] " se=" sqrt(V[`h',`h'])
}
postclose `P'
postclose `S'
use `coeffile', clear
format estimate standard_error %12.6g
export delimited using "`out'/coefficients.csv", replace
use `sumfile', clear
format ref_mean_outcome ref_mean_level control_mean_outcome headline_effect headline_se %12.6g
export delimited using "`out'/summary.csv", replace
timer off 1
quietly timer list 1
di "VALUE_V1_DONE seconds=" r(t1)
