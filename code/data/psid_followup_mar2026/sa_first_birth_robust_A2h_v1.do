* Rooms robustness on the A2h sample (Oct 19, 2026): the Table 8 rows rerun on the
* sample of the main text. Sample definition held fixed (treated adults reference
* person or spouse at -3/-2, never-parent controls, IW weights, entry at first
* adult observation); each variant changes one thing, as in the household-sample
* rows of sa_rooms_first_birth_v2.do (H21, Hb4, Hb6, Hshift, Hctrl):
*   base   A2h as in v2 (check: 1.465 at +3/+4)
*   r21    reference window -2/-1, reported window +2/+3
*   b4     reference window -5/-4, reported window +3/+4
*   b6     reference window -7/-6, reported window +3/+4
*   shift  harmonized room measure (shelf rooms aligned to the interview)
*   ctrl   last first-birth cohort as control cohort, never-parents retained
* Preparation copied from prepare_first_birth_outcomes_v3.do.
clear all
set more off
set processors 4
version 17.0
local out "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_first_birth_robust_A2h_v1"
capture mkdir "`out'"
* 1. Official year-specific family items, merged to the interview year.
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
use ID V603 V1274 V1979 V2577 V3110 V3524 V3941 V4452 V5366 V5866 V6484 V7089 V7700 V8369 V8999 V10447 V11628 V13037 V14140 V15148 V16649 V18087 V19387 V20687 V22441 ER2062 ER5061 ER7155 ER10072 ER13077 ER17088 ER21117 ER25098 ER36103 ER42132 ER47440 ER53140 ER60155 ER66156 ER72156 using "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta", clear
isid ID
rename V603 moved_raw1969
rename V1274 moved_raw1970
rename V1979 moved_raw1971
rename V2577 moved_raw1972
rename V3110 moved_raw1973
rename V3524 moved_raw1974
rename V3941 moved_raw1975
rename V4452 moved_raw1976
rename V5366 moved_raw1977
rename V5866 moved_raw1978
rename V6484 moved_raw1979
rename V7089 moved_raw1980
rename V7700 moved_raw1981
rename V8369 moved_raw1982
rename V8999 moved_raw1983
rename V10447 moved_raw1984
rename V11628 moved_raw1985
rename V13037 moved_raw1986
rename V14140 moved_raw1987
rename V15148 moved_raw1988
rename V16649 moved_raw1989
rename V18087 moved_raw1990
rename V19387 moved_raw1991
rename V20687 moved_raw1992
rename V22441 moved_raw1993
rename ER2062 moved_raw1994
rename ER5061 moved_raw1995
rename ER7155 moved_raw1996
rename ER10072 moved_raw1997
rename ER13077 moved_raw1999
rename ER17088 moved_raw2001
rename ER21117 moved_raw2003
rename ER25098 moved_raw2005
rename ER36103 moved_raw2007
rename ER42132 moved_raw2009
rename ER47440 moved_raw2011
rename ER53140 moved_raw2013
rename ER60155 moved_raw2015
rename ER66156 moved_raw2017
rename ER72156 moved_raw2019
reshape long moved_raw, i(ID) j(year)
drop if missing(moved_raw)
tempfile src_moved_raw
save `src_moved_raw'
use ID V604 V1275 V1980 V2578 V3111 V3525 V3943 V4454 V5368 V5868 V6486 V7091 V7702 V8370 V9001 V10449 V11630 V13039 V14142 V15150 V16651 V18089 V19389 V20689 V22444 ER2065 ER5064 ER7158 ER10075 ER13080 ER17091 ER21120 ER25101 ER36106 ER42135 ER47443 ER53143 ER60158 ER66159 ER72159 using "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta", clear
isid ID
rename V604 why_raw1969
rename V1275 why_raw1970
rename V1980 why_raw1971
rename V2578 why_raw1972
rename V3111 why_raw1973
rename V3525 why_raw1974
rename V3943 why_raw1975
rename V4454 why_raw1976
rename V5368 why_raw1977
rename V5868 why_raw1978
rename V6486 why_raw1979
rename V7091 why_raw1980
rename V7702 why_raw1981
rename V8370 why_raw1982
rename V9001 why_raw1983
rename V10449 why_raw1984
rename V11630 why_raw1985
rename V13039 why_raw1986
rename V14142 why_raw1987
rename V15150 why_raw1988
rename V16651 why_raw1989
rename V18089 why_raw1990
rename V19389 why_raw1991
rename V20689 why_raw1992
rename V22444 why_raw1993
rename ER2065 why_raw1994
rename ER5064 why_raw1995
rename ER7158 why_raw1996
rename ER10075 why_raw1997
rename ER13080 why_raw1999
rename ER17091 why_raw2001
rename ER21120 why_raw2003
rename ER25101 why_raw2005
rename ER36106 why_raw2007
rename ER42135 why_raw2009
rename ER47443 why_raw2011
rename ER53143 why_raw2013
rename ER60158 why_raw2015
rename ER66159 why_raw2017
rename ER72159 why_raw2019
reshape long why_raw, i(ID) j(year)
drop if missing(why_raw)
tempfile src_why_raw
save `src_why_raw'

* 2. Shelf columns.
use ID year AGEREP EDUYEAR SEX REL CURRENT HHID FID IW DEATHYEAR RELCHIREP RELCHI1BYEAR ACTUALROOMS_ HOMEOWN RELCHI1TYPE RELCHI1BYEAR RELCHI2TYPE RELCHI2BYEAR RELCHI3TYPE RELCHI3BYEAR RELCHI4TYPE RELCHI4BYEAR RELCHI5TYPE RELCHI5BYEAR RELCHI6TYPE RELCHI6BYEAR RELCHI7TYPE RELCHI7BYEAR RELCHI8TYPE RELCHI8BYEAR RELCHI9TYPE RELCHI9BYEAR RELCHI10TYPE RELCHI10BYEAR RELCHI11TYPE RELCHI11BYEAR RELCHI12TYPE RELCHI12BYEAR RELCHI13TYPE RELCHI13BYEAR RELCHI14TYPE RELCHI14BYEAR RELCHI15TYPE RELCHI15BYEAR RELCHI16TYPE RELCHI16BYEAR RELCHI17TYPE RELCHI17BYEAR RELCHI18TYPE RELCHI18BYEAR RELCHI19TYPE RELCHI19BYEAR RELCHI20TYPE RELCHI20BYEAR using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta", clear
sort ID year
by ID: gen double rooms_shift = ACTUALROOMS_[_n-1] if _n > 1 & inlist(year-year[_n-1],1,2)
foreach s in rooms own_raw moved_raw why_raw {
    merge 1:1 ID year using `src_`s'', keep(master match) nogenerate
}
* 3. Verification: rebuilt rooms vs shifted shelf column; own vs shelf HOMEOWN.
quietly count if CURRENT == 1 & !missing(rooms, rooms_shift)
local both = r(N)
quietly count if CURRENT == 1 & !missing(rooms, rooms_shift) & rooms != rooms_shift
local mismatch = r(N)
assert `mismatch' == 0
gen byte own = cond(own_raw == 1, 1, cond(own_raw == 5, 0, .)) if !missing(own_raw)
gen byte own_shelf = cond(HOMEOWN == 1, 1, cond(HOMEOWN == 2, 0, .)) if !missing(HOMEOWN)
quietly count if CURRENT == 1 & !missing(own, own_shelf)
local own_both = r(N)
quietly count if CURRENT == 1 & !missing(own, own_shelf) & own != own_shelf
local own_mismatch = r(N)
assert `own_mismatch' == 0
* 4. Outcome recodes (official codes: moved 1 yes / 5 no; why 3 = more space, 6 = neighbourhood; 9/98/99 DK-NA).
gen byte moved = cond(moved_raw == 1, 1, cond(moved_raw == 5, 0, .)) if !missing(moved_raw)
gen byte moved_space = .
replace moved_space = 0 if moved == 0
replace moved_space = 1 if moved == 1 & why_raw == 3
replace moved_space = 0 if moved == 1 & inlist(why_raw, 1, 2, 4, 5, 6, 7, 8)
gen byte moved_nbhd = .
replace moved_nbhd = 0 if moved == 0
replace moved_nbhd = 1 if moved == 1 & why_raw == 6
replace moved_nbhd = 0 if moved == 1 & inlist(why_raw, 1, 2, 3, 4, 5, 7, 8)
* Non-room codes to missing by interview year; zero retained.
quietly count if CURRENT == 1 & ((year <= 1984 & rooms == 9) | (inrange(year,1985,1993) & rooms == 99) | (year >= 1994 & inlist(rooms,98,99)))
local codes = r(N)
replace rooms = . if year <= 1984 & rooms == 9
replace rooms = . if inrange(year,1985,1993) & rooms == 99
replace rooms = . if year >= 1994 & inlist(rooms,98,99)
* 5. Person-level constructions before any restriction.
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
keep ID year AGEREP EDUYEAR rooms rooms_shift own moved moved_space moved_nbhd bio_first_year relchi1_year relchirep_max n_current_fids woman current rel hhid fid iw year_entry_adult
sort ID year
isid ID year
tempfile base
save `base'

tempname C S
postfile `C' str8 variant str8 window double estimate standard_error using "`out'/coef.dta", replace
postfile `S' str8 variant int baseline_lo int baseline_hi str8 headline double headline_effect headline_se from_m3m2 observations clusters treated_ids control_ids ref_mean seconds using "`out'/summary.dta", replace
foreach v in base r21 b4 b6 shift ctrl {
    timer clear 2
    timer on 2
    use `base', clear
    if "`v'" == "shift" {
        replace rooms = rooms_shift
        replace rooms = . if year <= 1984 & rooms == 9
        replace rooms = . if inrange(year,1985,1993) & rooms == 99
        replace rooms = . if year >= 1994 & inlist(rooms,98,99)
    }
    local blo = cond("`v'"=="r21", -2, cond("`v'"=="b4", -5, cond("`v'"=="b6", -7, -3)))
    local bhi = cond("`v'"=="r21", -1, cond("`v'"=="b4", -4, cond("`v'"=="b6", -6, -2)))
    keep if current & !missing(rooms, AGEREP, EDUYEAR) & AGEREP >= 18
    keep if !missing(iw) & iw > 0
    gen double f_c_y = bio_first_year
    drop if missing(f_c_y) & !(relchirep_max == 0)
    drop if !missing(f_c_y) & f_c_y < year_entry_adult
    gen byte control = missing(f_c_y) & relchirep_max == 0
    * A2h sample, always defined at -3/-2.
    gen byte hs_row = inlist(rel,1,2) & inrange(year - f_c_y, -3, -2)
    bysort ID: egen byte hs_base = max(hs_row)
    gen byte base_row = inrange(year - f_c_y, -3, -2)
    bysort ID: egen byte has_base = max(base_row)
    keep if control | (has_base & hs_base == 1)
    drop hs_row hs_base base_row has_base
    if "`v'" == "ctrl" {
        quietly summarize f_c_y
        replace control = f_c_y == r(max)
    }
    gen double K = year - f_c_y
    * Cohort support in the reference window used.
    gen byte in_reference = inrange(K, `blo', `bhi')
    bysort f_c_y: egen long reference_rows = total(in_reference)
    gen byte treated = !missing(f_c_y) & control == 0
    drop if treated & reference_rows == 0
    drop in_reference reference_rows treated
    local b = `bhi'
    gen byte Wleft = K <= `b'-6 & !missing(K)
    gen byte Wm2 = inrange(K,`b'-5,`b'-4)
    gen byte Wm1 = inrange(K,`b'-3,`b'-2)
    gen byte Wp1 = inrange(K,`b'+1,`b'+2)
    gen byte Wp2 = inrange(K,`b'+3,`b'+4)
    gen byte Wp3 = inrange(K,`b'+5,`b'+6)
    gen byte Wp4 = inrange(K,`b'+7,`b'+8)
    gen byte Wp5 = inrange(K,`b'+9,`b'+10)
    gen byte Wp6 = inrange(K,`b'+11,`b'+12)
    gen byte Wright = K >= `b'+13 & !missing(K)
    local dummies Wleft Wm2 Wm1 Wp1 Wp2 Wp3 Wp4 Wp5 Wp6 Wright
    local headline = cond("`v'"=="r21", "Wp2", "Wp" + string((4-`bhi')/2))
    * Window that covers -3/-2 when the reference is earlier.
    local m3 = cond("`v'"=="b4", "Wp1", cond("`v'"=="b6", "Wp2", ""))
    eventstudyinteract rooms `dummies' [pw=iw], vce(cluster ID) absorb(ID year) ///
        cohort(f_c_y) control_cohort(control) covariates(i.AGEREP i.EDUYEAR)
    matrix bb = e(b_iw)
    matrix V = e(V_iw)
    local n = e(N)
    local nc = e(N_clust)
    local names : colnames bb
    forvalues i = 1/`=colsof(bb)' {
        local w : word `i' of `names'
        post `C' ("`v'") ("`w'") (bb[1,`i']) (sqrt(V[`i',`i']))
    }
    local h = colnumb("bb","`headline'")
    local fm = .
    if "`m3'" != "" {
        local c3 = colnumb("bb","`m3'")
        local fm = bb[1,`h'] - bb[1,`c3']
    }
    quietly summarize rooms [aw=iw] if e(sample) & inrange(K,`blo',`bhi')
    local rm = r(mean)
    egen byte ttag = tag(ID) if e(sample) & !control & !missing(f_c_y)
    quietly count if ttag
    local tids = r(N)
    egen byte ctag = tag(ID) if e(sample) & control
    quietly count if ctag
    local cids = r(N)
    timer off 2
    quietly timer list 2
    post `S' ("`v'") (`blo') (`bhi') ("`headline'") (bb[1,`h']) (sqrt(V[`h',`h'])) (`fm') (`n') (`nc') (`tids') (`cids') (`rm') (r(t2))
    di "ROBUST_PASS `v' N=`n' `headline'=" bb[1,`h']
}
postclose `C'
postclose `S'
use "`out'/coef.dta", clear
format estimate standard_error %24.17g
export delimited using "`out'/coefficients.csv", replace
use "`out'/summary.dta", clear
format headline_effect headline_se from_m3m2 ref_mean %24.17g
export delimited using "`out'/summary.csv", replace
erase "`out'/coef.dta"
erase "`out'/summary.dta"
di "ROBUST_DONE"
