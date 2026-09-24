clear all
set more off
set processors 4
version 17.0
log using "/tmp/psid_rooms_v2_20260924/prepare.log", replace text
timer clear 1
timer on 1
* 1. Rooms from the official year-specific PSID variables (crosswalk of Sept 12).
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
gen double rooms_raw_source = rooms
tempfile rooms_source
save `rooms_source'
* 2. Shelf columns.
use ID year AGEREP EDUYEAR SEX REL CURRENT HHID FID IW DEATHYEAR RELCHIREP RELCHI1BYEAR ACTUALROOMS_ RELCHI1TYPE RELCHI1BYEAR RELCHI2TYPE RELCHI2BYEAR RELCHI3TYPE RELCHI3BYEAR RELCHI4TYPE RELCHI4BYEAR RELCHI5TYPE RELCHI5BYEAR RELCHI6TYPE RELCHI6BYEAR RELCHI7TYPE RELCHI7BYEAR RELCHI8TYPE RELCHI8BYEAR RELCHI9TYPE RELCHI9BYEAR RELCHI10TYPE RELCHI10BYEAR RELCHI11TYPE RELCHI11BYEAR RELCHI12TYPE RELCHI12BYEAR RELCHI13TYPE RELCHI13BYEAR RELCHI14TYPE RELCHI14BYEAR RELCHI15TYPE RELCHI15BYEAR RELCHI16TYPE RELCHI16BYEAR RELCHI17TYPE RELCHI17BYEAR RELCHI18TYPE RELCHI18BYEAR RELCHI19TYPE RELCHI19BYEAR RELCHI20TYPE RELCHI20BYEAR using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta", clear
sort ID year
by ID: gen double rooms_shift = ACTUALROOMS_[_n-1] if _n > 1 & inlist(year-year[_n-1],1,2)
merge 1:1 ID year using `rooms_source', keep(master match)
gen byte raw_matched = _merge == 3
drop _merge
* 3. Verification of the rebuilt rooms against the shifted shelf column.
quietly count if CURRENT == 1 & !missing(rooms, rooms_shift)
local both = r(N)
quietly count if CURRENT == 1 & !missing(rooms, rooms_shift) & rooms != rooms_shift
local mismatch = r(N)
quietly count if CURRENT == 1 & !missing(rooms) & missing(rooms_shift)
local raw_only = r(N)
quietly count if CURRENT == 1 & missing(rooms) & !missing(rooms_shift)
local shift_only = r(N)
assert `mismatch' == 0
* 4. Non-room codes to missing by interview year; zero (shared room) retained.
quietly count if CURRENT == 1 & ((year <= 1984 & rooms == 9) | (inrange(year,1985,1993) & rooms == 99) | (year >= 1994 & inlist(rooms,98,99)))
local codes = r(N)
replace rooms = . if year <= 1984 & rooms == 9
replace rooms = . if inrange(year,1985,1993) & rooms == 99
replace rooms = . if year >= 1994 & inlist(rooms,98,99)
quietly count if CURRENT == 1 & rooms == 0
local zeros = r(N)
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
* 6. Keep current person-years not after death.
keep if CURRENT == 1
drop if year > DEATHYEAR
keep ID year AGEREP EDUYEAR rooms rooms_shift raw_matched bio_first_year relchi1_year relchirep_max n_current_fids woman current rel hhid fid iw year_entry_adult
sort ID year
isid ID year
quietly count
local rows = r(N)
quietly count if !missing(rooms) & inlist(year,1969,1975,1976)
local recovered = r(N)
compress
save "/tmp/psid_rooms_v2_20260924/analysis_sample.dta", replace
timer off 1
quietly timer list 1
di "ROOMS_V2_PREPARATION_PASS rows=`rows' both=`both' mismatch=`mismatch' raw_only=`raw_only' shift_only=`shift_only' codes=`codes' zeros=`zeros' recovered_1969_1975_1976=`recovered' seconds=`=r(t1)'"
log close _all
