* 1968-2023 copy of sa_first_birth_hhcluster_prepare_v1.do (Oct 20 2026). Only changes: inputs are the 1968-2023 files from
* build_psid_a2h_inputs_2023.do (same items and SHELF definitions, 2021/2023 added), output paths, and
* why-moved codes 10-12 (COVID, evicted, foreclosed; new in 2021) counted as non-space, non-neighbourhood reasons.
* Copy of sa_own_first_birth_model_support_prepare.do for the household-clustered A2h event study (Oct 19 2026); only the scratch paths differ.
* Only changes: output paths (scratch, deleted after use) and four last-report fields used by first_birth_model_bridge.py.
clear all
set more off
set processors 4
version 17.0
log using "/private/tmp/claude-501/-Users-tommasodesanto-Desktop-Projects-Fertility-Fertility-Spring26/5d7c0112-2a2b-45f4-ab6e-9684dff40b7f/scratchpad/hhc2023/prepare.log", replace text
timer clear 1
timer on 1
* 1. Official year-specific family items, merged to the interview year.
use ID year rooms using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/extension_2023/items_1968_2023.dta" if !missing(rooms), clear
isid ID year
tempfile src_rooms
save `src_rooms'
use ID year own_raw using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/extension_2023/items_1968_2023.dta" if !missing(own_raw), clear
isid ID year
tempfile src_own_raw
save `src_own_raw'
use ID year moved_raw using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/extension_2023/items_1968_2023.dta" if !missing(moved_raw), clear
isid ID year
tempfile src_moved_raw
save `src_moved_raw'
use ID year why_raw using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/extension_2023/items_1968_2023.dta" if !missing(why_raw), clear
isid ID year
tempfile src_why_raw
save `src_why_raw'

* 2. Shelf columns.
use ID year AGEREP EDUYEAR SEX REL CURRENT HHID FID IW DEATHYEAR RELCHIREP RELCHI1BYEAR ACTUALROOMS_ HOMEOWN RELCHI1TYPE RELCHI1BYEAR RELCHI2TYPE RELCHI2BYEAR RELCHI3TYPE RELCHI3BYEAR RELCHI4TYPE RELCHI4BYEAR RELCHI5TYPE RELCHI5BYEAR RELCHI6TYPE RELCHI6BYEAR RELCHI7TYPE RELCHI7BYEAR RELCHI8TYPE RELCHI8BYEAR RELCHI9TYPE RELCHI9BYEAR RELCHI10TYPE RELCHI10BYEAR RELCHI11TYPE RELCHI11BYEAR RELCHI12TYPE RELCHI12BYEAR RELCHI13TYPE RELCHI13BYEAR RELCHI14TYPE RELCHI14BYEAR RELCHI15TYPE RELCHI15BYEAR RELCHI16TYPE RELCHI16BYEAR RELCHI17TYPE RELCHI17BYEAR RELCHI18TYPE RELCHI18BYEAR RELCHI19TYPE RELCHI19BYEAR RELCHI20TYPE RELCHI20BYEAR using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/extension_2023/shelf_1968_2023.dta", clear
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
replace moved_space = 0 if moved == 1 & inlist(why_raw, 1, 2, 4, 5, 6, 7, 8, 10, 11, 12)
gen byte moved_nbhd = .
replace moved_nbhd = 0 if moved == 0
replace moved_nbhd = 1 if moved == 1 & why_raw == 6
replace moved_nbhd = 0 if moved == 1 & inlist(why_raw, 1, 2, 3, 4, 5, 7, 8, 10, 11, 12)
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
* Added for the H5 model-support run (Oct 18 2026), before current-person filtering, as in the Oct 10 bridge
* panel: last fertility report = last year with CURRENT == 1 and nonmissing RELCHIREP.
gen double frep_year = year if CURRENT == 1 & !missing(RELCHIREP)
bysort ID: egen double last_fertility_report_year = max(frep_year)
gen double frep_age = AGEREP if year == last_fertility_report_year
bysort ID: egen double last_fertility_report_age = max(frep_age)
gen double cur_year = year if CURRENT == 1
bysort ID: egen double last_current_year = max(cur_year)
gen double cur_age = AGEREP if year == last_current_year
bysort ID: egen double last_current_age = max(cur_age)
drop frep_year frep_age cur_year cur_age
keep if CURRENT == 1
drop if year > DEATHYEAR
keep ID year AGEREP EDUYEAR rooms rooms_shift own moved moved_space moved_nbhd bio_first_year relchi1_year relchirep_max n_current_fids woman current rel hhid fid iw year_entry_adult last_fertility_report_year last_fertility_report_age last_current_year last_current_age
sort ID year
isid ID year
quietly count
local rows = r(N)
quietly count if !missing(own)
local n_own = r(N)
quietly count if !missing(moved)
local n_moved = r(N)
quietly count if !missing(moved_space)
local n_space = r(N)
quietly summarize moved
local mean_moved = r(mean)
quietly summarize moved_space
local mean_space = r(mean)
quietly summarize own
local mean_own = r(mean)
compress
save "/private/tmp/claude-501/-Users-tommasodesanto-Desktop-Projects-Fertility-Fertility-Spring26/5d7c0112-2a2b-45f4-ab6e-9684dff40b7f/scratchpad/hhc2023/analysis_sample.dta", replace
timer off 1
quietly timer list 1
di "ROOMS_V3_PREPARATION_PASS rows=`rows' rooms_both=`both' rooms_mismatch=`mismatch' own_both=`own_both' own_mismatch=`own_mismatch' codes=`codes' n_own=`n_own' n_moved=`n_moved' n_space=`n_space' mean_own=`mean_own' mean_moved=`mean_moved' mean_space=`mean_space' seconds=`=r(t1)'"
log close _all
