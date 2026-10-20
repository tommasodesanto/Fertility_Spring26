* 1968-2023 copy of sa_first_birth_descriptives_v1.do (Oct 20 2026). Only changes: inputs are the 1968-2023 files from
* build_psid_a2h_inputs_2023.do (same items and SHELF definitions, 2021/2023 added), output paths, and
* why-moved codes 10-12 (COVID, evicted, foreclosed; new in 2021) counted as non-space, non-neighbourhood reasons.
* Descriptive statistics for the A2h first-birth event-study sample (Oct 19, 2026).
* First-time parents: weighted (PSID longitudinal weight) means in the reference
* window -3/-2. Outcome rows use each outcome's own A2h sample (filter on that
* outcome's non-missing rows, as in the v3 event studies); demographic rows use
* the rooms A2h sample. Never-parents: all person-years in the same samples, and
* an age-matched version restricted to the 5th-95th percentile of first-time
* parents' reference-window age.
* Preparation copied from prepare_first_birth_outcomes_v3.do.
clear all
set more off
set processors 4
version 17.0
local out "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/sa_first_birth_descriptives_2023"
capture mkdir "`out'"
local cms "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta"
describe ACTUALROOMS_ RELCHI1TYPE RELCHI1BYEAR RELCHIREP HOMEOWN using "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/extension_2023/shelf_1968_2023.dta"
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
keep if CURRENT == 1
drop if year > DEATHYEAR
keep ID year AGEREP EDUYEAR rooms rooms_shift own moved moved_space moved_nbhd bio_first_year relchi1_year relchirep_max n_current_fids woman current rel hhid fid iw year_entry_adult
sort ID year
isid ID year
gen double age = AGEREP
gen double edu = EDUYEAR
gen double svy = year
tempfile base
save `base'

tempname P
tempfile meanfile
postfile `P' str12 variable str12 basis str20 grp double mean rows persons using `meanfile', replace
foreach basis in rooms own moved moved_space moved_nbhd {
    use `base', clear
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
    gen byte ref = !control & inrange(year - f_c_y, -3, -2)
    if "`basis'" == "rooms" {
        quietly count
        local n_all = r(N)
        egen byte pt = tag(ID)
        quietly count if pt
        local p_all = r(N)
        quietly count if pt & control
        local p_ctrl = r(N)
        quietly count if pt & !control
        local p_trt = r(N)
        quietly summarize year
        local y0 = r(min)
        local y1 = r(max)
        quietly summarize year if !control
        local ty0 = r(min)
        local ty1 = r(max)
        quietly summarize bio_first_year if !control
        local b0 = r(min)
        local b1 = r(max)
        di "SAMPLE rooms A2h: person-years=`n_all' persons=`p_all' treated=`p_trt' controls=`p_ctrl' years=`y0'-`y1' treated_years=`ty0'-`ty1' first_birth_years=`b0'-`b1'"
        drop pt
        _pctile age [aw=iw] if ref, p(5 95)
        local a5 = r(r1)
        local a95 = r(r2)
        di "AGE_RANGE treated ref window p5=`a5' p95=`a95'"
        local vars "rooms age edu woman svy"
    }
    else {
        local vars "`basis'"
    }
    foreach v of local vars {
        foreach g in parents_ref never_all never_agematch {
            if "`g'" == "parents_ref" local cond "ref"
            if "`g'" == "never_all" local cond "control"
            if "`g'" == "never_agematch" local cond "control & inrange(age, `a5', `a95')"
            quietly summarize `v' [aw=iw] if `cond'
            local m = r(mean)
            quietly count if `cond' & !missing(`v')
            local r = r(N)
            egen byte tg = tag(ID) if `cond' & !missing(`v')
            quietly count if tg == 1
            local p = r(N)
            drop tg
            post `P' ("`v'") ("`basis'") ("`g'") (`m') (`r') (`p')
        }
    }
    di "DESC_PASS `basis'"
}
postclose `P'
use `meanfile', clear
format mean %12.6g
export delimited using "`out'/descriptives.csv", replace
di "DESC_DONE"
