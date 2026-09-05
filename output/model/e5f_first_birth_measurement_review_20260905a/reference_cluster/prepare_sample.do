clear all
set more off
set processors 8
version 17.0

* Corrected first-birth housing-space target.
* Unit: one weighted woman (reference person or spouse/partner) per
* single-family-unit household-year.
* Outcome: ACTUALROOMS_ shifted forward by one observed interview within person,
* because the PSIDSHELF extraction attaches the next-wave rooms item to the
* preceding row. Event-study baseline is K=-2; controls are women whose full
* relationship history confirms no children.

local project "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"
local source  "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local outroot "`project'/code/data/psid_followup_mar2026/output"
local outdir  "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_first_birth_measurement_review_20260905a/reference_cluster"

cap mkdir "`outroot'"
cap mkdir "`outdir'"
capture confirm file "`outdir'/target_receipt.csv"
if !_rc {
    di as error "Refusing to overwrite completed target receipt: `outdir'/target_receipt.csv"
    exit 602
}

capture log close _all
log using "`outdir'/sa_rooms_first_birth_household_aligned_v1.log", replace text
timer clear 1
timer on 1

cap which eventstudyinteract
if _rc {
    di as error "eventstudyinteract is not installed. Aborting."
    exit 198
}
cap which svmat2
if _rc {
    di as error "svmat2 is not installed. Aborting."
    exit 198
}

di as text "Loading the PSID shelf and aligning ACTUALROOMS_ to its interview..."
use ID FID HHID year AGEREP EDUYEAR SEX DEATHYEAR REL CURRENT RELCHIREP ///
    RELCHI1TYPE RELCHI1BYEAR RELCHI2TYPE RELCHI2BYEAR ///
    RELCHI3TYPE RELCHI3BYEAR RELCHI4TYPE RELCHI4BYEAR ///
    RELCHI5TYPE RELCHI5BYEAR RELCHI6TYPE RELCHI6BYEAR ///
    RELCHI7TYPE RELCHI7BYEAR RELCHI8TYPE RELCHI8BYEAR ///
    RELCHI9TYPE RELCHI9BYEAR RELCHI10TYPE RELCHI10BYEAR ///
    RELCHI11TYPE RELCHI11BYEAR RELCHI12TYPE RELCHI12BYEAR ///
    RELCHI13TYPE RELCHI13BYEAR RELCHI14TYPE RELCHI14BYEAR ///
    RELCHI15TYPE RELCHI15BYEAR RELCHI16TYPE RELCHI16BYEAR ///
    RELCHI17TYPE RELCHI17BYEAR RELCHI18TYPE RELCHI18BYEAR ///
    RELCHI19TYPE RELCHI19BYEAR RELCHI20TYPE RELCHI20BYEAR ///
    ACTUALROOMS_ IW using "`source'", clear

sort ID year
by ID: gen double rooms = ACTUALROOMS_[_n-1] if _n > 1
by ID: gen double rooms_source_year = year[_n-1] if _n > 1

* Fertility history is a person-level object and is constructed before any
* sample restriction or household reporter selection. TYPE==1 is biological.
gen double first_biological_birth_candidate = .
forvalues child = 1/20 {
    replace first_biological_birth_candidate = RELCHI`child'BYEAR ///
        if RELCHI`child'TYPE == 1 & !missing(RELCHI`child'BYEAR) & ///
        (missing(first_biological_birth_candidate) | ///
         RELCHI`child'BYEAR < first_biological_birth_candidate)
}
bysort ID: egen double first_child_year = min(first_biological_birth_candidate)
bysort ID: egen double max_children_reported = max(RELCHIREP)
drop first_biological_birth_candidate

* HHID is a physical dwelling and can contain several PSID family units.
* Count current FIDs using every current member before selecting women.
egen byte current_fid_tag = tag(HHID year FID) ///
    if CURRENT == 1 & !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0
bysort HHID year: egen int n_current_fids = total(current_fid_tag)
drop current_fid_tag

* Missing/non-room codes follow the year to which the lagged value is aligned.
replace rooms = . if year <= 1984 & rooms == 9
replace rooms = . if inrange(year, 1985, 1993) & rooms == 99
replace rooms = . if year >= 1994 & inlist(rooms, 98, 99)

drop if year > DEATHYEAR
keep if CURRENT == 1
keep if !missing(AGEREP) & AGEREP >= 18
keep if !missing(IW) & IW > 0
keep if SEX == 2 & inlist(REL, 1, 2)
keep if !missing(HHID) & HHID > 0 & !missing(FID) & FID > 0

egen byte multi_fu_hhyear_tag = tag(HHID year) if n_current_fids > 1
quietly count if multi_fu_hhyear_tag
local excl_multi_hh = r(N)
quietly count if n_current_fids > 1
local excl_multi_women = r(N)
drop multi_fu_hhyear_tag
keep if n_current_fids == 1

* Household outcomes enter once per household-year. Reference women receive
* deterministic priority in the remaining within-family duplicates. Reporter
* selection does not condition on outcome or education availability.
gen byte household_priority = REL != 1
sort HHID year household_priority ID
by HHID year: gen int household_women_before_dedup = _N
by HHID year: keep if _n == 1
isid HHID year

gen double rooms_alignment_gap_years = year - rooms_source_year
quietly count if !missing(rooms) & !inlist(rooms_alignment_gap_years, 1, 2)
local excl_bad_room_gap = r(N)
replace rooms = . if !inlist(rooms_alignment_gap_years, 1, 2)
drop rooms_alignment_gap_years
drop if missing(rooms) | missing(EDUYEAR)

bysort ID: egen double year_entry = min(year)
drop if !missing(first_child_year) & first_child_year < year_entry

gen byte untimed_known_parent = ///
    missing(first_child_year) & !missing(max_children_reported) & ///
    max_children_reported > 0
gen byte unknown_child_hist = ///
    missing(first_child_year) & missing(max_children_reported)
egen byte excl_known_tag = tag(ID) if untimed_known_parent
quietly count if excl_known_tag
local excl_known_ids = r(N)
egen byte excl_unknown_tag = tag(ID) if unknown_child_hist
quietly count if excl_unknown_tag
local excl_unknown_ids = r(N)
drop excl_known_tag excl_unknown_tag
drop if untimed_known_parent | unknown_child_hist
gen byte never_treated = missing(first_child_year) & max_children_reported == 0
assert never_treated == 1 if missing(first_child_year)
gen double K = year - first_child_year

cap drop L*event F*event
forvalues k = 0/10 {
    gen byte L`k'event = K == `k'
}
gen byte L11event = K > 10 & !missing(K)
gen byte F1event = K == -1
gen byte F3event = K == -3
gen byte F4event = K == -4
gen byte F5event = K == -5
gen byte F6event = K == -6
gen byte F7event = K <= -7 & !missing(K)

quietly count
local input_obs = r(N)
egen byte input_id_tag = tag(ID)
quietly count if input_id_tag
local input_ids = r(N)
drop input_id_tag
quietly count if never_treated
local never_obs = r(N)
egen byte treated_id_tag = tag(ID) if !never_treated
quietly count if treated_id_tag
local treated_ids = r(N)
drop treated_id_tag
quietly summarize household_women_before_dedup
local max_women_before_dedup = r(max)


assert `input_obs'==49872
assert `input_ids'==4527
assert !missing(ID, year, AGEREP, EDUYEAR, IW, rooms)
keep ID year first_child_year never_treated K IW AGEREP EDUYEAR rooms L*event F*event
compress
save "/var/folders/xs/ttwplgz97gvdr_6gkg5r_r9h0000gn/T/e5f-reference-sample-2yjv_is7/analysis_sample.dta", replace
sysdir
foreach lib in lftools.mlib lreghdfe.mlib lmoremata.mlib {
    capture findfile `lib'
    if !_rc di as result "PINNED_LIBRARY `lib' `r(fn)'"
}
di as result "REFERENCE_SAMPLE_PREPARATION_PASS"
log close _all
