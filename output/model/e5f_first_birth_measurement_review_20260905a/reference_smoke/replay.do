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
local outdir  "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_first_birth_measurement_review_20260905a/reference_smoke"

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


do "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_first_birth_measurement_review_20260905a/reference_smoke/reference_check.do" toy "`outdir'"
