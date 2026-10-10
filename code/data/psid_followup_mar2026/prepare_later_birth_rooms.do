* Later-birth rooms event study: private panel preparation (October 10, 2026).
* Adds biological second- and third-birth years to the frozen first-birth
* bridge panel (A2h rules, rooms, weights and fertility-report fields unchanged).
* Birth order counts children: twins at order q give order q+1 the same year.
version 17.0
args bridge_panel outpanel outdir
set more off
local shelf "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local keepvars ID year RELCHIREP
forvalues c = 1/20 {
    local keepvars `keepvars' RELCHI`c'TYPE RELCHI`c'BYEAR
}
use `keepvars' using "`shelf'", clear
gen int nbio_row = 0
gen int nrec_row = 0
forvalues c = 1/20 {
    replace nbio_row = nbio_row + (RELCHI`c'TYPE == 1 & !missing(RELCHI`c'BYEAR))
    replace nrec_row = nrec_row + (!missing(RELCHI`c'TYPE))
}
bysort ID: egen int nbio_max = max(nbio_row)
bysort ID: egen int nbio_min = min(nbio_row)
quietly count if nbio_max != nbio_min
local varying_rows = r(N)
* Most complete history: the row with the most biological records, latest year on ties.
gsort ID -nbio_row -year
by ID: keep if _n == 1
keep ID nbio_row nrec_row RELCHI*TYPE RELCHI*BYEAR
reshape long RELCHI@TYPE RELCHI@BYEAR, i(ID) j(slot)
keep if RELCHITYPE == 1 & !missing(RELCHIBYEAR)
sort ID RELCHIBYEAR slot
by ID: gen int order = _n
keep if order <= 3
keep ID order RELCHIBYEAR
reshape wide RELCHIBYEAR, i(ID) j(order)
capture confirm variable RELCHIBYEAR3
rename RELCHIBYEAR1 hist_bio_b1
rename RELCHIBYEAR2 bio_second_year
rename RELCHIBYEAR3 bio_third_year
tempfile hist
save `hist'
* Child-record counts per person (all persons, including those without biological records).
use `keepvars' using "`shelf'", clear
gen int nbio_row = 0
gen int nrec_row = 0
forvalues c = 1/20 {
    replace nbio_row = nbio_row + (RELCHI`c'TYPE == 1 & !missing(RELCHI`c'BYEAR))
    replace nrec_row = nrec_row + (!missing(RELCHI`c'TYPE))
}
collapse (max) nbio_max=nbio_row nrec_max=nrec_row, by(ID)
merge 1:1 ID using `hist', nogenerate
tempfile persons
save `persons'
use "`bridge_panel'", clear
merge m:1 ID using `persons', keep(master match) nogenerate
* The first biological birth must agree with the frozen bridge construction.
quietly count if !missing(bio_first_year) & bio_first_year != hist_bio_b1
local b1_mismatch_rows = r(N)
quietly count if missing(bio_first_year) != missing(hist_bio_b1)
local b1_missing_mismatch_rows = r(N)
* Model age cells (four-year periods from age 18), same DOB proxy as the bridge.
gen double dob_c = year - AGEREP if inrange(AGEREP,18,100) & current == 1
bysort ID: egen double model_dob_proxy = median(dob_c)
drop dob_c
replace model_dob_proxy = . if !inrange(model_dob_proxy,1800,2082)
gen double model_age_index = floor((year - model_dob_proxy - 18)/4)
gen double b1_j = floor((bio_first_year - model_dob_proxy - 18)/4)
gen double b2_j = floor((bio_second_year - model_dob_proxy - 18)/4)
gen double b3_j = floor((bio_third_year - model_dob_proxy - 18)/4)
gen double censor_j = floor((last_fertility_report_age - 18)/4)
compress
save "`outpanel'", replace
clear
set obs 1
gen long varying_history_rows = `varying_rows'
gen long b1_mismatch_rows = `b1_mismatch_rows'
gen long b1_missing_mismatch_rows = `b1_missing_mismatch_rows'
export delimited using "`outdir'/preparation_receipt.csv", replace
di "LATER_BIRTH_PREP_PASS varying=`varying_rows' b1_mismatch=`b1_mismatch_rows' b1_missmismatch=`b1_missing_mismatch_rows'"
