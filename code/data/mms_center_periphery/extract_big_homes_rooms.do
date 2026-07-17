clear all
set more off
version 17.0

local dta "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/extract27.dta"
local outdir "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/mms_center_periphery/output_big_homes_rooms"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/extract_big_homes_rooms.log", replace text

* ==================================================================
* Rooms-based analog of the Redfin big-homes chart.
* Filter: owner-occupied, head-of-household, ACS 2022-2023, rooms >= 6.
* Cells: generation x household composition.
* This matches the model H_own grid (rooms units) at H >= 6.
* ==================================================================

use year age perwt gq pernum ownershp rooms serial nchild yngch ///
    using "`dta'", clear

keep if year >= 2022 & year <= 2023
keep if inlist(gq, 1, 2)

bysort year serial: gen numprec = _N

keep if pernum == 1
keep if ownershp == 1

* ACS rooms is the literal count; drop missing/N/A.
keep if !missing(rooms) & rooms > 0 & rooms < 99
keep if rooms >= 6

di _n "=== After heads + owner + rooms>=6 filter ==="
count

gen birth_year = year - age
gen str20 generation = ""
replace generation = "Silent"      if inrange(birth_year, 1928, 1945)
replace generation = "Boomer"      if inrange(birth_year, 1946, 1964)
replace generation = "Gen X"       if inrange(birth_year, 1965, 1980)
replace generation = "Millennial"  if inrange(birth_year, 1981, 1996)
replace generation = "Gen Z"       if inrange(birth_year, 1997, 2012)
drop if generation == ""

gen byte has_minor_children = (nchild > 0 & yngch < 18 & yngch != 99) if !missing(nchild)
gen str30 hh_type = ""
replace hh_type = "1-2 adults"             if has_minor_children == 0 & numprec <= 2
replace hh_type = "3+ people, no children" if has_minor_children == 0 & numprec >= 3
replace hh_type = "with own children"      if has_minor_children == 1
drop if hh_type == ""

di _n "=== Final sample ==="
count

preserve
    collapse (sum) cell_weight = perwt, by(generation hh_type)
    egen total_weight = total(cell_weight)
    gen share = cell_weight / total_weight
    gsort -share
    list, sepby(generation) noobs abbrev(25)
    export delimited using "`outdir'/big_homes_rooms_cells.csv", replace
restore

di _n "=== Done. Output: `outdir' ==="
log close
