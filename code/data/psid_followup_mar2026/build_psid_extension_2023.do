* Extend the PSID-SHELF columns used by the A2h first-birth event study to the 2021 and 2023 waves (Oct 20 2026).
* Applies PSID-SHELF's own construction rules (Construction_Files/Code, 1968-2019 beta) to the 2023 PSID release:
*   FAM2021ER, FAM2023ER, IND2023ER (cross-year individual 1968-2023), CAH85_23 (childbirth history 1985-2023).
* The same code builds 2017 and 2019 from the same raw files; those rows are compared with SHELF as a check.
* Writes (outside the repo, next to the PSID data):
*   relchi_cah85_23.dta      person-level child records (RELCHI1-20 TYPE/BYEAR, RELCHIREP), SHELF rules
*   shelf_rows_2021_2023.dta person-year SHELF columns for 2021 and 2023
*   raw_items_2021_2023.dta  person-year raw family items (rooms, own, moved, why moved) for 2021 and 2023
*   validation_2017_2019.csv agreement of the rebuilt 2017/2019 columns with SHELF
* Rule changes needed for the new release (disclosed):
*   CAH: SHELF's CAH106 and CAH108 are CAH116 and CAH118 in the 2023 release (same labels).
*   Deflator: SHELF's own index lookup (see section 3) extended to 2021 and 2023; no new price data needed.
version 17.0
clear all
set more off
local raw "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/raw_2021_2023"
local ext "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/extension_2023"
local shelf "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local cms "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta"
capture mkdir "`ext'"
log using "`ext'/build_psid_extension_2023.log", replace text

*------------------------------------------------------------------
* 1. Child records from CAH85_23 (SHELF "01 Collect relationship identifiers").
*------------------------------------------------------------------
infix CAH2 2-2 CAH3 3-6 CAH4 7-9 CAH15 30-33 CAH116 188-189 using "`raw'/cah85_23/CAH85_23.txt", clear
gen long fileorder = _n
gen double ID = CAH3*1000 + CAH4
bysort ID (fileorder): gen int k = _n
assert k <= 20
keep ID k CAH2 CAH15 CAH116
reshape wide CAH2 CAH15 CAH116, i(ID) j(k)
egen flagmax = rowmax(CAH116*)
gen tb = 0 if !missing(flagmax)
gen ta = 0 if !missing(flagmax)
forvalues i = 1/20 {
    capture confirm variable CAH2`i'
    if _rc {
        gen CAH2`i' = .
        gen CAH15`i' = .
        gen CAH116`i' = .
    }
    gen RELCHI`i'TYPE = CAH2`i' if inrange(CAH2`i',1,2)
    gen RELCHI`i'BYEAR = CAH15`i' if inrange(CAH15`i',1900,2099)
    gen rep`i' = CAH116`i' if inrange(CAH116`i',0,20)
    replace tb = rep`i' if CAH2`i'==1 & inrange(CAH116`i',1,20)
    replace tb = . if CAH2`i'==1 & missing(rep`i')
    replace ta = rep`i' if CAH2`i'==2 & inrange(CAH116`i',1,20)
    replace ta = . if CAH2`i'==2 & missing(rep`i')
}
gen RELCHIREP = tb + ta
keep ID RELCHIREP RELCHI*TYPE RELCHI*BYEAR
order ID RELCHIREP
compress
label data "Child records from PSID CAH 1985-2023, PSID-SHELF rules"
save "`ext'/relchi_cah85_23.dta", replace

*------------------------------------------------------------------
* 2. Individual file: person-level and 2017-2023 person-year items.
*------------------------------------------------------------------
infix ER30001 2-5 ER30002 6-8 ER32000 2057-2057 ER32050 2188-2191 ///
    fid2017 4870-4874 seq2017 4875-4876 relraw2017 4877-4878 age2017 4879-4881 edu2017 4980-4981 double wt2017 5238-5244 ///
    fid2019 5250-5254 seq2019 5255-5256 relraw2019 5257-5258 age2019 5259-5261 edu2019 5384-5385 double wt2019 5672-5678 ///
    fid2021 5684-5688 seq2021 5689-5690 relraw2021 5691-5692 age2021 5693-5695 edu2021 5818-5819 double wt2021 6107-6113 ///
    fid2023 6119-6123 seq2023 6124-6125 relraw2023 6126-6127 age2023 6128-6130 edu2023 6253-6254 double wt2023 6542-6548 ///
    using "`raw'/ind2023er/IND2023ER.txt", clear
gen double ID = ER30001*1000 + ER30002
isid ID
* sex (SHELF 01 demographic (1)).
gen SEX = cond(ER32000==1, 1, cond(ER32000==2, 2, .))
* death year (SHELF 01 demographic (4)-(6), 03 demographic (3)).
gen deathsy = ER32050 if inrange(ER32050, 1967, 2099)
gen drange = ER32050 if inrange(ER32050, 1, 1966) | inrange(ER32050, 2100, 9998)
gen f2 = floor(drange/100) if !missing(drange)
gen l2 = mod(drange, 100) if !missing(drange)
gen dll = cond(inrange(f2,67,99), f2+1900, cond(inrange(f2,1,66), f2+2000, .))
gen dul = cond(inrange(l2,67,99), l2+1900, cond(inrange(l2,1,66), l2+2000, .))
gen DEATHYEAR = deathsy
replace DEATHYEAR = round((dll+dul)*0.5, 1) if missing(deathsy) & !missing(dll, dul)
drop ER32000 ER32050 deathsy drange f2 l2 dll dul ER30001 ER30002
reshape long fid seq relraw age edu wt, i(ID) j(year)
* survey identifiers (SHELF 01/03 survey identifiers).
gen FID = fid if fid != 0
gen CURRENT = .
replace CURRENT = 0 if seq == 0 | inrange(seq, 21, 99)
replace CURRENT = 1 if inrange(seq, 1, 20)
gen REL = .
replace REL = 1 if inlist(relraw, 10)
replace REL = 2 if inlist(relraw, 20, 22, 90, 92)
replace REL = 3 if inlist(relraw, 30, 33, 35)
replace REL = 4 if inlist(relraw, 88, 38, 83, 37, 50, 57, 58, 60, 65, 66, 67, 68, 69, 40, 47, 48, 70, 71, 74, 75, 72, 73, 95, 96, 97, 98)
gen IW = wt if wt != 0
gen AGEREP = age if inrange(age, 1, 199)
gen eduyearrep = edu if inrange(edu, 1, 20)
tempfile ind
save `ind'

*------------------------------------------------------------------
* 3. Family items: 2021 and 2023 from the family files, 2017 and 2019 from the main-study file.
*------------------------------------------------------------------
infix fid 2-6 state 7-8 rooms 61-62 own_raw 63-63 double value 64-71 double mor1 107-113 double mor2 150-156 moved_raw 301-301 why_raw 308-309 hhidraw 7870-7874 using "`raw'/fam2021er/FAM2021ER.txt", clear
gen year = 2021
tempfile f21
save `f21'
infix fid 2-6 state 7-8 rooms 62-63 own_raw 64-64 double value 65-72 double mor1 107-113 double mor2 150-156 moved_raw 283-283 why_raw 290-291 hhidraw 7624-7628 using "`raw'/fam2023er/FAM2023ER.txt", clear
gen year = 2023
append using `f21'
isid year fid
tempfile fnew
save `fnew'
use `ind' if inlist(year, 2021, 2023), clear
merge m:1 year fid using `fnew', keep(master match) nogenerate
tempfile pnew
save `pnew'
use ID ER66003 ER66029 ER66030 ER66031 ER66051 ER66072 ER66156 ER66159 ER71560 ER72003 ER72029 ER72030 ER72031 ER72051 ER72072 ER72156 ER72159 ER77621 using "`cms'", clear
rename (ER66003 ER66029 ER66030 ER66031 ER66051 ER66072 ER66156 ER66159 ER71560) (state2017 rooms2017 own_raw2017 value2017 mor12017 mor22017 moved_raw2017 why_raw2017 hhidraw2017)
rename (ER72003 ER72029 ER72030 ER72031 ER72051 ER72072 ER72156 ER72159 ER77621) (state2019 rooms2019 own_raw2019 value2019 mor12019 mor22019 moved_raw2019 why_raw2019 hhidraw2019)
reshape long state rooms own_raw value mor1 mor2 moved_raw why_raw hhidraw, i(ID) j(year)
tempfile fold
save `fold'
use `ind' if inlist(year, 2017, 2019), clear
merge 1:1 ID year using `fold', keep(master match) nogenerate
append using `pnew'
sort ID year

* housing and state (SHELF 01 housing (1), 01 wealth (4)-(6), 03 wealth (1), 01 demographic (8)).
gen HHID = hhidraw if hhidraw != 0
gen HOMEOWN = cond(own_raw==1, 1, cond(own_raw==5, 2, cond(own_raw==8, 3, .)))
gen GEOSTATE = state if inrange(state,1,51) | inrange(state,61,66)
replace GEOSTATE = 60 if state == 0
gen homevaluen = value if inrange(value, 0, 9999997)
gen m1 = mor1 if inrange(mor1, 0, 9999997)
gen m2 = mor2 if inrange(mor2, 0, 9999997)
egen homemortotn = rowtotal(m1 m2) if !missing(m1) | (missing(m1) & m2 != 0 & !missing(m2)), missing
replace homemortotn = . if (missing(m1) & missing(m2)) | (missing(m1) & m2 == 0)
* PCE price index as SHELF 07 applies it: word (year-1966) and word (2022-1966) of a list that starts in 1959,
* i.e. index(year-8) and index(2014). Checked: this reproduces SHELF HOMEVALUER in 2017 and 2019 exactly.
* The "real 2022 USD" label in SHELF is therefore not what the numbers are; kept for comparability.
gen double pce = .
replace pce = 118.08571 if year == 2017
replace pce = 123.24161 if year == 2019
replace pce = 127.19513 if year == 2021
replace pce = 129.21168 if year == 2023
gen double HOMEVALUER = homevaluen * (128.9754/pce)
gen double HOMEMORTOTR = homemortotn * (128.9754/pce)
replace HOMEVALUER = 9999999 if homevaluen == 9999999
* family-level items missing for persons not in the family unit or an institution (SHELF 08).
foreach v of varlist HOMEOWN HOMEVALUER HOMEMORTOTR AGEREP {
    replace `v' = . if !inrange(seq, 1, 20) & !inrange(seq, 51, 59)
}
tempfile built
save `built'

*------------------------------------------------------------------
* 4. Education: report, else last reported value carried forward when FID and REL exist (SHELF 03 education (1)).
*------------------------------------------------------------------
use ID year EDUYEAR using "`shelf'" if year <= 2019 & !missing(EDUYEAR), clear
bysort ID (year): keep if _n == _N
keep ID year EDUYEAR
rename (year EDUYEAR) (edu_last_year edu_last)
tempfile edulast
save `edulast'
use ID year EDUYEAR using "`shelf'" if year <= 2015 & !missing(EDUYEAR), clear
bysort ID (year): keep if _n == _N
keep ID EDUYEAR
rename EDUYEAR edu_last_2015
tempfile edulast15
save `edulast15'
use `built', clear
merge m:1 ID using `edulast', keep(master match) nogenerate
merge m:1 ID using `edulast15', keep(master match) nogenerate
sort ID year
* carried value entering each year: 2017 from SHELF through 2015; 2019 from 2017; 2021 from SHELF through 2019; 2023 from 2021.
gen double carry = .
by ID: replace carry = edu_last_2015 if year == 2017
by ID: replace carry = cond(!missing(eduyearrep[_n-1]), eduyearrep[_n-1], carry[_n-1]) if year == 2019 & year[_n-1] == 2017
by ID: replace carry = edu_last if year == 2021
by ID: replace carry = cond(!missing(eduyearrep[_n-1]), eduyearrep[_n-1], carry[_n-1]) if year == 2023 & year[_n-1] == 2021
gen EDUYEAR = eduyearrep
replace EDUYEAR = carry if missing(EDUYEAR) & !missing(carry) & !missing(FID) & !missing(relraw) & relraw != 0
tempfile all
save `all'

*------------------------------------------------------------------
* 5. Check: rebuilt 2017 and 2019 columns against SHELF.
*------------------------------------------------------------------
use ID year AGEREP EDUYEAR SEX REL CURRENT HHID FID IW DEATHYEAR HOMEOWN HOMEVALUER HOMEMORTOTR GEOSTATE using "`shelf'" if inlist(year, 2017, 2019), clear
rename (AGEREP EDUYEAR SEX REL CURRENT HHID FID IW DEATHYEAR HOMEOWN HOMEVALUER HOMEMORTOTR GEOSTATE) s_=
merge 1:1 ID year using `all', keep(match) nogenerate
tempname V
postfile `V' str12 variable int year double both_missing agree differ only_shelf only_new using "`ext'/validation.dta", replace
foreach v in AGEREP EDUYEAR SEX REL CURRENT HHID FID IW DEATHYEAR HOMEOWN HOMEVALUER HOMEMORTOTR GEOSTATE {
    foreach y in 2017 2019 {
        local tol = cond(inlist("`v'", "HOMEVALUER", "HOMEMORTOTR"), 0.5, cond("`v'"=="IW", 0.006, 0))
        quietly count if year==`y' & missing(s_`v') & missing(`v')
        local bm = r(N)
        quietly count if year==`y' & !missing(s_`v', `v') & abs(s_`v' - `v') <= `tol'
        local ag = r(N)
        quietly count if year==`y' & !missing(s_`v', `v') & abs(s_`v' - `v') > `tol'
        local df = r(N)
        quietly count if year==`y' & !missing(s_`v') & missing(`v')
        local os = r(N)
        quietly count if year==`y' & missing(s_`v') & !missing(`v')
        local on = r(N)
        post `V' ("`v'") (`y') (`bm') (`ag') (`df') (`os') (`on')
    }
}
postclose `V'
* current persons only, the rows the analysis uses.
tempname W
postfile `W' str12 variable int year double agree differ only_shelf only_new using "`ext'/validation_current.dta", replace
foreach v in AGEREP EDUYEAR SEX REL HHID FID IW DEATHYEAR HOMEOWN HOMEVALUER HOMEMORTOTR GEOSTATE {
    foreach y in 2017 2019 {
        local tol = cond(inlist("`v'", "HOMEVALUER", "HOMEMORTOTR"), 0.5, cond("`v'"=="IW", 0.006, 0))
        quietly count if year==`y' & s_CURRENT==1 & !missing(s_`v', `v') & abs(s_`v' - `v') <= `tol'
        local ag = r(N)
        quietly count if year==`y' & s_CURRENT==1 & !missing(s_`v', `v') & abs(s_`v' - `v') > `tol'
        local df = r(N)
        quietly count if year==`y' & s_CURRENT==1 & !missing(s_`v') & missing(`v')
        local os = r(N)
        quietly count if year==`y' & s_CURRENT==1 & missing(s_`v') & !missing(`v')
        local on = r(N)
        post `W' ("`v'") (`y') (`ag') (`df') (`os') (`on')
    }
}
postclose `W'
gen double ratio_value = s_HOMEVALUER / HOMEVALUER if s_HOMEVALUER > 0 & HOMEVALUER > 0 & !missing(s_HOMEVALUER, HOMEVALUER)
tabstat ratio_value, by(year) statistics(min p50 max n) format(%12.6f)
di "why-moved codes among current movers, 2017 and 2019 (compare with 2021 and 2023 below)"
tab why_raw year if s_CURRENT == 1 & moved_raw == 1, m
use "`ext'/validation.dta", clear
export delimited using "`ext'/validation_2017_2019.csv", replace
list, noobs sep(0)
use "`ext'/validation_current.dta", clear
export delimited using "`ext'/validation_2017_2019_current.csv", replace
list, noobs sep(0)

*------------------------------------------------------------------
* 6. Outputs for 2021 and 2023.
*------------------------------------------------------------------
use `all' if inlist(year, 2021, 2023), clear
tab why_raw year if CURRENT == 1 & moved_raw == 1, m
preserve
    keep ID year rooms own_raw moved_raw why_raw
    keep if !missing(rooms) | !missing(own_raw) | !missing(moved_raw) | !missing(why_raw)
    compress
    label data "PSID raw family items 2021 and 2023 (A8 rooms, A19 own/rent, A49 moved, A50 why moved), by person"
    save "`ext'/raw_items_2021_2023.dta", replace
    tab year
    tab why_raw year, m
    tab moved_raw year, m
    tab own_raw year, m
    summarize rooms, detail
restore
keep ID year AGEREP EDUYEAR SEX REL CURRENT HHID FID IW DEATHYEAR HOMEOWN HOMEVALUER HOMEMORTOTR GEOSTATE
compress
label data "PSID-SHELF columns rebuilt for 2021 and 2023 with SHELF rules"
save "`ext'/shelf_rows_2021_2023.dta", replace
tab year CURRENT, m
summarize
log close
