clear all
set more off
version 17.0

* Synthetic fixture for the portable custom-driver smoke.  This deliberately
* avoids the real PSID shelf and exercises the dependency, cohort, covariance,
* and export path on a small balanced panel.
args driver fixture outroot ado_root
if "`driver'" == "" | "`fixture'" == "" | "`outroot'" == "" | "`ado_root'" == "" {
    di as error "usage: do test_psid_housing_fullsample_stata_smoke.do driver fixture outroot ado_root"
    exit 198
}
capture log close _all
log using "`outroot'/smoke_setup.log", replace text

* Keep package installation in the persistent task root so every fresh Stata
* process can reuse the same verified ado path.
local task_ado "`ado_root'"
cap mkdir "`task_ado'"
sysdir set PLUS "`task_ado'"
foreach pkg in require ftools reghdfe avar ivreg2 moremata {
    capture which `pkg'
    if _rc ssc install `pkg', replace
}
capture which svmat2
if _rc {
    copy "http://fmwww.bc.edu/repec/bocode/s/svmat2.ado" "`task_ado'/svmat2.ado", replace
    copy "http://fmwww.bc.edu/repec/bocode/s/svmat2.sthlp" "`task_ado'/svmat2.sthlp", replace
}
capture confirm file "`task_ado'/l/lmoremata.mlib"
if _rc ssc install moremata, replace
capture which eventstudyinteract
if _rc net install eventstudyinteract, from("https://raw.githubusercontent.com/lsun20/EventStudyInteract/main") replace
mata: mata mlib index
foreach pkg in require eventstudyinteract ftools reghdfe avar ivreg2 svmat2 {
    capture which `pkg'
    if _rc {
        di as error "required Stata package unavailable after setup: `pkg'"
        exit 199
    }
}
capture confirm file "`task_ado'/l/lmoremata.mlib"
if _rc {
    di as error "required moremata Mata library unavailable after setup"
    exit 199
}

set obs 1800
gen long ID = ceil(_n / 30)
gen long year = 1990 + mod(_n - 1, 30)
gen long HHID = ID
gen long FID = ID
gen byte CURRENT = 1
gen byte SEX = 2
gen byte REL = 1
gen int AGEREP = 25 + mod(ID, 12) + (year - 1990)
gen int EDUYEAR = 12 + mod(ID, 5)
gen int DEATHYEAR = 2100
gen double IW = 1 + mod(ID, 7) / 10
gen byte HOMEOWN = mod(ID + year, 2)
gen byte RELCHIREP = cond(mod(ID, 4) == 0, 0, 1)
gen double first_year = cond(mod(ID, 4) == 0, ., 2003 + mod(ID, 4))
gen double second_year = cond(mod(ID, 4) == 0, ., first_year + 4)

forvalues child = 1/20 {
    gen byte RELCHI`child'TYPE = .
    gen double RELCHI`child'BYEAR = .
}
replace RELCHI1TYPE = 1 if !missing(first_year)
replace RELCHI1BYEAR = first_year
replace RELCHI2TYPE = 1 if !missing(second_year)
replace RELCHI2BYEAR = second_year
replace RELCHIREP = 2 if !missing(second_year)
drop first_year second_year

sort ID year
save "`fixture'", replace

* The driver exits after the requested arm, so this fixture creation file is
* called once per arm by the surrounding smoke wrapper.
capture log close _all
do "`driver'" "`fixture'" "`outroot'" aligned_first_ownership all "`ado_root'"
