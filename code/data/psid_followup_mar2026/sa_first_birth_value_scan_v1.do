* Variable scan for the first-birth house-value event study (Oct 18, 2026).
* Lists geography-like shelf variables with years covered, and checks the
* timing of the shelf value/mortgage columns against official family items.
clear all
set more off
version 17.0
local shelf "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
describe using "`shelf'", varlist
local all `r(varlist)'
local geo
foreach v of local all {
    local u = upper("`v'")
    if strpos("`u'","GEO") | strpos("`u'","STATE") | strpos("`u'","COUNTY") | strpos("`u'","CNTY") | strpos("`u'","FIPS") | strpos("`u'","METRO") | strpos("`u'","MSA") | strpos("`u'","SMSA") | strpos("`u'","URBAN") | strpos("`u'","RURAL") | strpos("`u'","BEALE") | strpos("`u'","RUCC") | strpos("`u'","RUC") | strpos("`u'","CBSA") | strpos("`u'","REGION") | strpos("`u'","CITY") | strpos("`u'","PLACE") | strpos("`u'","TRACT") | strpos("`u'","ZIP") | strpos("`u'","PUMA") | strpos("`u'","SIZE") {
        local geo `geo' `v'
    }
}
di "GEO_CANDIDATES: `geo'"
use ID year `geo' HOMEVALUER HOMEMORTOTR HOMEOWN ACTUALROOMS_ CURRENT using "`shelf'", clear
describe `geo' HOMEVALUER HOMEMORTOTR HOMEOWN
foreach v in `geo' HOMEVALUER HOMEMORTOTR {
    quietly summarize year if !missing(`v')
    di "COVER `v' first=" r(min) " last=" r(max) " n=" r(N)
    quietly levelsof year if !missing(`v'), local(yrs)
    di "YEARS `v': `yrs'"
}
foreach v of local geo {
    capture noisily tab `v' if year==2001, missing
}
* Timing check against official house value items (1999, 2001, 2019).
local cms "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta"
describe ER13041 ER13047 ER17044 ER17052 ER72031 V5 V10438 using "`cms'"
sort ID year
tempfile shelf_small
save `shelf_small'
use ID ER13041 ER17044 using "`cms'", clear
rename ER13041 val1999
rename ER17044 val2001
merge 1:m ID using `shelf_small', keep(match) nogenerate
foreach y in 1999 2001 {
    di "YEAR `y': shelf HOMEVALUER at same year vs official"
    count if year==`y' & CURRENT==1 & !missing(HOMEVALUER) & HOMEVALUER>0 & val`y'>0 & val`y'<9999990
    corr HOMEVALUER val`y' if year==`y' & CURRENT==1 & HOMEVALUER>0 & val`y'>0 & val`y'<9999990
    corr HOMEVALUER val`y' if year==`y'-2 & CURRENT==1 & HOMEVALUER>0 & val`y'>0 & val`y'<9999990
    gen double ratio`y' = HOMEVALUER/val`y' if year==`y' & CURRENT==1 & HOMEVALUER>0 & val`y'>0 & val`y'<9999990
    summarize ratio`y', detail
}
