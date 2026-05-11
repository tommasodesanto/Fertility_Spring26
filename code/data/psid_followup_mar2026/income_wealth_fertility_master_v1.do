clear all
set more off
version 17.0

* ==================================================================
* Income / wealth / fertility master diagnostic
*
* Purpose:
*   Build one audited long-format table for fertility outcomes by
*   young-adult income and wealth bins. This consolidates the exploratory
*   fertility_wealth_v1-v12 outputs into a single reproducible object.
*
* Main output:
*   output/income_wealth_fertility_master_v1/income_wealth_fertility_master_v1.csv
*
* Notes:
*   - Wealth is measured at ages 25-30.
*   - Fertility outcomes use first-birth timing and children-count measures.
*   - ACS NCHILD-style current-household measures are not used here.
* ==================================================================

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local out_root "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output"
local outdir "`out_root'/income_wealth_fertility_master_v1"

cap mkdir "`out_root'"
cap mkdir "`outdir'"

capture log close _all
log using "`outdir'/income_wealth_fertility_master_v1.log", replace text

global weight IW

di as text "Loading PSID for income/wealth/fertility master diagnostic..."
use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR EARNINDR NETWORTHR NETWORTH2R NETWORTH3R HOMEEQUITYR ///
    HOMEOWN ${weight} using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

gen own = HOMEOWN
replace own = 0 if own == 2
replace own = . if own == 3

gen first_birth_year = RELCHI1BYEAR
gen n_children = RELCHIREP

* ------------------------------------------------------------------
* Entry income, wealth, and tenure at ages 25-30
* ------------------------------------------------------------------

preserve
    keep if inrange(AGEREP, 25, 30)
    keep if !missing(${weight})
    collapse ///
        (mean) sex_young = SEX ///
        (mean) total_nw = NETWORTHR ///
        (mean) liquid_nw = NETWORTH2R ///
        (mean) broad_nw = NETWORTH3R ///
        (mean) home_equity = HOMEEQUITYR ///
        (mean) family_income = INCFAMR ///
        (mean) individual_earnings = EARNINDR ///
        (mean) own_rate_2530 = own ///
        (mean) weight_young = ${weight}, by(ID)

    gen renter_2530 = (own_rate_2530 < 0.5) if !missing(own_rate_2530)
    gen owner_2530 = (own_rate_2530 >= 0.5) if !missing(own_rate_2530)

    gen total_nw_k = total_nw / 1000
    gen liquid_nw_k = liquid_nw / 1000
    gen broad_nw_k = broad_nw / 1000
    gen home_equity_k = home_equity / 1000
    gen income_k = family_income / 1000
    gen earnings_k = individual_earnings / 1000

    gen total_nw_to_inc = total_nw / family_income if family_income > 1000
    gen liquid_nw_to_inc = liquid_nw / family_income if family_income > 1000
    gen broad_nw_to_inc = broad_nw / family_income if family_income > 1000
    gen home_equity_to_inc = home_equity / family_income if family_income > 1000

    tempfile entry
    save `entry', replace
restore

* ------------------------------------------------------------------
* Ownership transition by 35
* ------------------------------------------------------------------

preserve
    keep if inrange(AGEREP, 31, 35)
    keep if !missing(own)
    collapse (mean) own_rate_3135 = own, by(ID)
    gen owner_by_35 = (own_rate_3135 >= 0.5) if !missing(own_rate_3135)
    keep ID owner_by_35 own_rate_3135
    tempfile tenure35
    save `tenure35', replace
restore

* ------------------------------------------------------------------
* Fertility timing from first-birth year
* ------------------------------------------------------------------

preserve
    collapse (min) min_year = year (min) min_age = AGEREP ///
        (first) first_birth_year_b = first_birth_year, by(ID)
    gen birth_year = min_year - min_age
    gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)

    foreach a in 30 35 40 45 {
        gen observed_by_`a' = (2019 - birth_year) >= `a'
        gen parent_by_`a' = (!missing(age_at_first_birth) & age_at_first_birth <= `a')
        gen childless_by_`a' = 1 - parent_by_`a'
    }

    keep ID birth_year age_at_first_birth observed_by_* parent_by_* childless_by_*
    tempfile fert_timing
    save `fert_timing', replace
restore

* ------------------------------------------------------------------
* Children-count outcomes near completed fertility
* ------------------------------------------------------------------

preserve
    keep if inrange(AGEREP, 35, 40)
    keep if !missing(n_children)
    collapse (max) children_35_40 = n_children, by(ID)
    tempfile children_3540
    save `children_3540', replace
restore

preserve
    keep if inrange(AGEREP, 43, 50)
    keep if !missing(n_children)
    collapse (max) children_43_50 = n_children, by(ID)
    tempfile children_4350
    save `children_4350', replace
restore

use `entry', clear
merge 1:1 ID using `tenure35', keep(match master) nogen
merge 1:1 ID using `fert_timing', keep(match) nogen
merge 1:1 ID using `children_3540', keep(match master) nogen
merge 1:1 ID using `children_4350', keep(match master) nogen

keep if !missing(weight_young) & weight_young > 0

save "`outdir'/income_wealth_fertility_master_sample_v1.dta", replace

* ------------------------------------------------------------------
* Long-format bin summaries
* ------------------------------------------------------------------

local py "`root'/Fertility_Spring26/code/model/.venv/bin/python"
local builder "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/build_income_wealth_fertility_master_v1.py"

di as text "Building long-format CSV from audited sample..."
shell "`py'" "`builder'"
capture confirm file "`outdir'/income_wealth_fertility_master_v1.csv"
if _rc {
    di as error "Python builder did not produce the income/wealth/fertility CSV."
    exit 601
}

di as result "Income/wealth/fertility master diagnostic complete."
di as result "Sample: `outdir'/income_wealth_fertility_master_sample_v1.dta"
di as result "CSV: `outdir'/income_wealth_fertility_master_v1.csv"

log close _all
