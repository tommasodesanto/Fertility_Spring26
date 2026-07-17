clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/hacamo_table8_replica_v1"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/hacamo_table8_replica_v1.log", replace text

* ==================================================================
* HACAMO TABLE 8 REPLICATION IN PSID
*
* Hacamo (2021) finds the credit-shock joint(buy AND newborn) effect
* is concentrated in LARGE home purchases (4+ rooms), not SMALL ones
* (<4 rooms). The mechanism is "access to space via ownership."
*
* Our PSID analog: among prime-age renters at 25-30, split joint(buy
* by 40 AND parent by 40) by the size of the home bought. Check
* whether DP-shortfall predicts joint_LARGE more than joint_SMALL.
* ==================================================================

* ------------------------------------------------------------------
* Build target P_target_t (median first-time-buyer home value)
* ------------------------------------------------------------------
use ID year AGEREP HOMEOWN HOMEVALUER HOMEMORTOTR YEARMOVED_ IW ///
    using "`dta'", clear
gen byte own = .
replace own = 1 if HOMEOWN == 1
replace own = 0 if HOMEOWN == 2
keep if inrange(year, 1984, 2019)
keep if !missing(AGEREP) & inrange(AGEREP, 18, 65)
sort ID year
by ID: gen any_renter_before = sum(own == 0)
by ID: gen byte first_purchase = own == 1 & any_renter_before[_n-1] >= 1 & ///
    sum(own == 1) == 1
keep if first_purchase == 1
keep if inrange(YEARMOVED_, year - 2, year)
keep if HOMEVALUER > 0 & !missing(HOMEVALUER)
gen year_bin = .
replace year_bin = 1 if inrange(year, 1984, 1994)
replace year_bin = 2 if inrange(year, 1995, 2004)
replace year_bin = 3 if inrange(year, 2005, 2010)
replace year_bin = 4 if inrange(year, 2011, 2019)
collapse (p50) p_target = HOMEVALUER [aweight=IW], by(year_bin)
sum p_target if !missing(p_target)
local fill = r(p50)
replace p_target = `fill' if missing(p_target) | p_target == 0
tempfile target
save `target', replace

* ------------------------------------------------------------------
* Get rooms at first owner spell for each ID
* (closest year of owner-state with ACTUALROOMS_ populated)
* ------------------------------------------------------------------
use ID year AGEREP HOMEOWN ACTUALROOMS_ using "`dta'", clear
drop if missing(AGEREP)
keep if inrange(year, 1984, 2019)
gen byte owner_year = HOMEOWN == 1
sort ID year
* First owner-year per ID
by ID: gen rank_own = sum(owner_year == 1) if owner_year == 1
keep if rank_own == 1
* Keep first such observation
keep if !missing(ACTUALROOMS_) & ACTUALROOMS_ > 0
collapse (first) first_owner_year = year first_owner_age = AGEREP ///
    rooms_at_first_own = ACTUALROOMS_, by(ID)
tempfile own_rooms
save `own_rooms', replace

* ------------------------------------------------------------------
* Built ownership-by-40 indicator (any owner spell at AGEREP <= 40)
* ------------------------------------------------------------------
use ID year AGEREP HOMEOWN using "`dta'", clear
drop if missing(AGEREP)
keep if inrange(year, 1984, 2019)
gen byte owner_year = HOMEOWN == 1
collapse (max) own_by_40 = owner_year if AGEREP <= 40, by(ID)
tempfile own40
save `own40', replace

* ------------------------------------------------------------------
* Build renter sample at 25-30 with wealth/income and DP shortfall
* ------------------------------------------------------------------
use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN IW using "`dta'", clear
drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)
gen byte renter = HOMEOWN == 2
keep if inrange(AGEREP, 25, 30)
keep if renter == 1
keep if !missing(NETWORTHR) & !missing(INCFAMR)
gen liquid_nw = cond(missing(NETWORTH2R), NETWORTHR, NETWORTH2R)
gen year_bin = .
replace year_bin = 1 if inrange(year, 1984, 1994)
replace year_bin = 2 if inrange(year, 1995, 2004)
replace year_bin = 3 if inrange(year, 2005, 2010)
replace year_bin = 4 if inrange(year, 2011, 2019)
merge m:1 year_bin using `target', keep(master match) nogen

gen byte sh_035 = liquid_nw < 0.035 * p_target if !missing(liquid_nw) & !missing(p_target)
gen byte sh_100 = liquid_nw < 0.10  * p_target if !missing(liquid_nw) & !missing(p_target)
gen byte sh_200 = liquid_nw < 0.20  * p_target if !missing(liquid_nw) & !missing(p_target)

collapse (mean) liquid_nw_mean = liquid_nw ///
         (mean) inc_mean = INCFAMR ///
         (mean) sh_035_m = sh_035 ///
         (mean) sh_100_m = sh_100 ///
         (mean) sh_200_m = sh_200 ///
         (mean) weight_young = IW, by(ID)
gen byte sh_035 = sh_035_m > 0.5 if !missing(sh_035_m)
gen byte sh_100 = sh_100_m > 0.5 if !missing(sh_100_m)
gen byte sh_200 = sh_200_m > 0.5 if !missing(sh_200_m)
tempfile entry
save `entry', replace

* ------------------------------------------------------------------
* Fertility outcomes
* ------------------------------------------------------------------
use ID year AGEREP RELCHI1BYEAR using "`dta'", clear
drop if missing(AGEREP)
collapse (min) min_year = year (min) min_age = AGEREP ///
    (first) first_birth_year_b = RELCHI1BYEAR, by(ID)
gen birth_year = min_year - min_age
gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
keep ID birth_year age_at_first_birth
tempfile fert
save `fert', replace

* Combine
use `entry', clear
merge 1:1 ID using `fert',       keep(master match) nogen
merge 1:1 ID using `own40',      keep(master match) nogen
merge 1:1 ID using `own_rooms',  keep(master match) nogen

gen byte parent_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
gen byte obs_by_40 = (2019 - birth_year) >= 40

* Define joint outcomes
* joint = owner by 40 AND parent by 40
* large home means rooms_at_first_own >= 4 (Hacamo's threshold)
gen byte joint_any   = own_by_40 == 1 & parent_by_40 == 1
gen byte joint_large = (own_by_40 == 1) & (rooms_at_first_own >= 4) & (parent_by_40 == 1)
gen byte joint_small = (own_by_40 == 1) & (rooms_at_first_own < 4) & (rooms_at_first_own >= 1) & (parent_by_40 == 1)

* Sample: observed by 40, with shortfall and income data
keep if obs_by_40 == 1
keep if !missing(sh_035) & !missing(inc_mean) & inc_mean > 0

di _n "=== Sample size ==="
count

di _n "=== Outcome means ==="
sum joint_any joint_large joint_small parent_by_40 own_by_40 [aweight=weight_young]

di _n "=== Rooms distribution conditional on having bought ==="
sum rooms_at_first_own [aweight=weight_young], detail

di _n "=== Joint outcome rates by shortfall status ==="
table sh_035 [aweight=weight_young], stat(mean joint_any joint_large joint_small) stat(frequency)

* ------------------------------------------------------------------
* Hacamo-style regressions: shortfall coefficient on joint_LARGE vs joint_SMALL
* ------------------------------------------------------------------
gen log_inc = log(inc_mean)
gen log_nw  = log(max(liquid_nw_mean, 1))

di _n "----- (1) joint_any (any home + baby) -----"
reg joint_any sh_035 log_inc log_nw [pweight=weight_young], robust

di _n "----- (2) joint_LARGE (4+ rooms + baby) -----"
reg joint_large sh_035 log_inc log_nw [pweight=weight_young], robust

di _n "----- (3) joint_SMALL (<4 rooms + baby) -----"
reg joint_small sh_035 log_inc log_nw [pweight=weight_young], robust

di _n "----- (4) joint_LARGE  with 20% shortfall instead of 3.5% -----"
reg joint_large sh_200 log_inc log_nw [pweight=weight_young], robust

di _n "----- (5) joint_SMALL  with 20% shortfall instead of 3.5% -----"
reg joint_small sh_200 log_inc log_nw [pweight=weight_young], robust

* Export summary table
preserve
    collapse (mean) joint_any joint_large joint_small (count) n_obs = ID ///
        [aweight=weight_young], by(sh_035)
    list
    export delimited using "`outdir'/joint_outcomes_by_shortfall_v1.csv", replace
restore

log close
