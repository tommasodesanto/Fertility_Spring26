clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/dp_joint_margin_v1"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/dp_joint_margin_v1.log", replace text

* ==================================================================
* JOINT MARGIN (HACAMO-STYLE OUTCOME) IN PSID
*
* Hacamo (2021) RFS estimates the joint probability:
*    P(buy a home AND have a first child)
* among young households exposed to mortgage-market deregulation.
* This is the right object for the user's model: the chain
*    credit -> DP -> ownership -> family-sized unit -> child
* happens inside one household. P(joint) is small if buying and
* birth are independent population responses; P(joint) is large
* if they are linked within-household.
*
* This do-file builds the cross-sectional analog in PSID. Sample:
* renters at 25-30 pre-first-birth, observed to 35 / 40. Outcomes:
* buy_by_X, birth_by_X, joint_by_X. Split by DP-shortfall (3.5% FHA
* threshold) x income tercile.
* ==================================================================

* ------------------------------------------------------------------
* STEP 1: target price by year-bin (re-uses dp_shortfall_v1 logic)
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
* STEP 2: build full ID-year panel with ownership and age info,
* to compute "buy by age X" outcomes
* ------------------------------------------------------------------
use ID year AGEREP HOMEOWN using "`dta'", clear
drop if missing(AGEREP)
keep if inrange(year, 1984, 2019)
gen byte owner_year = HOMEOWN == 1
* Birth year of the person
collapse (min) min_year = year (min) min_age = AGEREP ///
         (max) any_owner_by_35_raw = owner_year if AGEREP <= 35, by(ID)
gen birth_year = min_year - min_age
drop min_year min_age
rename any_owner_by_35_raw owner_by_35
tempfile ages
save `ages', replace

* Buy-by-40 and buy-by-45 separately
use ID year AGEREP HOMEOWN using "`dta'", clear
drop if missing(AGEREP)
keep if inrange(year, 1984, 2019)
gen byte owner_year = HOMEOWN == 1
collapse (max) owner_by_40 = owner_year if AGEREP <= 40, by(ID)
tempfile own40
save `own40', replace

use ID year AGEREP HOMEOWN using "`dta'", clear
drop if missing(AGEREP)
keep if inrange(year, 1984, 2019)
gen byte owner_year = HOMEOWN == 1
collapse (max) owner_by_45 = owner_year if AGEREP <= 45, by(ID)
tempfile own45
save `own45', replace

* ------------------------------------------------------------------
* STEP 3: build renter-entry sample at ages 25-30 with wealth/income
* ------------------------------------------------------------------
use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN IW ///
    using "`dta'", clear
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
gen byte sh_200 = liquid_nw < 0.20  * p_target if !missing(liquid_nw) & !missing(p_target)
collapse (mean) liquid_nw_mean = liquid_nw ///
         (mean) inc_mean = INCFAMR ///
         (mean) p_target_mean = p_target ///
         (mean) sh_035_m = sh_035 ///
         (mean) sh_200_m = sh_200 ///
         (mean) weight_young = IW, by(ID)
gen byte sh_035 = sh_035_m > 0.5
gen byte sh_200 = sh_200_m > 0.5
replace sh_035 = . if missing(sh_035_m)
replace sh_200 = . if missing(sh_200_m)
tempfile entry
save `entry', replace

* ------------------------------------------------------------------
* STEP 4: fertility outcomes + merge
* ------------------------------------------------------------------
use ID year AGEREP RELCHI1BYEAR using "`dta'", clear
drop if missing(AGEREP)
collapse (first) first_birth_year_b = RELCHI1BYEAR, by(ID)
tempfile fert
save `fert', replace

use `entry', clear
merge 1:1 ID using `fert',  keep(match master) nogen
merge 1:1 ID using `ages',  keep(match master) nogen
merge 1:1 ID using `own40', keep(match master) nogen
merge 1:1 ID using `own45', keep(match master) nogen

* Birth-by-X indicators
gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
gen byte birth_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)
gen byte birth_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
gen byte birth_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)

* Observation horizon: must be old enough by 2019
gen byte obs_by_35 = (2019 - birth_year) >= 35
gen byte obs_by_40 = (2019 - birth_year) >= 40
gen byte obs_by_45 = (2019 - birth_year) >= 45

* Joint margin
gen byte joint_35 = owner_by_35 == 1 & birth_by_35 == 1 if !missing(owner_by_35)
gen byte joint_40 = owner_by_40 == 1 & birth_by_40 == 1 if !missing(owner_by_40)
gen byte joint_45 = owner_by_45 == 1 & birth_by_45 == 1 if !missing(owner_by_45)

* ------------------------------------------------------------------
* HEADLINE: marginal vs joint by shortfall (3.5% FHA threshold)
* ------------------------------------------------------------------
preserve
    keep if obs_by_40 == 1
    keep if !missing(sh_035) & !missing(inc_mean)
    di _n "=== Sample size for joint-by-40 analysis ==="
    count

    xtile inc_t = inc_mean [pweight=weight_young], n(3)
    label define inct 1 "Bottom inc" 2 "Middle" 3 "Top"
    label values inc_t inct

    di _n "=== Joint, marginal own, marginal birth by 40 -- by shortfall x income tercile ==="
    table inc_t sh_035 [aweight=weight_young], ///
        stat(mean owner_by_40 birth_by_40 joint_40) stat(frequency)

    di _n "=== UNCONDITIONAL by shortfall (across all income terciles) ==="
    table sh_035 [aweight=weight_young], ///
        stat(mean owner_by_40 birth_by_40 joint_40) stat(frequency)

    collapse (mean) own_40 = owner_by_40 birth_40 = birth_by_40 joint_40 ///
             (count) n_obs = ID [aweight=weight_young], by(inc_t sh_035)
    list, sepby(inc_t)
    export delimited using "`outdir'/joint_margin_by40_by_shortfall_v1.csv", replace
restore

preserve
    keep if obs_by_45 == 1
    keep if !missing(sh_035) & !missing(inc_mean)
    di _n "=== Sample size for joint-by-45 analysis ==="
    count
    xtile inc_t = inc_mean [pweight=weight_young], n(3)
    label values inc_t inct
    collapse (mean) own_45 = owner_by_45 birth_45 = birth_by_45 joint_45 ///
             (count) n_obs = ID [aweight=weight_young], by(inc_t sh_035)
    list, sepby(inc_t)
    export delimited using "`outdir'/joint_margin_by45_by_shortfall_v1.csv", replace
restore

* ------------------------------------------------------------------
* REGRESSIONS: marginal outcomes vs joint outcome
* The model predicts the gap on JOINT is closer to the gap on EACH
* marginal than to their product (i.e., the within-household chain
* dominates the independent-effects story).
* ------------------------------------------------------------------
preserve
    keep if obs_by_40 == 1
    keep if !missing(sh_035) & !missing(inc_mean) & inc_mean > 0
    gen log_inc = log(inc_mean)
    gen log_nw  = log(max(liquid_nw_mean, 1))

    di _n "--- P(owner by 40) ~ sh_035 + log inc + log NW ---"
    reg owner_by_40 sh_035 log_inc log_nw [pweight=weight_young], robust

    di _n "--- P(birth by 40) ~ sh_035 + log inc + log NW ---"
    reg birth_by_40 sh_035 log_inc log_nw [pweight=weight_young], robust

    di _n "--- P(joint by 40) ~ sh_035 + log inc + log NW ---"
    reg joint_40    sh_035 log_inc log_nw [pweight=weight_young], robust
restore

preserve
    keep if obs_by_45 == 1
    keep if !missing(sh_035) & !missing(inc_mean) & inc_mean > 0
    gen log_inc = log(inc_mean)
    gen log_nw  = log(max(liquid_nw_mean, 1))

    di _n "--- P(owner by 45) ~ sh_035 + log inc + log NW ---"
    reg owner_by_45 sh_035 log_inc log_nw [pweight=weight_young], robust

    di _n "--- P(birth by 45) ~ sh_035 + log inc + log NW ---"
    reg birth_by_45 sh_035 log_inc log_nw [pweight=weight_young], robust

    di _n "--- P(joint by 45) ~ sh_035 + log inc + log NW ---"
    reg joint_45    sh_035 log_inc log_nw [pweight=weight_young], robust
restore

* ------------------------------------------------------------------
* FIGURE: bar of joint by 40 by shortfall x income tercile
* ------------------------------------------------------------------
preserve
    keep if obs_by_40 == 1
    keep if !missing(sh_035) & !missing(inc_mean)
    xtile inc_t = inc_mean [pweight=weight_young], n(3)
    label values inc_t inct

    collapse (mean) joint_40 (count) n_obs = ID [aweight=weight_young], ///
        by(inc_t sh_035)
    graph bar joint_40, ///
        over(sh_035, relabel(1 "DP suff." 2 "DP short")) ///
        over(inc_t) ///
        bar(1, fcolor(navy%70) lcolor(navy)) ///
        ytitle("P(own by 40 AND first birth by 40)", size(medium)) ///
        ylabel(0(0.1)0.7, angle(h)) ///
        blabel(bar, format(%4.2f) size(small)) ///
        title("Joint margin: ownership and first birth by 40", size(medlarge)) ///
        subtitle("PSID renters at 25-30, by FHA-floor shortfall and income tercile", size(small)) ///
        note("DP shortfall = liquid wealth < 3.5% of median first-time-buyer home value in same year-bin.", size(small)) ///
        graphregion(color(white)) scheme(s1color) legend(off)
    graph export "`outdir'/joint_by40_shortfall_inc_tercile_v1.png", replace width(1800)
restore

log close
