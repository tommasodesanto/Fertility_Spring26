clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local hpi  "`root'/Fertility_Spring26/code/data/fhfa_hpi/hpi_at_state.csv"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/dp_hpi_shock_v1"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/dp_hpi_shock_v1.log", replace text

* ==================================================================
* HPI SHOCK x WEALTH IN PSID
*
* Price-side analog to Hacamo's credit-side shock. The model
* constraint b >= 0.20 p_i H_k binds more tightly when p_i rises
* for the SAME renter. For each renter at 25-30, compute their
* state's HPI growth over the [t, t+5] window and interact with
* their pre-shock wealth-to-threshold ratio. The model predicts:
* HPI rises hurt the joint margin (own AND birth) MORE for
* renters closer to the threshold.
* ==================================================================

* ------------------------------------------------------------------
* STEP 0: prep FHFA state HPI
* ------------------------------------------------------------------
import delimited using "`hpi'", clear varnames(1)
keep state year hpi
rename state state_name
rename year hpi_year
rename hpi hpi_st
destring hpi_year hpi_st, replace force
keep if !missing(state_name) & !missing(hpi_year) & !missing(hpi_st)
* Lowercase normalized state name for safe merging
gen str40 state_name_lc = lower(strtrim(state_name))
drop state_name
rename hpi_year year
rename hpi_st   hpi
tempfile hpi_long
save `hpi_long', replace
di _n "FHFA state HPI panel:"
sum year hpi

* ------------------------------------------------------------------
* STEP 1: target price by year-bin (reuse)
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
* STEP 2: build full owner/age panel for outcomes
* ------------------------------------------------------------------
use ID year AGEREP HOMEOWN using "`dta'", clear
drop if missing(AGEREP)
keep if inrange(year, 1984, 2019)
gen byte owner_year = HOMEOWN == 1
collapse (min) min_year = year (min) min_age = AGEREP ///
         (max) own_by_35 = owner_year if AGEREP <= 35, by(ID)
gen birth_year = min_year - min_age
drop min_year min_age
tempfile ages35
save `ages35', replace

use ID year AGEREP HOMEOWN using "`dta'", clear
drop if missing(AGEREP)
keep if inrange(year, 1984, 2019)
gen byte owner_year = HOMEOWN == 1
collapse (max) own_by_40 = owner_year if AGEREP <= 40, by(ID)
tempfile own40
save `own40', replace

* ------------------------------------------------------------------
* STEP 3: renter sample at 25-30 with wealth/income AND state info
* ------------------------------------------------------------------
use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN IW GEOSTATE ///
    using "`dta'", clear
drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)
gen byte renter = HOMEOWN == 2
keep if inrange(AGEREP, 25, 30)
keep if renter == 1
keep if !missing(NETWORTHR) & !missing(INCFAMR) & !missing(GEOSTATE)
gen liquid_nw = cond(missing(NETWORTH2R), NETWORTHR, NETWORTH2R)
gen year_bin = .
replace year_bin = 1 if inrange(year, 1984, 1994)
replace year_bin = 2 if inrange(year, 1995, 2004)
replace year_bin = 3 if inrange(year, 2005, 2010)
replace year_bin = 4 if inrange(year, 2011, 2019)
merge m:1 year_bin using `target', keep(master match) nogen

* Decode GEOSTATE to string name and normalize for merging to FHFA
decode GEOSTATE, gen(state_name_raw)
gen str40 state_name_lc = lower(strtrim(state_name_raw))
drop state_name_raw

* Per-obs HPI in observation year AND 5 years later
rename year obs_year
preserve
    keep ID obs_year state_name_lc
    rename obs_year year
    merge m:1 state_name_lc year using `hpi_long', keep(master match) nogen
    rename hpi hpi_t
    rename year obs_year
    keep ID obs_year hpi_t
    tempfile hpi_t
    save `hpi_t', replace
restore
preserve
    gen year = obs_year + 5
    keep ID obs_year state_name_lc year
    merge m:1 state_name_lc year using `hpi_long', keep(master match) nogen
    rename hpi hpi_tp5
    drop year
    keep ID obs_year hpi_tp5
    tempfile hpi_tp5
    save `hpi_tp5', replace
restore

merge 1:1 ID obs_year using `hpi_t',   keep(master match) nogen
merge 1:1 ID obs_year using `hpi_tp5', keep(master match) nogen

gen dhpi_5 = log(hpi_tp5 / hpi_t) if !missing(hpi_t) & !missing(hpi_tp5) & hpi_t > 0

* Wealth-to-threshold ratio at the 3.5% line
gen rel_wealth_035 = liquid_nw / (0.035 * p_target) if !missing(liquid_nw) & !missing(p_target) & p_target > 0
gen byte sh_035 = liquid_nw < 0.035 * p_target if !missing(liquid_nw) & !missing(p_target)
gen byte sh_200 = liquid_nw < 0.20  * p_target if !missing(liquid_nw) & !missing(p_target)

* Take state in the first observed year of the 25-30 window
sort ID obs_year
by ID: gen byte first_in_id = _n == 1
preserve
    keep if first_in_id == 1
    keep ID state_name_lc
    rename state_name_lc state_name_first
    tempfile state_first
    save `state_first', replace
restore
* Collapse to one obs per ID
collapse (mean) liquid_nw_mean = liquid_nw ///
         (mean) inc_mean = INCFAMR ///
         (mean) p_target_mean = p_target ///
         (mean) sh_035_m = sh_035 ///
         (mean) sh_200_m = sh_200 ///
         (mean) dhpi_5_mean = dhpi_5 ///
         (mean) rel_wealth_035_mean = rel_wealth_035 ///
         (mean) weight_young = IW, by(ID)
merge 1:1 ID using `state_first', keep(master match) nogen
rename state_name_first state_name_lc
* Guard binary shortfall vars against Stata's missing-as-+inf rule
gen byte sh_035 = sh_035_m > 0.5 if !missing(sh_035_m)
gen byte sh_200 = sh_200_m > 0.5 if !missing(sh_200_m)
tempfile entry
save `entry', replace

* ------------------------------------------------------------------
* STEP 4: outcomes
* ------------------------------------------------------------------
use ID year AGEREP RELCHI1BYEAR using "`dta'", clear
drop if missing(AGEREP)
collapse (first) first_birth_year_b = RELCHI1BYEAR, by(ID)
tempfile fert
save `fert', replace

use `entry', clear
merge 1:1 ID using `fert',   keep(master match) nogen
merge 1:1 ID using `ages35', keep(master match) nogen
merge 1:1 ID using `own40',  keep(master match) nogen

gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
gen byte birth_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)
gen byte birth_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
gen byte obs_by_35 = (2019 - birth_year) >= 35
gen byte obs_by_40 = (2019 - birth_year) >= 40

gen byte joint_35 = own_by_35 == 1 & birth_by_35 == 1 if !missing(own_by_35)
gen byte joint_40 = own_by_40 == 1 & birth_by_40 == 1 if !missing(own_by_40)

* ------------------------------------------------------------------
* HEADLINE: HPI shock x wealth-to-threshold interaction on joint
* ------------------------------------------------------------------
preserve
    keep if obs_by_40 == 1
    keep if !missing(sh_035) & !missing(dhpi_5_mean) & !missing(inc_mean) & inc_mean > 0

    di _n "=== Sample size for HPI-shock joint-40 analysis ==="
    count

    sum dhpi_5_mean, detail
    sum rel_wealth_035_mean, detail

    gen log_inc = log(inc_mean)
    gen log_nw  = log(max(liquid_nw_mean, 1))

    di _n "--- P(joint by 40) ~ dHPI(5) + shortfall + interaction ---"
    reg joint_40 c.dhpi_5_mean##i.sh_035 log_inc log_nw [pweight=weight_young], robust

    di _n "--- P(own by 40) ~ dHPI(5) + shortfall + interaction ---"
    reg own_by_40 c.dhpi_5_mean##i.sh_035 log_inc log_nw [pweight=weight_young], robust

    di _n "--- P(birth by 40) ~ dHPI(5) + shortfall + interaction ---"
    reg birth_by_40 c.dhpi_5_mean##i.sh_035 log_inc log_nw [pweight=weight_young], robust

    di _n "--- Alternate: use rel_wealth_035 continuously ---"
    di "Higher rel_wealth = more surplus (lower DP-constraint pressure)"
    reg joint_40 c.dhpi_5_mean##c.rel_wealth_035_mean log_inc log_nw ///
        [pweight=weight_young], robust
restore

* ------------------------------------------------------------------
* BIN PLOT: joint by 40 as a function of dHPI(5) x shortfall
* ------------------------------------------------------------------
preserve
    keep if obs_by_40 == 1
    keep if !missing(sh_035) & !missing(dhpi_5_mean)
    xtile dhpi_q = dhpi_5_mean [pweight=weight_young], n(4)
    label define dq 1 "Q1 lowest dHPI" 2 "Q2" 3 "Q3" 4 "Q4 highest dHPI"
    label values dhpi_q dq

    di _n "=== joint_40 by dHPI quartile x shortfall ==="
    table dhpi_q sh_035 [aweight=weight_young], stat(mean joint_40) stat(frequency)

    collapse (mean) joint_40 (count) n_obs = ID [aweight=weight_young], ///
        by(dhpi_q sh_035)
    list, sepby(dhpi_q)
    export delimited using "`outdir'/joint_by40_dhpi_shortfall_v1.csv", replace

    graph bar joint_40, ///
        over(sh_035, relabel(1 "DP suff." 2 "DP short")) ///
        over(dhpi_q) ///
        bar(1, fcolor(navy%70) lcolor(navy)) ///
        ytitle("P(own by 40 AND first birth by 40)", size(medium)) ///
        ylabel(0(0.1)0.7, angle(h)) ///
        blabel(bar, format(%4.2f) size(small)) ///
        title("Joint margin by HPI shock and DP shortfall", size(medlarge)) ///
        subtitle("PSID renters at 25-30, by 5-yr state HPI growth quartile and FHA-shortfall", size(small)) ///
        note("FHFA state HPI; quartiles of log(HPI{sub:t+5}/HPI{sub:t}). Shortfall: liquid wealth < 3.5% of P_target.", size(small)) ///
        graphregion(color(white)) scheme(s1color) legend(off)
    graph export "`outdir'/joint_by40_dhpi_shortfall_v1.png", replace width(1800)
restore

log close
