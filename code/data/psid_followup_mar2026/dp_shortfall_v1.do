clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/dp_shortfall_v1"
capture mkdir "`outdir'"
capture log close _all
log using "`outdir'/dp_shortfall_v1.log", replace text

* ==================================================================
* DP-SHORTFALL VS CHILDLESSNESS (V1)
*
* Object: do prime-age renters who lack the down payment for a
* "typical" first-time-buyer home become disproportionately
* childless? This is the empirical analog of the model's collateral
* constraint b >= (1 - phi) * p_i * H_k.
*
* Step 1: build target price by year = median home value among
* observed first-time PSID buyers (national).
* Step 2: build a panel of pre-birth renters at ages 25-30.
* Step 3: compare childlessness by 45 across (DP-shortfall yes/no)
* within income terciles, so the comparison nets out the income
* channel the user already pushes on.
* ==================================================================

* ------------------------------------------------------------------
* STEP 1: target price by year-bin from first-time buyers
* (reuse the first-purchase logic from dp_at_first_purchase_v1)
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

* Year-bin median value (real); broader bins = more stable benchmark
gen year_bin = .
replace year_bin = 1 if inrange(year, 1984, 1994)
replace year_bin = 2 if inrange(year, 1995, 2004)
replace year_bin = 3 if inrange(year, 2005, 2010)
replace year_bin = 4 if inrange(year, 2011, 2019)

collapse (p50) p_target = HOMEVALUER ///
         (count) n_buyers = ID [aweight=IW], by(year_bin)
list
export delimited using "`outdir'/target_price_by_year_bin_v1.csv", replace

* Fill any empty year-bins with median across non-missing bins
sum p_target if !missing(p_target)
local fill = r(p50)
replace p_target = `fill' if missing(p_target) | p_target == 0

tempfile target
save `target', replace

* ------------------------------------------------------------------
* STEP 2: build pre-birth renter panel at ages 25-30
* ------------------------------------------------------------------
use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN IW ///
    using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

* Identify renter-year observations (HOMEOWN==2)
gen byte renter = HOMEOWN == 2

* Restrict to entry window 25-30 and observed-as-renter at that age
keep if inrange(AGEREP, 25, 30)
keep if renter == 1
keep if !missing(NETWORTHR) & !missing(INCFAMR)

* Liquid wealth proxy: use NETWORTH2R (PSID's net-worth-excluding-home).
* For renters NETWORTH2R == NETWORTHR mechanically, but using the
* exclude-home version is conservative and consistent with treating
* home equity as illiquid.
gen liquid_nw = cond(missing(NETWORTH2R), NETWORTHR, NETWORTH2R)

* Year-bin for merge with target price
gen year_bin = .
replace year_bin = 1 if inrange(year, 1984, 1994)
replace year_bin = 2 if inrange(year, 1995, 2004)
replace year_bin = 3 if inrange(year, 2005, 2010)
replace year_bin = 4 if inrange(year, 2011, 2019)

merge m:1 year_bin using `target', keep(master match) nogen

* Required DP at 20% threshold (the GSE conforming floor; the modal
* binding number both in the bunching plot and in the calibrated phi)
gen req_dp = 0.20 * p_target
gen byte shortfall = liquid_nw < req_dp if !missing(liquid_nw) & !missing(req_dp)
gen byte shortfall_strict = liquid_nw < 0.10 * p_target if !missing(liquid_nw) & !missing(req_dp)

* Collapse to one obs per ID: mean shortfall, wealth, income across the
* 25-30 window (the entry-period summary the v12 do-file already uses)
preserve
    collapse (mean) liquid_nw_mean = liquid_nw ///
             (mean) inc_mean = INCFAMR ///
             (mean) p_target_mean = p_target ///
             (mean) req_dp_mean = req_dp ///
             (mean) shortfall_mean = shortfall ///
             (mean) shortfall_strict_mean = shortfall_strict ///
             (mean) weight_young = IW, by(ID)
    gen byte shortfall_any = shortfall_mean > 0.5
    gen byte shortfall_strict_any = shortfall_strict_mean > 0.5
    tempfile entry
    save `entry', replace
restore

* ------------------------------------------------------------------
* STEP 3: bring in fertility outcomes and observation horizon
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

use `entry', clear
merge 1:1 ID using `fert', keep(match) nogen

gen byte observed_by_45 = (2019 - birth_year) >= 45
gen byte observed_by_40 = (2019 - birth_year) >= 40
gen byte observed_by_35 = (2019 - birth_year) >= 35

gen byte parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
gen byte parent_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
gen byte parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)
gen byte childless_by_45 = 1 - parent_by_45
gen byte childless_by_40 = 1 - parent_by_40
gen byte childless_by_35 = 1 - parent_by_35

di _n "=== Pre-birth renter sample at 25-30 with non-missing wealth/income ==="
count
di "Observed by 45:"
count if observed_by_45 == 1
di "Observed by 40:"
count if observed_by_40 == 1

* ------------------------------------------------------------------
* HEADLINE: childlessness by DP-shortfall status
* (unconditional and within income terciles)
* ------------------------------------------------------------------
preserve
    keep if observed_by_45 == 1
    keep if !missing(shortfall_any) & !missing(inc_mean)

    di _n "Sample size for headline:"
    count

    di _n "=== Share of pre-birth renters with DP shortfall ==="
    sum shortfall_any shortfall_strict_any [aweight=weight_young]

    di _n "=== Childless by 45 by shortfall status (unconditional) ==="
    table shortfall_any [aweight=weight_young], stat(mean childless_by_45 childless_by_40) ///
        stat(frequency)

    * Within income terciles
    xtile inc_tercile = inc_mean [pweight=weight_young], n(3)
    label define inct 1 "Bottom inc tercile" 2 "Middle" 3 "Top"
    label values inc_tercile inct

    di _n "=== Childless by 45 by shortfall x income tercile ==="
    table inc_tercile shortfall_any [aweight=weight_young], ///
        stat(mean childless_by_45) stat(frequency)

    collapse (mean) childless_by_45 childless_by_40 childless_by_35 ///
             (count) n_obs = ID [aweight=weight_young], by(inc_tercile shortfall_any)
    list, sepby(inc_tercile)
    export delimited using "`outdir'/childless_by_shortfall_inc_tercile_v1.csv", replace
restore

* ------------------------------------------------------------------
* Regression: childless_by_45 on shortfall, controlling for log income
* ------------------------------------------------------------------
preserve
    keep if observed_by_45 == 1
    keep if !missing(shortfall_any) & !missing(inc_mean) & inc_mean > 0

    gen log_inc = log(inc_mean)
    gen log_nw  = log(max(liquid_nw_mean, 1))

    di _n "--- Childless by 45 ~ shortfall ---"
    reg childless_by_45 shortfall_any [pweight=weight_young], robust

    di _n "--- Childless by 45 ~ shortfall + log income ---"
    reg childless_by_45 shortfall_any log_inc [pweight=weight_young], robust

    di _n "--- Childless by 45 ~ shortfall + log income + log NW ---"
    reg childless_by_45 shortfall_any log_inc log_nw [pweight=weight_young], robust

    di _n "--- Childless by 45 ~ STRICT shortfall (10% threshold) + log inc ---"
    reg childless_by_45 shortfall_strict_any log_inc [pweight=weight_young], robust
restore

* ------------------------------------------------------------------
* FIGURE: bar of childlessness by shortfall x income tercile
* ------------------------------------------------------------------
preserve
    keep if observed_by_45 == 1
    keep if !missing(shortfall_any) & !missing(inc_mean)

    xtile inc_tercile = inc_mean [pweight=weight_young], n(3)
    label define inct 1 "Bottom inc tercile" 2 "Middle" 3 "Top"
    label values inc_tercile inct

    collapse (mean) childless_by_45 (count) n_obs = ID ///
        [aweight=weight_young], by(inc_tercile shortfall_any)
    list
    export delimited using "`outdir'/childless_by_shortfall_inc_tercile_for_plot_v1.csv", replace

    graph bar childless_by_45, ///
        over(shortfall_any, relabel(1 "DP sufficient" 2 "DP shortfall")) ///
        over(inc_tercile) ///
        bar(1, fcolor(navy%70) lcolor(navy)) ///
        ytitle("Share childless by 45", size(medium)) ///
        ylabel(0(0.1)0.7, angle(h)) ///
        blabel(bar, format(%4.2f) size(small)) ///
        title("Renters who lack the 20% down payment remain childless more often") ///
        subtitle("PSID renters at 25-30, pre-first-birth, by DP shortfall and income tercile") ///
        note("DP shortfall = liquid wealth < 20% of median first-time-buyer home value in same year-bin.") ///
        graphregion(color(white)) scheme(s1color) legend(off)
    graph export "`outdir'/childless_by_shortfall_inc_tercile_v1.png", replace width(1800)
restore

log close
