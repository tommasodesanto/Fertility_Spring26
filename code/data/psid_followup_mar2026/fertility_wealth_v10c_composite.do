clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v10c.log", replace text

* ==================================================================
* V10c: COMPOSITE FIGURES
*
* The dose-response doesn't work for fertility directly because
* wealth has zero direct effect on fertility — the entire channel
* runs through ownership. So the right figure shows:
*   (a) Wealth determines who becomes an owner (strong gradient)
*   (b) Conditional on wealth, owners are much less childless (large gap)
*
* This IS the margin: wealth relaxes the ownership constraint,
* and the ownership constraint binds on fertility at every wealth level.
* ==================================================================

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen own = HOMEOWN

* --- Entry at 25-30 ---
preserve
    keep if inrange(AGEREP, 25, 30) & !missing(own) & !missing(NETWORTHR)
    collapse (mean) own_rate_2530 = own ///
        (mean) nw_dollars = NETWORTHR ///
        (mean) inc_dollars = INCFAMR ///
        (mean) weight_young = IW, by(ID)
    gen renter_2530 = (own_rate_2530 < 0.5)
    gen nw_k = nw_dollars / 1000
    tempfile entry
    save `entry', replace
restore

* --- Tenure at 31-35 ---
preserve
    keep if inrange(AGEREP, 31, 35) & !missing(own)
    collapse (mean) own_rate_3135 = own, by(ID)
    gen owner_by_35 = (own_rate_3135 >= 0.5) if !missing(own_rate_3135)
    keep ID owner_by_35
    tempfile t3135
    save `t3135', replace
restore

* --- Fertility ---
preserve
    collapse (min) min_year = year (min) min_age = AGEREP ///
        (first) first_birth_year_b = first_birth_year, by(ID)
    gen birth_year = min_year - min_age
    gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
    keep ID birth_year age_at_first_birth
    tempfile fert
    save `fert', replace
restore

use `entry', clear
merge 1:1 ID using `t3135', keep(match master) nogen
merge 1:1 ID using `fert', keep(match) nogen

gen observed_by_45 = (2019 - birth_year) >= 45
gen parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
gen childless_by_45 = 1 - parent_by_45

* Restrict to renters at 25-30, observed by 45
keep if renter_2530 == 1 & observed_by_45 == 1 & !missing(owner_by_35)

_pctile nw_dollars [pw=weight_young], p(1 99)
drop if nw_dollars < r(r1) | nw_dollars > r(r2)

di _n "=== Final sample ==="
count

* Dollar bins
gen nw_bin = .
replace nw_bin = 1 if nw_dollars < 0
replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

tab nw_bin

* ==================================================================
* 1. WITHIN-WEALTH CHILDLESSNESS GAP BY DOLLAR BINS
*    (Dollar-bin version of the v9 within_wealth_gap)
* ==================================================================
di _n "============================================="
di    "1. Within-wealth childlessness gap (dollar bins)"
di    "============================================="

preserve
    collapse (mean) childless_by_45 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin owner_by_35)
    gen transition = cond(owner_by_35 == 1, "Became owner", "Stayed renter")
    list, sep(2)
    export delimited using "`outdir'/within_wealth_gap_dollars_v10c.csv", replace

    twoway (connected childless_by_45 nw_bin if owner_by_35==1, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(vlarge) lwidth(medthick)) ///
           (connected childless_by_45 nw_bin if owner_by_35==0, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(vlarge) lwidth(medthick)) ///
        , legend(order(1 "Became owner by 35" 2 "Stayed renter by 35") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share childless by 45", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.8, angle(h)) ///
          title("Ownership is the binding constraint on fertility") ///
          subtitle("Childlessness by 45, among renters at 25-30 (PSID)") ///
          note("At every wealth level, those who transition to ownership" ///
               "are 20-35pp less likely to remain childless by 45.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/within_wealth_gap_dollars_v10c.png", replace width(1600)
restore

* ==================================================================
* 2. CONSTRAINT MARGIN: ownership rate + childlessness gap
*    Two-axis: bars = ownership rate, line = gap
* ==================================================================
di _n "============================================="
di    "2. Ownership rate + childlessness gap"
di    "============================================="

preserve
    * Compute childlessness by (bin, transition)
    collapse (mean) childless_by_45 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin owner_by_35)
    reshape wide childless_by_45 n_obs, i(nw_bin) j(owner_by_35)
    gen gap = childless_by_450 - childless_by_451
    gen own_rate = n_obs1 / (n_obs0 + n_obs1)
    rename childless_by_450 cless_renter
    rename childless_by_451 cless_owner

    list
    export delimited using "`outdir'/margin_gap_v10c.csv", replace

    twoway (bar own_rate nw_bin, barw(0.6) fcolor(navy%40) lcolor(navy)) ///
           (connected gap nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick) ///
            yaxis(2)) ///
        , legend(order(1 "Ownership transition rate" 2 "Childlessness gap (renter - owner)") ///
                 position(6) rows(1) size(small)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Ownership rate", size(medium) axis(1)) ///
          ytitle("Childlessness gap (pp)", size(medium) axis(2)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.7, angle(h) axis(1)) ///
          ylabel(0(0.1)0.5, angle(h) axis(2)) ///
          title("The ownership-fertility constraint margin") ///
          subtitle("Among renters at 25-30") ///
          note("Gap = renter childlessness minus owner childlessness." ///
               "Wealth raises ownership rate; at every level, owners have more children.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/margin_gap_v10c.png", replace width(1600)
restore

* ==================================================================
* 3. THREE-PANEL SUMMARY: ownership rate, cless if owner, cless if renter
*    All three series on one axis, by dollar bin
* ==================================================================
di _n "============================================="
di    "3. Three-series summary"
di    "============================================="

preserve
    collapse (mean) childless_by_45 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin owner_by_35)
    reshape wide childless_by_45 n_obs, i(nw_bin) j(owner_by_35)
    gen own_rate = n_obs1 / (n_obs0 + n_obs1)
    rename childless_by_450 cless_renter
    rename childless_by_451 cless_owner

    twoway (bar own_rate nw_bin, barw(0.6) fcolor(gs14) lcolor(gs8)) ///
           (connected cless_renter nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(vlarge) lwidth(thick)) ///
           (connected cless_owner nw_bin, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 "Ownership rate" 2 "Childless if renter" 3 "Childless if owner") ///
                 position(6) rows(1) size(small)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.8, angle(h)) ///
          title("Wealth relaxes the ownership constraint; ownership enables fertility") ///
          subtitle("Among renters at 25-30 (PSID 1984-2019)") ///
          note("Bars: share transitioning to ownership by 35. Lines: childlessness by 45." ///
               "Wealth raises ownership; at every wealth level, owners have more children.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/three_series_v10c.png", replace width(1600)
restore

* ==================================================================
* 4. SAMPLE SIZE TABLE
* ==================================================================
di _n "=== Cell sizes by bin and transition ==="
tab nw_bin owner_by_35

log close
