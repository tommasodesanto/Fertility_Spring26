clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v10.log", replace text

* ==================================================================
* V10: DOSE-RESPONSE OF WEALTH ON OWNERSHIP AND FERTILITY
*
* Goal: Show that wealth operates on fertility THROUGH the ownership
*       channel, concentrated at the constraint boundary.
*
* Four approaches:
*   A. Dollar-bin dose-response (absolute NW, down-payment thresholds)
*   B. Local polynomial smooth (lpoly)
*   C. Expanded sample (parenthood by 40, not 45)
*   D. All young adults (not just initial renters)
* ==================================================================

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR NETWORTH2R HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen own = HOMEOWN

* ==================================================================
* BUILD PERSON-LEVEL PANEL COMPONENTS
* ==================================================================

* --- Entry characteristics at 25-30 ---
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if !missing(own) & !missing(NETWORTHR)
    collapse (mean) own_rate_2530 = own ///
        (mean) nw_dollars = NETWORTHR ///
        (mean) liq_nw_dollars = NETWORTH2R ///
        (mean) inc_dollars = INCFAMR ///
        (mean) weight_young = IW, by(ID)
    gen renter_2530 = (own_rate_2530 < 0.5)
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

* --- Tenure at 31-40 (for expanded sample) ---
preserve
    keep if inrange(AGEREP, 31, 40) & !missing(own)
    collapse (mean) own_rate_3140 = own, by(ID)
    gen owner_by_40 = (own_rate_3140 >= 0.5) if !missing(own_rate_3140)
    keep ID owner_by_40
    tempfile t3140
    save `t3140', replace
restore

* --- Birth year and fertility ---
preserve
    collapse (min) min_year = year (min) min_age = AGEREP ///
        (first) first_birth_year_b = first_birth_year, by(ID)
    gen birth_year = min_year - min_age
    gen age_at_first_birth = first_birth_year_b - birth_year if !missing(first_birth_year_b)
    keep ID birth_year age_at_first_birth
    tempfile fert
    save `fert', replace
restore

* --- Merge ---
use `entry', clear
merge 1:1 ID using `t3135', keep(match master) nogen
merge 1:1 ID using `t3140', keep(match master) nogen
merge 1:1 ID using `fert', keep(match) nogen

gen observed_by_45 = (2019 - birth_year) >= 45
gen observed_by_40 = (2019 - birth_year) >= 40

gen parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
gen parent_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
gen parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)

* Real 2019 dollars, in thousands
gen nw_k = nw_dollars / 1000

di _n "=== FULL SAMPLE (all young adults at 25-30 with wealth data) ==="
count
tab renter_2530
tab observed_by_45
tab observed_by_40

* ==================================================================
* APPROACH A: DOLLAR-BIN DOSE-RESPONSE (RENTERS, OBSERVED BY 45)
* Bins map to down-payment thresholds for a median-priced home
* ==================================================================
di _n "============================================="
di    "APPROACH A: Dollar-bin dose-response"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_45 == 1 & !missing(owner_by_35)

    * Trim 1-99 on NW in dollars
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di _n "Sample size (renters at 25-30, observed by 45):"
    count
    sum nw_dollars [w=weight_young], d

    * Bins in absolute dollars (2019 real)
    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    tab nw_bin

    collapse (mean) owner_by_35 parent_by_45 ///
        (p50) median_nw_k = nw_k ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)
    list

    export delimited using "`outdir'/dose_response_dollars_v10.csv", replace

    twoway (bar owner_by_35 nw_bin, barw(0.6) fcolor(navy%50) lcolor(navy)) ///
           (connected parent_by_45 nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 "Became owner by 35" 2 "Became parent by 45") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.9, angle(h)) ///
          title("Wealth, ownership transition, and fertility") ///
          subtitle("Among renters at 25-30 (PSID 1984-2019)") ///
          note("N shown per bin. Bars = ownership transition rate; line = parenthood rate." ///
               "Bins correspond to down-payment thresholds for median-priced home.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/dose_response_dollars_v10.png", replace width(1600)
restore

* ==================================================================
* APPROACH B: LOCAL POLYNOMIAL SMOOTH (RENTERS, OBSERVED BY 45)
* ==================================================================
di _n "============================================="
di    "APPROACH B: Local polynomial smooth"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_45 == 1 & !missing(owner_by_35)

    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample for lpoly:"
    count

    * Cap x-axis at $150k for visual clarity
    gen nw_plot = min(nw_dollars/1000, 150)

    * lpoly does not accept pweights; use aweight instead
    twoway (lpoly owner_by_35 nw_plot [aw=weight_young], ///
            lcolor(navy) lwidth(thick) bwidth(20) degree(1)) ///
           (lpoly parent_by_45 nw_plot [aw=weight_young], ///
            lcolor(cranberry) lwidth(thick) lpattern(dash) bwidth(20) degree(1)) ///
        , legend(order(1 "Pr(own by 35)" 2 "Pr(parent by 45)") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 ($1000s, 2019 real)", size(medium)) ///
          ytitle("Probability", size(medium)) ///
          ylabel(0(0.1)0.9, angle(h)) ///
          xline(0, lcolor(gs10) lpattern(dot)) ///
          title("Ownership and parenthood rise together with wealth") ///
          subtitle("Local polynomial smooth, renters at 25-30") ///
          note("Epanechnikov kernel, bandwidth = $20k. Trimmed at p1/p99.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/lpoly_dual_v10.png", replace width(1600)
restore

* ==================================================================
* APPROACH C: EXPANDED SAMPLE — PARENTHOOD BY 40 (NOT 45)
* Roughly doubles/triples the sample
* ==================================================================
di _n "============================================="
di    "APPROACH C: Expanded sample (parenthood by 40)"
di    "============================================="

preserve
    keep if renter_2530 == 1 & observed_by_40 == 1 & !missing(owner_by_35)

    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample size (observed by 40):"
    count

    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    tab nw_bin

    collapse (mean) owner_by_35 parent_by_40 ///
        (p50) median_nw_k = nw_k ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)
    list

    export delimited using "`outdir'/dose_response_by40_v10.csv", replace

    twoway (bar owner_by_35 nw_bin, barw(0.6) fcolor(navy%50) lcolor(navy)) ///
           (connected parent_by_40 nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 "Became owner by 35" 2 "Became parent by 40") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.9, angle(h)) ///
          title("Wealth, ownership, and fertility (expanded sample)") ///
          subtitle("Among renters at 25-30, parenthood by 40") ///
          note("Larger sample from relaxing age-45 observability requirement.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/dose_response_by40_v10.png", replace width(1600)
restore

* ==================================================================
* APPROACH D: ALL YOUNG ADULTS (NOT JUST INITIAL RENTERS)
* Uses homeownership rate (not transition) as outcome
* ==================================================================
di _n "============================================="
di    "APPROACH D: All young adults at 25-30"
di    "============================================="

preserve
    keep if observed_by_45 == 1 & !missing(owner_by_35)

    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample size (all young adults, observed by 45):"
    count

    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    tab nw_bin

    collapse (mean) owner_by_35 parent_by_45 ///
        (p50) median_nw_k = nw_k ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)
    list

    export delimited using "`outdir'/dose_response_all_v10.csv", replace

    twoway (bar owner_by_35 nw_bin, barw(0.6) fcolor(navy%50) lcolor(navy)) ///
           (connected parent_by_45 nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 "Homeowner by 35" 2 "Parent by 45") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 (2019 $)", size(medium)) ///
          ytitle("Share", size(medium)) ///
          xlabel(1 `""< $0""' 2 `""$0-10k""' 3 `""$10-25k""' 4 `""$25-50k""' ///
                 5 `""$50-100k""' 6 `""$100k+""', angle(20)) ///
          ylabel(0(0.1)0.9, angle(h)) ///
          title("Wealth predicts both ownership and parenthood") ///
          subtitle("All young adults at 25-30 (PSID 1984-2019)") ///
          note("Includes both initial renters and owners at 25-30.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/dose_response_all_v10.png", replace width(1600)
restore

* ==================================================================
* APPROACH E: LPOLY ON ALL YOUNG ADULTS (largest sample)
* ==================================================================
di _n "============================================="
di    "APPROACH E: Local polynomial, all young adults"
di    "============================================="

preserve
    keep if observed_by_45 == 1 & !missing(owner_by_35)

    _pctile nw_dollars [pw=weight_young], p(2 98)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample for lpoly (all):"
    count

    gen nw_plot = min(nw_dollars/1000, 200)

    twoway (lpoly owner_by_35 nw_plot [aw=weight_young], ///
            lcolor(navy) lwidth(thick) bwidth(25) degree(1)) ///
           (lpoly parent_by_45 nw_plot [aw=weight_young], ///
            lcolor(cranberry) lwidth(thick) lpattern(dash) bwidth(25) degree(1)) ///
        , legend(order(1 "Pr(own by 35)" 2 "Pr(parent by 45)") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 ($1000s, 2019 real)", size(medium)) ///
          ytitle("Probability", size(medium)) ///
          ylabel(0(0.1)1.0, angle(h)) ///
          xline(0, lcolor(gs10) lpattern(dot)) ///
          title("Ownership and parenthood: nonparametric dose-response") ///
          subtitle("All young adults at 25-30") ///
          note("Local linear smooth, Epanechnikov kernel, bandwidth = $25k.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/lpoly_all_v10.png", replace width(1600)
restore

* ==================================================================
* SUMMARY TABLE: Sample sizes across approaches
* ==================================================================
di _n "============================================="
di    "SAMPLE SIZE COMPARISON"
di    "============================================="
di "Approach A (renters, obs by 45): " _continue
count if renter_2530==1 & observed_by_45==1 & !missing(owner_by_35)
di "Approach C (renters, obs by 40): " _continue
count if renter_2530==1 & observed_by_40==1 & !missing(owner_by_35)
di "Approach D (all, obs by 45):     " _continue
count if observed_by_45==1 & !missing(owner_by_35)
di "Approach D (all, obs by 40):     " _continue
count if observed_by_40==1 & !missing(owner_by_35)

log close
