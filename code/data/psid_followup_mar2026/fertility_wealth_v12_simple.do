clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v12.log", replace text

* ==================================================================
* V12: SIMPLE MOTIVATIONAL FACT — WEALTH PREDICTS FERTILITY
*
* No ownership decomposition. Just: wealthier young adults are more
* likely to have children by a given age. The model explains why.
*
* Key insight: the gradient is steep at 35 (timing), flat at 45
* (convergence). Use parenthood by 35 for the sharpest picture.
* ==================================================================

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

gen first_birth_year = RELCHI1BYEAR
gen own = HOMEOWN
replace own = 0 if own == 2
replace own = . if own == 3

* ==================================================================
* Entry wealth at 25-30, ALL young adults (not just renters)
* ==================================================================
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if !missing(NETWORTHR)
    collapse (mean) nw_dollars = NETWORTHR ///
        (mean) inc_dollars = INCFAMR ///
        (mean) weight_young = IW, by(ID)
    tempfile entry
    save `entry', replace
restore

* Fertility
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
merge 1:1 ID using `fert', keep(match) nogen

gen observed_by_45 = (2019 - birth_year) >= 45
gen observed_by_40 = (2019 - birth_year) >= 40
gen observed_by_35 = (2019 - birth_year) >= 35

gen parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
gen parent_by_40 = (!missing(age_at_first_birth) & age_at_first_birth <= 40)
gen parent_by_35 = (!missing(age_at_first_birth) & age_at_first_birth <= 35)
gen parent_by_30 = (!missing(age_at_first_birth) & age_at_first_birth <= 30)

gen childless_by_35 = 1 - parent_by_35
gen childless_by_40 = 1 - parent_by_40
gen childless_by_45 = 1 - parent_by_45

gen nw_k = nw_dollars / 1000

di _n "=== SAMPLE SIZES ==="
di "All at 25-30 with wealth:"
count
di "Observed by 35:"
count if observed_by_35 == 1
di "Observed by 40:"
count if observed_by_40 == 1
di "Observed by 45:"
count if observed_by_45 == 1

* ==================================================================
* FIGURE 1: Parenthood rate by wealth quintile, at ages 35, 40, 45
* Shows timing gradient: steep at 35, flat at 45
* ==================================================================
di _n "============================================="
di    "FIGURE 1: Parenthood by wealth quintile"
di    "============================================="

preserve
    keep if observed_by_35 == 1

    * Trim wealth at 1-99
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample (observed by 35):"
    count

    xtile wealth_q = nw_dollars [pw=weight_young], n(5)

    * Quintile medians for labels
    tabstat nw_k [aw=weight_young], by(wealth_q) stat(p50 n)

    collapse (mean) parent_by_35 parent_by_40 parent_by_45 ///
        (p50) median_nw_k = nw_k ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_q)
    list

    export delimited using "`outdir'/parenthood_by_quintile_v12.csv", replace

    twoway (connected parent_by_35 wealth_q, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick)) ///
           (connected parent_by_40 wealth_q, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(large) lwidth(medthick)) ///
           (connected parent_by_45 wealth_q, ///
            mcolor(gs8) lcolor(gs8) msymbol(S) msize(large) lwidth(medthick) lpattern(dash)) ///
        , legend(order(1 "Parent by 35" 2 "Parent by 40" 3 "Parent by 45") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth quintile at 25-30", size(medium)) ///
          ytitle("Share who are parents", size(medium)) ///
          xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5") ///
          ylabel(0.3(0.1)0.9, angle(h)) ///
          title("Wealth predicts fertility timing") ///
          subtitle("All young adults at 25-30 (PSID 1984-2019)") ///
          note("Quintiles of net worth at 25-30. N per quintile shown in output.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/parenthood_by_quintile_v12.png", replace width(1600)
restore

* ==================================================================
* FIGURE 2: Childlessness at 35 by dollar bins
* Down-payment-relevant thresholds
* ==================================================================
di _n "============================================="
di    "FIGURE 2: Childlessness at 35 by dollar bins"
di    "============================================="

preserve
    keep if observed_by_35 == 1
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    gen nw_bin = .
    replace nw_bin = 1 if nw_dollars < 0
    replace nw_bin = 2 if nw_dollars >= 0     & nw_dollars < 10000
    replace nw_bin = 3 if nw_dollars >= 10000  & nw_dollars < 25000
    replace nw_bin = 4 if nw_dollars >= 25000  & nw_dollars < 50000
    replace nw_bin = 5 if nw_dollars >= 50000  & nw_dollars < 100000
    replace nw_bin = 6 if nw_dollars >= 100000 & !missing(nw_dollars)

    tab nw_bin

    collapse (mean) childless_by_35 parent_by_35 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)
    list

    export delimited using "`outdir'/childless35_dollars_v12.csv", replace

    graph bar childless_by_35, ///
        over(nw_bin, relabel(1 `""< $0""' 2 `""$0-" "10k""' 3 `""$10-" "25k""' ///
             4 `""$25-" "50k""' 5 `""$50-" "100k""' 6 `""$100k+""')) ///
        bar(1, fcolor(cranberry%70) lcolor(cranberry)) ///
        ytitle("Share childless at 35", size(medium)) ///
        ylabel(0(0.1)0.7, angle(h)) ///
        blabel(bar, format(%4.2f) size(medlarge)) ///
        title("Wealth and childlessness at 35") ///
        subtitle("All young adults at 25-30 (PSID 1984-2019)") ///
        note("Net worth measured at 25-30 in 2019 dollars. Trimmed at p1/p99.") ///
        graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/childless35_dollars_v12.png", replace width(1600)
restore

* ==================================================================
* FIGURE 3: Local polynomial — parenthood by 35 as function of NW
* Smooth, no binning noise
* ==================================================================
di _n "============================================="
di    "FIGURE 3: lpoly parenthood by 35"
di    "============================================="

preserve
    keep if observed_by_35 == 1
    _pctile nw_dollars [pw=weight_young], p(2 98)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample for lpoly:"
    count

    gen nw_plot = min(nw_dollars/1000, 200)

    twoway (lpoly parent_by_35 nw_plot [aw=weight_young], ///
            lcolor(cranberry) lwidth(thick) bwidth(25) degree(1)) ///
           (lpoly parent_by_40 nw_plot [aw=weight_young], ///
            lcolor(navy) lwidth(thick) lpattern(dash) bwidth(25) degree(1)) ///
        , legend(order(1 "Parent by 35" 2 "Parent by 40") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Net worth at 25-30 ($1000s, 2019 real)", size(medium)) ///
          ytitle("Probability of parenthood", size(medium)) ///
          ylabel(0.3(0.1)0.9, angle(h)) ///
          xline(0, lcolor(gs10) lpattern(dot)) ///
          title("Wealth and fertility") ///
          subtitle("Local polynomial smooth, all young adults at 25-30") ///
          note("Local linear, Epanechnikov kernel, bw = $25k. Trimmed p2/p98.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/lpoly_parenthood_v12.png", replace width(1600)
restore

* ==================================================================
* FIGURE 4: Age at first birth by wealth tercile (survival curve)
* Cumulative parenthood by age — curves separate by wealth
* ==================================================================
di _n "============================================="
di    "FIGURE 4: Cumulative parenthood by age x wealth"
di    "============================================="

preserve
    keep if observed_by_45 == 1
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    di "Sample (observed by 45):"
    count

    xtile wealth_t = nw_dollars [pw=weight_young], n(3)
    label define wt 1 "Bottom third" 2 "Middle third" 3 "Top third"
    label values wealth_t wt

    * Expand to age grid
    gen parent_by_25 = (!missing(age_at_first_birth) & age_at_first_birth <= 25)
    gen parent_by_28 = (!missing(age_at_first_birth) & age_at_first_birth <= 28)
    gen parent_by_30 = (!missing(age_at_first_birth) & age_at_first_birth <= 30)
    gen parent_by_32 = (!missing(age_at_first_birth) & age_at_first_birth <= 32)
    gen parent_by_38 = (!missing(age_at_first_birth) & age_at_first_birth <= 38)

    collapse (mean) parent_by_25 parent_by_28 parent_by_30 parent_by_32 ///
        parent_by_35 parent_by_38 parent_by_40 parent_by_45 ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_t)
    list

    export delimited using "`outdir'/cum_parenthood_v12.csv", replace

    * Reshape to long for plotting
    reshape long parent_by_, i(wealth_t) j(age)
    rename parent_by_ cum_parent

    twoway (connected cum_parent age if wealth_t == 1, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(large) lwidth(thick)) ///
           (connected cum_parent age if wealth_t == 2, ///
            mcolor(gs6) lcolor(gs6) msymbol(O) msize(large) lwidth(medthick)) ///
           (connected cum_parent age if wealth_t == 3, ///
            mcolor(navy) lcolor(navy) msymbol(D) msize(large) lwidth(thick)) ///
        , legend(order(1 "Bottom third NW" 2 "Middle third" 3 "Top third NW") ///
                 position(6) rows(1) size(medium)) ///
          xtitle("Age", size(medium)) ///
          ytitle("Cumulative share who are parents", size(medium)) ///
          xlabel(25 28 30 32 35 38 40 45) ///
          ylabel(0(0.1)0.9, angle(h)) ///
          title("Wealth and fertility timing") ///
          subtitle("Cumulative parenthood by age and wealth tercile at 25-30") ///
          note("PSID 1984-2019. Terciles of net worth at 25-30.") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/cum_parenthood_v12.png", replace width(1600)
restore

* ==================================================================
* REGRESSIONS: simple reduced-form
* ==================================================================
di _n "============================================="
di    "REGRESSIONS"
di    "============================================="

preserve
    keep if observed_by_35 == 1
    _pctile nw_dollars [pw=weight_young], p(1 99)
    drop if nw_dollars < r(r1) | nw_dollars > r(r2)

    * Standardize NW for interpretability
    sum nw_dollars [aw=weight_young]
    gen nw_std = (nw_dollars - r(mean)) / r(sd)

    * NW in $10k units
    gen nw_10k = nw_dollars / 10000

    di _n "--- Parent by 35 on NW ($10k) ---"
    reg parent_by_35 nw_10k [pw=weight_young], robust

    di _n "--- Parent by 40 on NW ($10k) ---"
    reg parent_by_40 nw_10k [pw=weight_young], robust

    di _n "--- Parent by 45 on NW ($10k) ---"
    reg parent_by_45 nw_10k [pw=weight_young], robust

    di _n "--- Childless by 35 on NW ($10k) ---"
    reg childless_by_35 nw_10k [pw=weight_young], robust
restore

log close
