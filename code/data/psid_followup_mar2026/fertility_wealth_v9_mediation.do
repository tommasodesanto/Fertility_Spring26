clear all
set more off
version 17.0

local root "/Users/tommasodesanto/Desktop/Projects/Fertility"
local dta  "`root'/PSID/PSIDSHELF_MOBILITY.dta"
local outdir "`root'/Fertility_Spring26/code/data/psid_followup_mar2026/output/fertility_wealth_v1"
capture log close _all
log using "`outdir'/fertility_wealth_v9.log", replace text

use ID year AGEREP SEX DEATHYEAR RELCHIREP RELCHI1BYEAR ///
    INCFAMR NETWORTHR HOMEOWN IW using "`dta'", clear

drop if year > DEATHYEAR
drop if missing(AGEREP) | AGEREP < 18
keep if inrange(year, 1984, 2019)

replace HOMEOWN = 0 if HOMEOWN == 2
replace HOMEOWN = . if HOMEOWN == 3

gen first_birth_year = RELCHI1BYEAR
gen pre_first_birth = (missing(first_birth_year) | year < first_birth_year)
gen total_nw_to_inc = NETWORTHR / INCFAMR if INCFAMR > 1000 & !missing(NETWORTHR, INCFAMR)
gen own = HOMEOWN

* ------------------------------------------------------------------
* Sample: renters at 25-30, pre-birth
* ------------------------------------------------------------------
preserve
    keep if inrange(AGEREP, 25, 30)
    keep if pre_first_birth == 1
    keep if !missing(own)
    collapse (mean) own_rate_2530 = own ///
        (mean) total_nw_2530 = total_nw_to_inc ///
        (mean) nw_dollars = NETWORTHR ///
        (mean) inc_dollars = INCFAMR ///
        (mean) weight_young = IW, by(ID)
    gen renter_2530 = (own_rate_2530 < 0.5)
    keep if renter_2530 == 1
    keep if !missing(total_nw_2530)
    tempfile r2530
    save `r2530', replace
restore

* Tenure at 31-35
preserve
    keep if inrange(AGEREP, 31, 35)
    keep if !missing(own)
    collapse (mean) own_rate_3135 = own, by(ID)
    gen owner_by_35 = (own_rate_3135 >= 0.5) if !missing(own_rate_3135)
    keep ID owner_by_35
    tempfile t3135
    save `t3135', replace
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

use `r2530', clear
merge 1:1 ID using `t3135', keep(match master) nogen
merge 1:1 ID using `fert', keep(match) nogen

gen observed_by_45 = (2019 - birth_year) >= 45
keep if observed_by_45 == 1
keep if !missing(owner_by_35)

gen parent_by_45 = (!missing(age_at_first_birth) & age_at_first_birth <= 45)
gen childless_by_45 = 1 - parent_by_45

* Trim
_pctile total_nw_2530 [pw=weight_young], p(1 99)
drop if total_nw_2530 < r(r1) | total_nw_2530 > r(r2)

di _n "=== Sample size ==="
count

* ------------------------------------------------------------------
* APPROACH 1: Within-wealth conditional fertility gap
* For each wealth quintile: childlessness by transition status
* This shows: AT EVERY WEALTH LEVEL, the ownership constraint binds
* ------------------------------------------------------------------
xtile wealth_q = total_nw_2530 [pw=weight_young], n(5)

preserve
    collapse (mean) childless_by_45 ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_q owner_by_35)
    gen transition = cond(owner_by_35 == 1, "Became owner", "Stayed renter")
    list, sep(2)
    export delimited using "`outdir'/within_wealth_gap_v9.csv", replace

    twoway (connected childless_by_45 wealth_q if owner_by_35==1, ///
            mcolor(navy) lcolor(navy) msymbol(O) msize(large) lwidth(medthick)) ///
           (connected childless_by_45 wealth_q if owner_by_35==0, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(S) msize(large) lwidth(medthick)) ///
        , legend(order(1 "Became owner by 35" 2 "Stayed renter by 35") ///
                 position(6) rows(1)) ///
          xtitle("Entry wealth quintile at 25-30 (among renters)") ///
          ytitle("Share childless by 45") ///
          xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5") ///
          ylabel(0(0.1)0.8, angle(h)) ///
          title("Ownership gap in fertility persists across wealth levels") ///
          subtitle("Among initial renters at 25-30") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/within_wealth_gap_v9.png", replace width(1600)
restore

* ------------------------------------------------------------------
* APPROACH 2: Ownership transition rate + parenthood by quintile
* Two-panel: how wealth determines both
* ------------------------------------------------------------------
preserve
    collapse (mean) owner_by_35 parent_by_45 ///
        (p50) median_nw = nw_dollars ///
        (count) n_obs = ID [pweight=weight_young], by(wealth_q)
    list
    export delimited using "`outdir'/dual_rates_quintile_v9.csv", replace

    twoway (bar owner_by_35 wealth_q, barw(0.6) fcolor(navy%60) lcolor(navy)) ///
           (connected parent_by_45 wealth_q, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick) ///
            yaxis(1)) ///
        , legend(order(1 "Ownership rate by 35" 2 "Parenthood rate by 45") ///
                 position(6) rows(1)) ///
          xtitle("Entry wealth quintile at 25-30 (among renters)") ///
          ytitle("Share") ///
          xlabel(1 "Q1 (neg.)" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5 (high)") ///
          ylabel(0(0.1)0.8, angle(h)) ///
          title("Wealth predicts both ownership transition and parenthood") ///
          subtitle("Among renters at 25-30, by entry net worth / income") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/dual_rates_v9.png", replace width(1600)
restore

* ------------------------------------------------------------------
* APPROACH 3: Mediation analysis — does wealth affect fertility
* through ownership, or independently?
* ------------------------------------------------------------------
di _n "=== MEDIATION: Wealth → Ownership → Fertility ==="
di _n "--- (A) Total effect: Wealth → Fertility ---"
reg parent_by_45 total_nw_2530 [pw=weight_young], robust
est store total_effect

di _n "--- (B) First stage: Wealth → Ownership ---"
reg owner_by_35 total_nw_2530 [pw=weight_young], robust
est store first_stage

di _n "--- (C) Controlled: Wealth + Ownership → Fertility ---"
reg parent_by_45 total_nw_2530 owner_by_35 [pw=weight_young], robust
est store controlled

di _n "--- Summary ---"
est table total_effect first_stage controlled, b se

di _n "Interpretation:"
di "If the wealth coefficient shrinks substantially when controlling"
di "for ownership, the wealth→fertility channel operates THROUGH ownership."
di "That is the 'margin' — wealth relaxes the ownership constraint,"
di "and ownership enables fertility."

* ------------------------------------------------------------------
* APPROACH 4: Dose-response in meaningful wealth bins
* Use NW/income ratio bins that correspond to real thresholds
* (e.g., can you make a down payment?)
* ------------------------------------------------------------------
gen nw_bin = .
replace nw_bin = 1 if total_nw_2530 < -0.5
replace nw_bin = 2 if total_nw_2530 >= -0.5 & total_nw_2530 < 0
replace nw_bin = 3 if total_nw_2530 >= 0 & total_nw_2530 < 0.2
replace nw_bin = 4 if total_nw_2530 >= 0.2 & total_nw_2530 < 0.5
replace nw_bin = 5 if total_nw_2530 >= 0.5 & total_nw_2530 < 1.0
replace nw_bin = 6 if total_nw_2530 >= 1.0 & !missing(total_nw_2530)

label define nwbin 1 "< -0.5" 2 "-0.5 to 0" 3 "0 to 0.2" 4 "0.2 to 0.5" ///
    5 "0.5 to 1" 6 "> 1"
label values nw_bin nwbin

preserve
    collapse (mean) owner_by_35 parent_by_45 childless_by_45 ///
        (count) n_obs = ID [pweight=weight_young], by(nw_bin)
    list
    export delimited using "`outdir'/dose_response_v9.csv", replace

    twoway (bar owner_by_35 nw_bin, barw(0.6) fcolor(navy%50) lcolor(navy)) ///
           (connected parent_by_45 nw_bin, ///
            mcolor(cranberry) lcolor(cranberry) msymbol(D) msize(vlarge) lwidth(thick)) ///
        , legend(order(1 "Ownership rate by 35" 2 "Parenthood rate by 45") ///
                 position(6) rows(1)) ///
          xtitle("Net worth / income at 25-30") ///
          ytitle("Share") ///
          xlabel(1 "< -0.5" 2 "-0.5 to 0" 3 "0 to 0.2" 4 "0.2 to 0.5" ///
                 5 "0.5 to 1" 6 "> 1", angle(20)) ///
          ylabel(0(0.1)0.8, angle(h)) ///
          title("The down-payment margin: wealth, ownership, and fertility") ///
          subtitle("Among renters at 25-30") ///
          graphregion(color(white)) scheme(s1color)
    graph export "`outdir'/dose_response_v9.png", replace width(1600)
restore

log close
