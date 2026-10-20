* Aggregate wealth-to-earnings ratio with the 2021 and 2023 PSID waves (Oct 20 2026).
* Same statistic as code/data/starter_home_trap_20261009/08_wealth_earnings_recent.py:
*   reference persons (RELTOHEAD_ == 10) aged 18-85 in 1999+, PSID longitudinal weight > 0, net worth observed,
*   earnings observed and >= 0 for those aged 65 or less;
*   ratio_orig   = sum(w * NETWORTHR) / sum over age <= 65 of (w * EARNINDRRC)
*   ratio_narrow = same with net worth minus farm/business, other real estate and vehicles (missing parts = 0).
* 1999-2019 values come from PSIDSHELF_MOBILITY.dta. 2021 and 2023 are rebuilt with SHELF's rules
* (01 Collect wealth/income, 03 Generate wealth/income, 07 inflation lookup) from FAM2021ER/FAM2023ER;
* the same code rebuilds 2019 from the main-study file and is compared with SHELF as a check.
version 17.0
clear all
set more off
local raw "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/raw_2021_2023"
local shelf "/Users/tommasodesanto/Desktop/Projects/Fertility/PSID/PSIDSHELF_MOBILITY.dta"
local cms "/Users/tommasodesanto/Desktop/Projects/Datasets/PSID-SHELF/Construction_Files/Data/Users/DD/Dropbox (University of Michigan)/Data/PSID/PSID_CMS/PSID_COMPLETE_MAIN_STUDY_1968_2019.dta"
local out "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/output/wealth_earnings_2023"
capture mkdir "`out'"
log using "`out'/build_wealth_earnings_2023.log", replace text

* 1. Family items: RP wages, RP business labor, SP wages, SP business labor, farm/business assets and debts,
*    other real estate assets and debts, vehicles, wealth with equity (labels matched to the 2019 items).
infix fid 2-6 double wrp 6624-6630 double brp 6550-6556 double wsp 6733-6739 double bsp 6659-6665 ///
    double fa 7132-7140 double fd 7142-7150 double ra 7178-7186 double rd 7188-7196 double auto 7210-7218 double nw 7324-7332 ///
    using "`raw'/fam2021er/FAM2021ER.txt", clear
gen year = 2021
tempfile f21
save `f21'
infix fid 2-6 double wrp 6378-6384 double brp 6304-6310 double wsp 6487-6493 double bsp 6413-6419 ///
    double fa 6886-6894 double fd 6896-6904 double ra 6932-6940 double rd 6942-6950 double auto 6964-6972 double nw 7078-7086 ///
    using "`raw'/fam2023er/FAM2023ER.txt", clear
gen year = 2023
append using `f21'
isid year fid
tempfile fam
save `fam'

* 2. Persons: interview number, relation to reference person, age and weight from the 2023 individual file.
infix ER30001 2-5 ER30002 6-8 ///
    fid2019 5250-5254 seq2019 5255-5256 relraw2019 5257-5258 age2019 5259-5261 double wt2019 5672-5678 ///
    fid2021 5684-5688 seq2021 5689-5690 relraw2021 5691-5692 age2021 5693-5695 double wt2021 6107-6113 ///
    fid2023 6119-6123 seq2023 6124-6125 relraw2023 6126-6127 age2023 6128-6130 double wt2023 6542-6548 ///
    using "`raw'/ind2023er/IND2023ER.txt", clear
gen double ID = ER30001*1000 + ER30002
drop ER30001 ER30002
reshape long fid seq relraw age wt, i(ID) j(year)
keep if relraw == 10 & inrange(seq, 1, 20)
tempfile ind
save `ind'
keep if inlist(year, 2021, 2023)
merge m:1 year fid using `fam', keep(master match) nogenerate
tempfile pnew
save `pnew'
use ID ER77315 ER77296 ER77343 ER77324 ER77451 ER77453 ER77465 ER77467 ER77473 ER77511 using "`cms'", clear
rename (ER77315 ER77296 ER77343 ER77324 ER77451 ER77453 ER77465 ER77467 ER77473 ER77511) (wrp brp wsp bsp fa fd ra rd auto nw)
gen year = 2019
tempfile f19
save `f19'
use `ind' if year == 2019, clear
merge 1:1 ID year using `f19', keep(master match) nogenerate
append using `pnew'

* 3. SHELF rules. Wealth items 2011+: keep -99999997..999999997 (0 included). Earnings 2011+: keep 0..9999997.
foreach v in nw fa fd ra rd auto {
    gen double `v'_n = `v' if inrange(`v', -99999997, 999999997)
}
foreach v in wrp brp wsp bsp {
    gen double `v'_n = `v' if inrange(`v', 0, 9999997)
}
gen double earnrp = wrp_n + brp_n
gen double earnsp = wsp_n + bsp_n
gen double earnrc = earnrp + earnsp if !missing(earnrp, earnsp)
replace earnrc = earnrp if !missing(earnrp) & missing(earnsp)
* SHELF 07 index lookup: wealth at index(year-8), earnings (tax year) at index(year-9), radix index(2014).
gen double pw = cond(year==2019, 123.24161, cond(year==2021, 127.19513, 129.21168))
gen double pe = cond(year==2019, 120.20075, cond(year==2021, 125.5398, 128.9754))
gen double NETWORTHR = nw_n * 128.9754 / pw
gen double WLTHFBIZTOTR = (fa_n - fd_n) * 128.9754 / pw
gen double WLTHREALTOTR = (ra_n - rd_n) * 128.9754 / pw
gen double WLTHAUTOTOTR = auto_n * 128.9754 / pw
gen double EARNINDRRC = earnrc * 128.9754 / pe
gen double IW = wt if wt != 0
gen double AGEREP = age if inrange(age, 1, 199)
gen RELTOHEAD_ = relraw
keep ID year AGEREP RELTOHEAD_ IW NETWORTHR EARNINDRRC WLTHFBIZTOTR WLTHREALTOTR WLTHAUTOTOTR
tempfile rebuilt
save `rebuilt'

* 4. Check 2019 against SHELF (reference persons).
use ID year AGEREP RELTOHEAD_ IW NETWORTHR EARNINDRRC WLTHFBIZTOTR WLTHREALTOTR WLTHAUTOTOTR using "`shelf'" if year >= 1999, clear
tempfile shelfw
save `shelfw'
keep if year == 2019 & RELTOHEAD_ == 10
rename (AGEREP IW NETWORTHR EARNINDRRC WLTHFBIZTOTR WLTHREALTOTR WLTHAUTOTOTR) s_=
merge 1:1 ID year using `rebuilt', keep(master match)
tab _merge
foreach v in AGEREP IW NETWORTHR EARNINDRRC WLTHFBIZTOTR WLTHREALTOTR WLTHAUTOTOTR {
    quietly count if _merge == 3 & !missing(s_`v', `v')
    local n = r(N)
    quietly count if _merge == 3 & !missing(s_`v', `v') & abs(s_`v' - `v') > max(1, 0.0005*abs(s_`v'))
    local d = r(N)
    quietly count if _merge == 3 & missing(s_`v') != missing(`v')
    di "CHECK2019 `v' both=`n' differ=`d' missing_mismatch=`r(N)'"
}

* 5. Ratios by wave and pooled windows: SHELF through 2019, rebuilt 2021 and 2023; rebuilt 2019 shown as a check.
use `rebuilt' if inlist(year, 2021, 2023), clear
gen byte src = 1
append using `shelfw'
replace src = 0 if missing(src)
append using `rebuilt'
replace src = 2 if missing(src)
keep if RELTOHEAD_ == 10
keep if inrange(AGEREP, 18, 85) & !missing(IW) & IW > 0 & !missing(NETWORTHR) & (AGEREP > 65 | (!missing(EARNINDRRC) & EARNINDRRC >= 0))
foreach v in WLTHFBIZTOTR WLTHREALTOTR WLTHAUTOTOTR {
    replace `v' = 0 if missing(`v')
}
gen double nw_c = NETWORTHR - WLTHFBIZTOTR - WLTHREALTOTR - WLTHAUTOTOTR
* Sensitivity: SHELF's lookup puts wealth at the PCE index of year-8 and earnings at year-9. Undo it and
* convert nominal wealth (survey year) and nominal earnings (previous calendar year) at their own years' PCE.
* PCE list as in SHELF 07 (calendar years); 2023 = 2022 x 120.511/116.100 (FRED annual PCEPI).
gen double pidx_y = .
gen double pidx_y1 = .
gen double pidx_w = .
gen double pidx_e = .
local years 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
local vals 106.00946 109.06383 112.1349 115.01231 118.41583 118.08571 120.20075 123.24161 125.5398 127.19513 128.9754 129.21168 130.52007 132.79867 135.51683 137.46421 138.95344 144.75197 154.10264 159.95637
forvalues k = 1/20 {
    local yy : word `k' of `years'
    local vv : word `k' of `vals'
    replace pidx_y = `vv' if year == `yy'
    replace pidx_y1 = `vv' if year - 1 == `yy'
    replace pidx_w = `vv' if year - 8 == `yy'
    replace pidx_e = `vv' if year - 9 == `yy'
}
gen double cw = (pidx_w / 128.9754) * (154.10264 / pidx_y)
gen double ce = (pidx_e / 128.9754) * (154.10264 / pidx_y1)
gen double wnw = IW * NETWORTHR
gen double wnwc = IW * nw_c
gen double wearn = IW * EARNINDRRC if AGEREP <= 65
gen double wnw_pce = wnw * cw
gen double wnwc_pce = wnwc * cw
gen double wearn_pce = wearn * ce
tempname R
postfile `R' str24 window str10 source double N N_earn ratio_orig ratio_narrow ratio_orig_pce ratio_narrow_pce using "`out'/ratios.dta", replace
program define postratio
    args handle window source cond
    quietly count if `cond'
    local N = r(N)
    quietly count if `cond' & AGEREP <= 65
    local Ne = r(N)
    quietly summarize wnw if `cond'
    local a = r(sum)
    quietly summarize wnwc if `cond'
    local b = r(sum)
    quietly summarize wearn if `cond'
    local d = r(sum)
    quietly summarize wnw_pce if `cond'
    local a2 = r(sum)
    quietly summarize wnwc_pce if `cond'
    local b2 = r(sum)
    quietly summarize wearn_pce if `cond'
    local d2 = r(sum)
    post `handle' ("`window'") ("`source'") (`N') (`Ne') (`a'/`d') (`b'/`d') (`a2'/`d2') (`b2'/`d2')
end
foreach y in 2013 2015 2017 2019 {
    postratio `R' "`y'" "shelf" "src == 0 & year == `y'"
}
postratio `R' "2019" "rebuilt" "src == 2 & year == 2019"
foreach y in 2021 2023 {
    postratio `R' "`y'" "rebuilt" "src == 1 & year == `y'"
}
postratio `R' "pooled_2017_2019" "shelf" "src == 0 & inlist(year, 2017, 2019)"
postratio `R' "pooled_2015_2019" "shelf" "src == 0 & inlist(year, 2015, 2017, 2019)"
postratio `R' "pooled_2019_2021" "mixed" "(src == 0 & year == 2019) | (src == 1 & year == 2021)"
postratio `R' "pooled_2021_2023" "rebuilt" "src == 1 & inlist(year, 2021, 2023)"
postratio `R' "pooled_2019_2023" "mixed" "(src == 0 & year == 2019) | (src == 1 & inlist(year, 2021, 2023))"
postclose `R'
use "`out'/ratios.dta", clear
format ratio_orig ratio_narrow %9.4f
list, noobs sep(0)
export delimited using "`out'/wealth_earnings_by_wave_2023.csv", replace
log close
