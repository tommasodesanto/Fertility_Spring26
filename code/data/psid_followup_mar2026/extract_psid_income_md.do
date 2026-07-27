version 17
clear all
set more off

args source_path output_path
if `"`source_path'"' == "" | `"`output_path'"' == "" {
    display as error "Usage: do extract_psid_income_md.do <source.dta> <output.dta>"
    exit 198
}

use ID year RELTOHEAD_ AGEREP DEATHYEAR EARNINDRRC IW using `"`source_path'"', clear

keep if RELTOHEAD_ == 10
keep if inrange(year, 1984, 2019)
keep if missing(DEATHYEAR) | year <= DEATHYEAR
keep if inrange(AGEREP, 25, 60)
keep if EARNINDRRC > 0 & EARNINDRRC < .
keep if IW > 0 & IW < .

isid ID year
compress
save `"`output_path'"', replace
