clear all
set more off
set processors 1
input ID year
1 2000
1 2001
2 2000
2 2001
3 2002
4 2002
4 2003
end
local more = 1
local prune_iter = 0
local removed_person_ids = 0
local removed_years = 0
while `more' {
    sort ID year
    by ID: gen long audit_n_id = _N
    bysort year: gen long audit_n_year = _N
    quietly count if audit_n_id == 1
    local one_id_rows = r(N)
    quietly count if audit_n_year == 1
    local one_year_rows = r(N)
    if `one_id_rows' == 0 & `one_year_rows' == 0 {
        local more = 0
        drop audit_n_id audit_n_year
    }
    else {
        egen byte audit_id_tag = tag(ID) if audit_n_id == 1
        quietly count if audit_id_tag
        local one_ids = r(N)
        egen byte audit_year_tag = tag(year) if audit_n_year == 1
        quietly count if audit_year_tag
        local one_years = r(N)
        local removed_person_ids = `removed_person_ids' + `one_ids'
        local removed_years = `removed_years' + `one_years'
        local prune_iter = `prune_iter' + 1
        drop if audit_n_id == 1 | audit_n_year == 1
        drop audit_n_id audit_n_year audit_id_tag audit_year_tag
    }
}

assert _N == 4
assert ID <= 2
assert `prune_iter' == 2
assert `removed_person_ids' == 2
di "SINGLETON_SMOKE_PASSED"
