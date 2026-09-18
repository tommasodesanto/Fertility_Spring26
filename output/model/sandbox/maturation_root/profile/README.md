# Children at home by parent age: model vs. ACS

**Model**: spec `maturation_parent_age_psi_root` (retained E5F calibration, psi_child held fixed, one GE stationary evaluation, Nb=120/J=17). `m` = children currently at home, read off `sol.g[b,tenure,i,j,n,cs]` under `child_state_mode=independent_count` (solver.py `current_child_bin_dt`): m(n,cs) = cs if cs<=n else 0, where n = children ever born (parity axis) and cs = child-state axis. Model ages are 18 + j*4 for j=0..16 (P.age_start, P.da).

**Data**: ACS 2005-2006 pooled (extract27.dta), household heads only (gq in {1,2}, pernum==1, relate==1, hhwt>0, age 18-85), same head filter as output/model/e5f_matched_pf_20260909a/design_research/housing/inspect_early_housing.py minus its ownershp/rooms restriction (not relevant to dependents). "Children of any age at home" = IPUMS NCHILD directly. "Children under 18 (lower bound)" = NCHILD only for households where ELDCH<18 (so every counted child is confirmed a minor); mixed adult+minor-child households (YNGCH<18 and ELDCH>=18) are excluded from that lower bound rather than zeroed, since ACS gives no per-child age list to split them. 4-year age bins matching model ages.

**mu (implied per-period exit probability)**: A_m = 18 years (expected duration a child stays at home), period length = 4 years, so under the model's memoryless per-period exit process mu = period_years / A_m = 0.2222 per 4-year period (a per-child-period hazard, not an age-dependent schedule).

## Key numbers (mean children at home)
| Age | Model E[m] | ACS mean NCHILD (any age) | ACS mean <18 (lower bound) | Model share(m>0) | ACS share(NCHILD>0) |
|---:|---:|---:|---:|---:|---:|
| 26 | 0.897 | 0.896 | 0.896 | 0.603 | 0.476 |
| 30 | 1.180 | 1.176 | 1.169 | 0.685 | 0.582 |
| 42 | 1.464 | 1.291 | 1.081 | 0.751 | 0.647 |
| 58 | 0.062 | 0.298 | 0.061 | 0.061 | 0.220 |
| 66 | 0.000 | 0.169 | 0.012 | 0.000 | 0.139 |
| 74 | 0.000 | 0.139 | 0.003 | 0.000 | 0.119 |

## Reading
- Young (22-34): model 1.18 vs. ACS-under-18 1.17 at age 30 -- model overstates dependents at this margin.
- Middle (38-54): model 1.46 vs. ACS-under-18 1.08 at age 42 -- model overstates dependents at this margin.
- Old (58+): model 0.00 vs. ACS-under-18 0.01 at age 66 -- the model's memoryless mu-exit process (0.222/period) implies a geometric tail of children still coded 'at home' at old parent ages that ACS households do not show; model understates dependents at this margin.

ACS receipt: {"source": "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/Spatial_aggregate_withmicrodata/raw_data/extract27.dta", "years": [2005, 2006], "year_meta": [{"year": 2005, "raw_records": 2878380, "head_records": 1131894, "head_weight": 108794850.0}, {"year": 2006, "raw_records": 2969741, "head_records": 1135397, "head_weight": 109258245.0}], "n_heads_total": 2267291, "weight_total": 218053095.0}
