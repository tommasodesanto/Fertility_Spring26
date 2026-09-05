Verdict: **the claim is algebraically correct for this specification**, conditional on the stated cohort remaining in the regression sample. The primary contrast is not intrinsically normalization-invariant for cohorts with no observed \(K=-2\); the ado does not impose a restriction that fixes this.

Let \(x_{igt}=1\{G_i=g,K_{it}=t\}\). The builder passes every event dummy except \(K=-2\): \(K=-1\), \(0,\ldots,10\), \(K>10\), and both pre-event tails ([builder lines 137–147](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do:137)). If cohort \(g\) has no \(K=-2\) observation, then on every retained row for its members,
\[
\sum_{t\ne-2}x_{igt}=1.
\]
Thus, adding \(c_g\) to every supported \(g\times t\) coefficient and subtracting \(c_g\) from that cohort’s person fixed effects leaves fitted values unchanged. Year FE and age/education controls are unchanged. This is a genuine rank deficiency, not ordinary changing cohort composition.

The ado first calculates cohort-by-event shares before `reghdfe`: it marks `touse` ([ado lines 12–14](/Users/tommasodesanto/Library/Application%20Support/Stata/ado/plus/e/eventstudyinteract.ado:12)), estimates shares in weighted regressions ([lines 40–47](/Users/tommasodesanto/Library/Application%20Support/Stata/ado/plus/e/eventstudyinteract.ado:40)), and only later runs `reghdfe` ([line 85](/Users/tommasodesanto/Library/Application%20Support/Stata/ado/plus/e/eventstudyinteract.ado:85)). It forms \(b_{IW,t}=\sum_g w_{gt}\delta_{gt}\) ([lines 103–108](/Users/tommasodesanto/Library/Application%20Support/Stata/ado/plus/e/eventstudyinteract.ado:103)). Hence the null shift changes the contrast by
\[
(w_{g,+3}-w_{g,-1})c_g.
\]

The supplied input support confirms unequal endpoint weights among no-\(K=-2\) cohorts. Collectively, their endpoint shares are 0.419136 at \(-1\) and 0.460097 at \(+3\) ([summary lines 3 and 7](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/e5f_first_birth_measurement_review_20260905a/empirical_support_summary.csv:3)); e.g., cohort 1986 has positive unequal shares. A one-unit \(c_g\) shift is only an illustrative normalization change, not an economic effect.

Caveats: the CSV is explicitly pre-final-estimation input support. The builder already removes missing rooms, age, education, and nonpositive weights ([lines 86–115](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/code/data/psid_followup_mar2026/sa_rooms_first_birth_household_aligned_v1.do:86)), so those do not create an additional endpoint restriction. But actual `e(sample)`, singleton pruning, and `reghdfe`’s omitted-column normalization require the timed-out regression’s output. They determine the reported numerical coefficients, not the existence of this null direction for retained no-base cohorts.

Strongest defensible wording: “The aggregate input support establishes a potential cohort-specific fixed-effect normalization null direction; because the IW endpoint weights differ, the reported contrast requires verification against the realized estimation sample and omitted-column normalization before being treated as normalization-invariant.”