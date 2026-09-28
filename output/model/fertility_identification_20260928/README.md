# Current-point fertility identification and reoptimization

Author-authorized September28 follow-up. Smoke submitted as18715827; main search awaits lead acceptance.
Frozen reference is overnight selected block0506, primary loss19.581310760.
No model, earnings, entry distribution, target value, grid, bound, closure or
gate changes. Three untargeted checks remain visible. Child benefit is always
renormalized to completed fertility2.100, retaining the positive-benefit gate.

Plan: six exact full-loop reference smokes (two per free-weight lane) with lead acceptance; forty normalized
central differences (all10coordinates, full and half steps); six population
searches of at most80objectives each; two fresh repeats per selected lane and17
standard plots. Six-hour controller cap, search stops at five hours, repeats
stop at5h50, exports within six hours. Maximum538 objectives,24single-thread
Torch workers,192GiB;1800-second objective cap and unchanged23SS normalization
cap. Median prior objective769s,90th percentile1114s:538/24*769=4.8worker-wall
hours idealized, with six-hour hard cap rather than guaranteed completion.
No Mac model work. Checkpoint every case/at most5minutes;30min stale diagnosis.

Lanes: original primary weights;10x and100x earlyfertility weights with all10
coordinates free; fixed continuation-scale profiles at half,double,quadruple
reference with other9free, each using10x earlyweight. At least two starting
points (reference plus successful earlyfertility-oriented overnight point).
Compare all under common original primary weights, raw14moment gaps, and
other-moment primary loss excluding earlyfertility. Fixedparameter profiles
are not exact-target constraints. No target swapping or demotion is adopted.

Current Jacobian must report step-size sensitivity, derivative of normalized
child benefit, scaled10-by10 scored-moment SVD and remaining full14rows.
No global identification, unreachable-target or optimum claims from this finite
exercise. Candidate acceptance requires original scientific gates and exact
repeats, with all failures and unrun cases retained.
