# Fixed-beta calibration profiles

Array17403262 is running: task0 fixesannualbeta0.98; task1 fixes0.99. Both passed source/target/checkpoint preflight. Author approved these profiles on September11. Each profile reoptimizes
all eight other structural coordinates against the same12scoredmoments and
separate completed-fertility normalization2.1. No target, weight, original source,
observer, numerical tolerance or economic equation is changed. The unrestricted
selected calibration is initialization and a reference only: it cannot be
selected as a profile result. These profiles do not impose a permanent global
beta ceiling or establish identification.

Each of two independent Slurm tasks uses8CPUs/64GiB, with a3h cap,7800-second
search budget and3000seconds reserved for final verification. At most two rounds
of16coordinate derivatives and12joint proposals; up to3fixed-beta seed attempts
(the two H0±2%alternatives only if the central seed fails the exact known housing
gate), then two exact final repetitions. Maximum60case evaluations/61repetitions
and488stationarysolves perprofile. Prior case times300–490seconds imply roughly
70–120minutes excluding queue, but conservative per-stage guards may stop earlier.
Each two-wave search stage requires4200seconds remaining under the wrapper cap.
No time/case budget is an accuracy or global-optimality guarantee.

The controller pins the original wrapper/objective/source and the exact selected
seed receipts. It verifies every candidate's fixed beta and all other parameter
values. The unchanged raw scorer continues to describe its original nine-free
contract; separate profile tables/metadata explicitly mark eightfree parameters
and beta as an additional fixed restriction. Every target remains in the saved
full fit table, and selected states retain the stable17diagnostic graphs.
Known housing-equilibrium rejections do not stop unrelated valid candidates;
unexpected source, observer, accounting or code errors stop for review.

`run_profile.py` and `test_run_profile.py` implement/test the actual adaptive loop.
`plan_beta_098.json`, `plan_beta_099.json` and `submit.sh` define the independent
cluster tasks. The existing exact scored-model-loop smoke is pinned. A complete
stub-based adaptive-loop smoke tests both rounds, both beta profiles and final
repetition handling; actual fixed-beta seed solves must pass the complete
source/normalization/equilibrium/accounting/observer gates before any search.

Literature context: De Nardi, French and Jones (2010), *Why Do the Elderly Save?
The Role of Medical Expenses*, JPE118(1), sectionVIII.C,p69, reports beta0.99 in
the endogenous-medical-spending extension and0.97 in its benchmark. This provides
an example of0.99, not a universal upper bound or an identification argument for
this housing/fertility model. Primary fulltext:
https://users.nber.org/~denardim/research/De_Nardi_French_Jones_JPE_2010.pdf
