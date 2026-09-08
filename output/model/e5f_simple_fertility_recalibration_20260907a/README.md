# Bounded recalibration of simple fertility nests

**Latest recalibration check: improved candidate; final verification submitted.**
Job17145615 stopped after1h56m37s:36 valid full histories (2 exact anchor smokes,
23 coordinate,11 joint) and1 rejected joint candidate. Joint005 failed the
unchanged market gate3.588e-4>2e-4. Stage-stop correctly blocked final repeats;
no source, specification, bound, target, weight or numerical gate was relaxed.

Best valid candidate joint012 loss26.249682727266702, versus nested start36.3717
and sequential control30.4085 (13.6766% lower). It passes all recorded market,
measurement, mass, population and terminal budget/value checks. First-birth
housing response0.457034 vs target0.720246; ownership0.528988 vs0.575472.
First-child jump0.464931 below upper0.5 (7.0% of bound span remaining); theta1
retains generic near-lower-bound flag. Improvement is provisional, not certified
by final repeats or evidence of a global optimum/policy validity.

Lead independently checked108 collected summary/fit/parameter hashes and
recomputed all432 fit-row losses; rejected case has no complete fit. Original
remote collector checked all11 completed joint histories including checkpoint
hashes and updated cross-stage best. After reviewing the distinct inadmissible
proposal, the2 originally budgeted exact repeats of joint012 were submitted as
array17152974,1core/24GiB each,100min cap. This is verification after review,
not retrying the failed candidate or restarting search. Total attempted
histories remain<=39. Immutable repeatplanSHA
8ceee024ffba823327496d128d383fe7bed838cc7bd966e5ec6a7606a0b68079.
No monitor or policies. Full provisional fits/parameter tables and artifacts:
`output/model/e5f_simple_fertility_recalibration_20260907a/PROVISIONAL_RESULTS.md`.


Author authorized recalibration after the successful retained-parameter comparison.
Scientific model stays frozen at bundle4199e948c5f3625c4a2af106623344ddd8f0b032262f26a8d3973223f5bd63c8.
All12 targets, weights,11 estimated coordinates and original bounds are unchanged.
Housing kappa=.005, supply elasticity=.63 and first-child jump upper=.5.
The old-state fertility level remains separately normalized to2.1 at each case.

The existing direct-search planner controls proposals. Two exact-reference
full-history smoke repeats use the same process/collection loop as the search.
Only after both pass:23 cases (anchor and plus/minus.005 normalized-coordinate
steps for every estimated parameter); up to12 joint candidates combining
improving coordinate directions and child-space floor/jump reallocations;
then2 exact repetitions of the best across all stages. Maximum39 histories.
This is a bounded local recalibration attempt, not proof of a global optimum.

Single Torch job:23 CPUs,322GiB (observed peak about9.4GiB per new-model case),
8-hour wall cap,100-minute per-case cap. Four parallel waves at the observed
38-minute history time suggest3–4hours including variation/overhead; worst
case4waves at100minutes plus collection remains within8hours. No new stage
starts without enough remaining time for its full case cap and5min overhead.
Any failed case prevents subsequent stages; already-running cases finish within
their caps and their completed evidence is retained. No automatic retry.

One-minute controller/case heartbeats and per-case summaries are written.
`search/latest_completed_case.json`, `search/live_best.json`, and staged
collector `search/best_so_far.json` remain readable. Missing case heartbeat
for30minutes stops the run. Final exact-repeat receipt and full comparison
fits/parameter tables are saved next to the final receipt. Comparison uses the
verified sequential exhaustive-saving control at retained parameters, loss30.4085;
new model starting loss36.3717. Sequential interpolation/storage still differs,
as disclosed in the specification. No policies, plots or monitoring automation.

Controller:code/model/tools/run_e5f_simple_fertility_search.py.
Contract and frozen plans here contain source/input hashes, explicit limits and
reproduction commands. The complete final comparison must be collected and
reviewed before promoting any result. Production code remains unchanged.
