# Bounded recalibration of simple fertility nests

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
