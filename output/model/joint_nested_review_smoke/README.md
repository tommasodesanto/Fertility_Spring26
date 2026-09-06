**Latest author steering, September 6 at 22:30 UTC: review after two hours,
before the full overnight search.** Full search remains stopped until that
discussion. Smoke `17074777` completed all four certified histories, but the
policy finalizer rejected a missing original inherited population: its saved
population had already undergone feasibility projection of `2.32469e-15`.
Dependent long job `17075663` was automatically cancelled with zero runtime.
This was a checkpoint handoff failure, not evidence against equilibrium or
an unfit final calibration.

The isolated repair adds an owned copy of the input population to joint-mode
`PeriodEvaluation`, before any price-specific projection. Policy branches
start from that copy. The finalizer exactly replays the original feasibility
gate and checks the fitted population and projected mass before proceeding;
it does not replace the old zero-projection check with a looser tolerance.
Economic decisions, targets and numerical projection behavior are unchanged.
Eighteen local tests pass, including a deliberately nonzero projection that
checks raw and gated populations remain distinct and do not alias the input.

Bounded verification job `17076426` is running on two Torch CPUs, 64GB, with
an 80-minute scheduler cap. It reruns the default-off reference, two exact
histories, two all-coordinate probes and four two-date policy loops. No
long-search dependency is attached. New exclusive snapshot:
`/scratch/td2248/projects/Fertility_Spring26_joint_nested_full_20260906c`.
Scientific bundle `2c4adcdec5b85cff39b4c5d1466e6224db1dddf13052fe66c8a6be86e9a32951`;
contract SHA `416477db8ce66d22a6017aee24fd8a2a2d974c3fcf87bbed6bfe8f6f673c48ab`.
The local new contract is `review_smoke_contract.json` in
`output/model/e5f_joint_nested_full_20260906a/`; the earlier snapshot and
`frozen_contract.json` remain failure evidence. The finite monitor now follows
this verification and prepares the discussion packet by approximately
September 7 at 00:30 UTC (20:30 New York). It must not automatically launch or
release a full search; the author requested that discussion first.
