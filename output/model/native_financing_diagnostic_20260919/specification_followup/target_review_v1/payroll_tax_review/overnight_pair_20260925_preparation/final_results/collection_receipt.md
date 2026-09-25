# Paired overnight calibration — collection receipt (2026-09-25)

This is a bounded collection of the existing Torch run. No model job was added,
retried, or changed during collection. The two selected calibrations use the
same target, bounds, seed bank, native source and gates, with PAYGO rates 17.9%
and 8.751017424959717%. They are separately recalibrated selections, so their
loss difference is not a fixed-parameter tax effect.

| Arm | Selected loss | Objectives scored/planned | Native solves completed/started | Selected case |
| --- | ---: | ---: | ---: | --- |
| 17.9% | 291.73166273362773 | 364/364 | 2024/2024 | `worker_04/point_07` |
| 8.751% | 280.41079921120206 | 364/364 | 2217/2217 | `worker_05/point_04` |

Both arms have zero incomplete or unrun objectives. Across the pair, all 728
planned normalized objectives scored and all 4,241 native solves completed.
The two exact selected repetitions per arm passed at zero tolerance. The
original selected and both repetition checkpoints match in native price, value
and distribution arrays; evaluated policy value and distribution; stationary
pre-distribution; fertility scale; loss; and byte-identical complete target and
parameter CSVs. The detailed array shapes/hashes and every repeat comparison
are in `collection_scientific_repeat_audit.json`. The completed exporter also
compared `chain.extract_moments` at zero tolerance; it did not save separate
per-moment hashes. The inherited reference checkpoint SHA256 is
`83a28e46b36e2fbe30338d366611f3ec209f0c5a68309ee4ee9fa8523b66adee`,
distinct from the new selected-case SHA256s `de1da882335328f9c0a6ace673c85e10522c8e68c09555f6a7da1717d93c33d7`
and `d5ef71bdaf9960273035c722a2428a55f14bab160e0596881c8a981e67b8ead1`.

The four matched common-seed smoke losses were, respectively, 460.832391 and
412.742514 at 17.9%, and 422.547487 and 381.048574 at 8.751%, for the
September 23 B-floor and September 24 commute structural points. Fertility
scale was normalized separately in each arm; full 13-row fits and 25-row
parameter tables are retained in `matched_smokes/`.

Slurm reports every smoke, worker, repetition and export task completed with
exit code 0. Smokes ran 00:48:56–01:18:01 EDT; workers 01:18:05–06:12:39;
repetitions 06:12:45–06:25:16; exports 06:25:20–06:25:52. All work finished
before the authorized eight-hour 08:49:12 deadline; worker walltime and the
extended deadline did not truncate proposal coverage.

The saved fiscal accounts imply pension per gross working earnings of
51.15903474030076% and 25.010815891415753% in the respective selected arms.
All 728 scored receipts match the objective/source fingerprints and pass their
retained market, fiscal, birth-entry, purchase and value gates. Five high-tax
proposals have trace positive `budget_excess_mass`, at most 1.5773e-23; they
pass the actual frozen predicate in `run_e5f_matched_pf_smoke.py:175-183`,
which counts spending gaps above 1e-9 and raises when mass exceeds 2e-10.
Both selected cases have exactly zero budget-excess mass. The full gate audit
is `collection_gate_audit.json`.

Each selected PDF has 21 pages: four report/table pages and the unchanged 17
standard figure pages. All pages rendered on Torch; no blank pages, and every
target and parameter name is present. The original PDFs and hashes are in each
arm directory; all-page contact sheets and page-1–4 previews are in `pdf_qa/`.
At full resolution page 4 has visible separation between the longest parameter
names and the estimate column. Some inherited figure legends remain crowded.
The broad 1%-of-range `near_bound` rule flags both fertility κ estimates, but
neither equals its lower bound.

A separate report-only attempt to rebuild native runtime before hashing
checkpoints stopped at a purchase-income source patch anchor error. It made no
model solve and produced no scientific comparison. Direct read-only hashing of
the saved checkpoints then passed; a subsequent Torch check found all 644
frozen source files still equal to their manifest hashes. The numerical run
and its selected PDFs were already complete before this collection attempt.
