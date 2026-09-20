# Rental-wedge v2 collection

The smoke job `18079902` failed after 16 seconds while loading the first case `cap6zero`; production job `18079903` was cancelled by the submitted dependency. No retry, cancellation, or new job was initiated by collection.

The immediate failure is `ModuleNotFoundError: No module named \'run_e5f_perfect_foresight_person_demography\'` in `run_e5f_isolated_rental_wedge.py` while unpickling the checkpoint packet, before the first household solve. No `solve_receipt` exists, so the saved actual household-solve count is 0.

The pre-case source receipt passed: 163 hashes before and after, zero bad hashes, exact before/after equality, and 43 recorded import origins. Saving-gain, reproduction, budget/value/cohort/population/entry gates, plot manifests, and production source-after checks are unavailable because execution stopped before solving.

Collected artifacts: [collection_receipt.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/rental_wedge_v2/collection_receipt.json), `results/failed_receipt.json`, `results/failure_case_cap6zero_1789924067.json`, `results/verify_smoke_before.json`, and the immediate stderr tail.
