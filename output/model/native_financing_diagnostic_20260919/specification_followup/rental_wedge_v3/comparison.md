# Rental-wedge v3 collection

The corrected source-packaging run completed one control case and then stopped on the first positive-slope smoke case. Smoke job `18080236` failed after 1:25; production job `18080237` was cancelled by the dependency. No retry, new job, or cancellation was initiated by collection.

The completed `cap6zero` control recorded one household solve (`28.7358` seconds; case wall `52.5524` seconds). Its control reproduction passed at (10^{-10}), all reported control policy differences were zero, the saving audit passed with maximum value gain `4.44e-16` against tolerance (10^{-7}), maximum saving gap `8.88e-16`, budget excess mass zero, maximum occupied budget excess `7.11e-15`, population mass `1.000000000000023`, cohort initial mass 1, and 17 standard plots.

The next case, `cap10s02`, failed before solving with `NotImplementedError: Rental wedge requires the golden-section renter block`. Its failure receipt is retained. Thus no positive-slope result or strict gate is reported. The new pre-case verification passed checkpoint deserialization and recorded 537 source files, with 163 hash checks before and after and exact equality. No arrays were downloaded.

See [collection_receipt.json](/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26/output/model/native_financing_diagnostic_20260919/specification_followup/rental_wedge_v3/collection_receipt.json) and the retained `results/failure_case_cap10s02_1789924502.json`.
