# Superseded local fixed-price frontend snapshots

This directory is history only. It holds byte-for-byte copies of the two predeployment files supplied in `tmp/production_deployment_20261003/predeployment_frontend/` on October 3, 2026. The fixed-price `run_model.py` was superseded by the stationary-GE production runner at `code/model/run_model.py`. The saved storage-helper snapshot is retained here for historical reproduction. The current `code/model/tools/model_run_io.py` has since changed into a compatibility loader that defaults to the production cache; active plotters load production storage directly.

[`sha256_receipt.json`](sha256_receipt.json) records each archived file's SHA-256 and confirms it matched the predeployment snapshot at copy time. Do not use these archived files as production entry points.
