# Superseded local fixed-price frontend snapshots

This directory is history only. It holds byte-for-byte copies of the two predeployment files supplied in `tmp/production_deployment_20261003/predeployment_frontend/` on October 3, 2026. The fixed-price `run_model.py` was superseded by the stationary-GE production runner at `code/model/run_model.py`. The saved storage-helper snapshot is retained here for historical reproduction. The current `code/model/tools/model_run_io.py` has since changed into a compatibility loader that defaults to the production cache; active plotters load production storage directly.

[`sha256_receipt.json`](sha256_receipt.json) records each archived file's SHA-256 and confirms it matched the predeployment snapshot at copy time. Do not use these archived files as production entry points.

The pre-map [code/model README](README_code_model_before_map.md) preserves the complete former navigation and historical narrative before its October 3, 2026 replacement. Its SHA-256 is `aea46dc8c37697791d4c7f989640c2d9a9863acdda21c0603aaefe699d382ccd`. It is an archival record, not current workflow guidance; use [`code/model/README.md`](../../code/model/README.md) and [`CALIBRATION_STATUS.md`](../../CALIBRATION_STATUS.md) for navigation and live status.
