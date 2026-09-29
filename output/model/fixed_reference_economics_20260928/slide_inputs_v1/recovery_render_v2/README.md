# Housing-price response figure: corrected footer spacing

Reference: **2007 stationary reference — block0506, September 28 verified export**.

Torch **18817838** completed with exit 0:0 in seven seconds. This is a formatting-only successor to [v1](../recovery_render_v1/README.md): three footer coordinates and the plot-area bottom margin change. Validation, numerical inputs, wording and tables are unchanged. The lead reviewed the four-line diff and visually inspected the actual PNG; the footer is legible. No model imports or lifecycle solves were run.

The final [figure PDF](actual_output/price_response.pdf), [PNG](actual_output/price_response.png), [four-row table CSV](actual_output/local_elasticities.csv), [LaTeX table](actual_output/local_elasticities.tex) and [input manifest](actual_output/manifest.json) are ready. All existing 17-plot model diagnostic packets remain separate; this figure is supplemental. It compares prescribed price and mapped-rent changes, not GE or transition paths.

The frozen renderer is `/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v2/source_v2/build_three_price_inputs.py`, SHA-256 `6b6c4250326ccc590df969833bf874152f3c1eddd23f34586107870464e03a57`. Inputs are the passed `comparison.csv`, `elasticities.csv` and `completed.json` from `/scratch/td2248/projects/fixed_reference_elasticity_recovery_20260929/solve_v1/`. The renderer verifies the completion flags, exact three-price scope, input hashes and displayed elasticities against source levels before drawing.

Remote outputs: `/scratch/td2248/projects/fixed_reference_economic_figures_20260929/recovery_render_v2/actual_output/recovery_18815133_v2/`. [Remote hashes](remote_sha256.txt) pin the source and outputs; [render_actual.sbatch](render_actual.sbatch) contains the exact zero-solve command. Do not resubmit into the completed output directory. Prior schema/render tests are in [v1 validation](../recovery_render_v1/torch_validation/).
