# Intergenerational Housing Fertility Model

This folder will contain the quantitative implementation of the no-geography
intergenerational housing mismatch and fertility model.

The active calibrated center/periphery model remains under `code/model/dt_cp_model/`.
This implementation starts from its lifecycle, income, tenure, and distribution
patterns, but it does not modify that code path.

Before changing this folder, update or check `IMPLEMENTATION_STATUS.md`. That
file is the live record of intended model objects, simplifications, and missing
pieces. No simplification should be left implicit in code.

Run from `code/model`:

```bash
.venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
.venv/bin/python -m intergen_housing_fertility.cli solve --mode smoke --max-iter-eq 20 --quiet
.venv/bin/python -m intergen_housing_fertility.cli diagnostics --mode smoke --fixed-prices --outdir ../../output/model/intergen_housing_fertility_smoke_fixed --quiet
```

The `solve` command reports whether owner-size market clearing converged. With
the first-pass deterministic size choice it may return the best price iterate
rather than an exact clearing vector; this is recorded in
`IMPLEMENTATION_STATUS.md`.
