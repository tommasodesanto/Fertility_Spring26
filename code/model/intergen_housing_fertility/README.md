# Intergenerational Housing Fertility Model

This package is the new one-market quantitative implementation for the
intergenerational housing allocation and fertility project. It starts from the
active workhorse lifecycle code under `code/model/dt_cp_model/`, but keeps that
code path unchanged.

For the simplest current-model run, open and run:

```bash
code/model/run_intergen_model.py
```

That file lists the current source theta, target set, grid, and output folder
at the top.

Run from `code/model` when using the package CLI directly:

```bash
.venv/bin/python -m intergen_housing_fertility.cli smoke --quiet
.venv/bin/python -m intergen_housing_fertility.cli solve --max-iter-eq 20 --quiet
.venv/bin/python -m intergen_housing_fertility.cli diagnostics --fixed-prices --outdir ../../output/model/intergen_housing_fertility_smoke_fixed --quiet
.venv/bin/python tools/build_intergen_mechanics_packet.py --max-iter-eq 10
```

The current pass uses one aggregate housing-services market, a 4-year decision
period, one dependent-child stage, a persistent Markov income state, lifecycle
income by age, owner housing rungs, continuous renter housing, down-payment
constraints, transaction/sale wedges, and a payment-to-income screen for new
owner choices. It is not calibrated.

`IMPLEMENTATION_STATUS.md` is the live implementation record. Any future
simplification or deferred object should be added there in the same edit as the
code change.

`tools/build_intergen_mechanics_packet.py` is a non-production inspection
driver for the June 2026 one-market strand. It reads a saved theta, re-solves
the active intergen model, and writes standard diagnostics, target-fit tables,
room-bin/rung shares, owner-entry thresholds, and optional policy proof-of-
concept cases under `output/model/intergen_mechanics_packet_YYYYMMDD/`.
