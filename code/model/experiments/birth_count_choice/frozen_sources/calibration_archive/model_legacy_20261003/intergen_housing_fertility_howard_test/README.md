# Intergenerational Housing Fertility Model

## Howard Test Copy

This folder is an experimental copy for Markov-income Howard speed testing.
The production package remains `code/model/intergen_housing_fertility/`.
Do not treat this folder as a calibrated or production model unless the Howard
path is explicitly promoted back to the active package.

This package is the new one-market quantitative implementation for the
intergenerational housing allocation and fertility project. It starts from the
active workhorse lifecycle code under `code/model/dt_cp_model/`, but keeps that
code path unchanged.

Run from `code/model`:

```bash
.venv/bin/python -m intergen_housing_fertility_howard_test.cli smoke --quiet
.venv/bin/python -m intergen_housing_fertility_howard_test.cli solve --max-iter-eq 20 --quiet
.venv/bin/python -m intergen_housing_fertility_howard_test.howard_speed_audit --variant howard --refine-mode full --quiet
```

The current pass uses one aggregate housing-services market, a 4-year decision
period, one dependent-child stage, a persistent Markov income state, lifecycle
income by age, owner housing rungs, continuous renter housing, down-payment
constraints, transaction/sale wedges, and a payment-to-income screen for new
owner choices. It is not calibrated.

`IMPLEMENTATION_STATUS.md` is the live implementation record. Any future
simplification or deferred object should be added there in the same edit as the
code change.
