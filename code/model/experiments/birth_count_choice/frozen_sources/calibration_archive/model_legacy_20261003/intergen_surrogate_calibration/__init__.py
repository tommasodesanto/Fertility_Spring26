"""Surrogate-assisted calibration prototype for the one-market intergenerational
housing-fertility model.

This package trains a multi-output Gaussian-process emulator on existing model
evaluation logs (theta -> moment vector), then uses it for (a) smooth surrogate
identification diagnostics that mirror the finite-difference Jacobian audit, and
(b) Bayesian-optimization proposals for the next batch of model evaluations.

It is a diagnostic / methods prototype, not a production SMM calibration. See
README.md in this folder.
"""
