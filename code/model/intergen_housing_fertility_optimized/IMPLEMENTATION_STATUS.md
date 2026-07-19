# Optimized Package Implementation Status

Updated: 2026-07-19

This directory is an experimental physical refactor of
`code/model/intergen_housing_fertility/`. Production code and production
runners are unchanged.

## Accepted

- Complete package copy under a separate import namespace.
- Stable test entry point: 71 class-based and 10 top-level tests.
- Exact fixed-price M5 household-policy parity.
- Compiled-scatter distribution parity to machine precision.
- Confirmed audit defects C1--C5 corrected or hardened here.
- Complete 15-row M5 target system owned by the package.
- Exact package-owned M5 configuration.
- Direct one-market bracket/root search with directional expansion.
- Price evaluation cache.
- Accepted-price Bellman reuse for final statistics.
- Two strict optimized M5 repeats with identical price, moments, loss, and
  residual.
- Six-case local fixed-price parameter panel with exact policy/value parity and
  target-moment agreement within `4.09e-14`.

## Rejected and removed

- The first fused age kernel preserved mass and state distributions but was
  slower than the accepted compiled-scatter implementation. Its measurements
  remain in `REFACTOR_REPORT.md`; no dormant runtime branch remains.

## Not authorized for production

- No tool, calibration runner, cluster script, paper table, or diagnostic
  packet should import this package yet.
- No calibration loss from the optimized nearby equilibrium should replace the
  canonical M5 loss. The canonical comparison remains the saved-price mapping.
- Promotion requires extending the local panel to live parameter bounds and
  boundary stresses, plus strict Torch verification under the production CPU
  allocation.
