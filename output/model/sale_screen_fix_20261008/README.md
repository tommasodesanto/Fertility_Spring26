# Sale-to-rent screen fix (Oct 8 2026, author decision B)

**Production root (pin to this):** `tmp/sale_screen_fix_20261008/root` = byte copy of `tmp/rebate_entry_impl_20261006/root`
(the engine the 14.402 base imports; untracked) with one rule change. `tmp/rebate_entry_impl_20261006/root` is left untouched.
The tracked `code/model/experiments/birth_count_choice` engine is the older Oct 4 version (no net-estate floor, no rebate
closure changes) and is NOT the production engine; `tracked_oct4_engine.patch` shows the same edit against it, unapplied.

**Rule change.** Timing is interest -> income -> housing transaction for purchases, resizes and now sales. The pre-income
sale-to-rent screen R b + (1-psi) p H >= 0 (tenure kernels, active whenever a renter credit limit d_bar is set) is retired; a
sale into renting is feasible iff the renter budget after R b + S + y admits b' >= the renter floor (-d_bar, death floor),
c > 0 and the parent room requirement. Renter floors, purchase timing, interest, preferences, menus, grids, mortgage rules
unchanged. `P.legacy_pre_income_sale_screen = True` restores the old rule bit for bit. Patch: `production_root.patch`
(household.py one line + credit.py docstring). Hashes of every model file in the root: `production_root_model_sha256.txt`.

**Verification at the saved 14.402 parameters** (`run_v2.py`/`common.py` = copies of the 14.402 case driver pointed at the
new root; `compare_baseline.py` -> `baseline_verification.md`):
- Legacy switch on: all 96 solution arrays bitwise identical to the saved 14.402 packet; same for the residuals and 14-row fit.
- Fix, saved price 0.77941391535061: renewal residual 1.7e-5 (tolerance 1e-6), housing residual 0, rebate re-balanced
  (T 0.192173 -> 0.192147).
- Fix, price re-solved at the same parameters (`cases_fix/ge_baseline`, fixed_h0 closure): price 0.7794387507545272 (+0.003%),
  population scale 1.000177, renewal residual -1.2e-9. Loss 14.4024 -> 14.5717. Every row moves by <= 0.001 except
  recent_parent_ownership 0.12666 -> 0.12443 (target 0.12761), loss 0.024 -> 0.274 (weight 27,056).
- LTV 95 (stationary, fixed price, rebate held, `../origination_only_mortgage_20261008/tables_fix.md`): completed fertility
  -1.06% (was -1.35%), impact +0.41% (was +0.42%). Origination-only regime under the fix: impact -0.06%, completed -1.38%.

Old results produced with the screen are legacy-screen results.
