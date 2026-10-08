# Origination-only mortgage vs revolving collateral (Oct 8 2026, New York)

Question (Tommaso, Oct 8): does mortgage origination matter for fertility once incumbents can no longer borrow against the house?

Base 14.402 (Mac round-3 chain 0, promoted Oct 7), fixed price 0.77941391535061, rent per room = user cost x price, rebate T held at
0.19217. Engine: experiment copy `tmp/origination_only_20261008/root` = byte copy of `tmp/rental_menu_precaution_20261007/root` plus one
default-off switch `P.stayer_no_new_borrowing` (`kernels.py` full_owner_block_kernel arg `due_no_new_borrowing`; `household.py` call site).

| regime | stayer floor | purchase / resize | renters |
|---|---|---|---|
| revolving (current) | b' >= max(min(b, -phi p H), death floor) | R b + y >= (1-phi) p H, b' >= -phi p H | b' >= 0 |
| origination-only | b' >= max(min(b, 0), death floor) | unchanged (new mortgage on every purchase, old debt settled from sale) | b' >= 0 |

No mandatory amortization (debt can stay interest-only); one interest rate on both signs of b; death / net-estate floors unchanged.

- `run_cells.py R80,R95,O80,O95` solves (about 9 s each). `check_repro.py`: R80 bitwise identical to the saved 14.402 packet
  (`../jmp_draft_deck_20261005/refit_best_14p40_mac_r3_chain0/cases_v2/baseline`) and to rental-menu cell A6 on V, g, birth policies,
  pre-distribution, hR, b', c, tenure probabilities; R95 bitwise identical to rental-menu A6_LTV95.
- `analyze.py` -> `tables.md` (definitions in its docstring), `outcomes.json`.

Headline: LTV 80 -> 95 moves births +0.42% on impact / -1.35% completed under the revolving rule and +0.02% / -1.66% under
origination-only. Switching regime at 80% changes births -0.3% but cuts ownership 69 -> 58% (65-75: 84 -> 56%); the revolving rule
is used mainly as old-age equity release, not by would-be parents. Liquid resources at the fertility margin halve (0.38 -> 0.20) with
no birth response.

Not representable with one signed asset: positive liquid savings alongside mortgage debt; contractual amortization by loan age.

## Sale-to-rent screen (Oct 8, 15:45 ET)
- `sale_screen_count.py` -> `sale_screen_count.md`: the adopted screen R b + (1-psi) p H >= 0 (kernels.py:341/419) omits current income,
  unlike purchase/resize. At LTV 95 it bars 60.5% of childless owners 22-33 from selling into renting, nearly all budget-feasible with income.
- DIAGNOSTIC switch `P.sale_screen_includes_income` (copied engine only; R b + S + y >= 0). Cells R80I/R95I; `compare_screen.py` ->
  `sale_screen_income_diagnostic.md`. R80 rerun after the edit still bitwise identical. LTV 80->95 completed fertility -1.35% -> -1.06%;
  impact +0.42% -> +0.41%; screen is inert at LTV 80 (+0.00%).
