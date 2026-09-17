# Task for Muse (OpenCode): four default-off model switches for steady-state tests

Repository root is the current directory. Read `CLAUDE.md` first, then
`code/model/sandbox/README.md`. Python: `code/model/.venv/bin/python`. Run
tests with `code/model/.venv/bin/python -m pytest <file>`. Do the four parts
in order; after each part run the full package test suite
(`code/model/intergen_eqscale_seq_optimized/tests/`) and stop at the first
part whose tests do not pass, reporting what remains. Time limit: 5 hours.

Hard rules. Every switch is OFF by default and, when off, must leave every
existing array bitwise identical (`np.array_equal`, not `allclose`); write
that test first for each part. Do not change any default value, target,
calibration profile, driver under `code/model/tools/` (except where a part
says so), or any file under `output/`. Do not use git. Do not launch cluster
jobs. Do not run the sandbox. Never use the word "parity" in prose or comments
you write (existing identifiers may stay). Where a design choice is not fixed
below, choose the simplest correct one and state it in the receipt.

Notation used below: \(b\) liquid position at the start of the period (negative
is mortgage debt), \(b'\) the chosen position, \(P\) house price, \(h\) owner
size (rooms), \(h^R\) rented rooms, \(\phi=0.80\) the financed share, \(m\)
children at home, \(y^d\) after-tax earnings, \(r\) rent per room.

Code map (from the September 15 verification; line numbers are approximate,
grep for the identifiers): down-payment test on trades only
`kernels.py:640-662` (`tenure_choice_kernel`, `bar = bg_b+sp-hc; dpc = dpn-sp;
if bg_b < dpc or bar < bmn: NEG_INF`); collateral floor re-applied to stayers
`kernels.py:979-985` (`bf = bmo[...]`, `current_unsecured`, `total_floor`);
`dp_arr`/`bmo` built at `solver.py:2917-2922`; unsecured line `parameters.py:146`
`lambda_d`, taper `parameters.py:147-148, 643-647`, `build_debt_caps:677`;
owner services `kernels.py:957-970` (`ht_c = owner_service_premium*(hsv -
owner_h_bar_scale*hbc)`); renter rent in the budget: renter kernel in
`kernels.py` (search `hR`, `rent`) and `P.hR_max` `parameters.py:191`; user
cost `parameters.py:255` and `run_e5f_perfect_foresight_transition.py:233-256`;
income into the budget `solver.py:226-233` (`income_at_state`) and
`parameters.py:724-751` (`set_income_given_w_and_pension`, `resolve_pension_value`);
bequest wealth `solver.py:2452-2459` (`Vbq = bequest_utility_vec(b_grid + hv, ...)`,
`hv = p_hat*H_own` gross); bequest flow moment `solver.py:6001-6062`;
property-tax revenue and rebate `solver.py:7356-7371, 5054-5065`.

## Part 1. Child earnings penalty (switch `child_earnings_penalty`)

Economics: a per-period earnings cost of children at home, proportional to
earnings (a time cost). New parameter `child_earnings_penalty`: a list indexed
by children at home \(m=0,1,2,3+\), default `[0,0,0,0]`. When any entry is
nonzero, working-age after-tax earnings become \(y^d(a,z)\,(1-\tau_c(m))\) in
the household budget, for \(a<a_R\) only. Pensions, the payroll-tax base and
the PAYGO balance use unpenalized earnings (state this in the receipt and in a
comment). Implementation: the income array `P.income[i,j]` has no child
index; build a multiplier over \((j, m)\) and apply it where resources are
formed in the Bellman and in every forward/moment path that recomputes
resources (grep every use of `income_at_state` and `P.income`; cover each).
Tests: zeros bitwise identical; with `[0,0.2,0.2,0.2]` the resources at a
state with \(m\ge1\) equal 0.8 times the unpenalized value at working ages and
are unchanged at retirement ages; mass conservation in the forward step.
Spec key: `child_earnings_penalty: [0.0, 0.2, 0.2, 0.2]`.

## Part 2. Mortgage block (switches `mortgage_origination_only`, `mortgage_amortization`)

Two independent booleans/parameters, both default off.

(a) `mortgage_origination_only: false` (default). When true, the collateral
floor \(b'\ge-\phi P h\) is applied only on transactions (purchase or size
change), exactly where the down-payment test already fires
(`kernels.py:640-662`); for a stayer (`to == tn`) the floor is replaced by a
no-cash-out rule: if \(b<0\) then \(b'\ge b\) (debt may not increase); if
\(b\ge0\) then \(b'\ge 0\)... no: if \(b\ge0\) the stayer may still borrow up
to \(-\phi P h\) (a first mortgage on an owned house is an origination);
implement as: stayer floor \(=\min\{b,\,-\phi P h\}\) when \(b<0\), and
\(-\phi P h\) when \(b\ge0\). The underwater-rollover taper logic at
`kernels.py:979-985` must be bypassed when this switch is on (a stayer is never
forced below its own balance), leaving the code path untouched when off.

(b) `mortgage_amortization: 0.0` (default, off). When \(\alpha_m>0\), a
stayer with \(b<0\) must reduce debt by at least the share \(\alpha_m\) per
period: \(b'\ge b\,(1-\alpha_m)\) (note \(b<0\), so this raises the floor
toward zero). Combined with (a) the stayer floor is
\(\max\{\,b(1-\alpha_m),\ \text{floor from (a)}\}\). For a mover the origination
limit applies as today. A thirty-year mortgage at a 2 percent real rate
retires about 11 percent of principal in four years; the spec will use
`mortgage_amortization: 0.11`.

Tests: both off bitwise identical; with (a) on, on a tiny grid, a stayer with
\(b<0\) cannot choose \(b'<b\) and a buyer's feasible set is unchanged; with (b)
on, a stayer with \(b=-1\) cannot choose \(b'<-0.89\); market clearing and
mass conservation unaffected in structure (run the tiny smoke in both modes).

## Part 3. Size-dependent rental wedge (switch `rental_wedge`)

Economics: the landlord's operating cost rises with unit size, so rent per
room is \(r(h^R)=r+w(h^R)\) with \(w(h)=w_0+w_1\max\{0,h-h_k\}\). New
parameters `rental_wedge_intercept` \(w_0\) (default 0), `rental_wedge_slope`
\(w_1\) (default 0), `rental_wedge_knee` \(h_k\) (default 6.0). When both
\(w_0=w_1=0\) the code path is bitwise identical. The wedge is a real cost
paid to the outside landlord: it enters the renter's budget as
\(h^R\,[r+w(h^R)]\) and nowhere else (not the property-tax base, not the
rebate, not the household's wealth); state this. It does not change the user
cost \(r\) itself. Separately, `hR_max` is already a parameter; the spec will
set it to the owner maximum (11.0) so that the cap is replaced by the wedge.
The owner premium `chi` is already a parameter; the spec will set it to 1.0.
Tests: zeros bitwise identical; with \(w_0=0.02, w_1=0.05, h_k=6\) the renter
kernel's cost of 8 rooms equals \(8(r+0.02+0.10)\); the renter's optimal
\(h^R\) is weakly lower than without the wedge at every state.

## Part 4. Estates paid to households (switch `estate_receiver`)

Economics: today a dying household's estate \(w^e=b+P h\) (gross of the
selling cost) enters only the bequest utility and its wealth leaves the
economy. New switch `estate_receiver: none` (default) or `"ages_45_65"`. When
on: (i) the estate is valued net of the selling cost, \(w^e=b+(1-\psi^s)Ph\),
in both the bequest utility and the accounting (a second boolean
`bequest_net_of_selling_cost`, default false, controls the utility side alone
so the two can be separated); (ii) the aggregate estate flow of the period,
\(E=\sum \text{deaths}\times w^e\), is paid as an equal lump-sum transfer to
every household aged 45 to 65 inclusive (calendar ages, converted to age
indices as the package does), \(T_E=E/\mathcal H_{45\text{–}65}\), entering
income the way the property-tax rebate does (`income_at_state`), so it is in
the budget but not in the down-payment test. \(T_E\) is an equilibrium object:
add it to the stationary fixed point as an outer iteration (start at 0, solve,
compute \(E\), update \(T_E\) with damping 0.5, repeat until the relative
change is below 1e-6), not as a new Newton dimension; document where. Tests:
off bitwise identical; on, the estate flow paid equals the estate flow
generated to 1e-10 at the converged point; mass conservation; the tiny smoke
solves. Spec keys: `estate_receiver: ages_45_65`, `bequest_net_of_selling_cost: true`.

## Part 5. Two parameter-only specs (no code)

Write `code/model/sandbox/specs/credit_line_modest_psi_fixed.yaml` with
`lambda_d` equal to three quarters of one quarter of mean annual earnings in
the model's units (Kaplan and Violante 2014; compute the number from the
package's earnings normalization and state it), keeping the existing age
taper; and `code/model/sandbox/specs/earnings_e6b_psi_fixed.yaml` that
replaces the live Floden–Lindé AR(1) by the E6b PSID decomposition: annual
persistence 0.8863 and stationary log-variance 0.3319 for the persistent
component, with the same fixed-effect variance 0.3931 and the same after-tax
scaling applied to both components (find the override keys the calibration
layer uses in `local_panel.py`, `externals.py`, `e6b_profile.py`; convert to
four-year parameters the way `local_panel.py:1063-1071` does; state every
number). Both specs `psi_mode: fixed`.

## Deliverable

A receipt (Outcome / Verification / Artifacts / Unresolved / Reported cost)
listing, per part: files and line ranges changed, the exact spec override
keys with the values to use, the tests added and their results, and any
design choice you made where this file left one open. Also write
`code/model/sandbox/specs/all_switches_psi_fixed.yaml` combining every switch
above with `child_maturation_mode: parent_age` (mu_young 0.05, a_rise 34,
a_full 62), `child_benefit_form: log`, and `psi_mode: fixed`, and one spec per
single switch (`sw_penalty`, `sw_mortgage`, `sw_wedge`, `sw_estate`), each
`psi_fixed`.
