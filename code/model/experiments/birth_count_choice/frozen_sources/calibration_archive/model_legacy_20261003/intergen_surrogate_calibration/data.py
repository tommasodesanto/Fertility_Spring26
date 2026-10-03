"""Load and curate existing model-evaluation logs for surrogate training.

The single source of truth for parameter bounds, target sets, weights, and the
loss function is the model package itself
(``intergen_housing_fertility.calibration`` and ``...local_panel``). This module
imports those objects rather than re-encoding them, so the surrogate pipeline
cannot silently drift from the real calibration objective.

Key facts established by reconnaissance (2026-06-22):

* Tier-1 evaluation data (full 13-parameter theta + raw moment dict per record):
    - hown5 refinement, 7,901 records, logs all 13 ``candidate_no_timing_v0``
      moments.
    - replacement-wave (3 variants), 7,200 records, logs the replacement moment
      set (old nonhousing wealth, owner/renter mean rooms and >=6-room shares).
* Records store the *moment vector*, not just the loss, so the emulator predicts
  moments and any target-set loss is recomputed downstream with the model's own
  ``diagnostic_loss``.
* ``beta`` is stored in model-period units (beta_annual ** period_years); the
  search box is given in annual units, so the beta bound is transformed here.
"""

from __future__ import annotations

import glob
import json
import os
from dataclasses import dataclass

import numpy as np

from intergen_housing_fertility import calibration as C
from intergen_housing_fertility import local_panel as LP

REPO_ROOT = "/Users/tommasodesanto/Desktop/Projects/Fertility/Fertility_Spring26"

# Parameter order matches the finite-difference Jacobian audit
# (output/model/intergen_sensitivity_jacobian_20260618/audit_summary.json) so the
# surrogate Jacobian is directly comparable column-for-column.
THETA_PARAMS = [
    "beta", "alpha_cons", "b_entry_fixed", "c_bar_0", "c_bar_n",
    "h_bar_0", "h_bar_jump", "h_bar_n", "psi_child", "kappa_fert",
    "chi", "theta0", "theta_n",
]

# Tier-1 data sources (directories of task_*/cases.jsonl).
TIER1_GLOBS = [
    os.path.join(
        REPO_ROOT,
        "output/model/cluster_pulls/"
        "results_intergen_housing_fertility_intergen_fast_globalde_"
        "hown5_refinement_v1_2h_20260617/task_*/cases.jsonl",
    ),
    os.path.join(
        REPO_ROOT,
        "output/model/cluster_pulls/intergen_replacement_cluster_wave_20260618/"
        "results_intergen_housing_fertility_intergen_replacement_*/task_*/cases.jsonl",
    ),
]


def model_bounds() -> tuple[np.ndarray, np.ndarray]:
    """Lower/upper bounds per THETA_PARAMS, in MODEL units.

    GLOBAL_DE_BOUNDS gives ``beta_annual`` in annual units; the solver consumes
    ``beta = beta_annual ** period_years``. We transform that one dimension and
    pass the rest through unchanged.
    """
    period_years = float(getattr(C, "PERIOD_YEARS", 4.0))
    raw = {name: (lo, hi) for name, lo, hi in LP.GLOBAL_DE_BOUNDS}
    lo = np.zeros(len(THETA_PARAMS))
    hi = np.zeros(len(THETA_PARAMS))
    for i, p in enumerate(THETA_PARAMS):
        if p == "beta":
            blo, bhi = raw["beta_annual"]
            lo[i], hi[i] = blo ** period_years, bhi ** period_years
        else:
            lo[i], hi[i] = raw[p]
    return lo, hi


def target_set(name: str) -> tuple[dict, dict]:
    """Return (targets, weights) dicts for a registered target set."""
    targets, weights = C.TARGET_SETS[name]
    return dict(targets), dict(weights)


def load_records(globs: list[str] | None = None) -> list[dict]:
    """Load all Tier-1 evaluation records with a recoverable 13-param theta."""
    globs = globs or TIER1_GLOBS
    files: list[str] = []
    for g in globs:
        files.extend(sorted(glob.glob(g)))
    records: list[dict] = []
    for fp in files:
        with open(fp) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                theta = rec.get("theta")
                moments = rec.get("moments")
                if not isinstance(theta, dict) or not isinstance(moments, dict):
                    continue
                if not all(p in theta for p in THETA_PARAMS):
                    continue
                records.append(rec)
    return records


def diagnostic_loss_no_residual(moments: dict, targets: dict, weights: dict) -> float:
    """The model's own diagnostic_loss with a clean (zero) market residual.

    Surrogate-predicted points have no solved residual; setting it to 0 means the
    residual penalty never fires, so this returns exactly the weighted SSE part
    of ``calibration.diagnostic_loss``. Records that were actually solved use the
    full loss (see ``record_loss``).
    """
    m = dict(moments)
    m.setdefault("market_residual", 0.0)
    return float(C.diagnostic_loss(m, targets=targets, weights=weights))


def record_loss(rec: dict, targets: dict, weights: dict) -> float:
    """Full diagnostic loss for a solved record (uses its real market residual)."""
    return float(C.diagnostic_loss(rec["moments"], targets=targets, weights=weights))


@dataclass
class Dataset:
    target_set_name: str
    moment_names: list[str]
    targets: dict
    weights: dict
    theta_params: list[str]
    lo: np.ndarray
    hi: np.ndarray
    X: np.ndarray            # (N, 13) raw theta
    Xs: np.ndarray           # (N, 13) scaled to [0, 1] by bounds
    M: dict                  # moment_name -> (N,) array
    loss: np.ndarray         # (N,) full diagnostic loss under this target set
    source: np.ndarray       # (N,) source/target-set label per record

    def scale(self, theta_raw: np.ndarray) -> np.ndarray:
        return (np.asarray(theta_raw, float) - self.lo) / (self.hi - self.lo)

    def unscale(self, x_scaled: np.ndarray) -> np.ndarray:
        return self.lo + np.asarray(x_scaled, float) * (self.hi - self.lo)


def build_dataset(target_set_name: str, records: list[dict] | None = None,
                  require_finite_loss: bool = True,
                  max_residual: float = 5e-3) -> Dataset:
    """Assemble a training dataset for one target set.

    Keeps only records whose moment dict contains every target moment (finite)
    AND that cleared the housing market (``market_residual <= max_residual``).
    The convergence filter matters: ``diagnostic_loss`` adds a flat +100 penalty
    to non-converged solves, a discontinuity the moment emulator cannot
    represent, and those solves' moments are not trustworthy. After filtering,
    the training loss equals pure weighted SSE, consistent with the emulated
    loss. Only ~1% of records are dropped.
    """
    if records is None:
        records = load_records()
    targets, weights = target_set(target_set_name)
    moment_names = list(targets.keys())
    lo, hi = model_bounds()

    def _as_float(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return np.nan

    rows_X, rows_M, rows_loss, rows_src = [], [], [], []
    n_dropped_resid = 0
    for rec in records:
        moments = rec["moments"]
        resid = _as_float(moments.get("market_residual", np.nan))
        if not np.isfinite(resid) or resid > max_residual:
            n_dropped_resid += 1
            continue
        vals = [_as_float(moments.get(name, np.nan)) for name in moment_names]
        if not all(np.isfinite(v) for v in vals):
            continue
        loss = record_loss(rec, targets, weights)
        if require_finite_loss and not np.isfinite(loss):
            continue
        rows_X.append([float(rec["theta"][p]) for p in THETA_PARAMS])
        rows_M.append([float(v) for v in vals])
        rows_loss.append(loss)
        rows_src.append(str(rec.get("origin", rec.get("label", "?"))))

    X = np.asarray(rows_X, dtype=float)
    Mmat = np.asarray(rows_M, dtype=float)
    loss = np.asarray(rows_loss, dtype=float)
    Xs = (X - lo) / (hi - lo)
    M = {name: Mmat[:, j] for j, name in enumerate(moment_names)}
    return Dataset(
        target_set_name=target_set_name, moment_names=moment_names,
        targets=targets, weights=weights, theta_params=THETA_PARAMS,
        lo=lo, hi=hi, X=X, Xs=Xs, M=M, loss=loss,
        source=np.asarray(rows_src),
    )


def curate_indices(Xs: np.ndarray, loss: np.ndarray, n_keep: int,
                   frac_best: float = 0.5, seed: int = 0) -> np.ndarray:
    """Pick a GP-tractable training subset.

    Combines (a) the lowest-loss records -- the basin we care about for
    optimization and local identification -- with (b) a farthest-point
    space-filling sample of the remainder for global coverage. Both matter: the
    basin gives local accuracy near good fits, the space-filling part keeps the
    global response surface honest.
    """
    n = Xs.shape[0]
    if n <= n_keep:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    n_best = int(round(frac_best * n_keep))
    n_fill = n_keep - n_best

    order = np.argsort(loss)
    best_idx = order[:n_best]
    remaining = order[n_best:]

    # Farthest-point sampling over the remaining pool (greedy max-min).
    pool = remaining.copy()
    # Seed the fill set with a random point, then add farthest points greedily.
    chosen = [int(pool[rng.integers(len(pool))])]
    # Precompute distances incrementally.
    min_d = np.full(len(pool), np.inf)
    pool_X = Xs[pool]
    for _ in range(min(n_fill, len(pool)) - 1):
        last = Xs[chosen[-1]]
        d = np.sum((pool_X - last) ** 2, axis=1)
        min_d = np.minimum(min_d, d)
        # Avoid re-picking already chosen.
        for c in chosen:
            hit = np.where(pool == c)[0]
            if hit.size:
                min_d[hit[0]] = -1.0
        nxt = int(pool[np.argmax(min_d)])
        chosen.append(nxt)
    fill_idx = np.array(sorted(set(chosen)), dtype=int)
    keep = np.unique(np.concatenate([best_idx, fill_idx]))
    return keep


if __name__ == "__main__":
    recs = load_records()
    print(f"loaded {len(recs)} tier-1 records")
    ds = build_dataset("candidate_no_timing_v0", recs)
    print(f"candidate_no_timing_v0 usable rows: {ds.X.shape[0]}")
    print(f"moments ({len(ds.moment_names)}): {ds.moment_names}")
    print(f"loss: min={ds.loss.min():.3f} median={np.median(ds.loss):.3f} "
          f"max={ds.loss.max():.3f}")
    lo, hi = ds.lo, ds.hi
    print("bounds (model units):")
    for p, a, b in zip(ds.theta_params, lo, hi):
        print(f"  {p:16s} [{a:.4f}, {b:.4f}]")
    keep = curate_indices(ds.Xs, ds.loss, n_keep=1200)
    print(f"curated subsample: {keep.size} "
          f"(best-loss in subset: {ds.loss[keep].min():.3f})")
