"""Bayesian-optimization proposals and a frontier (Pareto) probe.

Acquisition is Monte-Carlo Expected Improvement on the *emulated loss*, where the
loss distribution is induced by the per-moment GP predictive uncertainty
(``MomentEmulator.predict_loss_mc``). EI is maximized over the parameter box by a
Latin-hypercube pool plus a short local hill-climb, then a diverse batch is
selected greedily by minimum scaled-distance.

The frontier probe uses the moment emulator to ask the question the SOTA index
poses directly: are the conflicting moments (e.g. prime-age ownership, old-age
ownership, housing user-cost share) jointly reachable, or is there a hard
trade-off? It maps the achievable Pareto front of per-moment absolute
deviations over a large emulated sample.
"""

from __future__ import annotations

import numpy as np

from .emulator import MomentEmulator


def latin_hypercube(n: int, d: int, rng: np.random.Generator) -> np.ndarray:
    """LHS sample in the unit cube [0,1]^d."""
    cut = np.linspace(0, 1, n + 1)
    u = rng.uniform(size=(n, d))
    a, b = cut[:n], cut[1:]
    pts = a[:, None] + u * (b - a)[:, None]
    for j in range(d):
        pts[:, j] = pts[rng.permutation(n), j]
    return pts


def expected_improvement(emu: MomentEmulator, Xs: np.ndarray, f_best: float,
                         n_samples: int = 64, seed: int = 0) -> np.ndarray:
    """MC expected improvement (minimization) at each row of Xs."""
    _, samples = emu.predict_loss_mc(Xs, n_samples=n_samples, seed=seed)
    improvement = np.maximum(f_best - samples, 0.0)
    return improvement.mean(axis=1)


def propose_batch(emu: MomentEmulator, f_best: float, q: int = 8,
                  pool_size: int = 6000, n_local: int = 6,
                  n_samples: int = 64, min_dist: float = 0.15,
                  seed: int = 0):
    """Propose q diverse high-EI parameter points (raw theta).

    Returns a list of dicts with raw theta, scaled point, EI, predicted loss, and
    predicted moments.
    """
    ds = emu.ds
    rng = np.random.default_rng(seed)
    d = len(ds.theta_params)

    pool = latin_hypercube(pool_size, d, rng)
    ei = expected_improvement(emu, pool, f_best, n_samples=n_samples, seed=seed)

    # Local hill-climb from the best pool points to sharpen maxima.
    top = np.argsort(ei)[::-1][: max(q * 4, 32)]
    cand = pool[top].copy()
    cand_ei = ei[top].copy()
    step = 0.06
    for _ in range(n_local):
        pert = cand + rng.normal(0, step, size=cand.shape)
        pert = np.clip(pert, 0.0, 1.0)
        pert_ei = expected_improvement(emu, pert, f_best,
                                       n_samples=n_samples, seed=seed + 1)
        better = pert_ei > cand_ei
        cand[better] = pert[better]
        cand_ei[better] = pert_ei[better]
        step *= 0.8

    # Greedy diverse selection by minimum scaled distance.
    order = np.argsort(cand_ei)[::-1]
    chosen: list[int] = []
    for idx in order:
        x = cand[idx]
        if all(np.linalg.norm(x - cand[c]) >= min_dist for c in chosen):
            chosen.append(idx)
        if len(chosen) >= q:
            break
    # If diversity filter is too strict, top up with next-best.
    for idx in order:
        if len(chosen) >= q:
            break
        if idx not in chosen:
            chosen.append(idx)

    means = emu.predict(cand[chosen])
    out = []
    for k, idx in enumerate(chosen):
        xs = cand[idx]
        theta_raw = ds.unscale(xs)
        pred_loss = float(emu.predict_loss(xs[None, :])[0])
        out.append({
            "rank": k,
            "theta": {p: float(v) for p, v in zip(ds.theta_params, theta_raw)},
            "scaled_point": xs.tolist(),
            "ei": float(cand_ei[idx]),
            "pred_loss": pred_loss,
            "pred_moments": {n: float(means[n][k]) for n in ds.moment_names},
        })
    return out


def propose_local(emu: MomentEmulator, train_idx: np.ndarray, q: int = 8,
                  top_k: int = 25, n_per: int = 600, sd: float = 0.04,
                  n_samples: int = 64, min_dist: float = 0.06, seed: int = 4):
    """Trust-region / local-refinement proposals around the best incumbents.

    Naive global EI in 13-D extrapolates the surrogate into sparsely-sampled
    regions where its mean is unreliable -- and where discrete moments (owner
    median rooms) snap to a different rung than the GP's smooth interpolant
    predicts -- so the real loss is far worse than predicted. This routine
    instead perturbs the lowest-loss training points by small Gaussian steps in
    scaled space, keeping proposals inside the data manifold so the surrogate is
    actually trustworthy and small steps do not flip owner rungs. This is the
    realistic use of the emulator: local screening with real-solver confirmation.
    """
    ds = emu.ds
    rng = np.random.default_rng(seed)
    Xtr = ds.Xs[train_idx]
    loss_tr = ds.loss[train_idx]
    f_best = float(loss_tr.min())
    centers = Xtr[np.argsort(loss_tr)[:top_k]]

    cand = []
    for c in centers:
        pert = c[None, :] + rng.normal(0, sd, size=(n_per, len(ds.theta_params)))
        cand.append(np.clip(pert, 0.0, 1.0))
    cand = np.vstack(cand)
    ei = expected_improvement(emu, cand, f_best, n_samples=n_samples, seed=seed)
    pred = emu.predict_loss(cand)

    # Rank by EI, then greedily pick diverse points.
    order = np.argsort(ei)[::-1]
    chosen: list[int] = []
    for idx in order:
        if all(np.linalg.norm(cand[idx] - cand[c]) >= min_dist for c in chosen):
            chosen.append(idx)
        if len(chosen) >= q:
            break
    for idx in order:
        if len(chosen) >= q:
            break
        if idx not in chosen:
            chosen.append(idx)

    means = emu.predict(cand[chosen])
    out = []
    for k, idx in enumerate(chosen):
        xs = cand[idx]
        theta_raw = ds.unscale(xs)
        out.append({
            "rank": k,
            "theta": {p: float(v) for p, v in zip(ds.theta_params, theta_raw)},
            "scaled_point": xs.tolist(),
            "ei": float(ei[idx]),
            "pred_loss": float(pred[idx]),
            "nn_dist_to_data": float(np.sqrt(((Xtr - xs) ** 2).sum(1)).min()),
            "pred_moments": {n: float(means[n][k]) for n in ds.moment_names},
        })
    return out


def propose_exploit(emu: MomentEmulator, pool_size: int = 15000,
                    n_local: int = 10, seed: int = 1):
    """Pure-exploitation point: argmin emulated loss (surrogate's best guess)."""
    ds = emu.ds
    rng = np.random.default_rng(seed)
    d = len(ds.theta_params)
    pool = latin_hypercube(pool_size, d, rng)
    loss = emu.predict_loss(pool)
    best = int(np.argmin(loss))
    x = pool[best].copy()
    best_loss = loss[best]
    step = 0.05
    for _ in range(n_local):
        pert = np.clip(x + rng.normal(0, step, size=(64, d)), 0.0, 1.0)
        pl = emu.predict_loss(pert)
        j = int(np.argmin(pl))
        if pl[j] < best_loss:
            x, best_loss = pert[j], pl[j]
        step *= 0.8
    theta_raw = ds.unscale(x)
    means = emu.predict(x[None, :])
    return {
        "theta": {p: float(v) for p, v in zip(ds.theta_params, theta_raw)},
        "scaled_point": x.tolist(),
        "pred_loss": float(best_loss),
        "pred_moments": {n: float(means[n][0]) for n in ds.moment_names},
    }


def pareto_front(emu: MomentEmulator, objective_moments: list[str],
                 pool_size: int = 25000, seed: int = 2):
    """Map the achievable Pareto front of per-moment |deviation| objectives.

    Each objective is |emulated_moment - target| for a moment in
    ``objective_moments``. Returns the nondominated emulated points, the
    achievable minimum of each objective marginally, and the best *joint* point
    (min over the sum of normalized objectives) -- the cleanest test of whether
    the conflicting targets are simultaneously reachable.
    """
    ds = emu.ds
    rng = np.random.default_rng(seed)
    d = len(ds.theta_params)
    pool = latin_hypercube(pool_size, d, rng)
    means = emu.predict(pool)

    objs = np.zeros((pool_size, len(objective_moments)))
    for j, name in enumerate(objective_moments):
        objs[:, j] = np.abs(means[name] - float(ds.targets[name]))

    # Nondominated set (minimization on all objectives). Computed on a
    # downsample for speed; the front is small and this is for visualization.
    nd_pool = min(pool_size, 8000)
    sub = rng.choice(pool_size, size=nd_pool, replace=False) if pool_size > nd_pool \
        else np.arange(pool_size)
    nd_sub = _nondominated(objs[sub])
    nd_idx = sub[nd_sub]

    marginal_min = {name: float(objs[:, j].min())
                    for j, name in enumerate(objective_moments)}
    # Best joint ("compromise") point: minimize the spread-normalized sum of
    # objectives, but restricted to the Pareto-nondominated set so the returned
    # point is guaranteed not to be dominated (a plain argmin over the whole pool
    # can return a dominated point). This is a scalarization compromise, not a
    # unique optimum; the marginal_min vs best_joint_objs comparison is what
    # answers "are these targets jointly reachable?".
    spread = objs.max(axis=0) - objs.min(axis=0)
    spread = np.where(spread > 0, spread, 1.0)
    joint = (objs / spread).sum(axis=1)
    # nd is computed on the downsample `sub`; map it back to full-pool indices.
    nd_full = np.zeros(pool_size, dtype=bool)
    nd_full[nd_idx] = True
    joint_masked = np.where(nd_full, joint, np.inf) if nd_idx.size else joint
    bj = int(np.argmin(joint_masked))
    theta_bj = ds.unscale(pool[bj])

    return {
        "objective_moments": objective_moments,
        "n_pool": pool_size,
        "n_nondominated": int(nd_idx.size),
        "nondominated_objs": objs[nd_idx],
        "nondominated_theta": ds.unscale(pool[nd_idx]),
        "marginal_min": marginal_min,
        "best_joint_theta": {p: float(v) for p, v in zip(ds.theta_params, theta_bj)},
        "best_joint_objs": {name: float(objs[bj, j])
                            for j, name in enumerate(objective_moments)},
        "best_joint_moments": {name: float(means[name][bj])
                               for name in ds.moment_names},
    }


def _nondominated(objs: np.ndarray, chunk: int = 1000) -> np.ndarray:
    """Boolean mask of Pareto-nondominated rows (minimization), vectorized.

    Point i is dominated iff some j has objs[j] <= objs[i] on all objectives and
    objs[j] < objs[i] on at least one. Computed in chunks of candidates against
    the full set via broadcasting, so cost is O(n^2 * d) in numpy rather than a
    Python double loop.
    """
    objs = np.asarray(objs, float)
    n = objs.shape[0]
    nd = np.ones(n, dtype=bool)
    P = objs  # (n, d)
    for s in range(0, n, chunk):
        e = min(s + chunk, n)
        C = objs[s:e]                              # (m, d)
        le = np.all(P[None, :, :] <= C[:, None, :], axis=2)   # (m, n)
        lt = np.any(P[None, :, :] < C[:, None, :], axis=2)    # (m, n)
        dominated = np.any(le & lt, axis=1)        # (m,)
        nd[s:e] = ~dominated
    return nd
