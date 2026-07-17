"""Multi-output moment emulator: theta -> moment vector.

One independent ARD-GP per target moment, trained on the shared (scaled) design
matrix. The per-dimension squared-difference tensor is computed once and reused
across moments. The emulator exposes:

* ``predict`` -> dict of moment means (and optional variances),
* ``predict_loss`` -> the model's own weighted-SSE loss recomputed from the
  emulated moment vector (so the surrogate optimizes exactly the real objective),
* ``predict_loss_mc`` -> Monte-Carlo expected loss / improvement using the GP
  predictive uncertainty across moments (used by the Bayesian-optimization
  acquisition),
* ``cross_validate`` -> K-fold held-out R^2 per moment and held-out loss
  accuracy. This is the decisive test: if the emulator cannot predict held-out
  moments, none of the downstream identification or BO claims hold.
"""

from __future__ import annotations

import pickle
import time

import numpy as np

from . import data as D
from .gp import ARDGaussianProcess, pairwise_sqdiffs


class MomentEmulator:
    def __init__(self, dataset: D.Dataset, n_restarts: int = 2, max_iter: int = 80,
                 seed: int = 0, verbose: bool = True):
        self.ds = dataset
        self.n_restarts = n_restarts
        self.max_iter = max_iter
        self.seed = seed
        self.verbose = verbose
        self.gps: dict[str, ARDGaussianProcess] = {}
        self.train_idx_: np.ndarray | None = None

    def fit(self, train_idx: np.ndarray):
        self.train_idx_ = np.asarray(train_idx)
        Xs = self.ds.Xs[self.train_idx_]
        D2 = pairwise_sqdiffs(Xs)  # shared across moments
        t0 = time.time()
        for j, name in enumerate(self.ds.moment_names):
            y = self.ds.M[name][self.train_idx_]
            gp = ARDGaussianProcess(n_restarts=self.n_restarts,
                                    max_iter=self.max_iter, seed=self.seed + j)
            gp.fit(Xs, y, D2=D2)
            self.gps[name] = gp
            if self.verbose:
                print(f"  fit {name:38s} lml={gp.lml_:8.1f} "
                      f"sn2={gp.sn2_:.2e}  ({time.time()-t0:5.1f}s)")
        return self

    # ---- persistence ----------------------------------------------------
    def save(self, path: str):
        """Pickle the fitted GPs + training indices (avoids re-fitting)."""
        with open(path, "wb") as fh:
            pickle.dump({"gps": self.gps, "train_idx": self.train_idx_,
                         "target_set": self.ds.target_set_name}, fh)

    @classmethod
    def load(cls, path: str, dataset: D.Dataset) -> "MomentEmulator":
        with open(path, "rb") as fh:
            blob = pickle.load(fh)
        emu = cls(dataset, verbose=False)
        emu.gps = blob["gps"]
        emu.train_idx_ = blob["train_idx"]
        return emu

    # ---- prediction -----------------------------------------------------
    def predict(self, Xs: np.ndarray, return_var: bool = False):
        Xs = np.atleast_2d(np.asarray(Xs, float))
        means, vars_ = {}, {}
        for name, gp in self.gps.items():
            if return_var:
                m, v = gp.predict(Xs, return_var=True)
                means[name], vars_[name] = m, v
            else:
                means[name] = gp.predict(Xs)
        return (means, vars_) if return_var else means

    def predict_loss(self, Xs: np.ndarray) -> np.ndarray:
        """Emulated weighted-SSE loss at each row of Xs (vectorized)."""
        Xs = np.atleast_2d(np.asarray(Xs, float))
        means = self.predict(Xs)
        loss = np.zeros(Xs.shape[0])
        for name in self.ds.moment_names:
            w = float(self.ds.weights.get(name, 1.0))
            t = float(self.ds.targets[name])
            loss += w * (means[name] - t) ** 2
        return loss

    def predict_loss_mc(self, Xs: np.ndarray, n_samples: int = 64,
                        seed: int = 0):
        """Monte-Carlo distribution of the loss using GP moment uncertainty.

        Returns (mean_loss, samples) where samples has shape (M, n_samples).
        Each moment is drawn from its independent predictive Gaussian and the
        loss is the model's weighted SSE; this captures that the loss is a sum of
        squared uncertain moments (a noncentral-chi-square-type object, not
        Gaussian), which a single scalar-loss GP would get wrong.

        Moments are treated as independent. This is exact for the emulator as
        built (one independent GP per moment, so there is no modeled
        cross-moment covariance), not a hidden approximation of a richer joint
        model. A multi-output GP with a coregionalization kernel would add that
        covariance; that is a documented scale-up path, not a correctness gap.
        """
        Xs = np.atleast_2d(np.asarray(Xs, float))
        means, vars_ = self.predict(Xs, return_var=True)
        rng = np.random.default_rng(seed)
        M = Xs.shape[0]
        samples = np.zeros((M, n_samples))
        for name in self.ds.moment_names:
            w = float(self.ds.weights.get(name, 1.0))
            t = float(self.ds.targets[name])
            mu = means[name][:, None]
            sd = np.sqrt(np.maximum(vars_[name], 0.0))[:, None]
            draws = mu + sd * rng.standard_normal((M, n_samples))
            samples += w * (draws - t) ** 2
        return samples.mean(axis=1), samples

    # ---- cross-validation ----------------------------------------------
    def cross_validate(self, idx: np.ndarray, k: int = 5, seed: int = 0):
        """K-fold CV on the curated subset; returns per-moment R^2 and loss fit.

        Trains a fresh emulator on k-1 folds and predicts the held-out fold,
        aggregating out-of-fold predictions so every point is predicted once by a
        model that did not see it.
        """
        idx = np.asarray(idx)
        n = idx.size
        rng = np.random.default_rng(seed)
        perm = rng.permutation(n)
        folds = np.array_split(perm, k)

        oof_mom = {name: np.full(n, np.nan) for name in self.ds.moment_names}
        oof_loss = np.full(n, np.nan)

        for fi, test_local in enumerate(folds):
            train_local = np.setdiff1d(perm, test_local)
            tr = idx[train_local]
            te = idx[test_local]
            Xtr = self.ds.Xs[tr]
            D2 = pairwise_sqdiffs(Xtr)
            gps = {}
            for j, name in enumerate(self.ds.moment_names):
                gp = ARDGaussianProcess(n_restarts=self.n_restarts,
                                        max_iter=self.max_iter, seed=self.seed + j)
                gp.fit(Xtr, self.ds.M[name][tr], D2=D2)
                gps[name] = gp
            Xte = self.ds.Xs[te]
            loss_pred = np.zeros(te.size)
            for name in self.ds.moment_names:
                m = gps[name].predict(Xte)
                oof_mom[name][test_local] = m
                w = float(self.ds.weights.get(name, 1.0))
                t = float(self.ds.targets[name])
                loss_pred += w * (m - t) ** 2
            oof_loss[test_local] = loss_pred
            if self.verbose:
                print(f"  CV fold {fi+1}/{k} done")

        # Per-moment out-of-fold R^2 and RMSE (oof arrays are in idx order).
        r2 = {}
        rmse = {}
        for name in self.ds.moment_names:
            y_true = self.ds.M[name][idx]
            y_pred = oof_mom[name]
            ss_res = np.nansum((y_true - y_pred) ** 2)
            ss_tot = np.nansum((y_true - np.nanmean(y_true)) ** 2)
            r2[name] = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
            rmse[name] = np.sqrt(np.nanmean((y_true - y_pred) ** 2))

        true_loss = self.ds.loss[idx]
        ss_res = np.nansum((true_loss - oof_loss) ** 2)
        ss_tot = np.nansum((true_loss - np.nanmean(true_loss)) ** 2)
        loss_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        # Spearman rank correlation (no scipy): rank then Pearson.
        loss_rank_corr = _rank_corr(true_loss, oof_loss)

        return {
            "moment_r2": r2,
            "moment_rmse": rmse,
            "loss_r2": float(loss_r2),
            "loss_rank_corr": float(loss_rank_corr),
            "oof_loss": oof_loss,
            "true_loss": true_loss,
            "idx": idx,
        }


def _rank_corr(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if a.size < 3:
        return np.nan
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    denom = np.sqrt(np.sum(ra ** 2) * np.sum(rb ** 2))
    return float(np.sum(ra * rb) / denom) if denom > 0 else np.nan
