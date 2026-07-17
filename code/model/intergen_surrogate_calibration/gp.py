"""Pure-numpy ARD Gaussian-process regression.

This module has no third-party dependencies beyond numpy. It exists because the
model virtualenv ships numpy/pandas/matplotlib/numba but NOT scipy or
scikit-learn, so the surrogate stack hand-rolls the GP it needs.

The GP is a zero-mean (on standardized targets) anisotropic squared-exponential
("ARD RBF") process. Hyperparameters -- one lengthscale per input dimension,
one signal variance, one observation-noise variance -- are fit by maximizing the
log marginal likelihood with Adam using analytic gradients.

Design choices relevant to the calibration use case:

* ARD lengthscales are reported back out. A long lengthscale on dimension d means
  the emulated moment barely moves with parameter d -- a direct, model-derived
  relevance / weak-identification signal that complements the structural
  Jacobian audit.
* ``predict_grad`` returns the analytic gradient of the posterior mean with
  respect to the (scaled) input, which the identification module turns into a
  smooth surrogate Jacobian ``d m / d theta``.
* ``pairwise_sqdiffs`` is exposed so the multi-output emulator can compute the
  per-dimension squared-difference tensor once on the shared design matrix and
  reuse it across every moment GP, instead of recomputing it per moment.

Run ``python -m intergen_surrogate_calibration.gp`` (or execute this file) to run
the self-test, which checks the analytic gradients against finite differences
and confirms the GP recovers a known smooth function out of sample.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def pairwise_sqdiffs(X: np.ndarray) -> np.ndarray:
    """Per-dimension squared differences.

    Returns an array ``D2`` of shape ``(D, N, N)`` with
    ``D2[d, i, j] = (X[i, d] - X[j, d]) ** 2``. Depends only on the design
    matrix, so it can be precomputed once and shared across output GPs.
    """
    X = np.asarray(X, dtype=float)
    # (N, 1, D) - (1, N, D) -> (N, N, D), then move D to front.
    diff = X[:, None, :] - X[None, :, :]
    return np.transpose(diff * diff, (2, 0, 1))


@dataclass
class ARDGaussianProcess:
    """Anisotropic-RBF GP regression in pure numpy.

    Inputs ``X`` are assumed already scaled by the caller (the emulator scales
    parameters into the unit cube). Targets ``y`` are standardized internally.
    """

    n_restarts: int = 2
    max_iter: int = 80
    patience: int = 12
    lr: float = 0.05
    jitter: float = 1e-6
    seed: int = 0
    # Hyperparameter box (log space) used to keep optimization well-behaved.
    log_ls_bounds: tuple[float, float] = (np.log(0.02), np.log(20.0))
    log_sf2_bounds: tuple[float, float] = (np.log(1e-3), np.log(1e3))
    log_sn2_bounds: tuple[float, float] = (np.log(1e-6), np.log(1.0))

    # Fitted state (filled by ``fit``).
    X_: np.ndarray = field(default=None, repr=False)
    ls_: np.ndarray = field(default=None, repr=False)
    sf2_: float = None
    sn2_: float = None
    y_mean_: float = None
    y_std_: float = None
    alpha_: np.ndarray = field(default=None, repr=False)
    Kinv_: np.ndarray = field(default=None, repr=False)
    lml_: float = None

    # ---- kernel ---------------------------------------------------------
    def _kernel_from_D2(self, D2: np.ndarray, ls: np.ndarray, sf2: float) -> np.ndarray:
        inv_l2 = 1.0 / (ls * ls)
        r2 = np.tensordot(inv_l2, D2, axes=([0], [0]))  # (N, N)
        return sf2 * np.exp(-0.5 * r2)

    # ---- objective: negative log marginal likelihood + grad -------------
    def _nlml_and_grad(self, params, D2, y, n):
        log_ls = params[:-2]
        log_sf2 = params[-2]
        log_sn2 = params[-1]
        ls = np.exp(log_ls)
        sf2 = np.exp(log_sf2)
        sn2 = np.exp(log_sn2)

        K = self._kernel_from_D2(D2, ls, sf2)
        Ky = K + (sn2 + self.jitter) * np.eye(n)
        # Cholesky for stable log-det; explicit inverse for the trace term.
        try:
            L = np.linalg.cholesky(Ky)
        except np.linalg.LinAlgError:
            Ky = Ky + 1e-3 * np.eye(n)
            L = np.linalg.cholesky(Ky)
        logdet = 2.0 * np.sum(np.log(np.diag(L)))
        Kinv = np.linalg.solve(Ky, np.eye(n))
        alpha = Kinv @ y

        # nlml is the NEGATIVE log marginal likelihood (the quantity Adam
        # minimizes); LML itself is reported as -nlml after the fit.
        nlml = 0.5 * float(y @ alpha) + 0.5 * logdet + 0.5 * n * np.log(2.0 * np.pi)

        # dLML/dp = 0.5 tr((alpha alpha^T - Kinv) dKy/dp); nlml grad is negative.
        A = np.outer(alpha, alpha) - Kinv  # (N, N)
        grad = np.zeros_like(params)
        inv_l2 = 1.0 / (ls * ls)
        for d in range(len(ls)):
            dK_d = K * (D2[d] * inv_l2[d])  # dK/d log_ls_d
            grad[d] = -0.5 * np.sum(A * dK_d)
        grad[-2] = -0.5 * np.sum(A * K)        # dK/d log_sf2 = K
        grad[-1] = -0.5 * np.sum(A * (sn2 * np.eye(n)))  # dKy/d log_sn2 = sn2 I
        return nlml, grad

    def _project(self, params):
        D = len(params) - 2
        params[:D] = np.clip(params[:D], *self.log_ls_bounds)
        params[D] = np.clip(params[D], *self.log_sf2_bounds)
        params[D + 1] = np.clip(params[D + 1], *self.log_sn2_bounds)
        return params

    def _adam(self, params0, D2, y, n):
        params = self._project(params0.copy())
        m = np.zeros_like(params)
        v = np.zeros_like(params)
        b1, b2, eps = 0.9, 0.999, 1e-8
        best_params, best_nlml = params.copy(), np.inf
        since_improve = 0
        for t in range(1, self.max_iter + 1):
            nlml, grad = self._nlml_and_grad(params, D2, y, n)
            if np.isfinite(nlml) and nlml < best_nlml - 1e-4:
                best_nlml, best_params = nlml, params.copy()
                since_improve = 0
            else:
                since_improve += 1
                if since_improve >= self.patience:
                    break  # marginal likelihood has plateaued
            m = b1 * m + (1 - b1) * grad
            v = b2 * v + (1 - b2) * (grad * grad)
            mhat = m / (1 - b1 ** t)
            vhat = v / (1 - b2 ** t)
            params = params - self.lr * mhat / (np.sqrt(vhat) + eps)
            params = self._project(params)
        # Evaluate the final point too.
        nlml, _ = self._nlml_and_grad(params, D2, y, n)
        if np.isfinite(nlml) and nlml < best_nlml:
            best_nlml, best_params = nlml, params.copy()
        return best_params, best_nlml

    def fit(self, X: np.ndarray, y: np.ndarray, D2: np.ndarray | None = None):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).reshape(-1)
        n, D = X.shape
        if D2 is None:
            D2 = pairwise_sqdiffs(X)

        self.y_mean_ = float(np.mean(y))
        # Floor the scale so a (near-)constant target does not blow up ys; a
        # constant moment carries no information and its GP collapses to a flat
        # mean, which is the correct behavior.
        self.y_std_ = max(float(np.std(y)), 1e-8)
        ys = (y - self.y_mean_) / self.y_std_

        rng = np.random.default_rng(self.seed)
        std_per_dim = np.std(X, axis=0)
        std_per_dim = np.where(std_per_dim > 1e-6, std_per_dim, 0.5)

        inits = []
        # Restart 0: fixed moderate lengthscale.
        inits.append(np.concatenate([np.log(np.full(D, 0.5)), [0.0, np.log(0.01)]]))
        # Restart 1: per-dimension std heuristic.
        inits.append(np.concatenate([np.log(std_per_dim), [0.0, np.log(0.01)]]))
        # Extra restarts: jittered.
        for _ in range(max(0, self.n_restarts - 2)):
            jit = rng.normal(0.0, 0.5, size=D)
            inits.append(np.concatenate([np.log(0.5) + jit, [0.0, np.log(0.05)]]))
        inits = inits[: max(1, self.n_restarts)]

        best_params, best_nlml = None, np.inf
        for p0 in inits:
            params, nlml = self._adam(np.asarray(p0, dtype=float), D2, ys, n)
            if nlml < best_nlml:
                best_nlml, best_params = nlml, params

        D_ = len(best_params) - 2
        self.ls_ = np.exp(best_params[:D_])
        self.sf2_ = float(np.exp(best_params[D_]))
        self.sn2_ = float(np.exp(best_params[D_ + 1]))
        self.lml_ = -best_nlml

        # Cache prediction objects.
        K = self._kernel_from_D2(D2, self.ls_, self.sf2_)
        Ky = K + (self.sn2_ + self.jitter) * np.eye(n)
        self.Kinv_ = np.linalg.solve(Ky, np.eye(n))
        self.alpha_ = self.Kinv_ @ ys
        self.X_ = X
        return self

    # ---- prediction -----------------------------------------------------
    def _cross_kernel(self, Xstar: np.ndarray) -> np.ndarray:
        """k(Xstar, X_train) -> (M, N), memory-efficient (no (M,N,D) tensor).

        Uses ``r2 = ||a||^2 + ||b||^2 - 2 a.b`` on inputs pre-scaled by 1/ls, so
        peak memory is O(M*N) rather than O(M*N*D).
        """
        inv_l = 1.0 / self.ls_
        A = Xstar * inv_l                 # (M, D)
        B = self.X_ * inv_l               # (N, D)
        a2 = np.sum(A * A, axis=1)        # (M,)
        b2 = np.sum(B * B, axis=1)        # (N,)
        r2 = a2[:, None] + b2[None, :] - 2.0 * (A @ B.T)
        np.maximum(r2, 0.0, out=r2)
        return self.sf2_ * np.exp(-0.5 * r2)

    def predict(self, Xstar: np.ndarray, return_var: bool = False,
                chunk: int = 4000):
        Xstar = np.asarray(Xstar, dtype=float)
        if Xstar.ndim == 1:
            Xstar = Xstar[None, :]
        M = Xstar.shape[0]
        mean = np.empty(M)
        var = np.empty(M) if return_var else None
        for s in range(0, M, chunk):
            e = min(s + chunk, M)
            Ks = self._cross_kernel(Xstar[s:e])          # (m, N)
            mean[s:e] = self.y_mean_ + self.y_std_ * (Ks @ self.alpha_)
            if return_var:
                KsKinv = Ks @ self.Kinv_                  # (m, N)
                v = np.sum(KsKinv * Ks, axis=1)           # (m,)
                var_s = np.maximum(self.sf2_ - v, 0.0)
                var[s:e] = (self.y_std_ ** 2) * var_s
        if not return_var:
            return mean
        return mean, var

    def predict_grad(self, xstar: np.ndarray) -> np.ndarray:
        """Gradient of the posterior mean wrt the (scaled) input, at one point.

        Returns a length-D vector ``d mean / d x`` in the scaled-input space.
        """
        xstar = np.asarray(xstar, dtype=float).reshape(1, -1)
        Ks = self._cross_kernel(xstar)[0]  # (N,)
        inv_l2 = 1.0 / (self.ls_ * self.ls_)
        diff = (xstar[0][None, :] - self.X_)  # (N, D)
        # d k(x*, x_i) / d x*_d = k * (-(x*_d - x_id) / ls_d^2)
        dk = -Ks[:, None] * diff * inv_l2[None, :]  # (N, D)
        grad_s = self.alpha_ @ dk  # (D,)
        return self.y_std_ * grad_s


# ----------------------------------------------------------------------------
# Self-test
# ----------------------------------------------------------------------------
def _self_test() -> None:
    rng = np.random.default_rng(0)
    D = 4
    n = 120

    def f(X):
        return (np.sin(3 * X[:, 0]) + 0.5 * X[:, 1] ** 2
                + 0.3 * X[:, 2] * X[:, 3] + 0.1 * X[:, 0] * X[:, 1])

    X = rng.uniform(0, 1, size=(n, D))
    y = f(X) + rng.normal(0, 0.01, size=n)

    gp = ARDGaussianProcess(n_restarts=2, max_iter=120, seed=1)
    gp.fit(X, y)

    # 1) Analytic vs finite-difference gradient of the NLML.
    D2 = pairwise_sqdiffs(X)
    ys = (y - gp.y_mean_) / gp.y_std_
    p = np.concatenate([np.log(np.full(D, 0.4)), [0.1, np.log(0.02)]])
    nlml0, grad = gp._nlml_and_grad(p, D2, ys, n)
    fd = np.zeros_like(p)
    eps = 1e-5
    for i in range(len(p)):
        pp = p.copy(); pp[i] += eps
        pm = p.copy(); pm[i] -= eps
        fd[i] = (gp._nlml_and_grad(pp, D2, ys, n)[0]
                 - gp._nlml_and_grad(pm, D2, ys, n)[0]) / (2 * eps)
    g_err = np.max(np.abs(grad - fd) / (np.abs(fd) + 1e-6))
    print(f"[self-test] NLML grad max rel error vs FD: {g_err:.2e}")
    assert g_err < 1e-3, "NLML gradient mismatch"

    # 2) Out-of-sample predictive accuracy.
    Xt = rng.uniform(0, 1, size=(400, D))
    yt = f(Xt)
    pred = gp.predict(Xt)
    ss_res = np.sum((yt - pred) ** 2)
    ss_tot = np.sum((yt - np.mean(yt)) ** 2)
    r2 = 1 - ss_res / ss_tot
    print(f"[self-test] out-of-sample R^2: {r2:.4f}")
    assert r2 > 0.95, "GP failed to recover smooth function"

    # 3) predict_grad vs finite difference of the posterior mean.
    x0 = rng.uniform(0.2, 0.8, size=D)
    g_an = gp.predict_grad(x0)
    g_fd = np.zeros(D)
    for i in range(D):
        xp = x0.copy(); xp[i] += eps
        xm = x0.copy(); xm[i] -= eps
        g_fd[i] = (gp.predict(xp)[0] - gp.predict(xm)[0]) / (2 * eps)
    pg_err = np.max(np.abs(g_an - g_fd) / (np.abs(g_fd) + 1e-4))
    print(f"[self-test] predict_grad max rel error vs FD: {pg_err:.2e}")
    assert pg_err < 1e-3, "predict_grad mismatch"

    # 4) Predictive variance is small in-sample, larger out in a far corner.
    _, var_near = gp.predict(X[:5], return_var=True)
    _, var_far = gp.predict(np.full((1, D), 5.0), return_var=True)
    print(f"[self-test] mean in-sample var {np.mean(var_near):.2e}, "
          f"far-corner var {var_far[0]:.2e}")
    assert var_far[0] > np.mean(var_near), "variance not larger away from data"

    print("[self-test] ARDGaussianProcess: ALL CHECKS PASSED")


if __name__ == "__main__":
    _self_test()
