"""Pure interface tests for the block-Toeplitz initial Jacobian; no native model."""
import unittest

import numpy as np

from e5f_ssj_toeplitz_jacobian import (assemble_jacobian, central_column, diagonal_default,
                                       lag_profiles, toeplitz_block)


def _toy_kernels(T):
    # D[i][j] is a dict lag -> coefficient, a stationary linear toy mapping.
    rng = np.random.default_rng(0)
    kernels = [[{k: float(rng.normal()) for k in range(-2, 3)} for _ in range(3)] for _ in range(3)]
    return kernels


def _toy_residual(u, kernels, T):
    # u is (3, T) in log coordinates; residual block i at t = sum_j sum_k D_ij[k] u_j[t-k].
    out = np.zeros((3, T))
    for i in range(3):
        for j in range(3):
            for k, coef in kernels[i][j].items():
                for t in range(T):
                    s = t - k
                    if 0 <= s < T:
                        out[i, t] += coef * u[j, s]
    return out.reshape(-1)


class TestToeplitzJacobian(unittest.TestCase):
    def test_recovers_exact_toeplitz_linear_map(self):
        T, s, h = 10, 5, 1e-4
        kernels = _toy_kernels(T)
        base = np.zeros((3, T))
        columns = []
        for j in range(3):
            up, dn = base.copy(), base.copy()
            up[j, s] += h; dn[j, s] -= h
            columns.append(central_column(_toy_residual(up, kernels, T), _toy_residual(dn, kernels, T), h, T))
        J, receipt = assemble_jacobian(columns, T, s)
        # Exact dense Jacobian of the linear toy map by columns.
        exact = np.zeros((3 * T, 3 * T))
        for col in range(3 * T):
            e = np.zeros(3 * T); e[col] = 1.0
            exact[:, col] = _toy_residual(e.reshape(3, T), kernels, T)
        np.testing.assert_allclose(J, exact, rtol=0, atol=1e-9)
        self.assertEqual(receipt["measured_lags"], list(range(-5, 5)))
        self.assertFalse(receipt["fake_news_derivatives_constructed"])

    def test_lag_profile_indexing_matches_column(self):
        T, s = 6, 2
        column = np.arange(18.0)
        lags, profiles = lag_profiles(column, T, s)
        self.assertEqual(lags.tolist(), [-2, -1, 0, 1, 2, 3])
        # rebate block at lag +1 is date s+1 = 3 of the third block.
        self.assertEqual(profiles[2, list(lags).index(1)], column[2 * T + 3])

    def test_toeplitz_block_zero_fills_unmeasured_lags(self):
        block = toeplitz_block(np.array([0, 1]), np.array([-2.0, -0.2]), 4)
        expected = np.array([[-2, 0, 0, 0], [-0.2, -2, 0, 0], [0, -0.2, -2, 0], [0, 0, -0.2, -2]], dtype=float)
        np.testing.assert_array_equal(block, expected)

    def test_diagonal_measurement_reproduces_default_layout(self):
        T = 3
        default = diagonal_default(T, 1.63)
        self.assertEqual(default.shape, (9, 9))
        np.testing.assert_array_equal(np.diag(default), [-1.63] * 3 + [-200.0] * 6)

    def test_rejects_bad_shapes_and_steps(self):
        with self.assertRaisesRegex(ValueError, "shape"):
            central_column(np.zeros(5), np.zeros(5), 1e-5, 2)
        with self.assertRaisesRegex(ValueError, "step"):
            central_column(np.zeros(6), np.zeros(6), 0.0, 2)
        with self.assertRaisesRegex(ValueError, "horizon"):
            lag_profiles(np.zeros(6), 2, 2)


if __name__ == "__main__":
    unittest.main()


class TestAssembleFromReceipt(unittest.TestCase):
    def test_extends_measured_profiles_to_longer_horizon(self):
        from e5f_ssj_toeplitz_jacobian import assemble_from_receipt
        T, s, h = 10, 5, 1e-4
        kernels = _toy_kernels(T)
        base = np.zeros((3, T))
        columns = []
        for j in range(3):
            up, dn = base.copy(), base.copy()
            up[j, s] += h; dn[j, s] -= h
            columns.append(central_column(_toy_residual(up, kernels, T), _toy_residual(dn, kernels, T), h, T))
        J10, receipt = assemble_jacobian(columns, T, s)
        J20, info = assemble_from_receipt(receipt, 20)
        self.assertEqual(J20.shape, (60, 60))
        # Same-horizon rebuild is exact; the longer horizon reproduces the toy's exact banded map.
        np.testing.assert_allclose(assemble_from_receipt(receipt, T)[0], J10, rtol=0, atol=1e-12)
        exact = np.zeros((60, 60))
        for col in range(60):
            e = np.zeros(60); e[col] = 1.0
            exact[:, col] = _toy_residual(e.reshape(3, 20), kernels, 20)
        np.testing.assert_allclose(J20, exact, rtol=0, atol=1e-9)
        self.assertEqual(info["source_horizon"], 10)
