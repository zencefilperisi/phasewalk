"""
Sanity tests for the analysis module.

Run with:
    pytest tests/test_analysis.py -v
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np

from classical.analysis import (
    classify_fixed_point,
    find_fixed_points,
    find_hopf_threshold,
    jacobian,
    sigmoid_derivative,
)
from classical.model import WCParams, sigmoid


def test_sigmoid_derivative_matches_numeric():
    """Closed-form derivative must match the central finite difference."""
    a, theta = 1.3, 4.0
    h = 1e-5
    for x in np.linspace(-2, 8, 50):
        analytic = sigmoid_derivative(x, a, theta)
        numeric = (sigmoid(x + h, a, theta) - sigmoid(x - h, a, theta)) / (2 * h)
        assert abs(analytic - numeric) < 1e-7, f"mismatch at x={x}"


def test_jacobian_shape():
    """Jacobian must be a 2x2 numpy array."""
    p = WCParams()
    J = jacobian([0.5, 0.5], p)
    assert J.shape == (2, 2)
    assert isinstance(J, np.ndarray)


def test_jacobian_matches_numeric():
    """Analytic Jacobian must agree with finite-difference Jacobian."""
    from classical.model import wilson_cowan_rhs

    p = WCParams()
    state = np.array([0.3, 0.4])
    h = 1e-6

    J_analytic = jacobian(state, p)

    J_numeric = np.zeros((2, 2))
    for k in range(2):
        plus = state.copy(); plus[k] += h
        minus = state.copy(); minus[k] -= h
        f_plus = np.array(wilson_cowan_rhs(0.0, plus, p))
        f_minus = np.array(wilson_cowan_rhs(0.0, minus, p))
        J_numeric[:, k] = (f_plus - f_minus) / (2 * h)

    assert np.allclose(J_analytic, J_numeric, atol=1e-6)


def test_classify_stable_focus():
    eigs = np.array([-1.0 + 2.0j, -1.0 - 2.0j])
    assert classify_fixed_point(eigs) == "stable focus"


def test_classify_unstable_focus():
    eigs = np.array([0.5 + 1.0j, 0.5 - 1.0j])
    assert classify_fixed_point(eigs) == "unstable focus"


def test_classify_saddle():
    eigs = np.array([1.0 + 0j, -1.5 + 0j])
    assert classify_fixed_point(eigs) == "saddle"


def test_classify_stable_node():
    eigs = np.array([-1.0 + 0j, -2.0 + 0j])
    assert classify_fixed_point(eigs) == "stable node"


def test_find_fixed_points_low_drive():
    """At very low drive there should be at least one stable fixed point
    near the origin (silent state)."""
    p = WCParams(P=0.0, Q=0.0)
    fps = find_fixed_points(p)
    assert len(fps) >= 1
    # At least one of them must be stable
    assert any(fp.is_stable for fp in fps)
    # And one must be near the origin
    assert any(fp.E < 0.2 and fp.I < 0.2 for fp in fps)


def test_find_fixed_points_oscillatory_regime():
    """At P=1.5 the unique fixed point should be unstable (limit cycle)."""
    p = WCParams(P=1.5)
    fps = find_fixed_points(p)
    assert len(fps) >= 1
    # No fixed point should be stable in the oscillatory regime
    assert not any(fp.is_stable for fp in fps)


def test_hopf_threshold_in_expected_range():
    """For default parameters the Hopf threshold (in P) should sit
    somewhere in (0.5, 1.5)."""
    base = WCParams()
    p_hopf = find_hopf_threshold(base, param_name="P", lo=0.0, hi=2.0)
    assert p_hopf is not None
    assert 0.5 < p_hopf < 1.5


if __name__ == "__main__":
    test_sigmoid_derivative_matches_numeric()
    test_jacobian_shape()
    test_jacobian_matches_numeric()
    test_classify_stable_focus()
    test_classify_unstable_focus()
    test_classify_saddle()
    test_classify_stable_node()
    test_find_fixed_points_low_drive()
    test_find_fixed_points_oscillatory_regime()
    test_hopf_threshold_in_expected_range()
    print("All analysis tests passed.")
