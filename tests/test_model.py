"""
Basic sanity tests for the Wilson-Cowan implementation.

Run with:
    pytest tests/

Or a single file:
    pytest tests/test_model.py -v
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np

from classical.model import WCParams, sigmoid, simulate, wilson_cowan_rhs


def test_sigmoid_bounds():
    """Sigmoid output must lie in [0, 1] for any real input."""
    x = np.linspace(-100, 100, 1000)
    y = sigmoid(x, a=1.0, theta=0.0)
    assert np.all(y >= 0.0)
    assert np.all(y <= 1.0)


def test_sigmoid_threshold():
    """S(theta) must equal 0.5 by construction."""
    assert abs(sigmoid(4.0, a=1.3, theta=4.0) - 0.5) < 1e-12


def test_rhs_returns_two_values():
    """The RHS must return [dE/dt, dI/dt]."""
    p = WCParams()
    out = wilson_cowan_rhs(0.0, [0.1, 0.1], p)
    assert len(out) == 2


def test_simulation_runs():
    """The integrator must run without raising."""
    p = WCParams()
    t, E, I = simulate(p, E0=0.1, I0=0.1, t_max=10, dt=0.1)
    assert len(t) == len(E) == len(I)
    assert len(t) > 0


def test_activities_in_unit_interval():
    """Activities must remain in [0, 1] (sigmoid is bounded)."""
    p = WCParams()
    _, E, I = simulate(p, E0=0.5, I0=0.5, t_max=50)
    # Allow small numerical tolerance
    assert np.all(E >= -1e-6) and np.all(E <= 1.0 + 1e-6)
    assert np.all(I >= -1e-6) and np.all(I <= 1.0 + 1e-6)


def test_low_input_settles_to_low_activity():
    """With no external drive, activities must decay to a low fixed point."""
    p = WCParams(P=0.0, Q=0.0)
    _, E, I = simulate(p, E0=0.05, I0=0.05, t_max=200)
    # Final tail must be close to zero
    assert E[-100:].mean() < 0.2
    assert I[-100:].mean() < 0.2


if __name__ == "__main__":
    # Allow execution without pytest as a smoke test
    test_sigmoid_bounds()
    test_sigmoid_threshold()
    test_rhs_returns_two_values()
    test_simulation_runs()
    test_activities_in_unit_interval()
    test_low_input_settles_to_low_activity()
    print("All tests passed.")
