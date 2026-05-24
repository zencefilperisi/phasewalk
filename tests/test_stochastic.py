"""
Sanity tests for the stochastic Wilson-Cowan implementation.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np

from classical.model import WCParams
from classical.stochastic import (
    StochasticParams,
    ensemble_simulate,
    lag1_autocorr_in_window,
    simulate_stochastic,
    variance_in_window,
)


def test_simulate_stochastic_runs():
    """Stochastic integrator must run without error."""
    p = WCParams(P=0.5)
    noise = StochasticParams(sigma_E=0.05, sigma_I=0.05)
    t, E, I = simulate_stochastic(p, noise, t_max=10.0, dt=0.01, seed=0)
    assert len(t) == len(E) == len(I)


def test_activities_remain_bounded():
    """Clipping must keep activities in [0, 1] even with large noise."""
    p = WCParams(P=1.5)
    noise = StochasticParams(sigma_E=0.5, sigma_I=0.5)  # deliberately large
    _, E, I = simulate_stochastic(p, noise, t_max=20.0, dt=0.01, seed=1)
    assert np.all(E >= 0.0) and np.all(E <= 1.0)
    assert np.all(I >= 0.0) and np.all(I <= 1.0)


def test_seed_reproducibility():
    """Same seed must produce identical trajectory."""
    p = WCParams(P=0.6)
    noise = StochasticParams(sigma_E=0.05, sigma_I=0.05)
    _, E1, _ = simulate_stochastic(p, noise, t_max=5.0, dt=0.01, seed=42)
    _, E2, _ = simulate_stochastic(p, noise, t_max=5.0, dt=0.01, seed=42)
    np.testing.assert_array_equal(E1, E2)


def test_different_seeds_differ():
    """Different seeds must produce different trajectories."""
    p = WCParams(P=0.6)
    noise = StochasticParams(sigma_E=0.05, sigma_I=0.05)
    _, E1, _ = simulate_stochastic(p, noise, t_max=5.0, dt=0.01, seed=1)
    _, E2, _ = simulate_stochastic(p, noise, t_max=5.0, dt=0.01, seed=2)
    assert not np.array_equal(E1, E2)


def test_ensemble_shape():
    """Ensemble output must have the right (n_trials, n_steps) shape."""
    p = WCParams(P=0.6)
    noise = StochasticParams(sigma_E=0.05, sigma_I=0.05)
    n_trials = 5
    t, E_all, I_all = ensemble_simulate(
        p, noise, n_trials=n_trials, t_max=5.0, dt=0.01, seed=0)
    assert E_all.shape == (n_trials, len(t))
    assert I_all.shape == (n_trials, len(t))


def test_ensemble_trajectories_distinct():
    """Different trials of an ensemble must be distinct."""
    p = WCParams(P=0.6)
    noise = StochasticParams(sigma_E=0.05, sigma_I=0.05)
    _, E_all, _ = ensemble_simulate(p, noise, n_trials=3,
                                     t_max=5.0, dt=0.01, seed=0)
    assert not np.array_equal(E_all[0], E_all[1])
    assert not np.array_equal(E_all[0], E_all[2])


def test_zero_noise_matches_deterministic_qualitatively():
    """With sigma=0 the stochastic solver should produce a smooth,
    deterministic-like trajectory that settles near a fixed point at
    low P. (We don't check bit-equality because the integrator is
    Euler-Maruyama vs RK45.)"""
    p = WCParams(P=0.0, Q=0.0)
    noise = StochasticParams(sigma_E=0.0, sigma_I=0.0)
    _, E, I = simulate_stochastic(p, noise, E0=0.05, I0=0.05,
                                   t_max=50.0, dt=0.01, seed=0)
    # Should decay toward low activity
    assert E[-500:].mean() < 0.2
    assert I[-500:].mean() < 0.2


def test_variance_in_window_basic():
    """Variance is 0 on a constant signal."""
    x = np.ones(2000)
    var = variance_in_window(x, window=500)
    valid = var[~np.isnan(var)]
    assert np.allclose(valid, 0.0)


def test_lag1_autocorr_basic():
    """Lag-1 autocorrelation of white noise is near zero;
    of a smooth signal it is near 1."""
    rng = np.random.default_rng(0)
    white = rng.standard_normal(2000)
    smooth = np.linspace(0, 1, 2000) + 0.01 * rng.standard_normal(2000)

    ac_white = lag1_autocorr_in_window(white, window=500)
    ac_smooth = lag1_autocorr_in_window(smooth, window=500)

    assert np.abs(np.nanmean(ac_white)) < 0.15
    assert np.nanmean(ac_smooth) > 0.85


if __name__ == "__main__":
    for name in dir():
        if name.startswith("test_"):
            globals()[name]()
    print("All stochastic tests passed.")
