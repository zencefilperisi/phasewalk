"""
Tests for the connectivity module.

Covers:
  compute_fcm (pearson):
    - Pearson matrix has the right shape and symmetry
    - diagonal is 1.0
    - values lie in [-1, 1]
    - recovers injected correlations
    - non-2-D signal is rejected
    - unknown measure is rejected

  compute_fcm (plv):
    - PLV matrix shape, symmetry
    - diagonal is 1.0
    - values lie in [0, 1]
    - detects phase-locking (constant phase offset -> PLV ~ 1)
    - phase-independent signals -> lower PLV

  fcm_to_graph:
    - diagonal always zeroed (both layers)
    - "positive" clips negatives; "signed" keeps them
    - symmetry preserved
    - input not mutated
    - non-square / non-2-D / unknown-layer rejected

  Integration: signal -> FCM -> graph
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pytest

from data.connectivity import compute_fcm, fcm_to_graph


# --------------------------------------------------------------------------- #
# compute_fcm — pearson
# --------------------------------------------------------------------------- #

def _synthetic_signal(n_channels=8, n_samples=2000, seed=0):
    rng = np.random.default_rng(seed)
    sig = rng.standard_normal((n_channels, n_samples))
    sig[1] = 0.9 * sig[0] + 0.1 * rng.standard_normal(n_samples)
    sig[3] = -0.9 * sig[2] + 0.1 * rng.standard_normal(n_samples)
    return sig


def test_pearson_shape_symmetry_diagonal():
    fcm = compute_fcm(_synthetic_signal(), measure="pearson")
    assert fcm.shape == (8, 8)
    assert np.allclose(fcm, fcm.T)
    assert np.allclose(np.diag(fcm), 1.0)


def test_pearson_values_in_unit_interval():
    fcm = compute_fcm(_synthetic_signal(), measure="pearson")
    assert np.all(fcm >= -1.0 - 1e-9)
    assert np.all(fcm <= +1.0 + 1e-9)


def test_pearson_recovers_injected_correlations():
    fcm = compute_fcm(_synthetic_signal(), measure="pearson")
    assert fcm[0, 1] > 0.8
    assert fcm[2, 3] < -0.8


def test_compute_fcm_rejects_1d():
    with pytest.raises(ValueError, match="2-D"):
        compute_fcm(np.zeros(100))


def test_compute_fcm_rejects_unknown_measure():
    with pytest.raises(ValueError, match="unknown measure"):
        compute_fcm(_synthetic_signal(), measure="coherence")


# --------------------------------------------------------------------------- #
# compute_fcm — plv
# --------------------------------------------------------------------------- #

def _phase_signals(fs=512, duration_s=10.0):
    """Three channels: two phase-locked at 10 Hz, one at 23 Hz (drifting)."""
    t = np.arange(0, duration_s, 1.0 / fs)
    x0 = np.sin(2 * np.pi * 10 * t)
    x1 = np.sin(2 * np.pi * 10 * t + 0.5)      # constant phase offset -> locked
    x2 = np.sin(2 * np.pi * 23 * t)            # different freq -> phase drifts
    return np.vstack([x0, x1, x2])


def test_plv_shape_symmetry_diagonal():
    plv = compute_fcm(_phase_signals(), measure="plv")
    assert plv.shape == (3, 3)
    assert np.allclose(plv, plv.T)
    assert np.allclose(np.diag(plv), 1.0)


def test_plv_values_in_unit_interval():
    plv = compute_fcm(_phase_signals(), measure="plv")
    assert np.all(plv >= -1e-9)
    assert np.all(plv <= 1.0 + 1e-9)


def test_plv_detects_phase_locking():
    plv = compute_fcm(_phase_signals(), measure="plv")
    # channels 0 and 1 share a constant phase offset -> near-perfect locking
    assert plv[0, 1] > 0.95


def test_plv_lower_for_independent_phases():
    plv = compute_fcm(_phase_signals(), measure="plv")
    # 0 vs 2 (different frequencies) is far less locked than 0 vs 1
    assert plv[0, 2] < plv[0, 1]
    assert plv[0, 2] < 0.5


def test_plv_is_real():
    plv = compute_fcm(_phase_signals(), measure="plv")
    assert np.isrealobj(plv)


# --------------------------------------------------------------------------- #
# fcm_to_graph — diagonal
# --------------------------------------------------------------------------- #

def _handmade_fcm():
    return np.array([
        [ 1.0,  0.5, -0.3,  0.0],
        [ 0.5,  1.0,  0.2, -0.7],
        [-0.3,  0.2,  1.0,  0.4],
        [ 0.0, -0.7,  0.4,  1.0],
    ])


def test_positive_layer_zeroes_diagonal():
    g = fcm_to_graph(_handmade_fcm(), layer="positive")
    assert np.all(np.diag(g) == 0.0)


def test_signed_layer_zeroes_diagonal():
    g = fcm_to_graph(_handmade_fcm(), layer="signed")
    assert np.all(np.diag(g) == 0.0)


# --------------------------------------------------------------------------- #
# fcm_to_graph — negatives handling
# --------------------------------------------------------------------------- #

def test_positive_layer_clips_negatives():
    g = fcm_to_graph(_handmade_fcm(), layer="positive")
    assert np.all(g >= 0.0)
    assert g[0, 1] == 0.5
    assert g[2, 3] == 0.4
    assert g[0, 2] == 0.0
    assert g[1, 3] == 0.0


def test_signed_layer_preserves_negatives():
    g = fcm_to_graph(_handmade_fcm(), layer="signed")
    assert g[0, 1] == 0.5
    assert g[0, 2] == -0.3
    assert g[1, 3] == -0.7
    assert g[2, 3] == 0.4


# --------------------------------------------------------------------------- #
# fcm_to_graph — invariants
# --------------------------------------------------------------------------- #

def test_symmetry_preserved_positive():
    g = fcm_to_graph(_handmade_fcm(), layer="positive")
    assert np.allclose(g, g.T)


def test_symmetry_preserved_signed():
    g = fcm_to_graph(_handmade_fcm(), layer="signed")
    assert np.allclose(g, g.T)


def test_input_fcm_not_mutated():
    fcm = _handmade_fcm()
    fcm_before = fcm.copy()
    _ = fcm_to_graph(fcm, layer="positive")
    _ = fcm_to_graph(fcm, layer="signed")
    assert np.array_equal(fcm, fcm_before)


# --------------------------------------------------------------------------- #
# fcm_to_graph — error paths
# --------------------------------------------------------------------------- #

def test_rejects_non_square_fcm():
    with pytest.raises(ValueError, match="square"):
        fcm_to_graph(np.zeros((4, 5)), layer="positive")


def test_rejects_non_2d_fcm():
    with pytest.raises(ValueError, match="square"):
        fcm_to_graph(np.zeros((4, 4, 4)), layer="positive")


def test_rejects_unknown_layer():
    with pytest.raises(ValueError, match="unknown layer"):
        fcm_to_graph(_handmade_fcm(), layer="magnitude")


# --------------------------------------------------------------------------- #
# Integration
# --------------------------------------------------------------------------- #

def test_pipeline_pearson_to_positive_graph():
    sig = _synthetic_signal(n_channels=10, n_samples=1000)
    g = fcm_to_graph(compute_fcm(sig), layer="positive")
    assert g.shape == (10, 10)
    assert np.all(np.diag(g) == 0.0)
    assert np.all(g >= 0.0)
    assert np.allclose(g, g.T)


def test_pipeline_plv_to_positive_graph():
    sig = _phase_signals()
    g = fcm_to_graph(compute_fcm(sig, measure="plv"), layer="positive")
    assert g.shape == (3, 3)
    assert np.all(np.diag(g) == 0.0)
    assert np.all(g >= 0.0)
    assert np.allclose(g, g.T)
