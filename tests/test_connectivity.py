"""
Tests for the connectivity module.

Covers:
  compute_fcm:
    - Pearson matrix has the right shape and symmetry
    - diagonal is 1.0
    - values lie in [-1, 1]
    - non-2-D signal is rejected
    - unknown measure is rejected

  fcm_to_graph:
    - diagonal is always zero after conversion (both layers)
    - "positive" clips negatives; "signed" keeps them
    - symmetry is preserved
    - original FCM is not mutated
    - non-square or non-2-D FCM is rejected
    - unknown layer is rejected
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pytest

from data.connectivity import compute_fcm, fcm_to_graph


# --------------------------------------------------------------------------- #
# compute_fcm
# --------------------------------------------------------------------------- #

def _synthetic_signal(n_channels=8, n_samples=2000, seed=0):
    """A random channel-time signal with a couple of injected correlations."""
    rng = np.random.default_rng(seed)
    sig = rng.standard_normal((n_channels, n_samples))
    # Inject a strong positive correlation between ch0 and ch1
    sig[1] = 0.9 * sig[0] + 0.1 * rng.standard_normal(n_samples)
    # And a strong negative correlation between ch2 and ch3
    sig[3] = -0.9 * sig[2] + 0.1 * rng.standard_normal(n_samples)
    return sig


def test_pearson_shape_symmetry_diagonal():
    sig = _synthetic_signal()
    fcm = compute_fcm(sig, measure="pearson")
    assert fcm.shape == (8, 8)
    assert np.allclose(fcm, fcm.T)
    assert np.allclose(np.diag(fcm), 1.0)


def test_pearson_values_in_unit_interval():
    sig = _synthetic_signal()
    fcm = compute_fcm(sig, measure="pearson")
    assert np.all(fcm >= -1.0 - 1e-9)
    assert np.all(fcm <= +1.0 + 1e-9)


def test_pearson_recovers_injected_correlations():
    sig = _synthetic_signal()
    fcm = compute_fcm(sig, measure="pearson")
    # ch0-ch1 injected as strongly positive (~+0.9)
    assert fcm[0, 1] > 0.8
    # ch2-ch3 injected as strongly negative (~-0.9)
    assert fcm[2, 3] < -0.8


def test_compute_fcm_rejects_1d():
    with pytest.raises(ValueError, match="2-D"):
        compute_fcm(np.zeros(100))


def test_compute_fcm_rejects_unknown_measure():
    sig = _synthetic_signal()
    with pytest.raises(ValueError, match="unknown measure"):
        compute_fcm(sig, measure="plv")


# --------------------------------------------------------------------------- #
# fcm_to_graph — diagonal
# --------------------------------------------------------------------------- #

def _handmade_fcm():
    """A small FCM with a mix of positive, negative, and zero entries."""
    return np.array([
        [ 1.0,  0.5, -0.3,  0.0],
        [ 0.5,  1.0,  0.2, -0.7],
        [-0.3,  0.2,  1.0,  0.4],
        [ 0.0, -0.7,  0.4,  1.0],
    ])


def test_positive_layer_zeroes_diagonal():
    fcm = _handmade_fcm()
    g = fcm_to_graph(fcm, layer="positive")
    assert np.all(np.diag(g) == 0.0)


def test_signed_layer_zeroes_diagonal():
    fcm = _handmade_fcm()
    g = fcm_to_graph(fcm, layer="signed")
    assert np.all(np.diag(g) == 0.0)


# --------------------------------------------------------------------------- #
# fcm_to_graph — negatives handling
# --------------------------------------------------------------------------- #

def test_positive_layer_clips_negatives():
    fcm = _handmade_fcm()
    g = fcm_to_graph(fcm, layer="positive")
    # No negative entries anywhere.
    assert np.all(g >= 0.0)
    # Entries that were positive are preserved (barring the diagonal).
    assert g[0, 1] == 0.5
    assert g[2, 3] == 0.4
    # Entries that were negative are now zero.
    assert g[0, 2] == 0.0
    assert g[1, 3] == 0.0


def test_signed_layer_preserves_negatives():
    fcm = _handmade_fcm()
    g = fcm_to_graph(fcm, layer="signed")
    # Off-diagonal entries pass through unchanged in sign and magnitude.
    assert g[0, 1] == 0.5
    assert g[0, 2] == -0.3
    assert g[1, 3] == -0.7
    assert g[2, 3] == 0.4


# --------------------------------------------------------------------------- #
# fcm_to_graph — invariants
# --------------------------------------------------------------------------- #

def test_symmetry_preserved_positive():
    fcm = _handmade_fcm()
    g = fcm_to_graph(fcm, layer="positive")
    assert np.allclose(g, g.T)


def test_symmetry_preserved_signed():
    fcm = _handmade_fcm()
    g = fcm_to_graph(fcm, layer="signed")
    assert np.allclose(g, g.T)


def test_input_fcm_not_mutated():
    fcm = _handmade_fcm()
    fcm_before = fcm.copy()
    _ = fcm_to_graph(fcm, layer="positive")
    _ = fcm_to_graph(fcm, layer="signed")
    assert np.array_equal(fcm, fcm_before), \
        "fcm_to_graph must not modify its input array"


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
    fcm = _handmade_fcm()
    with pytest.raises(ValueError, match="unknown layer"):
        fcm_to_graph(fcm, layer="magnitude")


# --------------------------------------------------------------------------- #
# Integration: signal -> FCM -> graph
# --------------------------------------------------------------------------- #

def test_pipeline_signal_to_positive_graph():
    """End-to-end: signal to graph without touching intermediate manually."""
    sig = _synthetic_signal(n_channels=10, n_samples=1000)
    fcm = compute_fcm(sig)
    g = fcm_to_graph(fcm, layer="positive")

    assert g.shape == (10, 10)
    assert np.all(np.diag(g) == 0.0)
    assert np.all(g >= 0.0)
    assert np.allclose(g, g.T)


def test_pipeline_signal_to_signed_graph():
    sig = _synthetic_signal(n_channels=10, n_samples=1000)
    fcm = compute_fcm(sig)
    g = fcm_to_graph(fcm, layer="signed")

    assert g.shape == (10, 10)
    assert np.all(np.diag(g) == 0.0)
    assert np.allclose(g, g.T)
    # For this synthetic signal we injected a negative correlation
    # between channels 2 and 3, so signed layer should keep it negative.
    assert g[2, 3] < 0.0
