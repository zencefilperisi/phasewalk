"""
Equivalence test for the batched time_averaged_distribution_all optimization.

The batched function must return the SAME numbers (to floating-point precision)
as looping the original per-node time_averaged_distribution over all start nodes.
This guards the ADR invariant that the eigendecomposition optimization is a speed
change only, never a methodology change.

Run:  py -m pytest tests/test_walks_fast.py -v
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pytest

from quantum_walks.walks import (
    time_averaged_distribution, time_averaged_distribution_all,
    participation_ratio,
)


def _layerP_graph(n=18, seed=0):
    """A symmetric, non-negative, zero-diagonal graph like an FCM Layer-P."""
    rng = np.random.default_rng(seed)
    M = rng.random((n, n))
    A = (M + M.T) / 2.0
    np.fill_diagonal(A, 0.0)
    return A


@pytest.mark.parametrize("quantum", [True, False])
def test_batched_matches_per_node(quantum):
    A = _layerP_graph(n=16, seed=1)
    times = np.linspace(0.1, 100.0, 60)
    alld = time_averaged_distribution_all(A, times, quantum=quantum)
    for s in range(A.shape[0]):
        ref = time_averaged_distribution(A, s, times, quantum=quantum)
        assert np.allclose(alld[:, s], ref, atol=1e-10), f"mismatch at start node {s}"


def test_pr_vector_matches():
    A = _layerP_graph(n=20, seed=2)
    times = np.linspace(0.1, 100.0, 80)
    alld = time_averaged_distribution_all(A, times, quantum=True)
    pr_fast = np.array([participation_ratio(alld[:, s]) for s in range(A.shape[0])])
    pr_ref = np.array([
        participation_ratio(time_averaged_distribution(A, s, times, quantum=True))
        for s in range(A.shape[0])
    ])
    assert np.allclose(pr_fast, pr_ref, atol=1e-9)


def test_quantum_distribution_normalized():
    A = _layerP_graph(n=12, seed=3)
    times = np.linspace(0.1, 50.0, 40)
    alld = time_averaged_distribution_all(A, times, quantum=True)
    # each start node's time-averaged distribution sums to 1
    assert np.allclose(alld.sum(axis=0), 1.0, atol=1e-9)