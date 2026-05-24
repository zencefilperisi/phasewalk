"""
Sanity tests for the quantum_walks module.

These tests encode the physical laws the walks must obey:
  - probability conservation (both classical and quantum),
  - correct initial condition,
  - Laplacian row-sum-zero property,
  - qualitative difference between quantum and classical return behaviour.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import networkx as nx

from quantum_walks.walks import (
    classical_diffusion,
    laplacian,
    localized_state,
    participation_ratio,
    quantum_walk,
    return_probability,
    time_averaged_distribution,
)


def _test_graph(n=20, k=4, p=0.3, seed=42):
    G = nx.watts_strogatz_graph(n, k, p, seed=seed)
    return nx.to_numpy_array(G)


def test_laplacian_rows_sum_to_zero():
    A = _test_graph()
    L = laplacian(A)
    assert np.allclose(L.sum(axis=1), 0.0)


def test_localized_state_real():
    s = localized_state(20, 3, complex_dtype=False)
    assert s.shape == (20,)
    assert s[3] == 1.0
    assert s.sum() == 1.0
    assert not np.iscomplexobj(s)


def test_localized_state_complex():
    s = localized_state(20, 3, complex_dtype=True)
    assert np.iscomplexobj(s)
    assert s[3] == 1.0


def test_classical_diffusion_conserves_probability():
    A = _test_graph()
    for t in [0.5, 2.0, 10.0, 50.0]:
        p = classical_diffusion(A, start_node=0, t=t)
        assert abs(p.sum() - 1.0) < 1e-9
        assert np.all(p >= -1e-12)


def test_quantum_walk_conserves_probability():
    A = _test_graph()
    for t in [0.5, 2.0, 10.0, 100.0]:
        p = quantum_walk(A, start_node=0, t=t)
        assert abs(p.sum() - 1.0) < 1e-9
        assert np.all(p >= -1e-12)


def test_walks_start_localized_at_t0():
    A = _test_graph()
    pc = classical_diffusion(A, start_node=5, t=0.0)
    pq = quantum_walk(A, start_node=5, t=0.0)
    assert abs(pc[5] - 1.0) < 1e-9
    assert abs(pq[5] - 1.0) < 1e-9


def test_classical_diffusion_spreads_to_uniform():
    A = _test_graph()
    n = A.shape[0]
    p_late = classical_diffusion(A, start_node=0, t=500.0)
    assert np.allclose(p_late, 1.0 / n, atol=1e-3)


def test_quantum_walk_does_not_settle_uniform():
    A = _test_graph()
    n = A.shape[0]
    p_late = quantum_walk(A, start_node=0, t=500.0)
    assert not np.allclose(p_late, 1.0 / n, atol=1e-3)


def test_participation_ratio_bounds():
    n = 20
    loc = np.zeros(n); loc[0] = 1.0
    assert abs(participation_ratio(loc) - 1.0) < 1e-12
    uni = np.full(n, 1.0 / n)
    assert abs(participation_ratio(uni) - n) < 1e-9


def test_hamiltonian_choice_changes_result():
    A = _test_graph()
    p_adj = quantum_walk(A, start_node=0, t=3.0, hamiltonian="adjacency")
    p_lap = quantum_walk(A, start_node=0, t=3.0, hamiltonian="laplacian")
    assert not np.allclose(p_adj, p_lap)


def test_invalid_hamiltonian_raises():
    A = _test_graph()
    try:
        quantum_walk(A, start_node=0, t=1.0, hamiltonian="nonsense")
        assert False, "should have raised ValueError"
    except ValueError:
        pass


def test_quantum_return_probability_oscillates():
    A = _test_graph()
    times = np.linspace(0.1, 30, 200)
    rp = return_probability(A, 0, times, quantum=True)
    later = rp[10:]
    assert later.max() - later.min() > 0.1


def test_classical_return_probability_monotone_ish():
    A = _test_graph()
    times = np.linspace(0.1, 30, 200)
    rp = return_probability(A, 0, times, quantum=False)
    later = rp[10:]
    assert later.max() - later.min() < 0.1


def test_time_averaged_distribution_normalized():
    A = _test_graph()
    times = np.linspace(0.1, 50, 100)
    avg = time_averaged_distribution(A, 0, times, quantum=True)
    assert abs(avg.sum() - 1.0) < 1e-9


if __name__ == "__main__":
    for name in list(globals()):
        if name.startswith("test_"):
            globals()[name]()
    print("All quantum-walk tests passed.")
