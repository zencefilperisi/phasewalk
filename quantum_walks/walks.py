"""
Continuous-time classical and quantum walks on graphs.

This module packages the dynamics explored interactively in the
micro-step scripts into reusable, tested functions.

Classical diffusion (heat equation on the graph):
    p(t) = expm(-L t) p(0)

Continuous-time quantum walk (Schrodinger equation, H = A):
    psi(t)        = expm(-i A t) psi(0)
    p_quantum(t)  = |psi(t)|^2

Hamiltonian convention: H = A (adjacency), the most common choice in
the quantum-walk literature (Mulken & Blumen, 2011). H = L (Laplacian)
is also used (Childs, 2010); switch via the `hamiltonian` argument.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm


def laplacian(adjacency):
    """Combinatorial graph Laplacian L = D - A. Each row sums to zero."""
    degree = np.diag(adjacency.sum(axis=1))
    return degree - adjacency


def localized_state(n, node, complex_dtype=False):
    """A state with all weight on one node.

    complex_dtype=True returns a complex vector (quantum walks);
    otherwise a real vector (classical diffusion).
    """
    state = np.zeros(n, dtype=complex if complex_dtype else float)
    state[node] = 1.0
    return state


def classical_diffusion(adjacency, start_node, t):
    """Probability distribution of classical diffusion at time t.

    Solves p(t) = expm(-L t) p0 with p0 localized on start_node.
    Returns a length-n probability vector.
    """
    n = adjacency.shape[0]
    L = laplacian(adjacency)
    p0 = localized_state(n, start_node, complex_dtype=False)
    return expm(-L * t) @ p0


def quantum_walk(adjacency, start_node, t, hamiltonian="adjacency"):
    """Probability distribution of a continuous-time quantum walk at time t.

    Evolves psi(t) = expm(-i H t) psi0 then returns |psi(t)|^2.
    hamiltonian: "adjacency" (default) or "laplacian".
    """
    n = adjacency.shape[0]
    if hamiltonian == "adjacency":
        H = adjacency
    elif hamiltonian == "laplacian":
        H = laplacian(adjacency)
    else:
        raise ValueError(
            f"hamiltonian must be 'adjacency' or 'laplacian', got {hamiltonian!r}")

    psi0 = localized_state(n, start_node, complex_dtype=True)
    psi_t = expm(-1j * H * t) @ psi0
    return np.abs(psi_t) ** 2


def participation_ratio(distribution):
    """PR = 1 / sum(p_i^2). PR=1 fully localized; PR=n fully uniform."""
    return 1.0 / np.sum(distribution ** 2)


def return_probability(adjacency, start_node, times, quantum=True,
                       hamiltonian="adjacency"):
    """Probability of being at the start node across a set of times.

    Quantum: oscillates (interference). Classical: decays monotonically.
    Returns a length-T array.
    """
    out = np.empty(len(times))
    for k, t in enumerate(times):
        if quantum:
            p = quantum_walk(adjacency, start_node, t, hamiltonian=hamiltonian)
        else:
            p = classical_diffusion(adjacency, start_node, t)
        out[k] = p[start_node]
    return out


def time_averaged_distribution(adjacency, start_node, times, quantum=True,
                               hamiltonian="adjacency"):
    """Average distribution over a set of times.

    For quantum walks the instantaneous distribution oscillates forever,
    so the time average is the physically meaningful long-run occupation.
    This is the principled fix for the noisy instantaneous PR curves.
    """
    n = adjacency.shape[0]
    acc = np.zeros(n)
    for t in times:
        if quantum:
            acc += quantum_walk(adjacency, start_node, t, hamiltonian=hamiltonian)
        else:
            acc += classical_diffusion(adjacency, start_node, t)
    return acc / len(times)
