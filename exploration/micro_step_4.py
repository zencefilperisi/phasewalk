"""
Micro-step 4: Continuous-time QUANTUM walk, compared with classical.

Goal: see how a quantum walk spreads differently from classical
diffusion on the SAME graph from the SAME starting node.

Classical:  p(t)        = expm(-L * t) @ p0           (probabilities)
Quantum:    psi(t)      = expm(-1j * H * t) @ psi0    (amplitudes)
            p_quantum(t) = |psi(t)|^2                  (probabilities)

We use H = A (the adjacency matrix) as the quantum Hamiltonian.
This is the most common convention in the quantum-walk literature
(Mulken & Blumen 2011).

Fill in the four TODO lines. Run with:  python micro_step_4.py
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import expm


# Rebuild graph and matrices (done for you).
G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
A = nx.to_numpy_array(G)
L = np.diag(A.sum(axis=1)) - A
n = G.number_of_nodes()

times = [0.5, 2.0, 10.0]

# Starting state: all probability/amplitude on node 0.
p0 = np.zeros(n)
p0[0] = 1.0


# ---- Classical diffusion (you already wrote this in micro-step 3) -------- #
def classical_distribution(t):
    return expm(-L * t) @ p0

classical = [classical_distribution(t) for t in times]


# --------------------------------------------------------------------------- #
# TODO 1: Set the quantum Hamiltonian H equal to the adjacency matrix A.
#
# Just one line: H = A
#
H = None   # <-- TODO 1
# --------------------------------------------------------------------------- #

if H is None:
    raise SystemExit("TODO 1 not done: H is still None.")


# --------------------------------------------------------------------------- #
# TODO 2: Set the initial quantum state psi0.
#
# It is the same shape as p0 (all weight on node 0), but it must be a
# COMPLEX vector because amplitudes can be complex.
#
# Use: psi0 = p0.astype(complex)
#
psi0 = None   # <-- TODO 2
# --------------------------------------------------------------------------- #

if psi0 is None:
    raise SystemExit("TODO 2 not done: psi0 is still None.")


# --------------------------------------------------------------------------- #
# TODO 3: Complete the quantum walk function.
#
# The amplitude vector at time t is:   expm(-1j * H * t) @ psi0
#   - 1j is Python's notation for the imaginary unit i.
#   - The probability at each node is the squared magnitude of the
#     amplitude: np.abs(psi)**2
#
# Complete the function (replace the `return None`):
#
def quantum_distribution(t):
    psi = expm(-1j * H * t) @ psi0        # complex amplitude vector
    return None   # <-- TODO 3  (return np.abs(psi)**2)
# --------------------------------------------------------------------------- #

quantum = [quantum_distribution(t) for t in times]

# Sanity check: both classical and quantum probabilities sum to ~1.
for t, pc, pq in zip(times, classical, quantum):
    if pq is None:
        raise SystemExit("TODO 3 not done: quantum_distribution returns None.")
    print(f"  t = {t:>4}:  classical sum = {pc.sum():.4f}   "
          f"quantum sum = {pq.sum():.4f}   (both should be ~1.0)")


# --------------------------------------------------------------------------- #
# TODO 4: Set the colormap name to "viridis" (same as before).
#
cmap_name = None   # <-- TODO 4
# --------------------------------------------------------------------------- #

if cmap_name is None:
    raise SystemExit("TODO 4 not done: cmap_name is still None.")


# Draw a 2x3 grid: top row classical, bottom row quantum (done for you).
pos = nx.spring_layout(G, seed=42)
vmax = max(max(d.max() for d in classical),
           max(d.max() for d in quantum))

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
cmap = plt.get_cmap(cmap_name)

for j, t in enumerate(times):
    nx.draw(G, pos, ax=axes[0, j], with_labels=True, node_size=380,
            node_color=classical[j], cmap=cmap, vmin=0, vmax=vmax,
            edge_color="lightgray", font_size=7)
    axes[0, j].set_title(f"Classical  t = {t}")

    nx.draw(G, pos, ax=axes[1, j], with_labels=True, node_size=380,
            node_color=quantum[j], cmap=cmap, vmin=0, vmax=vmax,
            edge_color="lightgray", font_size=7)
    axes[1, j].set_title(f"Quantum  t = {t}")

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=vmax))
fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02, label="probability")

fig.suptitle("Classical diffusion (top) vs quantum walk (bottom), "
             "same graph, same start (node 0)", fontsize=14)
plt.savefig("micro_step_4_comparison.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nDone. Figure saved as micro_step_4_comparison.png")
print("Compare the two rows: does the quantum walk spread, localize, or "
      "bounce back differently from classical diffusion?")
