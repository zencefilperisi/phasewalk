"""
Micro-step 5: Understand the quantum walk through three checks.

(A) Probability conservation: does the total probability stay 1?
(B) Starting node matters: hub vs peripheral start.
(C) Return probability: does the walker come back to where it started?

These three checks turn "I saw a nice picture" into "I understand
what the quantum walk does." Fill in the five TODOs.

Run with:  python micro_step_5.py
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

degrees = A.sum(axis=1)
hub_node = int(np.argmax(degrees))          # most-connected node
periph_node = int(np.argmin(degrees))       # least-connected node
print(f"Hub node: {hub_node} (degree {int(degrees[hub_node])})")
print(f"Peripheral node: {periph_node} (degree {int(degrees[periph_node])})")


def quantum_walk(start_node, t):
    """Return the probability distribution of a quantum walk started at
    `start_node`, evaluated at time `t`."""
    psi0 = np.zeros(n, dtype=complex)
    psi0[start_node] = 1.0
    # --------------------------------------------------------------------- #
    # TODO 1: compute psi(t) = expm(-1j * A * t) @ psi0
    #         (A is the adjacency matrix used as Hamiltonian)
    psi = None   # <-- TODO 1
    # --------------------------------------------------------------------- #
    if psi is None:
        raise SystemExit("TODO 1 not done.")
    return np.abs(psi) ** 2


# ----- Check A: probability conservation ---------------------------------- #
print("\n[Check A] Probability conservation over time:")
for t in [0.0, 1.0, 5.0, 20.0, 100.0]:
    p = quantum_walk(start_node=0, t=t)
    # --------------------------------------------------------------------- #
    # TODO 2: compute the total probability (sum of p) and store in `total`
    total = None   # <-- TODO 2  (use p.sum())
    # --------------------------------------------------------------------- #
    if total is None:
        raise SystemExit("TODO 2 not done.")
    print(f"  t = {t:>6}:  total probability = {total:.6f}")


# ----- Check B: hub vs peripheral start ----------------------------------- #
# We measure how "spread out" the distribution is using the participation
# ratio: PR = 1 / sum(p_i^2). PR = 1 means fully localized on one node;
# PR = n means perfectly spread over all n nodes.
def participation_ratio(p):
    # --------------------------------------------------------------------- #
    # TODO 3: return 1.0 / sum of (p**2)
    return None   # <-- TODO 3  (return 1.0 / np.sum(p**2))
    # --------------------------------------------------------------------- #


t_grid = np.linspace(0.1, 30, 200)
pr_hub = [participation_ratio(quantum_walk(hub_node, t)) for t in t_grid]
pr_periph = [participation_ratio(quantum_walk(periph_node, t)) for t in t_grid]

if pr_hub[0] is None:
    raise SystemExit("TODO 3 not done.")


# ----- Check C: return probability ---------------------------------------- #
# Probability of being back at the START node as a function of time.
return_from_0 = []
for t in t_grid:
    p = quantum_walk(start_node=0, t=t)
    # --------------------------------------------------------------------- #
    # TODO 4: append the probability at node 0, i.e. p[0], to the list
    # <-- TODO 4  (return_from_0.append(p[0]))
    # --------------------------------------------------------------------- #

if len(return_from_0) == 0:
    raise SystemExit("TODO 4 not done.")


# ----- Plot Checks B and C ------------------------------------------------ #
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(t_grid, pr_hub, label=f"start at hub (node {hub_node})")
axes[0].plot(t_grid, pr_periph, label=f"start at periphery (node {periph_node})")
axes[0].set_xlabel("time")
axes[0].set_ylabel("participation ratio (spread)")
axes[0].set_title("(B) How spread out: hub vs peripheral start")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].plot(t_grid, return_from_0, color="crimson")
axes[1].set_xlabel("time")
axes[1].set_ylabel("probability at node 0")
axes[1].set_title("(C) Return probability: walker coming back to node 0")
axes[1].grid(alpha=0.3)

# --------------------------------------------------------------------------- #
# TODO 5: give the whole figure a title with fig.suptitle(...)
#         Use any short descriptive string you like, e.g.:
#         "Quantum walk properties on the small-world graph"
# <-- TODO 5
# --------------------------------------------------------------------------- #

plt.tight_layout()
plt.savefig("micro_step_5_properties.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nDone. Figure saved as micro_step_5_properties.png")
