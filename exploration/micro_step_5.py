"""
Micro-step 5: Understand the quantum walk through three checks.
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import expm

G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
A = nx.to_numpy_array(G)
L = np.diag(A.sum(axis=1)) - A
n = G.number_of_nodes()

degrees = A.sum(axis=1)
hub_node = int(np.argmax(degrees))
periph_node = int(np.argmin(degrees))
print(f"Hub node: {hub_node} (degree {int(degrees[hub_node])})")
print(f"Peripheral node: {periph_node} (degree {int(degrees[periph_node])})")


def quantum_walk(start_node, t):
    psi0 = np.zeros(n, dtype=complex)
    psi0[start_node] = 1.0
    psi = expm(-1j * A * t) @ psi0
    return np.abs(psi) ** 2


print("\n[Check A] Probability conservation over time:")
for t in [0.0, 1.0, 5.0, 20.0, 100.0]:
    p = quantum_walk(start_node=0, t=t)
    total = p.sum()
    print(f"  t = {t:>6}:  total probability = {total:.6f}")


def participation_ratio(p):
    return 1.0 / np.sum(p ** 2)


t_grid = np.linspace(0.1, 30, 200)
pr_hub = [participation_ratio(quantum_walk(hub_node, t)) for t in t_grid]
pr_periph = [participation_ratio(quantum_walk(periph_node, t)) for t in t_grid]


return_from_0 = []
for t in t_grid:
    p = quantum_walk(start_node=0, t=t)
    return_from_0.append(p[0])


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

fig.suptitle("Quantum walk properties on the small-world graph")

plt.tight_layout()
plt.savefig("figures/micro_step_5_graph.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nDone. Figure saved as micro_step_5_properties.png")
