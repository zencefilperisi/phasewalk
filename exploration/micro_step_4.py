"""
Micro-step 4: Continuous-time QUANTUM walk, compared with classical.
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import expm

G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
A = nx.to_numpy_array(G)
L = np.diag(A.sum(axis=1)) - A
n = G.number_of_nodes()

times = [0.5, 2.0, 10.0]
p0 = np.zeros(n)
p0[0] = 1.0

def classical_distribution(t):
    return expm(-L * t) @ p0

classical = [classical_distribution(t) for t in times]

H = A
psi0 = p0.astype(complex)

def quantum_distribution(t):
    psi = expm(-1j * H * t) @ psi0
    return np.abs(psi) ** 2

quantum = [quantum_distribution(t) for t in times]

for t, pc, pq in zip(times, classical, quantum):
    print(f"  t = {t:>4}:  classical sum = {pc.sum():.4f}   "
          f"quantum sum = {pq.sum():.4f}   (both should be ~1.0)")

cmap_name = "viridis"

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
plt.savefig("figures/micro_step_4_graph.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nDone. Figure saved as micro_step_4_comparison.png")
