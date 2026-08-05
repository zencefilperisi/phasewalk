"""
Micro-step 3: Spread probability over the graph (classical diffusion).
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import expm

G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
A = nx.to_numpy_array(G)
L = np.diag(A.sum(axis=1)) - A
n = G.number_of_nodes()

p0 = np.zeros(n)
p0[0] = 1.0
print(f"Start: all probability on node 0. Sum of p0 = {p0.sum()} (should be 1.0)")

times = [0.5, 2.0, 10.0]

def distribution_at(t):
    return expm(-L * t) @ p0

distributions = [distribution_at(t) for t in times]

for t, p in zip(times, distributions):
    print(f"  t = {t:>4}:  sum of p(t) = {p.sum():.4f}   (should stay ~1.0)")

cmap_name = "viridis"

pos = nx.spring_layout(G, seed=42)
vmax = max(d.max() for d in distributions)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for ax, t, p in zip(axes, times, distributions):
    nx.draw(
        G, pos, ax=ax, with_labels=True, node_size=420,
        node_color=p, cmap=plt.get_cmap(cmap_name),
        vmin=0.0, vmax=vmax, edge_color="lightgray", font_size=8,
    )
    ax.set_title(f"Classical diffusion at t = {t}")

sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap_name),
                           norm=plt.Normalize(vmin=0, vmax=vmax))
fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02, label="probability")

fig.suptitle("Probability spreading from node 0 over the small-world graph",
             fontsize=14)
plt.savefig("figures/micro_step_3_graph.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nDone. Figure saved as micro_step_3_diffusion.png")
