"""
Micro-step 6: Resolve the hub-vs-periphery open question.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from quantum_walks.walks import (
    quantum_walk,
    participation_ratio,
    time_averaged_distribution,
)


G = nx.watts_strogatz_graph(20, 4, 0.3, seed=42)
A = nx.to_numpy_array(G)
n = A.shape[0]

degrees = A.sum(axis=1)
hub_node = int(np.argmax(degrees))
periph_node = int(np.argmin(degrees))
print(f"Hub node: {hub_node} (degree {int(degrees[hub_node])})")
print(f"Peripheral node: {periph_node} (degree {int(degrees[periph_node])})")

avg_times = np.linspace(0.1, 100.0, 500)

avg_hub = time_averaged_distribution(A, hub_node, avg_times, quantum=True)
avg_periph = time_averaged_distribution(A, periph_node, avg_times, quantum=True)

print(f"\nTime-averaged distribution sums: "
      f"hub = {avg_hub.sum():.4f}, periph = {avg_periph.sum():.4f}  "
      f"(both should be ~1.0)")

pr_hub = participation_ratio(avg_hub)
pr_periph = participation_ratio(avg_periph)

print(f"\nTime-averaged participation ratio:")
print(f"  hub start (node {hub_node}):        PR = {pr_hub:.3f}")
print(f"  peripheral start (node {periph_node}): PR = {pr_periph:.3f}")
print(f"  difference: {pr_hub - pr_periph:+.3f}")

pr_by_start = []
for start in range(n):
    avg = time_averaged_distribution(A, start, avg_times, quantum=True)
    pr_by_start.append(participation_ratio(avg))
pr_by_start = np.array(pr_by_start)

corr = np.corrcoef(degrees, pr_by_start)[0, 1]
print(f"\nAcross all {n} start nodes:")
print(f"  correlation(start-degree, time-averaged PR) = {corr:+.3f}")

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(degrees, pr_by_start, s=80, color="#4C72B0",
           edgecolor="white", linewidth=1.2, zorder=3)
for i in range(n):
    ax.annotate(str(i), (degrees[i], pr_by_start[i]),
                fontsize=7, ha="center", va="center", color="white")
ax.set_xlabel("degree of starting node")
ax.set_ylabel("time-averaged participation ratio (spread)")
ax.set_title(f"Does the starting node's degree predict quantum-walk spread?\n"
             f"correlation = {corr:+.3f}")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("figures/micro_step_6_graph.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_6_hub_vs_periphery.png")
