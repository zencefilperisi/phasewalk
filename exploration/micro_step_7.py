"""
Micro-step 7: Which network centrality predicts quantum-walk spread?
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from quantum_walks.walks import participation_ratio, time_averaged_distribution


G = nx.watts_strogatz_graph(20, 4, 0.3, seed=42)
A = nx.to_numpy_array(G)
n = A.shape[0]

avg_times = np.linspace(0.1, 100.0, 500)
pr_by_start = np.array([
    participation_ratio(time_averaged_distribution(A, s, avg_times, quantum=True))
    for s in range(n)
])

degree = np.array([dict(G.degree())[i] for i in range(n)])
betweenness = np.array([nx.betweenness_centrality(G)[i] for i in range(n)])
eigenvector = np.array([nx.eigenvector_centrality(G, max_iter=1000)[i] for i in range(n)])
clustering = np.array([nx.clustering(G)[i] for i in range(n)])

for name, vec in [("degree", degree), ("betweenness", betweenness),
                  ("eigenvector", eigenvector), ("clustering", clustering)]:
    if vec is None:
        raise SystemExit(f"{name} is still None.")

measures = {
    "degree": degree,
    "betweenness": betweenness,
    "eigenvector": eigenvector,
    "clustering": clustering,
}

print("Correlation of each centrality measure with time-averaged PR:")
correlations = {}
for name, vec in measures.items():
    c = np.corrcoef(vec, pr_by_start)[0, 1]
    correlations[name] = c
    print(f"  {name:12}: {c:+.3f}")

best = max(correlations, key=lambda k: abs(correlations[k]))
print(f"\nStrongest predictor (by |correlation|): {best} "
      f"({correlations[best]:+.3f})")

fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
for ax, (name, vec) in zip(axes, measures.items()):
    ax.scatter(vec, pr_by_start, s=70, color="#4C72B0",
               edgecolor="white", linewidth=1.0, zorder=3)
    ax.set_xlabel(name)
    ax.set_title(f"{name}\ncorr = {correlations[name]:+.3f}")
    ax.grid(alpha=0.3)
axes[0].set_ylabel("time-averaged participation ratio")

fig.suptitle("Which node property best predicts quantum-walk spread?",
             fontsize=14)
plt.tight_layout()
plt.savefig("figures/micro_step_7_graph.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_7_centrality.png")
