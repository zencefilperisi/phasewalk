"""
Micro-step 8: Is the clustering -> spread effect robust, or a fluke?
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from quantum_walks.walks import participation_ratio, time_averaged_distribution


N_GRAPHS = 100
N_NODES = 20
avg_times = np.linspace(0.1, 100.0, 300)

rng = np.random.default_rng(0)
correlations = []

print(f"Testing clustering-vs-spread correlation across {N_GRAPHS} graphs...")

for g in range(N_GRAPHS):
    seed = int(rng.integers(0, 10**6))
    G = nx.watts_strogatz_graph(N_NODES, 4, 0.3, seed=seed)
    A = nx.to_numpy_array(G)

    pr = np.array([
        participation_ratio(time_averaged_distribution(A, s, avg_times, quantum=True))
        for s in range(N_NODES)
    ])

    clustering = np.array([nx.clustering(G)[i] for i in range(N_NODES)])

    if clustering.std() == 0:
        continue

    c = np.corrcoef(clustering, pr)[0, 1]
    correlations.append(c)

correlations = np.array(correlations)

print(f"\nValid graphs: {len(correlations)}")
print(f"Mean correlation:      {correlations.mean():+.3f}")
print(f"Std dev:               {correlations.std():.3f}")
print(f"Fraction negative:     {(correlations < 0).mean():.1%}")
print(f"Range:                 [{correlations.min():+.3f}, {correlations.max():+.3f}]")

fig, ax = plt.subplots(figsize=(9, 5.5))
ax.hist(correlations, bins=20, color="#4C72B0",
        edgecolor="white", alpha=0.85)
ax.axvline(0, color="black", lw=1.2, ls="--", label="zero (no effect)")
ax.axvline(correlations.mean(), color="crimson", lw=2,
           label=f"mean = {correlations.mean():+.3f}")
ax.axvline(-0.49, color="green", lw=2, ls=":",
           label="single-graph value (-0.49)")
ax.set_xlabel("correlation(clustering, spread) per graph")
ax.set_ylabel("number of graphs")
ax.set_title(f"Distribution of the clustering-spread correlation "
             f"across {len(correlations)} random graphs")
ax.legend()
plt.tight_layout()
plt.savefig("figures/micro_step_8_graph.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_8_robustness.png")
