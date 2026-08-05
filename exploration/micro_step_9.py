"""
Micro-step 9: What really predicts quantum-walk spread?
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from quantum_walks.walks import participation_ratio, time_averaged_distribution


N_GRAPHS = 80
N_NODES = 20
avg_times = np.linspace(0.1, 100.0, 300)
rng = np.random.default_rng(0)

r2_degree, r2_clustering, r2_spectral = [], [], []

print(f"Comparing predictors across {N_GRAPHS} graphs...")

for g in range(N_GRAPHS):
    seed = int(rng.integers(0, 10**6))
    G = nx.watts_strogatz_graph(N_NODES, 4, 0.3, seed=seed)
    A = nx.to_numpy_array(G)

    pr = np.array([
        participation_ratio(time_averaged_distribution(A, s, avg_times, quantum=True))
        for s in range(N_NODES)
    ])

    degree = np.array([dict(G.degree())[i] for i in range(N_NODES)], float)
    clustering = np.array([nx.clustering(G)[i] for i in range(N_NODES)])

    evals, evecs = np.linalg.eigh(A)

    spectral = np.array([
        1.0 / np.sum((evecs[s, :] ** 2) ** 2) for s in range(N_NODES)
    ])

    def r2(predictor, target):
        if predictor.std() == 0 or target.std() == 0:
            return np.nan
        return np.corrcoef(predictor, target)[0, 1] ** 2

    r2_degree.append(r2(degree, pr))
    r2_clustering.append(r2(clustering, pr))
    r2_spectral.append(r2(spectral, pr))

r2_degree = np.array(r2_degree)
r2_clustering = np.array(r2_clustering)
r2_spectral = np.array(r2_spectral)

print("\nMean R^2 (variance of spread explained), across graphs:")
print(f"  degree                : {np.nanmean(r2_degree):.3f}")
print(f"  clustering            : {np.nanmean(r2_clustering):.3f}")
print(f"  spectral participation: {np.nanmean(r2_spectral):.3f}")
print(f"\nSpectral predictor R^2 > 0.5 in "
      f"{np.nanmean(r2_spectral > 0.5):.0%} of graphs.")

fig, ax = plt.subplots(figsize=(10, 5.5))
bins = np.linspace(0, 1, 21)
ax.hist(r2_degree, bins=bins, alpha=0.6, label="degree", color="#4C72B0")
ax.hist(r2_clustering, bins=bins, alpha=0.6, label="clustering", color="#DD8452")
ax.hist(r2_spectral, bins=bins, alpha=0.6, label="spectral participation",
        color="#55A467")
ax.set_xlabel("$R^2$  (fraction of spread variance explained)")
ax.set_ylabel("number of graphs")
ax.set_title("How well does each predictor explain quantum-walk spread?")
ax.legend()
plt.tight_layout()
plt.savefig("figures/micro_step_9_graph.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_9_predictors.png")
