"""
Micro-step 9: What really predicts quantum-walk spread?

Micro-steps 7-8 showed that classical network centralities (degree,
betweenness, eigenvector, clustering) do NOT reliably predict how far a
quantum walk spreads. This script tests a different idea, rooted in the
actual mathematics of the quantum walk.

A continuous-time quantum walk evolves as expm(-i A t). Diagonalizing
the adjacency matrix A = sum_k lambda_k |v_k><v_k|, the walk is just a
sum of modes |v_k>, each oscillating at its own frequency lambda_k.

How widely a start node spreads should therefore depend on how its
initial state distributes over these eigenvectors. We quantify that
with the SPECTRAL PARTICIPATION of the start node:

    weights_k = (component of start node on eigenvector k)^2
              = evecs[start, k]^2
    spectral_participation = 1 / sum_k weights_k^2

(the same participation-ratio formula as before, but over eigenvectors
instead of over nodes).

We compare, across many graphs:
  - how well classical centralities predict spread (R^2), and
  - how well spectral participation predicts spread (R^2).

R^2 is the fraction of variance explained: 0 = useless, 1 = perfect.

Fill in the TODOs. Run from inside phasewalk/:  python micro_step_9.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from quantum_walks.walks import participation_ratio, time_averaged_distribution


N_GRAPHS = 80
N_NODES = 20
avg_times = np.linspace(0.1, 100.0, 300)
rng = np.random.default_rng(0)

# We collect, per graph, the R^2 of each predictor against spread.
r2_degree, r2_clustering, r2_spectral = [], [], []

print(f"Comparing predictors across {N_GRAPHS} graphs...")

for g in range(N_GRAPHS):
    seed = int(rng.integers(0, 10**6))
    G = nx.watts_strogatz_graph(N_NODES, 4, 0.3, seed=seed)
    A = nx.to_numpy_array(G)

    # Target: time-averaged participation ratio (spread) per start node.
    pr = np.array([
        participation_ratio(time_averaged_distribution(A, s, avg_times, quantum=True))
        for s in range(N_NODES)
    ])

    # Classical predictors (done for you).
    degree = np.array([dict(G.degree())[i] for i in range(N_NODES)], float)
    clustering = np.array([nx.clustering(G)[i] for i in range(N_NODES)])

    # ------------------------------------------------------------------- #
    # TODO 1: eigen-decomposition of the adjacency matrix.
    #   np.linalg.eigh(A) returns (eigenvalues, eigenvectors).
    np.linalg.eigh(A)
    #   eigenvectors is a matrix whose COLUMN k is the k-th eigenvector,
    #   so evecs[node, k] is the component of `node` on eigenvector k.
    #
    #   evals, evecs = np.linalg.eigh(A)
    evals, evecs = np.linalg.eigh(A)   # <-- TODO 1
    # ------------------------------------------------------------------- #
    if evecs is None:
        raise SystemExit("TODO 1 not done.")

    # ------------------------------------------------------------------- #
    # TODO 2: spectral participation of each start node.
    #   For start node s:
    #     weights = evecs[s, :] ** 2          # length-n weights on modes
    #     spec = 1.0 / np.sum(weights ** 2)   # participation over modes
    #   Build an array `spectral` with one value per node.
    #
    #   Hint (one clean way):
    #     spectral = np.array([
    #         1.0 / np.sum((evecs[s, :] ** 2) ** 2) for s in range(N_NODES)
    #     ])
    spectral = np.array([
        1.0 / np.sum((evecs[s, :] ** 2) ** 2) for s in range(N_NODES)
    ])   # <-- TODO 2
    # ------------------------------------------------------------------- #
    if spectral is None:
        raise SystemExit("TODO 2 not done.")

    # R^2 = (correlation)^2 for a single predictor.
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


# --------------------------------------------------------------------------- #
# Summary (done for you).
# --------------------------------------------------------------------------- #
print("\nMean R^2 (variance of spread explained), across graphs:")
print(f"  degree                : {np.nanmean(r2_degree):.3f}")
print(f"  clustering            : {np.nanmean(r2_clustering):.3f}")
print(f"  spectral participation: {np.nanmean(r2_spectral):.3f}")
print(f"\nSpectral predictor R^2 > 0.5 in "
      f"{np.nanmean(r2_spectral > 0.5):.0%} of graphs.")


# --------------------------------------------------------------------------- #
# Plot: distribution of R^2 for each predictor (done for you).
# --------------------------------------------------------------------------- #
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
plt.savefig("micro_step_9_predictors.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_9_predictors.png")
print("\nThe question: does the spectral predictor (green) sit far to the")
print("right (high R^2, consistently) compared to the classical ones?")
