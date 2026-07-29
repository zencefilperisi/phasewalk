"""
Micro-step 7: Which network centrality predicts quantum-walk spread?

Micro-step 6 found that a node's DEGREE only weakly predicts how widely
a quantum walk started there spreads (correlation +0.25). The open
question: is there a better predictor among standard centrality measures?

We test four node properties as predictors of the time-averaged
participation ratio (PR):
  - degree         : number of neighbours
  - betweenness    : how often the node lies on shortest paths (a bridge)
  - eigenvector    : how connected the node is to other well-connected nodes
  - clustering     : how interconnected the node's own neighbours are

For each, we compute its correlation with PR across all 20 start nodes.
The measure with the strongest |correlation| is the best predictor.

Fill in the TODOs. Run from inside the phasewalk/ folder:
    python micro_step_7.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from quantum_walks.walks import participation_ratio, time_averaged_distribution


G = nx.watts_strogatz_graph(20, 4, 0.3, seed=42)
A = nx.to_numpy_array(G)
n = A.shape[0]

# Time-averaged participation ratio for each possible start node
# (this is the quantity we are trying to predict).
avg_times = np.linspace(0.1, 100.0, 500)
pr_by_start = np.array([
    participation_ratio(time_averaged_distribution(A, s, avg_times, quantum=True))
    for s in range(n)
])


# --------------------------------------------------------------------------- #
# Compute the four centrality measures. NetworkX returns each as a dict
# {node: value}; we convert to an array indexed by node number.
#
# TODO 1: degree centrality.
#   NetworkX: dict(G.degree()) gives {node: degree}.
#   Complete: degree = np.array([dict(G.degree())[i] for i in range(n)])
degree = np.array([dict(G.degree())[i] for i in range(n)])        # <-- TODO 1
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# TODO 2: betweenness centrality.
#   bc = nx.betweenness_centrality(G)   # returns {node: value}
bc = nx.betweenness_centrality(G)
#   betweenness = np.array([bc[i] for i in range(n)])
betweenness = np.array([bc[i] for i in range(n)])   # <-- TODO 2
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# TODO 3: eigenvector centrality.
ec = nx.eigenvector_centrality(G, max_iter=1000)
#   ec = nx.eigenvector_centrality(G, max_iter=1000)
#   eigenvector = np.array([ec[i] for i in range(n)])
eigenvector = np.array([ec[i] for i in range(n)])   # <-- TODO 3
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# TODO 4: clustering coefficient.
#   cl = nx.clustering(G)
cl = nx.clustering(G)
#   clustering = np.array([cl[i] for i in range(n)])
clustering = np.array([cl[i] for i in range(n)])    # <-- TODO 4
# --------------------------------------------------------------------------- #

for name, vec in [("degree", degree), ("betweenness", betweenness),
                  ("eigenvector", eigenvector), ("clustering", clustering)]:
    if vec is None:
        raise SystemExit(f"TODO not done: {name} is still None.")


# --------------------------------------------------------------------------- #
# Compute and print the correlation of each measure with PR.
# (done for you)
# --------------------------------------------------------------------------- #
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


# --------------------------------------------------------------------------- #
# Plot: four scatter panels, one per measure, PR on the y-axis.
# (done for you)
# --------------------------------------------------------------------------- #
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
plt.savefig("micro_step_7_centrality.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_7_centrality.png")
print("\nThink about: which measure has the strongest relationship?")
print("Is it positive or negative? What would a NEGATIVE correlation mean?")
