"""
Micro-step 8: Is the clustering -> spread effect robust, or a fluke?

Micro-step 7 found, on ONE graph (seed=42), that clustering correlates
negatively (-0.49) with quantum-walk spread. A single number from a
single graph proves nothing. The question now:

    Does this negative correlation hold ACROSS MANY random graphs,
    or was -0.49 just a property of that one particular graph?

We generate many independent Watts-Strogatz graphs, compute the
clustering-vs-PR correlation in each, and look at the DISTRIBUTION of
those correlations. If the effect is real, most graphs should show a
clearly negative correlation. If the correlations scatter around zero,
the single -0.49 was a fluke.

This is a "robustness check" — standard practice before believing any
result from a single sample.

Fill in the TODOs. Run from inside phasewalk/:  python micro_step_8.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from quantum_walks.walks import participation_ratio, time_averaged_distribution


N_GRAPHS = 100        # how many random graphs to test
N_NODES = 20
avg_times = np.linspace(0.1, 100.0, 300)

rng = np.random.default_rng(0)   # reproducible set of graph seeds
correlations = []

print(f"Testing clustering-vs-spread correlation across {N_GRAPHS} graphs...")

for g in range(N_GRAPHS):
    seed = int(rng.integers(0, 10**6))
    G = nx.watts_strogatz_graph(N_NODES, 4, 0.3, seed=seed)
    A = nx.to_numpy_array(G)

    # Spread (time-averaged PR) for each start node (done for you).
    pr = np.array([
        participation_ratio(time_averaged_distribution(A, s, avg_times, quantum=True))
        for s in range(N_NODES)
    ])

    # ------------------------------------------------------------------- #
    # TODO 1: clustering coefficient of each node as an array.
    #   cl = nx.clustering(G)          # {node: value}
    cl = nx.clustering(G)
    #   clustering = np.array([cl[i] for i in range(N_NODES)])
    clustering = np.array([cl[i] for i in range(N_NODES)])   # <-- TODO 1
    # ------------------------------------------------------------------- #
    if clustering is None:
        raise SystemExit("TODO 1 not done.")

    # Skip degenerate graphs where every node has identical clustering
    # (correlation undefined if a variable has zero variance).
    if clustering.std() == 0:
        continue

    # ------------------------------------------------------------------- #
    # TODO 2: correlation between clustering and pr.
    #   np.corrcoef(clustering, pr) returns a 2x2 matrix; the off-diagonal
    np.corrcoef(clustering, pr)
    #   entry [0, 1] is the correlation we want.
    c = (np.corrcoef(clustering, pr)[0, 1])   # <-- TODO 2  (np.corrcoef(clustering, pr)[0, 1])
    # ------------------------------------------------------------------- #
    if c is None:
        raise SystemExit("TODO 2 not done.")

    correlations.append(c)

correlations = np.array(correlations)


# --------------------------------------------------------------------------- #
# Summary statistics (done for you) — the decisive numbers.
# --------------------------------------------------------------------------- #
print(f"\nValid graphs: {len(correlations)}")
print(f"Mean correlation:      {correlations.mean():+.3f}")
print(f"Std dev:               {correlations.std():.3f}")
print(f"Fraction negative:     {(correlations < 0).mean():.1%}")
print(f"Range:                 [{correlations.min():+.3f}, {correlations.max():+.3f}]")


# --------------------------------------------------------------------------- #
# Histogram of correlations across graphs (done for you).
# --------------------------------------------------------------------------- #
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
plt.savefig("micro_step_8_robustness.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_8_robustness.png")
print("\nThe decisive question:")
print("  - If the histogram sits clearly LEFT of zero (mean strongly")
print("    negative, most graphs negative): the effect is ROBUST.")
print("  - If the histogram straddles zero (mean near 0, correlations")
print("    scatter both ways): the single -0.49 was NOT representative.")
