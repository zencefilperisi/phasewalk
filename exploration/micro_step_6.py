"""
Micro-step 6: Resolve the hub-vs-periphery open question.

In micro-step 5 the INSTANTANEOUS participation ratio oscillated so
wildly that we could not tell whether starting a quantum walk at a hub
node spreads differently from starting at a peripheral node.

The fix: compare TIME-AVERAGED quantities, which wash out the
oscillations and reveal the underlying long-run behaviour.

This script imports the tested functions from the quantum_walks module
(the production code), rather than re-implementing them. That is how
real analysis code is organized: the module holds the logic, the
analysis script just calls it.

Fill in the TODOs. Run with:  python micro_step_6.py
(Run it from inside the phasewalk/ project folder so the import works.)
"""
import sys
from pathlib import Path

# Make the project importable. This assumes the script sits in the
# project root (phasewalk/). If you put it elsewhere, adjust ROOT.
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# Import the PRODUCTION functions (not re-implementing them here).
from quantum_walks.walks import (
    quantum_walk,
    participation_ratio,
    time_averaged_distribution,
)


# Same synthetic graph as before, so results are comparable.
G = nx.watts_strogatz_graph(20, 4, 0.3, seed=42)
A = nx.to_numpy_array(G)
n = A.shape[0]

degrees = A.sum(axis=1)
hub_node = int(np.argmax(degrees))
periph_node = int(np.argmin(degrees))
print(f"Hub node: {hub_node} (degree {int(degrees[hub_node])})")
print(f"Peripheral node: {periph_node} (degree {int(degrees[periph_node])})")


# We average over a long, dense set of times to get a stable estimate.
avg_times = np.linspace(0.1, 100.0, 500)


# --------------------------------------------------------------------------- #
# TODO 1: Compute the time-averaged distribution for a hub start.
#
# Use time_averaged_distribution(A, hub_node, avg_times, quantum=True).
#
avg_hub = time_averaged_distribution(A, hub_node, avg_times, quantum=True)   # <-- TODO 1
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# TODO 2: Compute the time-averaged distribution for a peripheral start.
#
avg_periph = time_averaged_distribution(A, periph_node, avg_times, quantum=True)   # <-- TODO 2
# --------------------------------------------------------------------------- #

if avg_hub is None or avg_periph is None:
    raise SystemExit("TODO 1 or 2 not done.")

# Sanity: both are probability distributions.
print(f"\nTime-averaged distribution sums: "
      f"hub = {avg_hub.sum():.4f}, periph = {avg_periph.sum():.4f}  "
      f"(both should be ~1.0)")


# --------------------------------------------------------------------------- #
# TODO 3: Compute the participation ratio of each time-averaged
#         distribution. Use participation_ratio(...).
#
pr_hub = participation_ratio(avg_hub)       # <-- TODO 3a
pr_periph = participation_ratio(avg_periph)    # <-- TODO 3b
# --------------------------------------------------------------------------- #

if pr_hub is None or pr_periph is None:
    raise SystemExit("TODO 3 not done.")

print(f"\nTime-averaged participation ratio:")
print(f"  hub start (node {hub_node}):        PR = {pr_hub:.3f}")
print(f"  peripheral start (node {periph_node}): PR = {pr_periph:.3f}")
print(f"  difference: {pr_hub - pr_periph:+.3f}")


# --------------------------------------------------------------------------- #
# Repeat the analysis for ALL nodes as start nodes, to see whether PR
# correlates with the starting node's degree across the whole graph.
# (This is done for you — it is the real test of the open question.)
# --------------------------------------------------------------------------- #
pr_by_start = []
for start in range(n):
    avg = time_averaged_distribution(A, start, avg_times, quantum=True)
    pr_by_start.append(participation_ratio(avg))
pr_by_start = np.array(pr_by_start)

# Correlation between starting-node degree and resulting spread.
corr = np.corrcoef(degrees, pr_by_start)[0, 1]
print(f"\nAcross all {n} start nodes:")
print(f"  correlation(start-degree, time-averaged PR) = {corr:+.3f}")


# --------------------------------------------------------------------------- #
# Plot: PR vs starting-node degree (the decisive figure).
# --------------------------------------------------------------------------- #
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
plt.savefig("micro_step_6_hub_vs_periphery.png", dpi=150, bbox_inches="tight")
plt.show()

print("\nDone. Figure saved as micro_step_6_hub_vs_periphery.png")
print("\nInterpretation guide:")
print("  - If correlation is strongly POSITIVE: higher-degree starts spread")
print("    more -> degree matters (instantaneous noise was hiding it).")
print("  - If correlation is near ZERO: starting degree does NOT predict")
print("    spread -> a genuine quantum effect worth reporting.")
