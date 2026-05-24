"""
Micro-step 3: Spread probability over the graph (classical diffusion).

Goal: see your first DYNAMICS. We drop all the probability on node 0
at time 0, then watch it spread across the network as time increases.

The governing equation is the heat / diffusion equation on a graph:

    p(t) = expm(-L * t) @ p0

where L is the Laplacian, p0 is the starting distribution, and
expm is the matrix exponential.

Fill in the three TODO lines. Run with:  python micro_step_3.py
"""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import expm


# Rebuild graph and matrices (done for you, from micro-steps 1-2).
G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
A = nx.to_numpy_array(G)
L = np.diag(A.sum(axis=1)) - A
n = G.number_of_nodes()


# --------------------------------------------------------------------------- #
# TODO 1: Build the starting probability vector p0.
#
# It should be a length-20 vector of zeros, except a 1.0 at node 0.
#   - np.zeros(n) makes a length-n vector of zeros.
#   - Then set the first entry to 1.0.
#
# Two lines:
#     p0 = np.zeros(n)
#     p0[0] = 1.0
#
p0 = None   # <-- TODO 1 (replace with the two lines above)
# --------------------------------------------------------------------------- #


if p0 is None:
    raise SystemExit("TODO 1 not done: p0 is still None.")
print(f"Start: all probability on node 0. Sum of p0 = {p0.sum()} (should be 1.0)")


# The three times at which we observe the spreading.
times = [0.5, 2.0, 10.0]


# --------------------------------------------------------------------------- #
# TODO 2: Compute p(t) = expm(-L * t) @ p0 for one time t.
#
# Write a function that takes a time value and returns the distribution.
# The matrix exponential is expm(-L * t). Matrix-vector multiply with @.
#
# Complete the function body (replace the `return None` line):
#
def distribution_at(t):
    # expm(-L * t) gives the propagator; multiply by p0 with the @ operator
    return None   # <-- TODO 2  (return expm(-L * t) @ p0)
# --------------------------------------------------------------------------- #


# Compute the distribution at each time.
distributions = [distribution_at(t) for t in times]

# Sanity check: probability must stay normalized (sum ~ 1) at all times.
for t, p in zip(times, distributions):
    if p is None:
        raise SystemExit("TODO 2 not done: distribution_at returns None.")
    print(f"  t = {t:>4}:  sum of p(t) = {p.sum():.4f}   (should stay ~1.0)")


# --------------------------------------------------------------------------- #
# TODO 3: Pick the colormap name for the node colours.
#
# Use the string "viridis" (a standard perceptually-uniform colormap).
# Just assign the string to the variable cmap_name.
#
cmap_name = None   # <-- TODO 3  (set to "viridis")
# --------------------------------------------------------------------------- #


if cmap_name is None:
    raise SystemExit("TODO 3 not done: cmap_name is still None.")


# Draw three panels: probability distribution at each time (done for you).
pos = nx.spring_layout(G, seed=42)
vmax = max(d.max() for d in distributions)   # shared colour scale

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
for ax, t, p in zip(axes, times, distributions):
    nx.draw(
        G, pos, ax=ax, with_labels=True, node_size=420,
        node_color=p, cmap=plt.get_cmap(cmap_name),
        vmin=0.0, vmax=vmax, edge_color="lightgray", font_size=8,
    )
    ax.set_title(f"Classical diffusion at t = {t}")

# A single shared colourbar.
sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(cmap_name),
                           norm=plt.Normalize(vmin=0, vmax=vmax))
fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02, label="probability")

fig.suptitle("Probability spreading from node 0 over the small-world graph",
             fontsize=14)
plt.savefig("micro_step_3_diffusion.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nDone. Figure saved as micro_step_3_diffusion.png")
print("Watch how probability concentrated on node 0 spreads outward as t grows.")
