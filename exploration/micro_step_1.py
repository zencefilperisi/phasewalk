"""
Micro-step 1: Create a small-world network and draw it.

Goal: learn how to build a graph with NetworkX and visualize it.
Nothing else. No walks, no matrices, no dynamics yet.

Fill in the three TODO lines. Everything else is done for you.
Run with:  python micro_step_1.py
"""
import matplotlib.pyplot as plt
import networkx as nx


# --------------------------------------------------------------------------- #
# TODO 1: Create a Watts-Strogatz small-world graph.
#
# Use the function nx.watts_strogatz_graph(n, k, p, seed=...).
#   n    = number of nodes        -> use 20
#   k    = neighbours per node    -> use 4
#   p    = rewiring probability   -> use 0.3
#   seed = random seed            -> use 42 (so it is reproducible)
#
# Assign the result to a variable called G.
#
# Replace the line below:
G = None   # <-- TODO 1
# --------------------------------------------------------------------------- #


# This block just checks you did TODO 1 correctly. Don't change it.
if G is None:
    raise SystemExit("TODO 1 is not done yet: G is still None.")
print(f"Graph created: {G.number_of_nodes()} nodes, "
      f"{G.number_of_edges()} edges.")


# --------------------------------------------------------------------------- #
# TODO 2: Compute a layout (node positions) for drawing.
#
# Use nx.spring_layout(G, seed=42).
# A "layout" is just a dictionary mapping each node to an (x, y) position
# so the drawing looks tidy. The seed makes the picture reproducible.
#
# Assign the result to a variable called pos.
#
# Replace the line below:
pos = None   # <-- TODO 2
# --------------------------------------------------------------------------- #


if pos is None:
    raise SystemExit("TODO 2 is not done yet: pos is still None.")


# --------------------------------------------------------------------------- #
# TODO 3: Draw the graph.
#
# Use nx.draw(G, pos, with_labels=True, node_color="skyblue",
#             node_size=500, edge_color="gray").
# This draws nodes at the positions in `pos`, labels them with their
# numbers, and colours them.
#
# Just write that one line below (no variable assignment needed):
# <-- TODO 3
# --------------------------------------------------------------------------- #


plt.title("Watts-Strogatz small-world graph (n=20, k=4, p=0.3)")
plt.savefig("micro_step_1_graph.png", dpi=150, bbox_inches="tight")
plt.show()
print("Done. A figure was saved as micro_step_1_graph.png")
