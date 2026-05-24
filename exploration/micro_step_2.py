"""
Micro-step 2: Turn the graph into matrices.

Goal: learn how to extract the adjacency matrix and the Laplacian
matrix from a NetworkX graph. These matrices are what the walks
(classical and quantum) actually operate on.

Fill in the two TODO lines. Everything else is done for you.
Run with:  python micro_step_2.py
"""
import numpy as np
import networkx as nx


# Recreate the same graph as in micro-step 1 (this part is done for you).
G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")


# --------------------------------------------------------------------------- #
# TODO 1: Get the adjacency matrix as a NumPy array.
#
# Use nx.to_numpy_array(G).
#   - It takes the graph G and returns a 20x20 NumPy array.
#   - Entry [i, j] is 1.0 if nodes i and j are connected, else 0.0.
#   - NOTE: no `seed` argument here! The matrix is fully determined by
#     the graph; there is nothing random to seed. (This is the mistake
#     from last time — now you know why it failed.)
#
# Assign the result to a variable called A.
#
A = None   # <-- TODO 1
# --------------------------------------------------------------------------- #


if A is None:
    raise SystemExit("TODO 1 not done: A is still None.")

print("\nAdjacency matrix A:")
print(f"  shape: {A.shape}        (should be (20, 20))")
print(f"  A is symmetric: {np.allclose(A, A.T)}   (should be True)")
print(f"  total number of 1s: {int(A.sum())}   "
      f"(should be 2 x edges = {2 * G.number_of_edges()})")


# --------------------------------------------------------------------------- #
# TODO 2: Build the Laplacian matrix L = D - A.
#
# The degree matrix D is a diagonal matrix whose [i, i] entry is the
# number of neighbours of node i. That number equals the sum of row i
# of A, i.e. A[i, :].sum(), or for all rows at once: A.sum(axis=1).
#
# NumPy can build a diagonal matrix from a vector with np.diag(vector).
#
# So:  D = np.diag(A.sum(axis=1))
#      L = D - A
#
# Write those two lines below (replace the single placeholder line):
L = None   # <-- TODO 2  (you may use one or two lines here)
# --------------------------------------------------------------------------- #


if L is None:
    raise SystemExit("TODO 2 not done: L is still None.")

print("\nLaplacian matrix L = D - A:")
print(f"  shape: {L.shape}        (should be (20, 20))")
print(f"  each row of L sums to ~0: {np.allclose(L.sum(axis=1), 0)}   "
      f"(should be True)")
print(f"  diagonal of L = node degrees: {np.diag(L).astype(int)}")

print("\nDone. You now have A (adjacency) and L (Laplacian).")
print("Next micro-step will use L to make probability spread over the graph.")
