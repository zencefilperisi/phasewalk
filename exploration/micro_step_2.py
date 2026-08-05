"""
Micro-step 2: Turn the graph into matrices.
"""
import numpy as np
import networkx as nx

G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
print(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

A = nx.to_numpy_array(G)

print("\nAdjacency matrix A:")
print(f"  shape: {A.shape}        (should be (20, 20))")
print(f"  A is symmetric: {np.allclose(A, A.T)}   (should be True)")
print(f"  total number of 1s: {int(A.sum())}   "
      f"(should be 2 x edges = {2 * G.number_of_edges()})")

D = np.diag(A.sum(axis=1))
L = D - A

print("\nLaplacian matrix L = D - A:")
print(f"  shape: {L.shape}        (should be (20, 20))")
print(f"  each row of L sums to ~0: {np.allclose(L.sum(axis=1), 0)}   "
      f"(should be True)")
print(f"  diagonal of L = node degrees: {np.diag(L).astype(int)}")

print("\nDone. You now have A (adjacency) and L (Laplacian).")
