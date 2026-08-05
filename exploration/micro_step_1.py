"""
Micro-step 1: Create a small-world network and draw it.
"""
import matplotlib.pyplot as plt
import networkx as nx

G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
print(f"Graph created: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color="skyblue",
        node_size=500, edge_color="gray")

plt.title("Watts-Strogatz small-world graph (n=20, k=4, p=0.3)")
plt.savefig("figures/micro_step_1_graph.png", dpi=150, bbox_inches="tight")
plt.show()
print("Done. A figure was saved as micro_step_1_graph.png")
