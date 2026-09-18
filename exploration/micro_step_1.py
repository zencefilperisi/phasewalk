"""
Micro-step 1: Create a small-world network and draw it.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx

FIG_DIR = Path(__file__).resolve().parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

G = nx.watts_strogatz_graph(n=20, k=4, p=0.3, seed=42)
print(f"Graph created: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")

pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color="skyblue",
        node_size=500, edge_color="gray")

plt.title("Watts-Strogatz small-world graph (n=20, k=4, p=0.3)")
plt.savefig(FIG_DIR / "micro_step_1_graph.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"Done. Figure saved as {FIG_DIR / 'micro_step_1_graph.png'}")
