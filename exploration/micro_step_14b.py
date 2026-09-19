"""
Micro-step 14b: Reframed CTQW spread analysis on a real seizure network.

After micro-step 14 exposed the saturation issue (classical PR pinned
at ~n regardless of start node), ADR 0009 reframed the Phase 2C
central analysis:
  - CTQW time-averaged PR is now the primary quantity of interest.
  - Classical diffusion is kept only as a sanity baseline (should be
    near n; if it is not, the graph is disconnected).
  - Rankings only, no thresholds.

This script produces three panels on ID1/Sz1 (Layer P):
  1. Quantum PR sorted from lowest to highest, with electrode labels
     visible — identifies the most-localized and most-spread electrodes.
  2. Classical PR across the same electrodes — a flat line near n,
     confirming the reframe motivation.
  3. Quantum PR vs weighted node degree, with Spearman rho reported.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
FIG_DIR = Path(__file__).resolve().parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from data.loader import load_seizure
from data.connectivity import compute_fcm, fcm_to_graph
from quantum_walks.walks import (
    participation_ratio,
    time_averaged_distribution,
)


MAT_PATH = r"C:\Users\User\Desktop\SWEC\ID1\Sz1.mat"

sz = load_seizure(MAT_PATH)
ictal = sz.segment("ictal")
fcm = compute_fcm(ictal)
A = fcm_to_graph(fcm, layer="positive")
n = A.shape[0]

print(f"Loaded {sz.patient_id}/{sz.seizure_id}, Layer P graph {A.shape}")


# Weighted degree: sum of edge weights incident to each node.
# For a weighted graph this is the natural analogue of degree.
weighted_degree = A.sum(axis=1)


# Time-averaged PR per start node.
t_grid = np.linspace(0.1, 100.0, 300)
pr_classical = np.empty(n)
pr_quantum = np.empty(n)

print("Computing time-averaged spread per start node...")
for s in range(n):
    avg_c = time_averaged_distribution(A, s, t_grid, quantum=False)
    avg_q = time_averaged_distribution(A, s, t_grid, quantum=True)
    pr_classical[s] = participation_ratio(avg_c)
    pr_quantum[s] = participation_ratio(avg_q)


# --- Reporting: rankings, not thresholds --------------------------------- #
order_q = np.argsort(pr_quantum)          # low to high
k = 5

lowest = order_q[:k]
highest = order_q[-k:][::-1]

print(f"\nQuantum PR summary:")
print(f"  range: {pr_quantum.min():.2f} to {pr_quantum.max():.2f}  (n = {n})")
print(f"  mean:  {pr_quantum.mean():.2f}")

print(f"\n{k} most LOCALIZED electrodes (lowest quantum PR):")
for rank, idx in enumerate(lowest, start=1):
    print(f"  #{rank}: electrode {idx:>3}  PR = {pr_quantum[idx]:6.2f}  "
          f"(weighted deg = {weighted_degree[idx]:.2f})")

print(f"\n{k} most SPREAD electrodes (highest quantum PR):")
for rank, idx in enumerate(highest, start=1):
    print(f"  #{rank}: electrode {idx:>3}  PR = {pr_quantum[idx]:6.2f}  "
          f"(weighted deg = {weighted_degree[idx]:.2f})")


# Classical baseline: should be near n and near-constant.
print(f"\nClassical PR baseline (sanity check):")
print(f"  range: {pr_classical.min():.2f} to {pr_classical.max():.2f}  "
      f"(should be near n = {n})")
print(f"  standard deviation: {pr_classical.std():.4f}  "
      f"(should be near 0 for a connected graph)")


# Quantum PR vs weighted degree.
rho_deg, pval_deg = spearmanr(weighted_degree, pr_quantum)
print(f"\nSpearman rank correlation between weighted degree and "
      f"quantum PR: rho = {rho_deg:+.3f}  (p = {pval_deg:.3g})")


# --- Figure --------------------------------------------------------------- #
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Panel 1: quantum PR sorted from lowest to highest.
ax = axes[0]
sorted_pr = pr_quantum[order_q]
positions = np.arange(n)
colors = ["#4C72B0"] * n
for idx in lowest:
    colors[np.where(order_q == idx)[0][0]] = "#DD8452"  # highlight lowest
for idx in highest:
    colors[np.where(order_q == idx)[0][0]] = "#55A467"  # highlight highest
ax.bar(positions, sorted_pr, color=colors, edgecolor="white", linewidth=0.4)
ax.set_xticks(positions)
ax.set_xticklabels([str(i) for i in order_q], rotation=90, fontsize=7)
ax.set_xlabel("electrode (sorted by quantum PR)")
ax.set_ylabel("quantum time-averaged PR")
ax.set_title(f"Quantum spread per electrode (sorted)\n"
             f"orange = lowest {k}, green = highest {k}")
ax.grid(axis="y", alpha=0.3)

# Panel 2: classical baseline in the same electrode order.
ax = axes[1]
ax.bar(positions, pr_classical[order_q], color="#8899AA",
       edgecolor="white", linewidth=0.4)
ax.axhline(n, color="black", lw=1, ls="--",
           label=f"ceiling = n = {n}")
ax.set_xticks(positions)
ax.set_xticklabels([str(i) for i in order_q], rotation=90, fontsize=7)
ax.set_xlabel("electrode (same order as panel 1)")
ax.set_ylabel("classical time-averaged PR")
ax.set_title("Classical baseline (should be flat near n)")
ax.legend(loc="lower right")
ax.grid(axis="y", alpha=0.3)

# Panel 3: quantum PR vs weighted degree.
ax = axes[2]
ax.scatter(weighted_degree, pr_quantum, s=70, color="#4C72B0",
           edgecolor="white", linewidth=1.0, zorder=3)
for i in range(n):
    ax.annotate(str(i), (weighted_degree[i], pr_quantum[i]),
                fontsize=6, ha="center", va="center", color="white")
ax.set_xlabel("weighted degree of start electrode")
ax.set_ylabel("quantum time-averaged PR")
ax.set_title(f"Quantum spread vs weighted degree\n"
             f"Spearman rho = {rho_deg:+.3f}")
ax.grid(alpha=0.3)

fig.suptitle(f"CTQW spread analysis on real ictal network "
             f"({sz.patient_id}/{sz.seizure_id}, Layer P)  |  "
             f"reframed per ADR 0009", fontsize=12)
plt.tight_layout()
out = FIG_DIR / "micro_step_14b_ctqw_spread.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out}")