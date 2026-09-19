"""
Micro-step 14: Classical vs quantum walk on a real seizure network.

Central question of Phase 2C: on the ictal FCM of ID1/Sz1 (Layer P,
per ADR 0008), does a continuous-time quantum walk assign different
importance to electrodes than a classical diffusion does?

Method:
  1. Load ID1/Sz1 and extract the ictal segment.
  2. Compute the Pearson FCM (ADR 0007, primary measure).
  3. Convert to Layer P graph (ADR 0008: diagonal zeroed, negatives clipped).
  4. For each electrode as start node, run the time-averaged
     distribution for classical diffusion and CTQW; take its
     participation ratio (spread).
  5. Compare the two spread rankings across electrodes:
     Spearman rank correlation and a scatter plot.

Output: one figure with two panels
  - left  : per-electrode spread, classical vs quantum, as a scatter
  - right : per-electrode difference (quantum PR - classical PR),
            ordered by electrode index

This is a Layer P (positive-only) analysis. Layer S (signed graph,
CTQW only) is deferred to a later micro-step per ADR 0008.
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
    classical_diffusion,
    quantum_walk,
    participation_ratio,
    time_averaged_distribution,
)


MAT_PATH = r"C:\Users\User\Desktop\SWEC\ID1\Sz1.mat"

# Load and build the Layer P graph.
sz = load_seizure(MAT_PATH)
ictal = sz.segment("ictal")
fcm = compute_fcm(ictal)
A = fcm_to_graph(fcm, layer="positive")

n = A.shape[0]
print(f"Loaded {sz.patient_id}/{sz.seizure_id}")
print(f"  ictal duration: {sz.ictal_duration_s:.1f} s")
print(f"  graph shape: {A.shape} (Layer P)")
print(f"  edge weight range: [{A[A>0].min():.3f}, {A.max():.3f}]")
print(f"  fraction of non-zero off-diagonal entries: "
      f"{(A > 0).sum() / (n * (n - 1)):.2%}")


# Time grid for the time-averaged distribution. The graph has weights
# in [0, 1] like our synthetic ones, so the same grid works.
# Using 300 points 0.1..100 as in earlier micro-steps.
t_grid = np.linspace(0.1, 100.0, 300)


# For each start node, compute time-averaged spread under both walks.
print("\nComputing time-averaged spread per start node...")
pr_classical = np.empty(n)
pr_quantum = np.empty(n)

for s in range(n):
    avg_c = time_averaged_distribution(A, s, t_grid, quantum=False)
    avg_q = time_averaged_distribution(A, s, t_grid, quantum=True)
    pr_classical[s] = participation_ratio(avg_c)
    pr_quantum[s] = participation_ratio(avg_q)

# Sanity: both are per-node spread values in [1, n].
print(f"\nClassical PR: min={pr_classical.min():.2f}, "
      f"max={pr_classical.max():.2f}, mean={pr_classical.mean():.2f}")
print(f"Quantum   PR: min={pr_quantum.min():.2f}, "
      f"max={pr_quantum.max():.2f}, mean={pr_quantum.mean():.2f}")


# Rank agreement between the two rankings.
rho, pval = spearmanr(pr_classical, pr_quantum)
print(f"\nSpearman rank correlation between classical and quantum "
      f"spread rankings: rho = {rho:+.3f}  (p = {pval:.3g})")

# Also Pearson (linear) correlation as a secondary summary.
r_pearson = np.corrcoef(pr_classical, pr_quantum)[0, 1]
print(f"Pearson correlation of the values (not ranks): r = {r_pearson:+.3f}")


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Left: per-electrode scatter, classical vs quantum PR.
ax = axes[0]
ax.scatter(pr_classical, pr_quantum, s=70, color="#4C72B0",
           edgecolor="white", linewidth=1.0, zorder=3)
for i in range(n):
    ax.annotate(str(i), (pr_classical[i], pr_quantum[i]),
                fontsize=6, ha="center", va="center", color="white")

# Diagonal reference line (equal spread).
lo = min(pr_classical.min(), pr_quantum.min())
hi = max(pr_classical.max(), pr_quantum.max())
pad = 0.05 * (hi - lo)
ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
        color="black", lw=1, ls="--", label="equal spread")
ax.set_xlim(lo - pad, hi + pad)
ax.set_ylim(lo - pad, hi + pad)
ax.set_xlabel("classical spread (time-averaged PR)")
ax.set_ylabel("quantum spread (time-averaged PR)")
ax.set_title(f"Per-electrode spread: classical vs quantum\n"
             f"Spearman rho = {rho:+.3f}   Pearson r = {r_pearson:+.3f}")
ax.legend(loc="upper left")
ax.grid(alpha=0.3)

# Right: per-electrode difference, quantum minus classical.
ax = axes[1]
diff = pr_quantum - pr_classical
colors = ["#DD8452" if d > 0 else "#4C72B0" for d in diff]
ax.bar(range(n), diff, color=colors, edgecolor="white", linewidth=0.5)
ax.axhline(0, color="black", lw=1)
ax.set_xlabel("electrode index")
ax.set_ylabel("quantum PR - classical PR")
ax.set_title("Where does quantum spread more (orange) or less (blue) "
             "than classical?")
ax.grid(axis="y", alpha=0.3)

fig.suptitle(f"Classical vs quantum walk on the real ictal network "
             f"({sz.patient_id}/{sz.seizure_id}, Layer P)", fontsize=13)
plt.tight_layout()

out_path = FIG_DIR / "micro_step_14_classical_vs_quantum.png"
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out_path}")