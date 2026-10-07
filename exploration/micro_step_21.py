"""
Micro-step 21: Layer S (signed graph) robustness check for CTQW.

ADR 0008 committed to a signed-graph robustness check: Layer P (primary)
clips negative FCM entries to zero; Layer S keeps them. The question:
does discarding negative correlations change the CTQW per-electrode
ranking? If Layer S ~= Layer P, discarding negatives is harmless and
Layer P stands. If they diverge, negative correlations carry structure
the CTQW ranking responds to — a finding, and a reason to reconsider
Layer S for Phase 3.

Layer S is CTQW-only: classical diffusion is undefined on a signed
graph (ADR 0008), so there is no classical arm here. We compare Layer S
CTQW rankings against Layer P CTQW rankings (loaded from the saved npz).

Two measurements per patient:
  (1) Per-seizure agreement: Spearman(Layer P ranking, Layer S ranking)
      for each seizure — does keeping negatives change the within-seizure
      ranking?
  (2) Within-cohort consistency of Layer S (mean pairwise Spearman)
      compared with Layer P's (+0.410 / +0.630 / +0.535) — is Layer S
      also reproducible?

Interpretation set BEFORE computing:
  - high per-seizure agreement (>0.7): Layer P robust, negatives don't
    change CTQW much.
  - low/moderate agreement: negatives carry CTQW-relevant structure.
  - Layer S consistency ~ Layer P consistency: both stable reps.

Seizure order taken from each npz's sz_labels. Layer S CTQW is computed
fresh (walks re-run); Layer P loaded from npz.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
FIG_DIR = Path(__file__).resolve().parent / "figures"

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from data.loader import load_seizure
from data.connectivity import compute_fcm, fcm_to_graph
from quantum_walks.walks import participation_ratio, time_averaged_distribution


SWEC = Path(r"C:\Users\User\Desktop\SWEC")
t_grid = np.linspace(0.1, 100.0, 300)

patients = [
    ("ID1", "micro_step_16_L1_vs_L2.npz", +0.410),
    ("ID2", "micro_step_17_ID2.npz",      +0.630),
    ("ID3", "micro_step_18_ID3.npz",      +0.535),
]


def mean_pairwise_spearman(value_matrix):
    n = value_matrix.shape[0]
    if n < 2:
        return np.nan, np.array([])
    rhos = []
    for i in range(n):
        for j in range(i + 1, n):
            rhos.append(spearmanr(value_matrix[i], value_matrix[j])[0])
    return np.mean(rhos), np.array(rhos)


results = {}
t_start = time.time()
for pid, npz_name, lp_consistency in patients:
    d = np.load(FIG_DIR / npz_name, allow_pickle=True)
    pr_LP = d["pr_L1"]
    sz_labels = [str(x) for x in d["sz_labels"]]
    n_seizures, n_electrodes = pr_LP.shape

    pr_LS = np.empty((n_seizures, n_electrodes))
    for si, label in enumerate(sz_labels):
        sz = load_seizure(str(SWEC / pid / f"{label}.mat"))
        ictal = sz.segment("ictal")
        fcm = compute_fcm(ictal)
        A_S = fcm_to_graph(fcm, layer="signed")   # keep negatives
        for s in range(n_electrodes):
            avg_q = time_averaged_distribution(A_S, s, t_grid, quantum=True)
            pr_LS[si, s] = participation_ratio(avg_q)
        print(f"  {pid}/{label}: Layer S done "
              f"(elapsed {(time.time()-t_start)/60:.1f}m)")

    # (1) Per-seizure agreement between Layer P and Layer S rankings
    per_sz = np.array([
        spearmanr(pr_LP[si], pr_LS[si])[0] for si in range(n_seizures)
    ])

    # (2) Within-cohort consistency of Layer S
    ls_consistency, ls_pairs = mean_pairwise_spearman(pr_LS)

    results[pid] = dict(
        n_seizures=n_seizures, n_electrodes=n_electrodes,
        sz_labels=sz_labels, pr_LP=pr_LP, pr_LS=pr_LS,
        per_sz=per_sz, lp_consistency=lp_consistency,
        ls_consistency=ls_consistency,
    )

    print("=" * 70)
    print(f"{pid}: {n_seizures} seizures, {n_electrodes} electrodes")
    print("=" * 70)
    print("  (1) Per-seizure Layer P vs Layer S ranking agreement:")
    for si, label in enumerate(sz_labels):
        print(f"        {label:<8}  Spearman = {per_sz[si]:+.3f}")
    print(f"        mean = {per_sz.mean():+.3f}")
    print("  (2) Within-cohort consistency:")
    print(f"        Layer P: {lp_consistency:+.3f}")
    print(f"        Layer S: {ls_consistency:+.3f}")
    print()


# Save
npz_out = FIG_DIR / "micro_step_21_layerS.npz"
np.savez(
    npz_out,
    **{f"{pid}_pr_LS": results[pid]["pr_LS"] for pid, _, _ in patients},
    **{f"{pid}_per_sz": results[pid]["per_sz"] for pid, _, _ in patients},
)
print(f"Data saved to {npz_out}")


# --- Figure --------------------------------------------------------- #
fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

# Panel 1: per-seizure P-vs-S agreement, all seizures, colored by patient
ax = axes[0]
colors = {"ID1": "#4C72B0", "ID2": "#DD8452", "ID3": "#55A467"}
xpos = 0
xticks, xticklabels = [], []
for pid, _, _ in patients:
    per_sz = results[pid]["per_sz"]
    labels = results[pid]["sz_labels"]
    for si in range(len(per_sz)):
        ax.bar(xpos, per_sz[si], color=colors[pid], edgecolor="white", linewidth=0.5)
        xticks.append(xpos)
        xticklabels.append(f"{pid}/{labels[si]}")
        xpos += 1
    xpos += 0.8   # gap between patients
ax.axhline(0.7, color="grey", lw=0.8, ls=":", alpha=0.7, label="0.7 (high agreement)")
ax.axhline(0, color="black", lw=1)
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, fontsize=6, rotation=90)
ax.set_ylabel("Spearman (Layer P vs Layer S ranking)")
ax.set_title("(1) Per-seizure: does keeping negatives change the CTQW ranking?")
ax.set_ylim(-0.1, 1.05)
ax.legend(fontsize=8)
ax.grid(axis="y", alpha=0.3)

# Panel 2: within-cohort consistency, Layer P vs Layer S
ax = axes[1]
pids = [p for p, _, _ in patients]
x = np.arange(len(pids))
w = 0.35
ax.bar(x - w/2, [results[p]["lp_consistency"] for p in pids], w,
       label="Layer P", color="#4C72B0")
ax.bar(x + w/2, [results[p]["ls_consistency"] for p in pids], w,
       label="Layer S", color="#C44E52")
ax.set_xticks(x)
ax.set_xticklabels(pids)
ax.set_ylabel("within-cohort consistency (mean pairwise Spearman)")
ax.set_title("(2) Is Layer S as reproducible as Layer P?")
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)

fig.suptitle("Layer S (signed) robustness check for CTQW — ADR 0008 commitment",
             fontsize=12)
plt.tight_layout()
out = FIG_DIR / "micro_step_21_layerS.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out}")

# Verdict
print("\n" + "=" * 70)
print("READING:")
print("=" * 70)
for pid, _, _ in patients:
    r = results[pid]
    m = r["per_sz"].mean()
    verdict = ("Layer P robust (negatives barely change ranking)" if m > 0.7
               else "negatives DO change ranking (moderate)" if m > 0.4
               else "negatives strongly change ranking")
    print(f"  {pid}: per-seizure P-vs-S mean {m:+.3f} -> {verdict}")
    print(f"        Layer S consistency {r['ls_consistency']:+.3f} "
          f"vs Layer P {r['lp_consistency']:+.3f}")