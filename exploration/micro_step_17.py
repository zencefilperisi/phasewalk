"""
Micro-step 17: ID2 replication of micro-steps 15 and 16.

ADR 0006 says ID2/ID3 are replication patients. Micro-steps 15 and 16
established on ID1 that:
  (A) CTQW per-electrode PR rankings have modest cross-seizure
      consistency (mean pairwise Spearman +0.410 on long seizures).
  (B) A candidate patient-signature electrode (28) stays stable across
      both Layer 1 and Layer 2.
  (C) Layer 2 (first 30s) is LESS consistent than Layer 1 — the onset
      network is per-seizure heterogeneous, full ictal more stereotyped.

This script re-runs the same pipeline on ID2 to see whether (A)-(C)
replicate, or whether ID1 was an unusual case.

ID2 structure:
  - 4 seizures (Sz1 96s, Sz2 260s, Sz3 236s, Sz4 301s)
  - 42 electrodes
  - ALL seizures >=30 s, so no short-seizure confounder (unlike ID1)

Statistical note: 4 seizures => only 6 pairs for within-cohort rho.
Confidence intervals are wide. The question we can answer is
"does the signal point the same way as ID1?", not "is the mean rho
significantly different from zero?".
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
FIG_DIR = Path(__file__).resolve().parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from data.loader import load_seizure
from data.connectivity import compute_fcm, fcm_to_graph
from quantum_walks.walks import participation_ratio, time_averaged_distribution


PATIENT_ID = "ID2"
PATIENT_FOLDER = Path(r"C:\Users\User\Desktop\SWEC") / PATIENT_ID
t_grid = np.linspace(0.1, 100.0, 300)
LAYER2_WINDOW_S = 30.0


def sz_num(p: Path) -> int:
    return int(p.stem[2:])

sz_files = sorted(PATIENT_FOLDER.glob("Sz*.mat"), key=sz_num)
print(f"Found {len(sz_files)} seizure files for {PATIENT_ID}")


# Load once to fix electrode count
first_sz = load_seizure(str(sz_files[0]))
n_electrodes = first_sz.eeg.shape[0]
n_seizures = len(sz_files)
window_samples = int(LAYER2_WINDOW_S * first_sz.fs)

print(f"Patient {PATIENT_ID}: {n_seizures} seizures, {n_electrodes} electrodes.\n")


pr_L1 = np.empty((n_seizures, n_electrodes))
pr_L2 = np.empty((n_seizures, n_electrodes))
sz_labels = []
sz_durations = []

t_start = time.time()
for si, sz_path in enumerate(sz_files):
    sz = load_seizure(str(sz_path))
    sz_labels.append(sz.seizure_id)
    sz_durations.append(sz.ictal_duration_s)

    ictal_full = sz.segment("ictal")
    ictal_first30 = ictal_full[:, :window_samples]

    # Layer 1: full ictal
    fcm_L1 = compute_fcm(ictal_full)
    A_L1 = fcm_to_graph(fcm_L1, layer="positive")
    for s in range(n_electrodes):
        avg_q = time_averaged_distribution(A_L1, s, t_grid, quantum=True)
        pr_L1[si, s] = participation_ratio(avg_q)

    # Layer 2: first 30s
    fcm_L2 = compute_fcm(ictal_first30)
    A_L2 = fcm_to_graph(fcm_L2, layer="positive")
    for s in range(n_electrodes):
        avg_q = time_averaged_distribution(A_L2, s, t_grid, quantum=True)
        pr_L2[si, s] = participation_ratio(avg_q)

    elapsed = time.time() - t_start
    print(f"  [{si+1}/{n_seizures}] {sz.seizure_id:<8} ictal {sz.ictal_duration_s:5.0f}s  "
          f"L1 PR [{pr_L1[si].min():5.2f}, {pr_L1[si].max():5.2f}]  "
          f"L2 PR [{pr_L2[si].min():5.2f}, {pr_L2[si].max():5.2f}]  "
          f"| elapsed {elapsed/60:.2f}m")


# --- Analysis ------------------------------------------------------- #

ranks_L1 = np.empty_like(pr_L1)
ranks_L2 = np.empty_like(pr_L2)
for si in range(n_seizures):
    ranks_L1[si] = pr_L1[si].argsort().argsort() + 1
    ranks_L2[si] = pr_L2[si].argsort().argsort() + 1

mean_rank_L1 = ranks_L1.mean(axis=0)
std_rank_L1 = ranks_L1.std(axis=0)
order_L1 = mean_rank_L1.argsort()

# Within-cohort pairs
pairs_L1 = []
pairs_L2 = []
for i in range(n_seizures):
    for j in range(i + 1, n_seizures):
        pairs_L1.append(spearmanr(pr_L1[i], pr_L1[j])[0])
        pairs_L2.append(spearmanr(pr_L2[i], pr_L2[j])[0])
pairs_L1 = np.array(pairs_L1)
pairs_L2 = np.array(pairs_L2)

# L1 vs L2 same-seizure
per_sz_L1_vs_L2 = np.array([
    spearmanr(pr_L1[si], pr_L2[si])[0] for si in range(n_seizures)
])

print("\n" + "=" * 70)
print(f"ID2 Within-Layer 1  ({len(pairs_L1)} pairs, full ictal):")
print("=" * 70)
print(f"  mean Spearman:   {pairs_L1.mean():+.3f}")
print(f"  median:          {np.median(pairs_L1):+.3f}")
print(f"  individual rhos: {np.round(pairs_L1, 3).tolist()}")

print(f"\nID2 Within-Layer 2  ({len(pairs_L2)} pairs, first 30s):")
print("=" * 70)
print(f"  mean Spearman:   {pairs_L2.mean():+.3f}")
print(f"  median:          {np.median(pairs_L2):+.3f}")
print(f"  individual rhos: {np.round(pairs_L2, 3).tolist()}")

delta = pairs_L2.mean() - pairs_L1.mean()
print(f"\n  Layer 2 - Layer 1 = {delta:+.3f}")
print(f"  ID1 reference    = -0.169 (Layer 2 less consistent)")
if delta < -0.05:
    print(f"  --> SAME direction as ID1 (onset less consistent than full ictal)")
elif delta > 0.05:
    print(f"  --> OPPOSITE direction (Layer 2 more consistent — unexpected)")
else:
    print(f"  --> No clear direction, essentially flat")

print(f"\nID2 per-seizure L1 vs L2 Spearman:")
for si in range(n_seizures):
    print(f"  {sz_labels[si]:<8}  rho = {per_sz_L1_vs_L2[si]:+.3f}  ({sz_durations[si]:5.0f}s)")
print(f"  mean: {per_sz_L1_vs_L2.mean():+.3f}  (ID1 reference: +0.478)")


# Top-5 most localized / spread
print("\n" + "=" * 70)
print(f"ID2 top-5 most-consistently LOCALIZED electrodes (Layer 1):")
print("=" * 70)
for i, idx in enumerate(order_L1[:5], 1):
    r = ranks_L1[:, idx].astype(int)
    print(f"  #{i}: electrode {idx:>3}  "
          f"mean rank = {mean_rank_L1[idx]:5.1f} ± {std_rank_L1[idx]:4.1f}  "
          f"(across seizures: {r.tolist()})")

print(f"\nID2 top-5 most-consistently SPREAD electrodes (Layer 1):")
for i, idx in enumerate(order_L1[-5:][::-1], 1):
    r = ranks_L1[:, idx].astype(int)
    print(f"  #{i}: electrode {idx:>3}  "
          f"mean rank = {mean_rank_L1[idx]:5.1f} ± {std_rank_L1[idx]:4.1f}  "
          f"(across seizures: {r.tolist()})")


# Save
npz_path = FIG_DIR / "micro_step_17_ID2.npz"
np.savez(
    npz_path,
    pr_L1=pr_L1, pr_L2=pr_L2,
    ranks_L1=ranks_L1, ranks_L2=ranks_L2,
    sz_labels=np.array(sz_labels),
    sz_durations=np.array(sz_durations),
    pairs_L1=pairs_L1, pairs_L2=pairs_L2,
    per_sz_L1_vs_L2=per_sz_L1_vs_L2,
    t_grid=t_grid,
    patient_id=PATIENT_ID,
    n_electrodes=n_electrodes,
)
print(f"\nData saved to {npz_path}")


# --- Figure --------------------------------------------------------- #
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5),
                          gridspec_kw={"width_ratios": [1.3, 1.0, 1.0]})

# Panel 1: Layer 1 heatmap (4 seizures × 42 electrodes, sorted by mean rank)
ax = axes[0]
im = ax.imshow(pr_L1[:, order_L1], aspect="auto", cmap="viridis")
ax.set_yticks(range(n_seizures))
ax.set_yticklabels(sz_labels, fontsize=9)
ax.set_xticks(range(n_electrodes))
ax.set_xticklabels([str(e) for e in order_L1], fontsize=6, rotation=90)
ax.set_xlabel("electrode (sorted by mean rank, Layer 1)")
ax.set_ylabel("seizure")
ax.set_title(f"ID2 Layer 1: Quantum PR across {n_seizures} seizures × {n_electrodes} electrodes")
plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="quantum PR")

# Panel 2: L1 vs L2 pairwise Spearman comparison
ax = axes[1]
bins = np.linspace(-0.2, 1.0, 15)
ax.hist(pairs_L1, bins=bins, alpha=0.6, color="#4C72B0",
        label=f"Layer 1 (mean {pairs_L1.mean():+.3f})",
        edgecolor="white")
ax.hist(pairs_L2, bins=bins, alpha=0.6, color="#DD8452",
        label=f"Layer 2 (mean {pairs_L2.mean():+.3f})",
        edgecolor="white")
ax.axvline(0, color="black", lw=1, ls="--")
ax.axvline(pairs_L1.mean(), color="#4C72B0", lw=1.5)
ax.axvline(pairs_L2.mean(), color="#DD8452", lw=1.5)
ax.set_xlabel("pairwise Spearman rho")
ax.set_ylabel("count")
ax.set_title(f"ID2 within-cohort consistency\n({len(pairs_L1)} pairs each)")
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)

# Panel 3: Per-seizure L1 vs L2 (SAME seizure, two windows)
ax = axes[2]
xs = np.arange(n_seizures)
colors = ["#DD8452" if r < 0.3 else "#55A467" if r > 0.7 else "#4C72B0"
          for r in per_sz_L1_vs_L2]
ax.bar(xs, per_sz_L1_vs_L2, color=colors, edgecolor="white", linewidth=0.5)
ax.axhline(per_sz_L1_vs_L2.mean(), color="black", lw=1.2, ls="--",
           label=f"mean = {per_sz_L1_vs_L2.mean():+.3f}")
ax.axhline(0.7, color="grey", lw=0.8, ls=":", alpha=0.7)
ax.axhline(0.3, color="grey", lw=0.8, ls=":", alpha=0.7)
ax.set_xticks(xs)
ax.set_xticklabels(sz_labels, fontsize=9, rotation=45)
ax.set_ylabel("Spearman rho (L1 vs L2 ranking, same seizure)")
ax.set_title("Per-seizure: does first 30s match full ictal?\ngreen >0.7, orange <0.3")
ax.set_ylim(-0.1, 1.05)
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)

fig.suptitle(f"ID2 replication of micro-steps 15 & 16  "
             f"({n_seizures} seizures, all >=30s, Layer P)",
             fontsize=12)
plt.tight_layout()

out = FIG_DIR / "micro_step_17_ID2.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out}")