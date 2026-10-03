"""
Micro-step 18: ID3 replication — completing the ADR 0006 cohort.

ID3 is the last replication patient. Its purpose is NOT statistical
power (only 2 seizures => 1 pair), but:
  - completing the ADR 0006 cohort (all three patients examined)
  - checking the pipeline runs on a larger, different implant
    (98 electrodes vs 47/42) with no code changes

ID3 structure:
  - 2 seizures (Sz1 125s, Sz2 73s), both >= 30s
  - 98 electrodes

With 1 pair, we report the single Layer-1 and Layer-2 Spearman rho and
compare qualitatively with ID1 (+0.410) and ID2 (+0.630). We also list
the candidate consistent electrodes (mean rank over the 2 seizures),
though "consistency over 2 seizures" is a weak notion — an electrode
low in both is suggestive, nothing more.

The figure uses a Sz1-vs-Sz2 scatter (each electrode a point) rather
than a histogram, because a 1-pair histogram is meaningless.
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


PATIENT_ID = "ID3"
PATIENT_FOLDER = Path(r"C:\Users\User\Desktop\SWEC") / PATIENT_ID
t_grid = np.linspace(0.1, 100.0, 300)
LAYER2_WINDOW_S = 30.0


def sz_num(p: Path) -> int:
    return int(p.stem[2:])

sz_files = sorted(PATIENT_FOLDER.glob("Sz*.mat"), key=sz_num)
print(f"Found {len(sz_files)} seizure files for {PATIENT_ID}")

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

    fcm_L1 = compute_fcm(ictal_full)
    A_L1 = fcm_to_graph(fcm_L1, layer="positive")
    for s in range(n_electrodes):
        avg_q = time_averaged_distribution(A_L1, s, t_grid, quantum=True)
        pr_L1[si, s] = participation_ratio(avg_q)

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


# --- Analysis (1 pair) --------------------------------------------- #
rho_L1 = spearmanr(pr_L1[0], pr_L1[1])[0]
rho_L2 = spearmanr(pr_L2[0], pr_L2[1])[0]

ranks_L1 = np.empty_like(pr_L1)
ranks_L2 = np.empty_like(pr_L2)
for si in range(n_seizures):
    ranks_L1[si] = pr_L1[si].argsort().argsort() + 1
    ranks_L2[si] = pr_L2[si].argsort().argsort() + 1

mean_rank_L1 = ranks_L1.mean(axis=0)
order_L1 = mean_rank_L1.argsort()

per_sz_L1_vs_L2 = np.array([
    spearmanr(pr_L1[si], pr_L2[si])[0] for si in range(n_seizures)
])

print("\n" + "=" * 70)
print("ID3 single-pair Spearman (2 seizures => 1 pair):")
print("=" * 70)
print(f"  Layer 1 (full ictal):  rho = {rho_L1:+.3f}")
print(f"  Layer 2 (first 30s):   rho = {rho_L2:+.3f}")
print(f"  Layer 2 - Layer 1    = {rho_L2 - rho_L1:+.3f}")
print(f"\n  Cross-patient Layer 1 comparison:")
print(f"    ID1 (9 long sz):  +0.410")
print(f"    ID2 (4 sz):       +0.630")
print(f"    ID3 (2 sz):       {rho_L1:+.3f}")

print(f"\n  Onset-heterogeneity check (L2 - L1):")
print(f"    ID1:  -0.169 (onset less consistent)")
print(f"    ID2:  +0.004 (no difference)")
print(f"    ID3:  {rho_L2 - rho_L1:+.3f}")

print(f"\nID3 per-seizure L1 vs L2 Spearman:")
for si in range(n_seizures):
    print(f"  {sz_labels[si]:<8}  rho = {per_sz_L1_vs_L2[si]:+.3f}  ({sz_durations[si]:5.0f}s)")


# Candidate electrodes (low in both seizures)
print("\n" + "=" * 70)
print("ID3 top-5 electrodes by mean rank (localized end) — weak, 2 seizures:")
print("=" * 70)
for i, idx in enumerate(order_L1[:5], 1):
    r = ranks_L1[:, idx].astype(int)
    print(f"  #{i}: electrode {idx:>3}  mean rank {mean_rank_L1[idx]:5.1f}  "
          f"(Sz1 rank {r[0]}, Sz2 rank {r[1]})")

print(f"\nID3 top-5 electrodes by mean rank (spread end):")
for i, idx in enumerate(order_L1[-5:][::-1], 1):
    r = ranks_L1[:, idx].astype(int)
    print(f"  #{i}: electrode {idx:>3}  mean rank {mean_rank_L1[idx]:5.1f}  "
          f"(Sz1 rank {r[0]}, Sz2 rank {r[1]})")


npz_path = FIG_DIR / "micro_step_18_ID3.npz"
np.savez(
    npz_path,
    pr_L1=pr_L1, pr_L2=pr_L2, ranks_L1=ranks_L1, ranks_L2=ranks_L2,
    sz_labels=np.array(sz_labels), sz_durations=np.array(sz_durations),
    rho_L1=rho_L1, rho_L2=rho_L2, per_sz_L1_vs_L2=per_sz_L1_vs_L2,
    t_grid=t_grid, patient_id=PATIENT_ID, n_electrodes=n_electrodes,
)
print(f"\nData saved to {npz_path}")


# --- Figure --------------------------------------------------------- #
fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                          gridspec_kw={"width_ratios": [1.3, 1.0]})

# Panel 1: heatmap (2 seizures × 98 electrodes, sorted by mean rank)
ax = axes[0]
im = ax.imshow(pr_L1[:, order_L1], aspect="auto", cmap="viridis")
ax.set_yticks(range(n_seizures))
ax.set_yticklabels(sz_labels, fontsize=10)
ax.set_xticks(range(0, n_electrodes, 3))
ax.set_xticklabels([str(order_L1[i]) for i in range(0, n_electrodes, 3)],
                   fontsize=6, rotation=90)
ax.set_xlabel("electrode (sorted by mean rank, Layer 1) — every 3rd labelled")
ax.set_ylabel("seizure")
ax.set_title(f"ID3 Layer 1: Quantum PR, {n_seizures} seizures × {n_electrodes} electrodes")
plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="quantum PR")

# Panel 2: Sz1 vs Sz2 scatter (the single pair, honest viz)
ax = axes[1]
ax.scatter(pr_L1[0], pr_L1[1], s=45, color="#4C72B0",
           edgecolor="white", linewidth=0.8, zorder=3, alpha=0.85)
lo = min(pr_L1.min(), pr_L1.min())
hi = max(pr_L1.max(), pr_L1.max())
pad = 0.05 * (hi - lo)
ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
        color="black", lw=1, ls="--", label="equal")
ax.set_xlim(lo - pad, hi + pad)
ax.set_ylim(lo - pad, hi + pad)
ax.set_xlabel(f"{sz_labels[0]} quantum PR")
ax.set_ylabel(f"{sz_labels[1]} quantum PR")
ax.set_title(f"ID3: the single seizure pair (Layer 1)\n"
             f"Spearman rho = {rho_L1:+.3f}")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
ax.set_aspect("equal")

fig.suptitle(f"ID3 replication ({n_seizures} seizures, {n_electrodes} electrodes, Layer P) — "
             f"completes ADR 0006 cohort", fontsize=12)
plt.tight_layout()

out = FIG_DIR / "micro_step_18_ID3.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out}")