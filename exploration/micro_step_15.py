"""
Micro-step 15: Cross-seizure consistency of CTQW localization on ID1.

Per ADR 0009's Open Questions: does the ranking of quantum PR change
substantially across different seizures of the same patient? Micro-step
14b showed the pattern on ID1/Sz1; here we test all 13 seizures of ID1.

Interpretation (ADR 0006: ID1 is the development patient):
  - High agreement across seizures  -> CTQW captures patient-level
    structure. "Electrode X is localized" is a patient-level statement.
  - Low agreement                    -> CTQW captures seizure-specific
    structure. "Electrode X is localized in seizure Y" only.
  - Middle                           -> partial consistency; some
    electrodes stable, others not.

Layer P per ADR 0008 (positive-only). t_grid identical to micro-step
14b (0.1..100, 300 points) — no ID1-tuned parameter change.

Output:
  - exploration/figures/micro_step_15_cross_seizure.png (3 panels)
  - exploration/figures/micro_step_15_pr_matrix.npz (raw PR matrix
    plus metadata, so downstream plots don't have to recompute)
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


PATIENT_ID = "ID1"
PATIENT_FOLDER = Path(r"C:\Users\User\Desktop\SWEC") / PATIENT_ID
t_grid = np.linspace(0.1, 100.0, 300)


# Sort seizure files by number (not lexicographically, so Sz2 comes
# before Sz10 — otherwise reporting order gets confusing).
def sz_num(p: Path) -> int:
    return int(p.stem[2:])   # "Sz10" -> 10

sz_files = sorted(PATIENT_FOLDER.glob("Sz*.mat"), key=sz_num)
if not sz_files:
    raise SystemExit(f"No Sz*.mat files found in {PATIENT_FOLDER}")

print(f"Found {len(sz_files)} seizure files in {PATIENT_FOLDER.name}")


# Load once to fix the electrode count (should be 47 for ID1 per ADR 0006).
first = load_seizure(str(sz_files[0]))
n_electrodes = first.eeg.shape[0]
n_seizures = len(sz_files)

print(f"Patient {PATIENT_ID}: {n_seizures} seizures, {n_electrodes} electrodes.\n")


# Results matrix: rows = seizures, cols = electrodes.
pr_matrix = np.empty((n_seizures, n_electrodes))
sz_labels = []
sz_durations = []

t_start = time.time()
for si, sz_path in enumerate(sz_files):
    sz = load_seizure(str(sz_path))
    sz_labels.append(sz.seizure_id)
    sz_durations.append(sz.ictal_duration_s)

    ictal = sz.segment("ictal")
    fcm = compute_fcm(ictal)
    A = fcm_to_graph(fcm, layer="positive")

    for s in range(n_electrodes):
        avg_q = time_averaged_distribution(A, s, t_grid, quantum=True)
        pr_matrix[si, s] = participation_ratio(avg_q)

    elapsed = time.time() - t_start
    eta = elapsed / (si + 1) * (n_seizures - si - 1)
    print(f"  [{si+1}/{n_seizures}] {sz.seizure_id:<8} "
          f"ictal {sz.ictal_duration_s:5.0f}s  "
          f"quantum PR range [{pr_matrix[si].min():5.2f}, {pr_matrix[si].max():5.2f}]  "
          f"| elapsed {elapsed/60:.1f}m, ETA {eta/60:.1f}m")


# --- Rank consistency analysis ------------------------------------------ #

# Per-seizure ranks: 1 = most localized (lowest PR), n_electrodes = most spread.
ranks_matrix = np.empty_like(pr_matrix)
for si in range(n_seizures):
    ranks_matrix[si] = pr_matrix[si].argsort().argsort() + 1

mean_rank = ranks_matrix.mean(axis=0)
std_rank = ranks_matrix.std(axis=0)

# Order electrodes by mean rank across seizures.
order = mean_rank.argsort()

print("\n" + "=" * 70)
print(f"5 most consistently LOCALIZED electrodes across {n_seizures} seizures:")
print("=" * 70)
for i, idx in enumerate(order[:5], 1):
    r = ranks_matrix[:, idx].astype(int)
    print(f"  #{i}: electrode {idx:>3}  "
          f"mean rank = {mean_rank[idx]:5.1f} ± {std_rank[idx]:4.1f}  "
          f"(range across seizures: {r.min()}–{r.max()})")

print(f"\n5 most consistently SPREAD electrodes:")
for i, idx in enumerate(order[-5:][::-1], 1):
    r = ranks_matrix[:, idx].astype(int)
    print(f"  #{i}: electrode {idx:>3}  "
          f"mean rank = {mean_rank[idx]:5.1f} ± {std_rank[idx]:4.1f}  "
          f"(range across seizures: {r.min()}–{r.max()})")


# Pairwise Spearman rank correlations between seizures.
pairs = []
for i in range(n_seizures):
    for j in range(i + 1, n_seizures):
        rho, _ = spearmanr(pr_matrix[i], pr_matrix[j])
        pairs.append(rho)
pairs = np.array(pairs)

print("\n" + "=" * 70)
print(f"Pairwise Spearman rank correlation ({len(pairs)} seizure pairs):")
print("=" * 70)
print(f"  mean:   {pairs.mean():+.3f}")
print(f"  median: {np.median(pairs):+.3f}")
print(f"  min:    {pairs.min():+.3f}")
print(f"  max:    {pairs.max():+.3f}")
print(f"  fraction with rho > 0.5: {(pairs > 0.5).mean():.1%}")
print(f"  fraction with rho < 0.0: {(pairs < 0).mean():.1%}")


# Follow-up on electrode 3 (the star of Sz1 per micro-step 14b).
print("\n" + "=" * 70)
print("Follow-up on electrode 3 (was #1 most localized in Sz1):")
print("=" * 70)
for si in range(n_seizures):
    r = int(ranks_matrix[si, 3])
    marker = "  <-- localized" if r <= 5 else ("  <-- spread" if r >= n_electrodes - 4 else "")
    print(f"  {sz_labels[si]:<8}  PR = {pr_matrix[si, 3]:6.2f}  "
          f"rank = {r:>2}/{n_electrodes}{marker}")


# Save raw data for future plots without re-computing.
npz_path = FIG_DIR / "micro_step_15_pr_matrix.npz"
np.savez(
    npz_path,
    pr_matrix=pr_matrix,
    ranks_matrix=ranks_matrix,
    sz_labels=np.array(sz_labels),
    sz_durations=np.array(sz_durations),
    t_grid=t_grid,
    patient_id=PATIENT_ID,
    n_electrodes=n_electrodes,
)
print(f"\nRaw PR matrix saved to {npz_path}")


# --- Figure: 3 panels --------------------------------------------------- #
fig, axes = plt.subplots(1, 3, figsize=(20, 6),
                          gridspec_kw={"width_ratios": [1.4, 1.1, 0.8]})

# Panel 1: heatmap (seizures × electrodes, columns ordered by mean rank).
ax = axes[0]
im = ax.imshow(pr_matrix[:, order], aspect="auto", cmap="viridis")
ax.set_yticks(range(n_seizures))
ax.set_yticklabels(sz_labels, fontsize=8)
ax.set_xticks(range(n_electrodes))
ax.set_xticklabels([str(e) for e in order], fontsize=6, rotation=90)
ax.set_xlabel("electrode (sorted by mean rank across seizures)")
ax.set_ylabel("seizure")
ax.set_title(f"Quantum PR across {n_seizures} seizures × {n_electrodes} electrodes\n"
             "consistent COLUMN colour → patient-level structure")
plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="quantum PR")

# Panel 2: mean rank ± std per electrode (ordered).
ax = axes[1]
x = np.arange(n_electrodes)
colors = []
for idx in order:
    if idx == 3:
        colors.append("#DD8452")   # highlight electrode 3
    else:
        colors.append("#4C72B0")
ax.errorbar(x, mean_rank[order], yerr=std_rank[order],
            fmt="o", markersize=4, capsize=2,
            ecolor="#4C72B0", elinewidth=0.8, alpha=0.85,
            markerfacecolor="none")
# Overlay coloured markers on top
for xi, idx in enumerate(order):
    ax.plot(xi, mean_rank[idx], "o", markersize=5,
            color=colors[xi], zorder=3)
ax.set_xticks(x)
ax.set_xticklabels([str(e) for e in order], fontsize=6, rotation=90)
ax.set_xlabel("electrode (sorted)")
ax.set_ylabel("mean rank across seizures  (± std)")
ax.set_title("Per-electrode rank consistency\n"
             "narrow error bars → stable ranking (orange = electrode 3)")
ax.grid(alpha=0.3)

# Panel 3: pairwise Spearman histogram.
ax = axes[2]
ax.hist(pairs, bins=15, color="#4C72B0", edgecolor="white", alpha=0.85)
ax.axvline(0, color="black", lw=1, ls="--")
ax.axvline(pairs.mean(), color="crimson", lw=1.5,
           label=f"mean = {pairs.mean():+.3f}")
ax.set_xlabel("Spearman rho between seizure pairs")
ax.set_ylabel("count")
ax.set_title(f"Ranking agreement\n({len(pairs)} seizure pairs)")
ax.legend(fontsize=8)
ax.grid(axis="y", alpha=0.3)

fig.suptitle(f"Cross-seizure consistency of CTQW localization "
             f"({PATIENT_ID}, {n_seizures} seizures, Layer P) — ADR 0006/0009",
             fontsize=12)
plt.tight_layout()

out = FIG_DIR / "micro_step_15_cross_seizure.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out}")