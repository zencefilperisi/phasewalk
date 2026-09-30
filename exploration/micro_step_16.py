"""
Micro-step 16: Layer 1 vs Layer 2 comparison on ID1.

Micro-step 15 gave a mean pairwise Spearman of +0.31 across 13 seizures
(Layer P, full ictal). Question: how much of that "low consistency" is
the duration confounder — 4 of the 13 seizures are under 30 s and may
give noisy FCMs from too little data?

Per ADR 0007's Comparison Protocol, we test this by restricting to
seizures with ictal duration >= 30 s (9 seizures for ID1) and computing
two versions of the analysis:

  Layer 1 subset:  full ictal, same 9 seizures     (what micro-step 15
                                                    did, but filtered)
  Layer 2:         first 30 s of ictal, same 9    (ADR 0007 sub-4)

Three comparisons:
  (i)   within-Layer 1 subset consistency  — pairwise Spearman across
                                              the 9 seizures' rankings
  (ii)  within-Layer 2 consistency          — same, on the first-30 s
                                              rankings
  (iii) Layer 1 vs Layer 2, per seizure    — how similar are the two
                                              rankings for the SAME
                                              seizure?

Interpretations, chosen BEFORE computing (per ADR 0006):
  - If (ii) >> (i): duration was the confounder, Layer 2 is cleaner.
  - If (ii) ~= (i): duration wasn't the issue; the +0.31 is genuine
    noise. The "modest consistency" story stands.
  - If (iii) is high (~0.7+): first-30 s captures the full-ictal
    ranking well; downstream analysis can use Layer 2 as a proxy.
  - If (iii) is low: first-30 s and full ictal describe different
    dynamics; both must be reported.

No thresholding, no free parameter tuning. 30-s window is fixed by
ADR 0007 sub-decision 4 (Bastos & Schoffelen 2016).
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
LAYER2_WINDOW_S = 30.0     # fixed by ADR 0007 sub-4


def sz_num(p: Path) -> int:
    return int(p.stem[2:])

all_sz_files = sorted(PATIENT_FOLDER.glob("Sz*.mat"), key=sz_num)
print(f"Found {len(all_sz_files)} seizure files.")


# --- First pass: filter to seizures with ictal >= 30 s ---------------- #
qualifying = []
for p in all_sz_files:
    sz = load_seizure(str(p))
    if sz.ictal_duration_s >= LAYER2_WINDOW_S:
        qualifying.append((p, sz))
    else:
        print(f"  Excluded {sz.seizure_id}: ictal {sz.ictal_duration_s:.0f}s < {LAYER2_WINDOW_S:.0f}s")

n_seizures = len(qualifying)
if n_seizures == 0:
    raise SystemExit("No seizures qualify for Layer 2.")

first_sz = qualifying[0][1]
n_electrodes = first_sz.eeg.shape[0]
window_samples = int(LAYER2_WINDOW_S * first_sz.fs)

print(f"\nQualifying seizures: {n_seizures}, electrodes: {n_electrodes}, "
      f"window: {window_samples} samples ({LAYER2_WINDOW_S:.0f}s at {first_sz.fs}Hz)\n")


# --- Compute Layer 1 subset and Layer 2 PR matrices ------------------- #
pr_L1 = np.empty((n_seizures, n_electrodes))
pr_L2 = np.empty((n_seizures, n_electrodes))
sz_labels = []
sz_durations = []

t_start = time.time()
for si, (path, sz) in enumerate(qualifying):
    sz_labels.append(sz.seizure_id)
    sz_durations.append(sz.ictal_duration_s)

    ictal_full = sz.segment("ictal")
    ictal_first30 = ictal_full[:, :window_samples]

    # Layer 1 subset: full ictal
    fcm_L1 = compute_fcm(ictal_full)
    A_L1 = fcm_to_graph(fcm_L1, layer="positive")
    for s in range(n_electrodes):
        avg_q = time_averaged_distribution(A_L1, s, t_grid, quantum=True)
        pr_L1[si, s] = participation_ratio(avg_q)

    # Layer 2: first 30 s
    fcm_L2 = compute_fcm(ictal_first30)
    A_L2 = fcm_to_graph(fcm_L2, layer="positive")
    for s in range(n_electrodes):
        avg_q = time_averaged_distribution(A_L2, s, t_grid, quantum=True)
        pr_L2[si, s] = participation_ratio(avg_q)

    elapsed = time.time() - t_start
    eta = elapsed / (si + 1) * (n_seizures - si - 1)
    print(f"  [{si+1}/{n_seizures}] {sz.seizure_id:<8} ictal {sz.ictal_duration_s:5.0f}s  "
          f"L1 PR [{pr_L1[si].min():5.2f}, {pr_L1[si].max():5.2f}]  "
          f"L2 PR [{pr_L2[si].min():5.2f}, {pr_L2[si].max():5.2f}]  "
          f"| elapsed {elapsed/60:.1f}m, ETA {eta/60:.1f}m")


# --- Analysis --------------------------------------------------------- #

# Comparison (i): within Layer 1 subset — pairwise Spearman across the 9
pairs_L1 = []
for i in range(n_seizures):
    for j in range(i + 1, n_seizures):
        rho, _ = spearmanr(pr_L1[i], pr_L1[j])
        pairs_L1.append(rho)
pairs_L1 = np.array(pairs_L1)

# Comparison (ii): within Layer 2 — pairwise Spearman across the 9
pairs_L2 = []
for i in range(n_seizures):
    for j in range(i + 1, n_seizures):
        rho, _ = spearmanr(pr_L2[i], pr_L2[j])
        pairs_L2.append(rho)
pairs_L2 = np.array(pairs_L2)

# Comparison (iii): Layer 1 vs Layer 2 for the SAME seizure
per_sz_L1_vs_L2 = []
for si in range(n_seizures):
    rho, _ = spearmanr(pr_L1[si], pr_L2[si])
    per_sz_L1_vs_L2.append(rho)
per_sz_L1_vs_L2 = np.array(per_sz_L1_vs_L2)


print("\n" + "=" * 70)
print("Within-Layer 1 subset  (9 seizures, full ictal, 36 pairs):")
print("=" * 70)
print(f"  mean Spearman:   {pairs_L1.mean():+.3f}")
print(f"  median:          {np.median(pairs_L1):+.3f}")
print(f"  fraction > 0.5:  {(pairs_L1 > 0.5).mean():.1%}")

print("\n" + "=" * 70)
print("Within-Layer 2         (9 seizures, first 30 s, 36 pairs):")
print("=" * 70)
print(f"  mean Spearman:   {pairs_L2.mean():+.3f}")
print(f"  median:          {np.median(pairs_L2):+.3f}")
print(f"  fraction > 0.5:  {(pairs_L2 > 0.5).mean():.1%}")

delta = pairs_L2.mean() - pairs_L1.mean()
print(f"\n  --> Layer 2 mean vs Layer 1 subset mean: {delta:+.3f}")
if delta > 0.10:
    print("      Layer 2 substantially more consistent — duration was a confounder.")
elif delta > 0.03:
    print("      Layer 2 slightly more consistent — mild duration effect.")
elif abs(delta) <= 0.03:
    print("      No meaningful change — duration was NOT the main confounder.")
else:
    print("      Layer 2 less consistent — restricting to onset hurts, unexpected.")


print("\n" + "=" * 70)
print("Per-seizure  Layer 1 vs Layer 2  (same seizure, two windows):")
print("=" * 70)
for si in range(n_seizures):
    r = per_sz_L1_vs_L2[si]
    marker = " <-- diverges" if r < 0.3 else (" <-- close match" if r > 0.7 else "")
    print(f"  {sz_labels[si]:<8}  rho = {r:+.3f}  ({sz_durations[si]:5.0f}s ictal){marker}")
print(f"\n  mean: {per_sz_L1_vs_L2.mean():+.3f}")


# Electrode 28 follow-up (the winner from micro-step 15)
ranks_L1 = np.empty_like(pr_L1)
ranks_L2 = np.empty_like(pr_L2)
for si in range(n_seizures):
    ranks_L1[si] = pr_L1[si].argsort().argsort() + 1
    ranks_L2[si] = pr_L2[si].argsort().argsort() + 1

print("\n" + "=" * 70)
print("Electrode 28 (micro-step 15's most-consistently-localized) across L1 / L2:")
print("=" * 70)
for si in range(n_seizures):
    r1 = int(ranks_L1[si, 28])
    r2 = int(ranks_L2[si, 28])
    print(f"  {sz_labels[si]:<8}  L1 rank {r1:>2}/{n_electrodes}  |  L2 rank {r2:>2}/{n_electrodes}")


# Save data for future plots
npz_path = FIG_DIR / "micro_step_16_L1_vs_L2.npz"
np.savez(
    npz_path,
    pr_L1=pr_L1, pr_L2=pr_L2,
    ranks_L1=ranks_L1, ranks_L2=ranks_L2,
    sz_labels=np.array(sz_labels),
    sz_durations=np.array(sz_durations),
    pairs_L1=pairs_L1, pairs_L2=pairs_L2,
    per_sz_L1_vs_L2=per_sz_L1_vs_L2,
    t_grid=t_grid,
    layer2_window_s=LAYER2_WINDOW_S,
)
print(f"\nRaw data saved to {npz_path}")


# --- Figure: 3 panels ------------------------------------------------- #
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Panel 1: within-cohort Spearman histograms, L1 vs L2 overlaid
ax = axes[0]
bins = np.linspace(-0.2, 1.0, 25)
ax.hist(pairs_L1, bins=bins, alpha=0.55, color="#4C72B0",
        label=f"Layer 1 subset (mean {pairs_L1.mean():+.3f})",
        edgecolor="white")
ax.hist(pairs_L2, bins=bins, alpha=0.55, color="#DD8452",
        label=f"Layer 2 first 30s (mean {pairs_L2.mean():+.3f})",
        edgecolor="white")
ax.axvline(0, color="black", lw=1, ls="--")
ax.axvline(pairs_L1.mean(), color="#4C72B0", lw=1.5)
ax.axvline(pairs_L2.mean(), color="#DD8452", lw=1.5)
ax.set_xlabel("Spearman rho between seizure pairs")
ax.set_ylabel("count")
ax.set_title(f"Within-cohort consistency\n(same 9 seizures, 36 pairs each)")
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)

# Panel 2: Layer 1 vs Layer 2 per-seizure agreement
ax = axes[1]
xs = np.arange(n_seizures)
colors = ["#DD8452" if r < 0.3 else "#55A467" if r > 0.7 else "#4C72B0"
          for r in per_sz_L1_vs_L2]
ax.bar(xs, per_sz_L1_vs_L2, color=colors, edgecolor="white", linewidth=0.5)
ax.axhline(per_sz_L1_vs_L2.mean(), color="black", lw=1.2, ls="--",
           label=f"mean = {per_sz_L1_vs_L2.mean():+.3f}")
ax.axhline(0.7, color="grey", lw=0.8, ls=":", alpha=0.7)
ax.axhline(0.3, color="grey", lw=0.8, ls=":", alpha=0.7)
ax.set_xticks(xs)
ax.set_xticklabels(sz_labels, fontsize=8, rotation=45)
ax.set_ylabel("Spearman rho  (L1 vs L2 ranking, same seizure)")
ax.set_title("How well does first 30 s stand for the whole ictal?\n"
             "green > 0.7 close match, orange < 0.3 diverges")
ax.set_ylim(-0.1, 1.05)
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)

# Panel 3: mean rank comparison per electrode (Layer 1 vs Layer 2)
ax = axes[2]
mean_rank_L1 = ranks_L1.mean(axis=0)
mean_rank_L2 = ranks_L2.mean(axis=0)
electrodes = np.arange(n_electrodes)
ax.scatter(mean_rank_L1, mean_rank_L2, s=60, color="#4C72B0",
           edgecolor="white", linewidth=1.0, zorder=3, alpha=0.85)
# Highlight electrode 28
ax.scatter(mean_rank_L1[28], mean_rank_L2[28], s=140, facecolor="#DD8452",
           edgecolor="black", linewidth=1.2, zorder=4, label="electrode 28")
# Equal-line
lo, hi = 0, n_electrodes + 1
ax.plot([lo, hi], [lo, hi], color="black", lw=1, ls="--", label="equal")
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_xlabel("mean rank across seizures — Layer 1")
ax.set_ylabel("mean rank across seizures — Layer 2")
ax.set_title("Per-electrode: does its mean rank move\nwhen restricting to first 30 s?")
ax.legend(fontsize=9, loc="upper left")
ax.grid(alpha=0.3)
ax.set_aspect("equal")

fig.suptitle(f"Layer 1 vs Layer 2 comparison "
             f"({PATIENT_ID}, {n_seizures} seizures ≥ 30s) — ADR 0007",
             fontsize=12)
plt.tight_layout()

out = FIG_DIR / "micro_step_16_L1_vs_L2.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out}")