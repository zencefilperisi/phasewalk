"""
Micro-step 22: PLV robustness check — do the CTQW findings hold under a
phase-based connectivity measure, or are they specific to Pearson?

ADR 0007 names PLV as the robustness measure. ADR 0007 note (2026-10-07)
fixed it to broadband PLV. This step re-runs the CTQW pipeline with PLV
FCMs (Layer P) and compares against the Pearson-based CTQW rankings.

Two measurements per patient:
  (1) Per-seizure agreement: Spearman between the Pearson-based CTQW
      ranking and the PLV-based CTQW ranking for each seizure — does
      switching the connectivity measure change the ranking?
  (2) Within-cohort consistency of PLV-based CTQW (mean pairwise
      Spearman) vs Pearson's (+0.410 / +0.630 / +0.535) — is the CTQW
      ranking reproducible under PLV too?

Interpretation set BEFORE computing:
  - high per-seizure agreement + similar within-cohort consistency:
    findings are robust to the connectivity measure; CTQW
    reproducibility is not a Pearson artefact.
  - low agreement or very different consistency: findings are
    measure-dependent.

Pearson CTQW PR loaded from npz (not recomputed). PLV CTQW computed
fresh (walks re-run; slow). Seizure order from each npz's sz_labels.
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
        return np.nan
    rhos = []
    for i in range(n):
        for j in range(i + 1, n):
            rhos.append(spearmanr(value_matrix[i], value_matrix[j])[0])
    return np.mean(rhos)


results = {}
t_start = time.time()
for pid, npz_name, pearson_consistency in patients:
    d = np.load(FIG_DIR / npz_name, allow_pickle=True)
    pr_pearson = d["pr_L1"]
    sz_labels = [str(x) for x in d["sz_labels"]]
    n_seizures, n_electrodes = pr_pearson.shape

    pr_plv = np.empty((n_seizures, n_electrodes))
    for si, label in enumerate(sz_labels):
        sz = load_seizure(str(SWEC / pid / f"{label}.mat"))
        ictal = sz.segment("ictal")
        plv_fcm = compute_fcm(ictal, measure="plv")
        A = fcm_to_graph(plv_fcm, layer="positive")   # PLV has no negatives
        for s in range(n_electrodes):
            avg_q = time_averaged_distribution(A, s, t_grid, quantum=True)
            pr_plv[si, s] = participation_ratio(avg_q)
        print(f"  {pid}/{label}: PLV CTQW done "
              f"(elapsed {(time.time()-t_start)/60:.1f}m)")

    # (1) Per-seizure Pearson vs PLV agreement
    per_sz = np.array([
        spearmanr(pr_pearson[si], pr_plv[si])[0] for si in range(n_seizures)
    ])
    # (2) Within-cohort consistency of PLV-based CTQW
    plv_consistency = mean_pairwise_spearman(pr_plv)

    results[pid] = dict(
        n_seizures=n_seizures, n_electrodes=n_electrodes, sz_labels=sz_labels,
        pr_pearson=pr_pearson, pr_plv=pr_plv, per_sz=per_sz,
        pearson_consistency=pearson_consistency, plv_consistency=plv_consistency,
    )

    print("=" * 70)
    print(f"{pid}: {n_seizures} seizures, {n_electrodes} electrodes")
    print("=" * 70)
    print("  (1) Per-seizure Pearson vs PLV CTQW ranking agreement:")
    for si, label in enumerate(sz_labels):
        print(f"        {label:<8}  Spearman = {per_sz[si]:+.3f}")
    print(f"        mean = {per_sz.mean():+.3f}")
    print("  (2) Within-cohort consistency:")
    print(f"        Pearson: {pearson_consistency:+.3f}")
    print(f"        PLV:     {plv_consistency:+.3f}")
    print()


# Save
npz_out = FIG_DIR / "micro_step_22_plv.npz"
np.savez(
    npz_out,
    **{f"{pid}_pr_plv": results[pid]["pr_plv"] for pid, _, _ in patients},
    **{f"{pid}_per_sz": results[pid]["per_sz"] for pid, _, _ in patients},
)
print(f"Data saved to {npz_out}")


# --- Figure --------------------------------------------------------- #
fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

# Panel 1: per-seizure Pearson-vs-PLV agreement
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
    xpos += 0.8
ax.axhline(0.7, color="grey", lw=0.8, ls=":", alpha=0.7, label="0.7 (high agreement)")
ax.axhline(0, color="black", lw=1)
ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels, fontsize=6, rotation=90)
ax.set_ylabel("Spearman (Pearson vs PLV ranking)")
ax.set_title("(1) Per-seizure: does switching to PLV change the CTQW ranking?")
ax.set_ylim(-0.3, 1.05)
ax.legend(fontsize=8)
ax.grid(axis="y", alpha=0.3)

# Panel 2: within-cohort consistency, Pearson vs PLV
ax = axes[1]
pids = [p for p, _, _ in patients]
x = np.arange(len(pids))
w = 0.35
ax.bar(x - w/2, [results[p]["pearson_consistency"] for p in pids], w,
       label="Pearson", color="#4C72B0")
ax.bar(x + w/2, [results[p]["plv_consistency"] for p in pids], w,
       label="PLV", color="#8172B3")
ax.set_xticks(x)
ax.set_xticklabels(pids)
ax.set_ylabel("within-cohort consistency (mean pairwise Spearman)")
ax.set_title("(2) Is CTQW reproducible under PLV as under Pearson?")
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)

fig.suptitle("PLV robustness check for CTQW — does the connectivity measure matter? "
             "(ADR 0007)", fontsize=12)
plt.tight_layout()
out = FIG_DIR / "micro_step_22_plv.png"
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
    agree = ("robust (ranking barely changes)" if m > 0.7
             else "moderate change" if m > 0.4
             else "measure-dependent (ranking changes a lot)")
    print(f"  {pid}: per-seizure Pearson-vs-PLV mean {m:+.3f} -> {agree}")
    print(f"        consistency: Pearson {r['pearson_consistency']:+.3f}, "
          f"PLV {r['plv_consistency']:+.3f}")