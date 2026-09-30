"""
Micro-step 15b: Diagnostic — is the +0.31 -> +0.41 shift a genuine
duration confounder, or just a sampling artefact of dropping n from
13 to 9?

Answered by splitting the 78 seizure pairs from micro-step 15 into
two groups:

  Group A: at least one of the pair is a short seizure (< 30s)
           i.e. Sz4 (10s), Sz9 (11s), Sz10 (12s), Sz11 (20s)
  Group B: both seizures are long (>= 30s)

If duration is a real confounder:
  Group A mean rho  <  Group B mean rho
  (pairs involving short seizures are systematically less consistent)

If the +0.10 shift is just from dropping n:
  Group A mean rho  ~=  Group B mean rho
  (the individual pair rhos are similar; the mean only moved because
   we removed some of the pairs)

This reads the saved .npz from micro-step 15 — no recomputation.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
FIG_DIR = Path(__file__).resolve().parent / "figures"

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu


data = np.load(FIG_DIR / "micro_step_15_pr_matrix.npz", allow_pickle=True)
pr_matrix = data["pr_matrix"]
sz_labels = data["sz_labels"]
sz_durations = data["sz_durations"]
n_seizures = len(sz_labels)

DURATION_THRESHOLD = 30.0
short_mask = sz_durations < DURATION_THRESHOLD
short_idx = np.where(short_mask)[0]
long_idx = np.where(~short_mask)[0]

print(f"{'seizure':<8} {'ictal(s)':>10}  {'group':<7}")
print("-" * 30)
for i in range(n_seizures):
    g = "short" if short_mask[i] else "long"
    print(f"  {sz_labels[i]:<6} {sz_durations[i]:>10.0f}  {g}")

print(f"\n{len(short_idx)} short, {len(long_idx)} long")


# Split the 78 pairs by group
pairs_A = []   # at least one short
pairs_B = []   # both long
pair_info = []

for i in range(n_seizures):
    for j in range(i + 1, n_seizures):
        rho, _ = spearmanr(pr_matrix[i], pr_matrix[j])
        i_short = short_mask[i]
        j_short = short_mask[j]

        if i_short or j_short:
            pairs_A.append(rho)
            grp = "A (>=1 short)"
        else:
            pairs_B.append(rho)
            grp = "B (both long)"

        pair_info.append((sz_labels[i], sz_labels[j], rho, grp))

pairs_A = np.array(pairs_A)
pairs_B = np.array(pairs_B)

# Summaries
print("\n" + "=" * 70)
print(f"Group A ({len(pairs_A)} pairs, at least one seizure < 30s):")
print("=" * 70)
print(f"  mean rho:   {pairs_A.mean():+.3f}")
print(f"  median:     {np.median(pairs_A):+.3f}")
print(f"  min / max:  {pairs_A.min():+.3f} / {pairs_A.max():+.3f}")

print(f"\nGroup B ({len(pairs_B)} pairs, both seizures >= 30s):")
print("=" * 70)
print(f"  mean rho:   {pairs_B.mean():+.3f}")
print(f"  median:     {np.median(pairs_B):+.3f}")
print(f"  min / max:  {pairs_B.min():+.3f} / {pairs_B.max():+.3f}")

delta = pairs_B.mean() - pairs_A.mean()
print(f"\n  Group B mean - Group A mean = {delta:+.3f}")
print(f"  (positive = pairs involving short seizures are LESS consistent)")

# Mann-Whitney U test — non-parametric, doesn't assume normality
u_stat, p_value = mannwhitneyu(pairs_B, pairs_A, alternative="greater")
print(f"\nMann-Whitney U test (H1: Group B > Group A):")
print(f"  U = {u_stat:.1f}   p = {p_value:.4f}")
if p_value < 0.05:
    verdict = ("Group B pairs are significantly more consistent — "
               "duration IS a real confounder.")
elif p_value < 0.15:
    verdict = ("Weak evidence for a duration effect. "
               "Not conclusive at conventional thresholds.")
else:
    verdict = ("No significant difference between groups. "
               "The +0.10 shift is likely a sampling artefact, not a duration effect.")
print(f"\n  --> {verdict}")


# Verify: mean of Group B should approximate the +0.41 from micro-step 16
print(f"\nSanity check:")
print(f"  Group B (both-long) mean from THIS test:      {pairs_B.mean():+.3f}")
print(f"  Layer 1 subset mean from micro-step 16:       +0.410  (expected match)")
print(f"  All-13-seizures mean from micro-step 15:      +0.310")


# Figure: two overlapping histograms
fig, ax = plt.subplots(figsize=(10, 5.5))
bins = np.linspace(-0.15, 0.85, 22)
ax.hist(pairs_A, bins=bins, alpha=0.6, color="#DD8452",
        label=f"Group A ({len(pairs_A)} pairs, >=1 short seizure)\n"
              f"mean = {pairs_A.mean():+.3f}",
        edgecolor="white")
ax.hist(pairs_B, bins=bins, alpha=0.6, color="#4C72B0",
        label=f"Group B ({len(pairs_B)} pairs, both long)\n"
              f"mean = {pairs_B.mean():+.3f}",
        edgecolor="white")
ax.axvline(0, color="black", lw=1, ls="--")
ax.axvline(pairs_A.mean(), color="#DD8452", lw=1.5)
ax.axvline(pairs_B.mean(), color="#4C72B0", lw=1.5)
ax.set_xlabel("Spearman rho between seizure pairs")
ax.set_ylabel("count")
ax.set_title(f"Diagnostic: are pairs involving short seizures less consistent?\n"
             f"Mann-Whitney U p = {p_value:.4f}")
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
out = FIG_DIR / "micro_step_15b_diagnostic.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out}")