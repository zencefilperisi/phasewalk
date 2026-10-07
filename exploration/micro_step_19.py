"""
Micro-step 19: Null-model test for the cross-seizure consistency and
signature-electrode phenomena.

The Phase 2C interim summary flagged (section 4.3, section 6) that the
"signature electrode" observation rests on low within-patient rank
variance but has NOT been tested against a null model. This step runs
that test, converting the observation into either a significant result
or a null one.

Two tests per patient, by permutation:

  Test 1 (overall consistency):
    Observed statistic = mean pairwise Spearman of CTQW PR rankings
    across the patient's seizures.
    Null = each seizure's ranking is an independent uniform random
    permutation of the electrodes. We draw N_PERM null cohorts and
    recompute the statistic.
    p = fraction of null >= observed.
    => Is the ranking reproducible beyond chance?

  Test 2 (signature electrodes):
    Observed statistic = the SMALLEST per-electrode rank std across
    electrodes (the single most stable electrode).
    Null = same permutation scheme; recompute the minimum std.
    p = fraction of null <= observed.
    => Is the most stable electrode more stable than full randomness
       would produce? (Note: this does not separate "this electrode is
       special" from "there is general consistency" — it confirms the
       extreme is real, given consistency.)

Data read from the saved .npz of micro-steps 16 (ID1, 9 long sz),
17 (ID2), 18 (ID3). No recomputation of walks.

For ID1 we use the 9-long-seizure subset (micro-step 16) to match the
+0.410 figure and avoid the short-seizure noise confound (micro-step 15b).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
FIG_DIR = Path(__file__).resolve().parent / "figures"

import numpy as np
import matplotlib.pyplot as plt


N_PERM = 20000
rng = np.random.default_rng(0)


def load_ranks(npz_name, key="ranks_L1"):
    """Return the ranks matrix (n_seizures x n_electrodes) from a saved npz."""
    d = np.load(FIG_DIR / npz_name, allow_pickle=True)
    return d[key]


def mean_pairwise_spearman_from_ranks(ranks):
    """Spearman == Pearson on ranks. Mean over upper-triangle pairs."""
    n_seizures = ranks.shape[0]
    if n_seizures < 2:
        raise ValueError("need >= 2 seizures")
    C = np.corrcoef(ranks)
    iu = np.triu_indices(n_seizures, k=1)
    return C[iu].mean()


def permutation_nulls(n_seizures, n_electrodes, n_perm, rng):
    """Draw null distributions of (mean pairwise Spearman, min rank std)."""
    null_rho = np.empty(n_perm)
    null_min_std = np.empty(n_perm)
    for k in range(n_perm):
        perms = np.empty((n_seizures, n_electrodes))
        for s in range(n_seizures):
            perms[s] = rng.permutation(n_electrodes) + 1
        C = np.corrcoef(perms)
        iu = np.triu_indices(n_seizures, k=1)
        null_rho[k] = C[iu].mean()
        null_min_std[k] = perms.std(axis=0).min()
    return null_rho, null_min_std


patients = [
    ("ID1", "micro_step_16_L1_vs_L2.npz"),
    ("ID2", "micro_step_17_ID2.npz"),
    ("ID3", "micro_step_18_ID3.npz"),
]

results = {}
for pid, npz_name in patients:
    ranks = load_ranks(npz_name)
    n_seizures, n_electrodes = ranks.shape

    obs_rho = mean_pairwise_spearman_from_ranks(ranks)
    obs_min_std = ranks.std(axis=0).min()

    null_rho, null_min_std = permutation_nulls(n_seizures, n_electrodes, N_PERM, rng)

    # Test 1: is observed consistency HIGHER than null? (one-sided, upper)
    p_rho = (np.sum(null_rho >= obs_rho) + 1) / (N_PERM + 1)
    z_rho = (obs_rho - null_rho.mean()) / null_rho.std()

    # Test 2: is observed min-std LOWER than null? (one-sided, lower)
    p_std = (np.sum(null_min_std <= obs_min_std) + 1) / (N_PERM + 1)
    z_std = (obs_min_std - null_min_std.mean()) / null_min_std.std()

    results[pid] = dict(
        n_seizures=n_seizures, n_electrodes=n_electrodes,
        obs_rho=obs_rho, null_rho=null_rho, p_rho=p_rho, z_rho=z_rho,
        obs_min_std=obs_min_std, null_min_std=null_min_std,
        p_std=p_std, z_std=z_std,
    )

    print("=" * 70)
    print(f"{pid}: {n_seizures} seizures, {n_electrodes} electrodes, "
          f"{N_PERM} permutations")
    print("=" * 70)
    print(f"  TEST 1 — cross-seizure consistency (mean pairwise Spearman)")
    print(f"    observed:          {obs_rho:+.3f}")
    print(f"    null mean ± std:   {null_rho.mean():+.3f} ± {null_rho.std():.3f}")
    print(f"    null 95th pct:     {np.percentile(null_rho, 95):+.3f}")
    print(f"    z-score:           {z_rho:+.1f}")
    print(f"    p-value:           {p_rho:.5f}  "
          f"{'*** significant' if p_rho < 0.05 else 'n.s.'}")
    print(f"  TEST 2 — most stable electrode (min per-electrode rank std)")
    print(f"    observed min std:  {obs_min_std:.3f}")
    print(f"    null mean ± std:   {null_min_std.mean():.3f} ± {null_min_std.std():.3f}")
    print(f"    null 5th pct:      {np.percentile(null_min_std, 5):.3f}")
    print(f"    z-score:           {z_std:+.1f}")
    print(f"    p-value:           {p_std:.5f}  "
          f"{'*** significant' if p_std < 0.05 else 'n.s.'}")
    print()


# --- Figure: 2 rows (Test 1, Test 2) x 3 cols (patients) ------------- #
fig, axes = plt.subplots(2, 3, figsize=(18, 9))

for col, (pid, _) in enumerate(patients):
    r = results[pid]

    # Row 0: Test 1 (consistency)
    ax = axes[0, col]
    ax.hist(r["null_rho"], bins=50, color="#8899AA", alpha=0.8, edgecolor="none")
    ax.axvline(r["obs_rho"], color="crimson", lw=2,
               label=f"observed {r['obs_rho']:+.3f}")
    ax.axvline(np.percentile(r["null_rho"], 95), color="black", lw=1, ls=":",
               label="null 95th pct")
    ax.set_title(f"{pid} — Test 1: consistency\n"
                 f"p = {r['p_rho']:.5f}, z = {r['z_rho']:+.1f}")
    ax.set_xlabel("mean pairwise Spearman")
    ax.set_ylabel("null count")
    ax.legend(fontsize=8)

    # Row 1: Test 2 (signature electrode)
    ax = axes[1, col]
    ax.hist(r["null_min_std"], bins=50, color="#8899AA", alpha=0.8, edgecolor="none")
    ax.axvline(r["obs_min_std"], color="crimson", lw=2,
               label=f"observed {r['obs_min_std']:.2f}")
    ax.axvline(np.percentile(r["null_min_std"], 5), color="black", lw=1, ls=":",
               label="null 5th pct")
    ax.set_title(f"{pid} — Test 2: most stable electrode\n"
                 f"p = {r['p_std']:.5f}, z = {r['z_std']:+.1f}")
    ax.set_xlabel("minimum per-electrode rank std")
    ax.set_ylabel("null count")
    ax.legend(fontsize=8)

fig.suptitle(f"Null-model tests: is CTQW ranking consistency beyond chance? "
             f"({N_PERM} permutations per patient)", fontsize=13)
plt.tight_layout()
out = FIG_DIR / "micro_step_19_null_model.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"Done. Figure saved as {out}")

# Summary line
print("\n" + "=" * 70)
print("SUMMARY (p < 0.05 = beyond chance):")
print("=" * 70)
for pid, _ in patients:
    r = results[pid]
    t1 = "YES" if r["p_rho"] < 0.05 else "no"
    t2 = "YES" if r["p_std"] < 0.05 else "no"
    print(f"  {pid}:  consistency beyond chance? {t1:<4}  "
          f"signature electrode beyond chance? {t2}")