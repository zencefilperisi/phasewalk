"""
Micro-step 20: Is CTQW reproducibility better than / different from
classical centrality reproducibility on the same FCMs?

Micro-step 19 showed CTQW rankings are reproducible beyond chance. That
does not show CTQW beats a classical centrality — if eigenvector
centrality were equally reproducible, the reproducibility would just
reflect stable FCM structure, not anything quantum-specific. This step
runs the comparison the project's core thesis rests on.

For each patient, on the SAME seizures and SAME Layer P graphs used for
the CTQW analysis, we compute the cross-seizure consistency (mean
pairwise Spearman) of three per-electrode measures:

  - CTQW time-averaged PR      (loaded from the saved npz; the slow part
                                is not recomputed)
  - eigenvector centrality     (principal eigenvector of A — the natural
                                classical spectral analogue)
  - weighted degree / strength (A.sum(axis=1) — simplest classical)

Two questions:
  (1) Reproducibility: is CTQW consistency higher than, similar to, or
      lower than the classical measures' consistency?
  (2) Distinctness: how different is the CTQW ranking from each classical
      ranking (mean within-seizure Spearman between them)? CTQW being
      reproducible AND distinct from classical is the interesting case.

Interpretation set BEFORE computing:
  - CTQW consistency >> classical  => quantum captures a more
    reproducible structure (strong).
  - CTQW ~= classical, but CTQW distinct from classical => both stable,
    but CTQW tracks a different, reproducible feature (Phase 3 decides
    if it is clinically better).
  - CTQW << classical => classical more stable; weak for the thesis.

Seizure order is taken from each npz's sz_labels so the recomputed
classical measures line up with the loaded CTQW PR.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
FIG_DIR = Path(__file__).resolve().parent / "figures"

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from data.loader import load_seizure
from data.connectivity import compute_fcm, fcm_to_graph


SWEC = Path(r"C:\Users\User\Desktop\SWEC")
LAYER2_WINDOW_S = 30.0   # not used here (full ictal), kept for reference

patients = [
    ("ID1", "micro_step_16_L1_vs_L2.npz"),
    ("ID2", "micro_step_17_ID2.npz"),
    ("ID3", "micro_step_18_ID3.npz"),
]


def mean_pairwise_spearman(value_matrix):
    """value_matrix: (n_seizures, n_electrodes). Mean pairwise Spearman."""
    n = value_matrix.shape[0]
    rhos = []
    for i in range(n):
        for j in range(i + 1, n):
            rhos.append(spearmanr(value_matrix[i], value_matrix[j])[0])
    return np.mean(rhos), np.array(rhos)


def classical_measures(A):
    """Return dict of per-electrode classical centralities on graph A."""
    strength = A.sum(axis=1)                       # weighted degree
    evals, evecs = np.linalg.eigh(A)               # A is symmetric
    ev_cent = np.abs(evecs[:, -1])                 # principal eigenvector
    return {"eigenvector": ev_cent, "weighted_degree": strength}


results = {}
for pid, npz_name in patients:
    d = np.load(FIG_DIR / npz_name, allow_pickle=True)
    pr_ctqw = d["pr_L1"]                      # (n_seizures, n_electrodes)
    sz_labels = [str(x) for x in d["sz_labels"]]
    n_seizures, n_electrodes = pr_ctqw.shape

    # Recompute classical measures on the same seizures / graphs.
    ev_mat = np.empty((n_seizures, n_electrodes))
    deg_mat = np.empty((n_seizures, n_electrodes))
    for si, label in enumerate(sz_labels):
        sz = load_seizure(str(SWEC / pid / f"{label}.mat"))
        ictal = sz.segment("ictal")
        fcm = compute_fcm(ictal)
        A = fcm_to_graph(fcm, layer="positive")
        cm = classical_measures(A)
        ev_mat[si] = cm["eigenvector"]
        deg_mat[si] = cm["weighted_degree"]

    # (1) Reproducibility of each measure
    c_ctqw, _ = mean_pairwise_spearman(pr_ctqw)
    c_ev, _ = mean_pairwise_spearman(ev_mat)
    c_deg, _ = mean_pairwise_spearman(deg_mat)

    # (2) Distinctness: within-seizure Spearman between CTQW and classical
    dist_ev = np.mean([spearmanr(pr_ctqw[s], ev_mat[s])[0] for s in range(n_seizures)])
    dist_deg = np.mean([spearmanr(pr_ctqw[s], deg_mat[s])[0] for s in range(n_seizures)])

    results[pid] = dict(
        n_seizures=n_seizures, n_electrodes=n_electrodes,
        c_ctqw=c_ctqw, c_ev=c_ev, c_deg=c_deg,
        dist_ev=dist_ev, dist_deg=dist_deg,
    )

    print("=" * 70)
    print(f"{pid}: {n_seizures} seizures, {n_electrodes} electrodes")
    print("=" * 70)
    print("  (1) Cross-seizure consistency (mean pairwise Spearman):")
    print(f"        CTQW PR            : {c_ctqw:+.3f}")
    print(f"        eigenvector cent.  : {c_ev:+.3f}")
    print(f"        weighted degree    : {c_deg:+.3f}")
    print("  (2) Distinctness — CTQW vs classical (within-seizure Spearman):")
    print(f"        CTQW vs eigenvector: {dist_ev:+.3f}")
    print(f"        CTQW vs degree     : {dist_deg:+.3f}")
    print()


# --- Verdicts ------------------------------------------------------- #
print("=" * 70)
print("READING (per patient):")
print("=" * 70)
for pid, _ in patients:
    r = results[pid]
    # reproducibility verdict
    best = max([("CTQW", r["c_ctqw"]), ("eigenvector", r["c_ev"]),
                ("degree", r["c_deg"])], key=lambda t: t[1])
    repro = (f"CTQW {'>' if r['c_ctqw'] > max(r['c_ev'], r['c_deg']) else '<='} "
             f"classical (most consistent: {best[0]})")
    # distinctness verdict
    distinct = ("distinct from eigenvector"
                if abs(r["dist_ev"]) < 0.5 else "similar to eigenvector")
    print(f"  {pid}: {repro};  CTQW ranking is {distinct} "
          f"(rho {r['dist_ev']:+.2f})")


# --- Figure --------------------------------------------------------- #
fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

# Panel 1: reproducibility bars, grouped by patient
ax = axes[0]
pids = [p for p, _ in patients]
x = np.arange(len(pids))
w = 0.26
ax.bar(x - w, [results[p]["c_ctqw"] for p in pids], w,
       label="CTQW PR", color="#4C72B0")
ax.bar(x, [results[p]["c_ev"] for p in pids], w,
       label="eigenvector centrality", color="#DD8452")
ax.bar(x + w, [results[p]["c_deg"] for p in pids], w,
       label="weighted degree", color="#55A467")
ax.set_xticks(x)
ax.set_xticklabels(pids)
ax.set_ylabel("cross-seizure consistency (mean pairwise Spearman)")
ax.set_title("(1) Reproducibility: CTQW vs classical measures")
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)

# Panel 2: distinctness — CTQW vs classical rankings
ax = axes[1]
ax.bar(x - w/2, [results[p]["dist_ev"] for p in pids], w,
       label="CTQW vs eigenvector", color="#DD8452")
ax.bar(x + w/2, [results[p]["dist_deg"] for p in pids], w,
       label="CTQW vs degree", color="#55A467")
ax.axhline(0, color="black", lw=1)
ax.axhline(0.5, color="grey", lw=0.8, ls=":", alpha=0.7)
ax.axhline(-0.5, color="grey", lw=0.8, ls=":", alpha=0.7)
ax.set_xticks(x)
ax.set_xticklabels(pids)
ax.set_ylabel("within-seizure Spearman (CTQW vs classical)")
ax.set_title("(2) Distinctness: is CTQW ranking different from classical?\n"
             "near 0 = distinct, near ±1 = redundant")
ax.set_ylim(-1.05, 1.05)
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)

fig.suptitle("CTQW vs classical centrality — reproducibility and distinctness "
             "(Layer 1, same seizures as micro-steps 16-18)", fontsize=12)
plt.tight_layout()
out = FIG_DIR / "micro_step_20_ctqw_vs_classical.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nDone. Figure saved as {out}")