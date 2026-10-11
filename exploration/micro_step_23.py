"""
Micro-step 23: Characterizing the CTQW ranking — which interpretable
graph property (if any) explains it on real ictal networks?

Micro-step 20 showed the CTQW per-electrode ranking is reproducible and
distinct from eigenvector centrality and degree. This step asks what the
ranking actually tracks, using interpretable, independently-meaningful
node properties on the real Layer P graphs:

    weighted degree (strength), eigenvector centrality,
    betweenness, weighted clustering coefficient.

The naive characterization — "spectral participation predicts CTQW" —
is partly tautological (Ziegler 2024: the time-averaged return
probability is definitionally tied to spectral participation), so it is
computed only as a REFERENCE CEILING, clearly flagged, not as a finding.

For each measure we compute the Spearman correlation with CTQW PR per
seizure and average across each patient's seizures.

Expected honest outcome (set before running): the interpretable
centralities will NOT strongly predict CTQW PR (consistent with the
synthetic result, micro-step 7), while spectral participation will —
i.e. CTQW captures a whole-spectrum property not reducible to a single
named centrality. If some interpretable measure DOES predict it well,
that would be a (better) interpretable characterization.

No walks here — CTQW PR loaded from npz; only graph measures recomputed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
FIG_DIR = Path(__file__).resolve().parent / "figures"

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from data.loader import load_seizure
from data.connectivity import compute_fcm, fcm_to_graph


SWEC = Path(r"C:\Users\User\Desktop\SWEC")

patients = [
    ("ID1", "micro_step_16_L1_vs_L2.npz"),
    ("ID2", "micro_step_17_ID2.npz"),
    ("ID3", "micro_step_18_ID3.npz"),
]

INTERPRETABLE = ["weighted_degree", "eigenvector", "betweenness", "clustering"]
REFERENCE = "spectral_participation"   # tautological ceiling


def node_measures(A):
    """Compute interpretable centralities + spectral participation on A."""
    n = A.shape[0]
    out = {}

    # weighted degree (node strength)
    out["weighted_degree"] = A.sum(axis=1)

    # eigenvector centrality = principal eigenvector (A symmetric, >=0)
    evals, evecs = np.linalg.eigh(A)
    out["eigenvector"] = np.abs(evecs[:, -1])

    # weighted clustering coefficient
    G = nx.from_numpy_array(A)
    cl = nx.clustering(G, weight="weight")
    out["clustering"] = np.array([cl[i] for i in range(n)])

    # betweenness with distance = 1 / weight (stronger edge = shorter path)
    Gd = nx.Graph()
    Gd.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            w = A[i, j]
            if w > 0:
                Gd.add_edge(i, j, distance=1.0 / w)
    bc = nx.betweenness_centrality(Gd, weight="distance")
    out["betweenness"] = np.array([bc.get(i, 0.0) for i in range(n)])

    # spectral participation (reference ceiling, partly tautological):
    # for node s, how spread its initial state is over the eigenbasis.
    out["spectral_participation"] = np.array([
        1.0 / np.sum((evecs[s, :] ** 2) ** 2) for s in range(n)
    ])

    return out


results = {}
for pid, npz_name in patients:
    d = np.load(FIG_DIR / npz_name, allow_pickle=True)
    pr_ctqw = d["pr_L1"]
    sz_labels = [str(x) for x in d["sz_labels"]]
    n_seizures, n_electrodes = pr_ctqw.shape

    measures = INTERPRETABLE + [REFERENCE]
    # per-seizure Spearman(measure, CTQW PR)
    rho_per_sz = {m: [] for m in measures}

    for si, label in enumerate(sz_labels):
        sz = load_seizure(str(SWEC / pid / f"{label}.mat"))
        ictal = sz.segment("ictal")
        A = fcm_to_graph(compute_fcm(ictal), layer="positive")
        nm = node_measures(A)
        for m in measures:
            rho = spearmanr(nm[m], pr_ctqw[si])[0]
            rho_per_sz[m].append(rho)

    mean_rho = {m: float(np.mean(rho_per_sz[m])) for m in measures}
    results[pid] = dict(n_seizures=n_seizures, mean_rho=mean_rho,
                        rho_per_sz=rho_per_sz)

    print("=" * 70)
    print(f"{pid}: mean Spearman(measure, CTQW PR) across {n_seizures} seizures")
    print("=" * 70)
    for m in INTERPRETABLE:
        print(f"  {m:<22}: {mean_rho[m]:+.3f}")
    print(f"  {'-'*40}")
    print(f"  {REFERENCE:<22}: {mean_rho[REFERENCE]:+.3f}  (tautological ceiling)")
    best_interp = max(INTERPRETABLE, key=lambda m: abs(mean_rho[m]))
    print(f"  best interpretable: {best_interp} ({mean_rho[best_interp]:+.3f})")
    print()


# --- Figure --------------------------------------------------------- #
fig, ax = plt.subplots(figsize=(13, 6))
pids = [p for p, _ in patients]
all_measures = INTERPRETABLE + [REFERENCE]
x = np.arange(len(all_measures))
w = 0.25
colors = {"ID1": "#4C72B0", "ID2": "#DD8452", "ID3": "#55A467"}
for k, pid in enumerate(pids):
    vals = [abs(results[pid]["mean_rho"][m]) for m in all_measures]
    ax.bar(x + (k - 1) * w, vals, w, label=pid, color=colors[pid])
ax.axvline(len(INTERPRETABLE) - 0.5, color="black", lw=1, ls="--")
ax.text(len(INTERPRETABLE) - 0.5, ax.get_ylim()[1]*0.95, "  tautological ->",
        fontsize=8, va="top")
ax.set_xticks(x)
ax.set_xticklabels([m.replace("_", "\n") for m in all_measures], fontsize=9)
ax.set_ylabel("|mean Spearman| with CTQW PR")
ax.set_title("Which graph property explains the CTQW ranking on real ictal networks?\n"
             "interpretable centralities (left) vs spectral participation "
             "(right, tautological reference)")
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
out = FIG_DIR / "micro_step_23_characterization.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.show()
print(f"Done. Figure saved as {out}")


# --- Verdict -------------------------------------------------------- #
print("\n" + "=" * 70)
print("READING:")
print("=" * 70)
for pid, _ in patients:
    mr = results[pid]["mean_rho"]
    best_interp = max(INTERPRETABLE, key=lambda m: abs(mr[m]))
    interp_strength = abs(mr[best_interp])
    taut = abs(mr[REFERENCE])
    if interp_strength < 0.4:
        line = (f"no interpretable centrality explains CTQW well "
                f"(best {best_interp} |rho|={interp_strength:.2f}); "
                f"spectral participation {taut:.2f} -> CTQW is spectrally "
                f"distinct, not reducible to a named centrality")
    else:
        line = (f"{best_interp} explains CTQW moderately "
                f"(|rho|={interp_strength:.2f}) -> partial interpretable "
                f"characterization")
    print(f"  {pid}: {line}")