"""
ds003029 - Cohort analysis (executes the ADR 0013 pre-registered plan).

Reads cohort_frozen.json (produced by ds003029_cohort_build.py), computes, per
subject, the across-seizure mean CTQW PR and the classical centralities over the
frozen common channels, then the per-subject AUC of each (score oriented as
pre-specified: CTQW = -PR; classical = +centrality). Runs the primary test
(Engel-I, one-sided Wilcoxon vs 0.5), the label-permutation null, the classical
baseline comparison, and the short-seizure sensitivity. Saves per-subject AUCs
to cohort_results.csv.

This is run AFTER the cohort is frozen; nothing here can change the cohort or the
SOZ sets.

Requires: mne, pandas, numpy, scipy, networkx, openpyxl
"""
import sys
import json
from pathlib import Path


def _find_root():
    for p in list(Path(__file__).resolve().parents) + [Path.cwd().resolve()]:
        if (p / "data" / "loader_bids.py").exists():
            return p
    return None


_root = _find_root()
if _root is None:
    sys.exit("Could not find repo root (folder with data/loader_bids.py).")
sys.path.insert(0, str(_root))

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import rankdata, wilcoxon

from data.loader_bids import discover_ictal_runs, load_ictal_window
from data.connectivity import compute_fcm, fcm_to_graph
from quantum_walks.walks import time_averaged_distribution_all, participation_ratio

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")
HERE = Path(__file__).resolve().parent
COHORT_JSON = HERE / "cohort_frozen.json"
T_GRID = np.linspace(0.1, 100.0, 300)
N_PERM = 5000
RNG = np.random.default_rng(0)

# pre-specified score orientations (ADR 0013): higher score = more SOZ-predicted
METHODS = ["ctqw_pr", "weighted_degree", "eigenvector", "betweenness"]
SIGN = {"ctqw_pr": -1.0, "weighted_degree": +1.0,
        "eigenvector": +1.0, "betweenness": +1.0}


def node_measures(A):
    """Raw per-node measures on a Layer-P graph A (symmetric, >=0, zero diag)."""
    out = {}
    avg = time_averaged_distribution_all(A, T_GRID, quantum=True)
    out["ctqw_pr"] = np.array([participation_ratio(avg[:, s]) for s in range(A.shape[0])])
    out["weighted_degree"] = A.sum(axis=1)
    evals, evecs = np.linalg.eigh(A)
    out["eigenvector"] = np.abs(evecs[:, -1])
    Gd = nx.Graph()
    Gd.add_nodes_from(range(A.shape[0]))
    for i in range(A.shape[0]):
        for j in range(i + 1, A.shape[0]):
            if A[i, j] > 0:
                Gd.add_edge(i, j, distance=1.0 / A[i, j])
    bc = nx.betweenness_centrality(Gd, weight="distance")
    out["betweenness"] = np.array([bc.get(i, 0.0) for i in range(A.shape[0])])
    return out


def auc(score, soz_mask):
    soz_mask = np.asarray(soz_mask, dtype=bool)
    n_pos, n_neg = int(soz_mask.sum()), int((~soz_mask).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    r = rankdata(score)
    return (r[soz_mask].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


# --------------------------------------------------------------------------- #
# Per-subject computation
# --------------------------------------------------------------------------- #
cohort = json.loads(COHORT_JSON.read_text())
print(f"frozen cohort: {len(cohort)} subjects\n")

subj_rows = []     # dicts with label, engel, aucs, scores, soz_mask, short
for label, info in cohort.items():
    common = info["common_channels"]
    soz_set = {s.upper() for s in info["soz"]}
    soz_mask = np.array([n.upper() in soz_set for n in common], dtype=bool)

    runs = discover_ictal_runs(DATA_ROOT / f"sub-{label}")
    per_method_runs = {k: [] for k in METHODS}
    short_any = False
    n_used = 0
    for rp in runs:
        ir = load_ictal_window(rp)
        if ir.dropped_reason:
            continue
        short_any = short_any or ir.window.short
        A = fcm_to_graph(compute_fcm(ir.signal, measure="pearson"), layer="positive")
        nm = node_measures(A)
        for k in METHODS:
            per_method_runs[k].append(dict(zip(ir.good_names, nm[k])))
        n_used += 1

    # across-seizure mean over the frozen common channels, oriented score
    scores = {}
    aucs = {}
    for k in METHODS:
        vals = np.array([np.mean([d[n] for d in per_method_runs[k] if n in d])
                         for n in common])
        scores[k] = SIGN[k] * vals
        aucs[k] = auc(scores[k], soz_mask)
    subj_rows.append(dict(label=label, engel=info["engel"], n_used=n_used,
                          n_good=len(common), n_soz=int(soz_mask.sum()),
                          soz_frac=info["soz_fraction"], short=short_any,
                          soz_mask=soz_mask, scores=scores, aucs=aucs))
    print(f"  {label:<9} engel={str(info['engel']):>3} runs={n_used} "
          f"good={len(common)} soz={int(soz_mask.sum())} "
          f"AUC(ctqw)={aucs['ctqw_pr']:.3f}")

# save per-subject AUCs
df = pd.DataFrame([{**{"label": r["label"], "engel": r["engel"], "n_used": r["n_used"],
                       "n_good": r["n_good"], "n_soz": r["n_soz"],
                       "soz_frac": r["soz_frac"], "short": r["short"]},
                    **{f"auc_{k}": r["aucs"][k] for k in METHODS}}
                   for r in subj_rows])
df.to_csv(HERE / "cohort_results.csv", index=False)


def group_auc(rows, method="ctqw_pr"):
    return np.array([r["aucs"][method] for r in rows], dtype=float)


def report_primary(rows, title):
    a = group_auc(rows, "ctqw_pr")
    a = a[~np.isnan(a)]
    if len(a) < 3:
        print(f"\n{title}: n={len(a)} too few"); return
    stat, p = wilcoxon(a - 0.5, alternative="greater")
    print(f"\n{title}  (n={len(a)})")
    print(f"  CTQW AUC: median={np.median(a):.3f}  mean={np.mean(a):.3f}  "
          f"IQR=[{np.percentile(a,25):.3f},{np.percentile(a,75):.3f}]")
    print(f"  one-sided Wilcoxon (AUC>0.5): W={stat:.1f}  p={p:.4f}")


# --------------------------------------------------------------------------- #
# ADR 0013 analyses
# --------------------------------------------------------------------------- #
print("\n" + "=" * 70)
print("PRIMARY + SECONDARY (ADR 0013)")
print("=" * 70)

engel_i = [r for r in subj_rows if r["engel"] == 1]
engel_hi = [r for r in subj_rows if r["engel"] not in (1, None)]

report_primary(engel_i, "PRIMARY: Engel-I")
report_primary(subj_rows, "SECONDARY: all outcomes")
report_primary(engel_hi, "SECONDARY: Engel > I (expected weaker)")

# short-seizure sensitivity (primary population)
report_primary([r for r in engel_i if not r["short"]],
               "SENSITIVITY: Engel-I, excluding short-seizure subjects")

# classical baseline comparison (Engel-I, paired)
print("\n" + "-" * 70)
print("CLASSICAL BASELINE (Engel-I, paired Wilcoxon vs CTQW)")
ctqw = group_auc(engel_i, "ctqw_pr")
for k in ["weighted_degree", "eigenvector", "betweenness"]:
    other = group_auc(engel_i, k)
    ok = ~(np.isnan(ctqw) | np.isnan(other))
    med = np.median(other[ok])
    try:
        _, p = wilcoxon(ctqw[ok] - other[ok])
    except ValueError:
        p = np.nan
    print(f"  {k:<16} AUC median={med:.3f}   CTQW median={np.median(ctqw[ok]):.3f}"
          f"   paired p={p:.3f}")

# label-permutation null (Engel-I group mean AUC)
print("\n" + "-" * 70)
print("LABEL-PERMUTATION NULL (Engel-I, CTQW)")
obs = np.nanmean(group_auc(engel_i, "ctqw_pr"))
null = np.empty(N_PERM)
for b in range(N_PERM):
    vals = []
    for r in engel_i:
        mask = r["soz_mask"].copy()
        RNG.shuffle(mask)
        vals.append(auc(r["scores"]["ctqw_pr"], mask))
    null[b] = np.nanmean(vals)
p_perm = (np.sum(null >= obs) + 1) / (N_PERM + 1)
print(f"  observed group-mean AUC = {obs:.3f}")
print(f"  null mean={null.mean():.3f} sd={null.std():.3f}  p={p_perm:.4f}")

# descriptive: AUC vs SOZ fraction (no filtering)
print("\n" + "-" * 70)
print("DESCRIPTIVE: AUC vs SOZ fraction (Engel-I)")
for r in sorted(engel_i, key=lambda x: x["soz_frac"]):
    print(f"  {r['label']:<9} soz_frac={r['soz_frac']:.2f}  AUC={r['aucs']['ctqw_pr']:.3f}")

print("\n" + "=" * 70)
print(f"Per-subject AUCs saved -> {HERE / 'cohort_results.csv'}")
print("Analysis done. Paste the full output back.")
print("=" * 70)