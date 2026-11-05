"""
ds003029 - confound check (exploratory, on the DISCOVERY data).

Phase 3 left an exploratory lead: in Engel-I patients, SOZ electrodes tend to
have LOWER weighted degree / eigenvector centrality (peripheral) than non-SOZ
electrodes (reverse-direction AUC ~0.76-0.78). Before investing in a held-out
dataset and a pre-registration, we check the most obvious confound: electrode
type (depth SEEG vs surface ECoG). Depth electrodes can be intrinsically less
connected; if SOZ sits on depths, "peripheral = SOZ" could be an electrode-
placement artefact, not a network property.

Key structural point: a subject whose electrodes are ALL one type has no
within-subject type confound (its AUC compares SOZ vs non-SOZ among same-type
electrodes). The confound only bites for mixed-montage subjects. So the decisive
question is: does the low-degree -> SOZ effect hold within PURE-ECoG subjects?

This recomputes per-electrode weighted degree + eigenvector centrality (Pearson
FCM -> Layer P), averaged across each subject's seizures over the frozen common
channels; NO CTQW, so it is fast. Reads exploration/ds003029/cohort_frozen.json.
Exploratory only (discovery data) -- NOT a confirmatory test.

Requires: mne, pandas, numpy, scipy, openpyxl
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
from scipy.stats import rankdata, wilcoxon

from data.loader_bids import discover_ictal_runs, load_ictal_window
from data.connectivity import compute_fcm, fcm_to_graph

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")
HERE = Path(__file__).resolve().parent
COHORT_JSON = HERE / "cohort_frozen.json"


def auc(score, mask):
    mask = np.asarray(mask, bool)
    npos, nneg = int(mask.sum()), int((~mask).sum())
    if npos == 0 or nneg == 0:
        return np.nan
    r = rankdata(score)
    return (r[mask].sum() - npos * (npos + 1) / 2) / (npos * nneg)


def name_to_type(channels_tsv):
    c = pd.read_csv(channels_tsv, sep="\t")
    return {str(r["name"]): str(r["type"]) for _, r in c.iterrows()}


cohort = json.loads(COHORT_JSON.read_text())
engel_i = {k: v for k, v in cohort.items() if v.get("engel") == 1}
print(f"Engel-I subjects: {len(engel_i)}\n")

rows = []           # per-subject summary
elec_records = []   # per-electrode (for pooled type x SOZ)
for label, info in engel_i.items():
    common = info["common_channels"]
    soz_set = {s.upper() for s in info["soz"]}
    runs = discover_ictal_runs(DATA_ROOT / f"sub-{label}")

    deg_runs, eig_runs = [], []
    typ_map = None
    for rp in runs:
        ir = load_ictal_window(rp)
        if ir.dropped_reason:
            continue
        A = fcm_to_graph(compute_fcm(ir.signal, measure="pearson"), layer="positive")
        deg = A.sum(axis=1)
        evals, evecs = np.linalg.eigh(A)
        eig = np.abs(evecs[:, -1])
        deg_runs.append(dict(zip(ir.good_names, deg)))
        eig_runs.append(dict(zip(ir.good_names, eig)))
        if typ_map is None:
            typ_map = name_to_type(rp.channels_tsv)
    if not deg_runs:
        continue

    deg = np.array([np.mean([d[n] for d in deg_runs if n in d]) for n in common])
    eig = np.array([np.mean([d[n] for d in eig_runs if n in d]) for n in common])
    typ = np.array([typ_map.get(n, "?") for n in common])
    soz = np.array([n.upper() in soz_set for n in common], dtype=bool)

    types_present = sorted(set(typ))
    montage = types_present[0] if len(types_present) == 1 else "MIXED"

    # AUC in the lead's direction: low degree / low eigenvector -> SOZ
    a_deg = auc(-deg, soz)
    a_eig = auc(-eig, soz)
    rows.append(dict(label=label, n=len(common), montage=montage,
                     n_soz=int(soz.sum()), n_soz_seeg=int(soz[typ == "SEEG"].sum()),
                     n_soz_ecog=int(soz[typ == "ECOG"].sum()),
                     auc_deg=a_deg, auc_eig=a_eig,
                     deg=deg, eig=eig, typ=typ, soz=soz))
    for n, t, s in zip(common, typ, soz):
        elec_records.append(dict(subject=label, type=t, soz=bool(s)))

df = pd.DataFrame([{k: r[k] for k in
                    ["label", "n", "montage", "n_soz", "n_soz_seeg", "n_soz_ecog",
                     "auc_deg", "auc_eig"]} for r in rows])

# ------------------------------------------------------------------ #
print(f"{'subject':<10}{'n':>4}{'montage':>8}{'nSOZ':>5}{'SEEG':>5}{'ECoG':>5}"
      f"{'AUC_deg':>9}{'AUC_eig':>9}")
print("-" * 60)
for _, r in df.iterrows():
    print(f"{r['label']:<10}{r['n']:>4}{r['montage']:>8}{r['n_soz']:>5}"
          f"{r['n_soz_seeg']:>5}{r['n_soz_ecog']:>5}{r['auc_deg']:>9.3f}{r['auc_eig']:>9.3f}")

print("\n" + "=" * 66)
print("CONFOUND CHECK")
print("=" * 66)

# 1. montage composition
comp = df["montage"].value_counts().to_dict()
print(f"montage composition (subjects): {comp}")

# 2. electrode-level type x SOZ (pooled)
e = pd.DataFrame(elec_records)
print("\npooled electrode type x SOZ:")
for t in sorted(e["type"].unique()):
    sub = e[e["type"] == t]
    print(f"  {t:<6}: {len(sub):>5} electrodes, SOZ rate {sub['soz'].mean():.3f}")

# 3. the decisive split: effect within pure-ECoG subjects (type-confound-free)
def summ(sub, name):
    a = sub["auc_deg"].dropna().values
    if len(a) < 3:
        print(f"  {name}: n={len(a)} too few"); return
    stat, p = wilcoxon(a - 0.5, alternative="greater")
    print(f"  {name:<28} n={len(a):>2}  median AUC_deg={np.median(a):.3f}  "
          f"Wilcoxon(>0.5) p={p:.4f}")

print("\neffect (low degree -> SOZ) by subject montage:")
summ(df, "ALL Engel-I")
summ(df[df["montage"] == "ECOG"], "pure-ECoG subjects")
summ(df[df["montage"] == "SEEG"], "pure-SEEG subjects")
summ(df[df["montage"] == "MIXED"], "mixed-montage subjects")

# 4. within-type AUC for mixed subjects (removes the type confound directly)
mixed = [r for r in rows if r["montage"] == "MIXED"]
if mixed:
    print("\nmixed-montage subjects, AUC within each electrode type:")
    for r in mixed:
        for t in ["ECOG", "SEEG"]:
            m = r["typ"] == t
            if m.sum() >= 3 and 0 < r["soz"][m].sum() < m.sum():
                print(f"  {r['label']:<9} {t:<5}: AUC_deg={auc(-r['deg'][m], r['soz'][m]):.3f}"
                      f"  ({int(r['soz'][m].sum())}/{int(m.sum())} SOZ)")
else:
    print("\n(no mixed-montage subjects -> the within-subject effect is already "
          "type-confound-free by construction)")

df.to_csv(HERE / "confound_type_results.csv", index=False)
print("\n" + "=" * 66)
print("READING: if the effect holds in pure-ECoG subjects (median AUC_deg well")
print("above 0.5, significant), it is NOT a depth-electrode artefact and the lead")
print("survives -> worth pre-registering + validating on held-out data.")
print(f"\nsaved -> {HERE / 'confound_type_results.csv'}")
print("Paste the full output back.")