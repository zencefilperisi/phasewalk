"""
ds003029 - mechanism check (exploratory, discovery data; ADR 0014 item iii).

The Phase-3 lead is that SOZ electrodes have LOW functional degree (peripheral).
Proposed mechanism: the seizure-generating core decouples from propagating
activity DURING the seizure, lowering its ictal degree. This script tests
whether that peripheral position is ictal-specific or already present
pre-ictally, by computing, per subject, the AUC(-degree -> SOZ) in a pre-ictal
window and in the ictal window, and comparing them.

  ictal AUC >> pre-ictal AUC  -> SOZ becomes peripheral during the seizure
                                 (supports the decoupling-core mechanism)
  ictal AUC ~= pre-ictal AUC  -> peripheral position is static (would also mean
                                 an interictal validation could work)

Windows (relative to electrographic onset): pre-ictal [onset-35s, onset-5s],
ictal [onset, onset+30s]. Pearson FCM -> Layer P -> weighted degree, no CTQW.
Engel-I subjects from cohort_frozen.json. Exploratory only.

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
import mne
from scipy.stats import rankdata, wilcoxon

from data.loader_bids import discover_ictal_runs, parse_seizure_window
from data.connectivity import compute_fcm, fcm_to_graph

mne.set_log_level("ERROR")

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")
HERE = Path(__file__).resolve().parent
COHORT_JSON = HERE / "cohort_frozen.json"
PRE = (-35.0, -5.0)      # pre-ictal window relative to onset
ICTAL_LEN = 30.0


def auc(score, mask):
    mask = np.asarray(mask, bool)
    npos, nneg = int(mask.sum()), int((~mask).sum())
    if npos == 0 or nneg == 0:
        return np.nan
    r = rankdata(score)
    return (r[mask].sum() - npos * (npos + 1) / 2) / (npos * nneg)


def good_names(channels_tsv):
    c = pd.read_csv(channels_tsv, sep="\t")
    g = c[(c["type"].isin(["ECOG", "SEEG"])) & (c["status"] == "good")]
    return g["name"].astype(str).tolist()


def two_window_degrees(rp):
    """Return (good_names, preictal_degree, ictal_degree) or None."""
    ev = pd.read_csv(rp.events_tsv, sep="\t")
    w = parse_seizure_window(ev, ICTAL_LEN)
    if w is None:
        return None
    onset = w.onset
    if onset + PRE[0] < 0:                      # not enough pre-ictal signal
        return None
    raw = mne.io.read_raw_brainvision(str(rp.vhdr), preload=True, verbose="ERROR")
    sf = float(raw.info["sfreq"])
    dur = raw.n_times / sf
    i0, i1 = onset, min(w.win_end, dur)
    p0, p1 = onset + PRE[0], onset + PRE[1]
    if (i1 - i0) < 1.0 or (p1 - p0) < 1.0:
        return None
    gn = [n for n in good_names(rp.channels_tsv) if n in raw.ch_names]
    idx = [raw.ch_names.index(n) for n in gn]
    data = raw.get_data()[idx]

    def deg(a, b):
        A = fcm_to_graph(compute_fcm(data[:, int(a * sf):int(b * sf)],
                                     measure="pearson"), layer="positive")
        return A.sum(axis=1)

    return gn, deg(p0, p1), deg(i0, i1)


cohort = json.loads(COHORT_JSON.read_text())
engel_i = {k: v for k, v in cohort.items() if v.get("engel") == 1}
print(f"Engel-I subjects: {len(engel_i)}\n")

rows = []
for label, info in engel_i.items():
    common = info["common_channels"]
    soz_set = {s.upper() for s in info["soz"]}
    pre_runs, ict_runs = [], []
    for rp in discover_ictal_runs(DATA_ROOT / f"sub-{label}"):
        res = two_window_degrees(rp)
        if res is None:
            continue
        gn, pre_d, ict_d = res
        pre_runs.append(dict(zip(gn, pre_d)))
        ict_runs.append(dict(zip(gn, ict_d)))
    if not pre_runs:
        print(f"  {label}: no usable run (pre-ictal too short?) -> skipped")
        continue
    pre = np.array([np.mean([d[n] for d in pre_runs if n in d]) for n in common])
    ict = np.array([np.mean([d[n] for d in ict_runs if n in d]) for n in common])
    soz = np.array([n.upper() in soz_set for n in common], dtype=bool)
    a_pre, a_ict = auc(-pre, soz), auc(-ict, soz)
    rows.append(dict(label=label, n_runs=len(pre_runs), n_soz=int(soz.sum()),
                     auc_pre=a_pre, auc_ict=a_ict))
    print(f"  {label:<9} runs={len(pre_runs)} soz={int(soz.sum()):>2}  "
          f"AUC pre-ictal={a_pre:.3f}  ictal={a_ict:.3f}  Δ={a_ict - a_pre:+.3f}")

df = pd.DataFrame(rows)
df.to_csv(HERE / "mechanism_ictal_vs_preictal.csv", index=False)

print("\n" + "=" * 60)
print("MECHANISM: is the peripheral-SOZ position ictal-specific?")
print("=" * 60)
pre = df["auc_pre"].dropna().values
ict = df["auc_ict"].dropna().values
both = df.dropna(subset=["auc_pre", "auc_ict"])
print(f"subjects: {len(df)}  (with both windows: {len(both)})")
print(f"  pre-ictal AUC: median={np.median(pre):.3f}")
print(f"  ictal     AUC: median={np.median(ict):.3f}")
if len(both) >= 3:
    stat, p = wilcoxon(both["auc_ict"] - both["auc_pre"], alternative="greater")
    print(f"  paired Wilcoxon (ictal > pre-ictal): p={p:.4f}")
print("\nreading:")
print("  ictal >> pre-ictal (p<0.05)  -> SOZ decouples during the seizure")
print("                                  (ictal-specific; interictal unlikely)")
print("  ictal ~= pre-ictal           -> peripheral position is static")
print("                                  (interictal validation could work)")
print(f"\nsaved -> {HERE / 'mechanism_ictal_vs_preictal.csv'}")
print("Paste the full output back.")