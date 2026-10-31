"""
ds003029 - Step 5: single-subject end-to-end integration check.

Wires the new BIDS loader to the frozen CTQW pipeline on ONE real subject,
before scaling to the whole cohort. This is an integration check, NOT a result:
one subject's number proves nothing scientifically -- it proves the chain runs
and produces sane values.

Chain:  load_ictal_window (ADR 0011/0012)
        -> compute_fcm (Pearson) -> fcm_to_graph (Layer P)
        -> time_averaged_distribution_all (CTQW, fast) -> per-electrode PR
        -> compare PR ranking to the clinician SOZ

Default subject is sub-pt2 (Engel I, has signal + SOZ). sub-pt1 and sub-jh107
are currently missing their .eeg signal files in the local download.

Place this anywhere (repo root, or exploration/ds003029/); it locates the repo
root by finding the folder that contains data/ and quantum_walks/.

Requires: mne, pandas, numpy, scipy, openpyxl
"""
import sys
from pathlib import Path

# Locate the repo root by finding the folder that actually contains
# data/loader_bids.py, searching this script's ancestors and the current
# working directory. Robust to wherever the script is placed.
def _find_repo_root():
    for p in list(Path(__file__).resolve().parents) + [Path.cwd().resolve()]:
        if (p / "data" / "loader_bids.py").exists():
            return p
    return None

_root = _find_repo_root()
if _root is None:
    sys.exit("Could not find repo root (a folder containing data/loader_bids.py).\n"
             "Run this from inside the phasewalk repo, with loader_bids.py in data/.")
sys.path.insert(0, str(_root))
print(f"[repo root: {_root}]")

import numpy as np
from scipy.stats import rankdata

from data.loader_bids import (
    load_clinical_table, discover_ictal_runs, load_ictal_window,
    match_soz, normalize_id,
)
from data.connectivity import compute_fcm, fcm_to_graph
from quantum_walks.walks import time_averaged_distribution_all, participation_ratio

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")
SUBJECT = "sub-jh105"   # Engel I, 5 seizures, 26 SOZ; fully downloaded
T_GRID = np.linspace(0.1, 100.0, 300)   # same grid as Phase 2C / SWEC


def auc_minus_pr(pr, soz_mask):
    """AUC of score=-pr predicting SOZ membership. >0.5 => SOZ more localized."""
    soz_mask = np.asarray(soz_mask, dtype=bool)
    n_pos = int(soz_mask.sum())
    n_neg = int((~soz_mask).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    r = rankdata(-pr)                       # lower PR -> higher rank
    return (r[soz_mask].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


print("=" * 70)
print(f"ds003029 integration -- {SUBJECT}")
print("=" * 70)

clinical = load_clinical_table(DATA_ROOT)
row = clinical[normalize_id(SUBJECT.replace("sub-", ""))]
soz_raw = row.get("soz_contacts")
print(f"outcome={row.get('outcome')}  engel={row.get('engel_score')}  "
      f"soz_contacts={soz_raw!r}")

runs = discover_ictal_runs(DATA_ROOT / SUBJECT)
print(f"ictal runs discovered: {len(runs)}")

per_run = []          # (run, {name: pr}, auc)
for rp in runs:
    ir = load_ictal_window(rp)
    if ir.dropped_reason:
        print(f"  {rp.run}: DROPPED ({ir.dropped_reason})")
        continue
    A = fcm_to_graph(compute_fcm(ir.signal, measure="pearson"), layer="positive")
    avg = time_averaged_distribution_all(A, T_GRID, quantum=True)
    pr = np.array([participation_ratio(avg[:, s]) for s in range(A.shape[0])])
    names = ir.good_names
    m = match_soz(soz_raw, names, ir.all_names)
    soz_set = {x.upper() for x in m.matched}
    soz_mask = np.array([n.upper() in soz_set for n in names], dtype=bool)
    auc = auc_minus_pr(pr, soz_mask)
    per_run.append((rp.run, dict(zip(names, pr)), auc))
    print(f"  {rp.run}: {len(names)} good ch, window {ir.window.duration:.1f}s"
          f"{' [short]' if ir.window.short else ''}, "
          f"SOZ {soz_mask.sum()}/{len(m.matched)+len(m.in_bad)+len(m.unmatched)} matched, "
          f"PR[min,med,max]=[{pr.min():.1f},{np.median(pr):.1f},{pr.max():.1f}], "
          f"AUC={auc:.3f}")

if not per_run:
    print("\nNo usable runs (missing signal?). Nothing to integrate.")
    sys.exit(0)

# --- across-seizure mean PR ------------------------------------------------ #
common = sorted(set.intersection(*[set(d) for _, d, _ in per_run]))
mean_pr = np.array([np.mean([d[n] for _, d, _ in per_run]) for n in common])
m = match_soz(soz_raw, common, common)
soz_set = {x.upper() for x in m.matched}
soz_mask = np.array([n.upper() in soz_set for n in common], dtype=bool)

print("\n" + "-" * 70)
print(f"across {len(per_run)} seizures, {len(common)} common channels, "
      f"{int(soz_mask.sum())} SOZ")
print(f"  per-seizure AUC : {[round(a, 3) for _, _, a in per_run]}")
print(f"  mean-PR AUC     : {auc_minus_pr(mean_pr, soz_mask):.3f}  "
      f"(>0.5 => SOZ electrodes more localized)")
print(f"  mean PR  SOZ    : {mean_pr[soz_mask].mean():.2f}")
print(f"  mean PR non-SOZ : {mean_pr[~soz_mask].mean():.2f}")
order = np.argsort(mean_pr)
rank_of = {i: r + 1 for r, i in enumerate(order)}
soz_ranks = sorted(rank_of[i] for i in np.where(soz_mask)[0])
print(f"  SOZ localization ranks (1=most localized of {len(common)}): {soz_ranks}")

print("\n" + "=" * 70)
print("Integration check done. Paste the full output back.")
print("=" * 70)