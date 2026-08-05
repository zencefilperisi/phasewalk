"""
Micro-step 11: Scan all seizure files across ID1, ID2, ID3.
"""
from pathlib import Path

import numpy as np
import scipy.io as sio


BASE = Path(r"C:\Users\User\Desktop\SWEC")

PATIENTS = ["ID1", "ID2", "ID3"]
FS = 512


def load_and_summarize(mat_path):
    try:
        mat = sio.loadmat(str(mat_path))
    except Exception as e:
        return {"error": f"loadmat failed: {e}"}

    if "EEG" not in mat:
        return {"error": "no 'EEG' variable in file"}

    EEG = mat["EEG"]
    if EEG.ndim != 2:
        return {"error": f"EEG is not 2-D: ndim={EEG.ndim}"}

    n_rows, n_cols = EEG.shape

    if n_rows > n_cols:
        n_time, n_electrodes = n_rows, n_cols
        orientation = "T x M (as documented)"
    else:
        n_time, n_electrodes = n_cols, n_rows
        orientation = "M x T (transposed — swapped)"

    duration_s = n_time / FS
    onset_sample = FS * 3 * 60
    end_sample = n_time - FS * 3 * 60
    seizure_len_s = (end_sample - onset_sample) / FS

    return {
        "shape": EEG.shape,
        "orientation": orientation,
        "n_time": n_time,
        "n_electrodes": n_electrodes,
        "duration_s": duration_s,
        "seizure_len_s": seizure_len_s,
        "dtype": str(EEG.dtype),
    }


rows = []
for pid in PATIENTS:
    folder = BASE / pid
    if not folder.is_dir():
        print(f"WARNING: folder not found: {folder}")
        continue
    sz_files = sorted(folder.glob("Sz*.mat"))
    if not sz_files:
        print(f"WARNING: no Sz*.mat files in {folder}")
        continue
    for f in sz_files:
        summary = load_and_summarize(f)
        rows.append({"patient": pid, "file": f.name, **summary})


print()
print(f"{'patient':<8} {'file':<10} {'shape':<20} {'n_elec':>7} "
      f"{'dur (s)':>9} {'seiz (s)':>9}  {'notes'}")
print("-" * 90)

for r in rows:
    if "error" in r:
        print(f"{r['patient']:<8} {r['file']:<10} ERROR: {r['error']}")
        continue
    notes = "" if "as documented" in r["orientation"] else "TRANSPOSED"
    print(f"{r['patient']:<8} {r['file']:<10} "
          f"{str(r['shape']):<20} {r['n_electrodes']:>7} "
          f"{r['duration_s']:>9.1f} {r['seizure_len_s']:>9.1f}  {notes}")


print()
print("Per-patient consistency:")
for pid in PATIENTS:
    p_rows = [r for r in rows if r.get("patient") == pid and "error" not in r]
    if not p_rows:
        continue
    n_el = set(r["n_electrodes"] for r in p_rows)
    print(f"  {pid}: {len(p_rows)} seizure(s); electrode counts = {n_el}")
    if len(n_el) > 1:
        print(f"    WARNING: electrode count varies within {pid} — investigate.")
