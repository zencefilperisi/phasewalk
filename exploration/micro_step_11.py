"""
Micro-step 11: Scan all seizure files across ID1, ID2, ID3.

Micro-step 10 opened ONE file and looked at it visually. That is the
right first step, but before writing analysis code we need to know
whether the other 18 files (13+4+2) are structurally similar or hide
surprises.

This script visits every Sz*.mat file in the three patient folders
and reports, in a compact table:
    patient, seizure file, shape, orientation OK?, duration, n_electrodes

We are looking for:
  - Do all patients have the SAME format we saw for ID1/Sz1?
  - Are electrode counts consistent WITHIN a patient?
  - Are electrode counts DIFFERENT ACROSS patients? (they will be,
    but by how much?)
  - Any file that fails to load? Any weird shape? Any missing file?

Fill in the TODOs and run.
"""
from pathlib import Path

import numpy as np
import scipy.io as sio


# --------------------------------------------------------------------------- #
# TODO 1: set BASE to the folder that contains ID1, ID2, ID3 as subfolders.
#
# Example (Windows) — the parent folder of ID1, ID2, ID3:
#   BASE = Path(r"C:\Users\User\Desktop\SWEC")
#
# NOT the ID1 folder itself. The folder ABOVE the patient folders.
#
BASE = Path(r"C:\Users\User\Desktop\SWEC")   # <-- TODO 1
# --------------------------------------------------------------------------- #

if BASE is None:
    raise SystemExit("TODO 1 not done: set BASE.")

PATIENTS = ["ID1", "ID2", "ID3"]
FS = 512  # sampling frequency, Hz


def load_and_summarize(mat_path):
    """Return a dict summarizing one .mat file (or an error field)."""
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

    # Short-term docs: T x M with T=time, M=electrodes.
    # Electrode count is 36-100; time is hundreds of thousands.
    # So bigger axis is time.
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


# --------------------------------------------------------------------------- #
# Walk the three patient folders and collect summaries.
# --------------------------------------------------------------------------- #
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


# --------------------------------------------------------------------------- #
# Print a compact, aligned table.
# --------------------------------------------------------------------------- #
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


# --------------------------------------------------------------------------- #
# Per-patient consistency check (done for you).
# --------------------------------------------------------------------------- #
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