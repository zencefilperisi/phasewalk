"""
ds003029 - Cohort build (the freeze step for the ADR 0013 pre-registration).

Enumerates the usable cohort and freezes each subject's common good-channel set
and SOZ reference to cohort_frozen.json, BEFORE any AUC is computed. Run once,
inspect the printed table, and commit cohort_frozen.json: the analysis script
reads this file, so the cohort and SOZ sets cannot drift after results are seen.

Usable (ADR 0011/0012): a matched clinical row, >= 1 COMPLETE ictal run
(header + .eeg + channels + events), and >= 1 SOZ contact matched to a good
channel. Excludes pt01 (duplicate). Subjects without a clinical row / signal /
matched SOZ are reported but omitted.

Requires: pandas, numpy  (no mne, no CTQW -- fast)
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

import mne
import pandas as pd
from data.loader_bids import (
    load_clinical_table, discover_ictal_runs, match_soz, normalize_id,
    parse_seizure_window, WINDOW_SECONDS, EXCLUDE_SUBJECTS,
)

mne.set_log_level("ERROR")


def window_fits_signal(rp):
    """True if the ADR-0012 ictal window falls within the published signal.
    Excludes runs (e.g. umf002-005) whose events annotate a seizure beyond the
    end of the trimmed BrainVision file. Reads header only (preload=False)."""
    try:
        ev = pd.read_csv(rp.events_tsv, sep="\t")
        w = parse_seizure_window(ev, WINDOW_SECONDS)
        if w is None:
            return False
        raw = mne.io.read_raw_brainvision(str(rp.vhdr), preload=False, verbose="ERROR")
        dur = raw.n_times / raw.info["sfreq"]
        return (w.win_start + 1.0) <= dur
    except Exception:
        return False

DATA_ROOT = Path(r"C:\Users\User\Desktop\ds003029")
OUT = Path(__file__).resolve().parent / "cohort_frozen.json"


def engel_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def good_and_all(channels_tsv):
    c = pd.read_csv(channels_tsv, sep="\t")
    good = c[(c["type"].isin(["ECOG", "SEEG"])) & (c["status"] == "good")]
    return good["name"].astype(str).tolist(), c["name"].astype(str).tolist()


clinical = load_clinical_table(DATA_ROOT)
subjects = sorted(p for p in DATA_ROOT.glob("sub-*") if p.is_dir())

cohort = {}
skipped = []
for sdir in subjects:
    label = sdir.name.replace("sub-", "")
    if label in EXCLUDE_SUBJECTS:
        skipped.append((label, "excluded-duplicate"))
        continue
    row = clinical.get(normalize_id(label))
    if row is None:
        skipped.append((label, "no-clinical-row"))
        continue
    runs = discover_ictal_runs(sdir)
    if not runs:
        skipped.append((label, "no-complete-ictal-signal"))
        continue
    runs = [rp for rp in runs if window_fits_signal(rp)]
    if not runs:
        skipped.append((label, "signal-shorter-than-seizure"))
        continue
    per_run_good, all_names = [], set()
    for rp in runs:
        g, a = good_and_all(rp.channels_tsv)
        per_run_good.append(set(g))
        all_names.update(a)
    common = sorted(set.intersection(*per_run_good))
    m = match_soz(row.get("soz_contacts"), common, sorted(all_names))
    if not m.matched:
        skipped.append((label, "no-soz-in-good"))
        continue
    cohort[label] = dict(
        center=str(row.get("clinical_center")),
        outcome=str(row.get("outcome")),
        engel=engel_int(row.get("engel_score")),
        runs=[rp.run for rp in runs],
        n_runs=len(runs),
        common_channels=common,
        n_good=len(common),
        soz=list(m.matched),
        n_soz=len(m.matched),
        soz_in_bad=len(m.in_bad),
        soz_unmatched=len(m.unmatched),
        soz_fraction=round(len(m.matched) / len(common), 3),
    )

OUT.write_text(json.dumps(cohort, indent=2))

# ---- report --------------------------------------------------------------- #
print(f"{'subject':<10}{'center':>7}{'eng':>4}{'runs':>5}{'good':>5}"
      f"{'nSOZ':>5}{'bad':>4}{'frac':>7}")
print("-" * 50)
for lab, info in cohort.items():
    print(f"{lab:<10}{info['center']:>7}{str(info['engel']):>4}{info['n_runs']:>5}"
          f"{info['n_good']:>5}{info['n_soz']:>5}{info['soz_in_bad']:>4}"
          f"{info['soz_fraction']:>7.2f}")

engel_i = [l for l, i in cohort.items() if i["engel"] == 1]
print("\n" + "=" * 60)
print(f"USABLE cohort: {len(cohort)} subjects")
print(f"  Engel I (seizure-free): {len(engel_i)}  -> {engel_i}")
from collections import Counter
print(f"  Engel distribution: {dict(Counter(i['engel'] for i in cohort.values()))}")
print(f"\nskipped ({len(skipped)}):")
for lab, why in skipped:
    print(f"  {lab}: {why}")
print(f"\nfrozen -> {OUT}")
print("Inspect the table, then run the analysis script.")